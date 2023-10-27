## Imports
import os
from sys import argv
import string
import random
import pandas as pd
import numpy as np
from tqdm import tqdm
import plotly.express as px
from draw_architecture import main_draw_controller

## Variables

big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv",sep='\t')
domdf = pd.read_pickle("data/datasets/huge_proteins_pfam.pkl")

def expand_list_values_to_df(pdseries):
    s = set()
    n = len(pdseries)
    for x in pdseries.iloc[:n]:
        [s.add(z) for z in x]
        
    return s

def count_bin_number_of_domains(maindf,coluna='pfam_acession'):
    l_of_pids = maindf.prot_id.to_list()
    columns = expand_list_values_to_df(maindf[coluna])
    shape=(len(l_of_pids),len(columns))

    z_df= pd.DataFrame(np.zeros(shape),
                   columns=columns,
                   index=l_of_pids)

    t_df = maindf.copy()

    t_df.set_index('prot_id',inplace=True)
    for x in tqdm(l_of_pids):
        for hk in t_df.at[x,coluna]:
            z_df.at[x,hk] += 1
    
    return (z_df,pd.DataFrame(z_df.sum()).rename(columns={0:'total_found'}).sort_values('total_found',ascending=False))

def count_full_number_of_domains(maindf,coluna_main='pfam_acession',ocurrence='ocurrence',coords='coord_on_prot'):
    l_of_pids = maindf.prot_id.to_list()
    columns = expand_list_values_to_df(maindf[coluna_main])
    shape=(len(l_of_pids),len(columns))

    z_df= pd.DataFrame(np.zeros(shape),
                   columns=columns,
                   index=l_of_pids)
    

    t_df = maindf.copy()

    t_df.set_index('prot_id',inplace=True)
    for x in tqdm(l_of_pids):
        for hk in t_df.at[x,coluna_main]:
            hk_index = t_df.at[x,coluna_main].index(hk)
            dom_occurrence = t_df.at[x,ocurrence][hk_index]
            z_df.at[x,hk] += dom_occurrence
    
    return (z_df,pd.DataFrame(z_df.sum()).rename(columns={0:'total_found'}).sort_values('total_found',ascending=False))

def deduce_dompattern(tempdf,acess='pfam_acession',occur='ocurrence',coord='coord_on_prot'):
    alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z','δ','λ','μ','π', 'Σ','τ', 'Φ', 'Ψ', 'ω','ζ','1','2','3','4','5','6','7','8','9']
    alphabet = [x.lower() for x in alphabet]
    
    x,y,z = tempdf[[acess,occur,coord]]
    dom2occur = list(zip(x,y))
        
    pat_keys = alphabet[:len(x)+1]
    l = []
    for i in dom2occur:
        for _ in range(i[1]):
            l.append(pat_keys[x.index(i[0])])
    
    final_pattern = [(x,y[0],y[1]) for x,y in list(zip(l,z))]
    final_pattern = sorted(final_pattern, key=lambda x: (x[1], -x[2]))
    final_pattern = "".join([letter[0] for letter in final_pattern])
    
    pat2dom = dict(zip(pat_keys,x))
    
    return pat2dom,final_pattern

def iter_dompattern(maindf):
    l_of_patterns = []
    for p in maindf.index:
         l_of_patterns.append(deduce_dompattern(maindf.iloc[p])[1])
            
    return tuple(l_of_patterns)

def compare_archs(df:pd.DataFrame):
    archs = df.copy()[['comp','pattern']].reset_index()

    archs['comp'] = archs['comp'].apply(tuple)
    groups = archs.groupby('comp')
    
    d = {}
    c = 0
    # loop through each row of the dataframe
    for _, group in tqdm(groups):
        patterns = group.groupby('pattern')
        for _, pattern in patterns:
            same_archs = pattern['prot_id'].to_list()
            d[c] = same_archs
            c += 1
        
    compare_archdf = pd.DataFrame(pd.Series(d),columns=['Same_Arch'])
    compare_archdf['count_same'] = compare_archdf.Same_Arch.apply(lambda x: len(x))
    compare_archdf = compare_archdf.sort_values(by='count_same',ascending=False)
    return compare_archdf

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def get_arch_code(protein_id, df_):
    if isinstance(protein_id,list):
        protein_id = protein_id[0]    
    return df_[df_.prot_id == protein_id]['arch_id'].iloc[0]

def flatten(l, repr=False):
    if not repr:
        return [item for sublist in l for item in sublist]
    else:
        return [sublist[0] for sublist in l]

## Processing/Analysis

def analyse(analysis="all",level='superkingdom'):
    cp_ar_df_save_to = f"data/results/architectures/huge_proteins_architectures_{analysis}.tsv"
    if not os.path.exists(cp_ar_df_save_to):
        
        archdf = pd.DataFrame(columns= ['comp','counts','pattern'], index=domdf.prot_id) ## Defined empty df
        archdf['comp'] = list(domdf.pipe(count_bin_number_of_domains)[0].to_numpy())
        archdf['counts'] = list(domdf.pipe(count_full_number_of_domains)[0].to_numpy())
        archdf['pattern'] = domdf.pipe(iter_dompattern)

        simpledf = big_df[['org_id','superkingdom','prot_id','org_id']]
        arch_org_df = archdf.merge(simpledf,left_index=True,right_on='prot_id').set_index('prot_id')

        if analysis in big_df[level].unique():
            ### Filter the working dataset to only contain the wanted organisms
            list_of_prots = big_df.query(f"{level} == '{analysis}'").prot_id.to_list()
            to_use_df = arch_org_df.filter(items=list_of_prots, axis=0)
        
        elif analysis == "all":
            to_use_df = arch_org_df

        ### Run the comparison of protein architechtures // carefull it takes a long time...
        compare_archdf = compare_archs(to_use_df)
        
        ## Assign the correct mapping of arch_code // logic to make consistent across analysis
        if analysis == "all":
            if not os.path.exists("data/outputs/mapping/code2arch.tsv"):
                ### Add and save the information about the "Arch ID code" created
                compare_archdf['arch_id'] = compare_archdf.Same_Arch.apply(lambda x:"".join(x)).apply(lambda x: id_generator(6))
                compare_archdf.set_index('arch_id').explode('Same_Arch').rename(columns={"Same_Arch":"prot_id"}).to_csv("data/outputs/mapping/code2arch.tsv", sep='\t')
            else:
                arch_codes_df = pd.read_csv("data/outputs/mapping/code2arch.tsv",sep='\t', usecols=['arch_id','prot_id'])
                compare_archdf['arch_id'] = compare_archdf.Same_Arch.apply(lambda x: get_arch_code(x,arch_codes_df))
        
        else:
            arch_codes_df = pd.read_csv("data/outputs/mapping/code2arch.tsv",sep='\t', usecols=['arch_id','prot_id'])
            compare_archdf['arch_id'] = compare_archdf.Same_Arch.apply(lambda x: get_arch_code(x,arch_codes_df))

        compare_archdf.to_pickle(cp_ar_df_save_to)
    else:
        compare_archdf = pd.read_pickle(cp_ar_df_save_to)

    ## Draw representative architectures
    print(f"Runned arch comparison for {analysis}, drawing representatives ...")
    list_of_architectures = flatten(compare_archdf.Same_Arch.to_list(), repr=True)
    main_draw_controller(list_of_architectures,draw_archs=True,save_dir=f"data/results/architectures/{analysis}")

    ## Plot
    print(f"Ploting {analysis}")
    fig = px.bar(compare_archdf,x='arch_id',y='count_same',
             range_x=(-1,20.5),
             template='ggplot2',
             labels={'arch_id':'Domain Architecture','count_same':'Frequency'})

    fig.update_layout(
        font=dict(
            size=13
        ),
        title=f"Big Proteins architecture in {analysis}",
        width=1000,
        height=700
    )

    fig.write_html(f"data/results/architectures/plots/domain_arch_most_common_{analysis}.html")
    
def main(analysis_type):
    print(f"Starting {analysis_type}")
    if analysis_type == 'all':
        print(f"Doing {analysis_type} ...")
        analyse()
    elif analysis_type == 'euka':
        print(f"Doing {analysis_type} ...")
        analyse("Eukaryota")
    elif analysis_type == 'bac':
        print(f"Doing {analysis_type} ...")
        analyse("Bacteria")
    elif analysis_type == 'arc':
        print(f"Doing {analysis_type} ...")
        analyse("Archaea")
    elif analysis_type == 'pvc':
        print(f"Doing {analysis_type} ...")
        analyse(analysis="planctomycetes", level="phylum")
        analyse(analysis="verrucomicrobia", level="phylum")
        analyse(analysis="chlamydiae", level="phylum")
    elif analysis_type == 'int_bac':
        print(f"Doing {analysis_type} ...")
        analyse(analysis="actinobacteria", level="phylum")
        analyse(analysis="proteobacteria", level="phylum")
        analyse(analysis="bacteroidetes", level="phylum")
        analyse(analysis="firmicutes", level="phylum")

if __name__ == "__main__":
    main('euka')
    main('bac')
    main('arc')
