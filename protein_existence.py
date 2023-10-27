from concurrent.futures import ThreadPoolExecutor
import os
import requests
from tqdm import tqdm
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from time import sleep


big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv", sep='\t',
                     usecols=['prot_id','org_id','prot_exist','phylum','superkingdom','upid'])
proteome_df = pd.read_csv("data/datasets/uniprotkb-proteomes-busco.tsv",sep='\t',
                          usecols=['Proteome Id','Protein count']).rename(columns={'Proteome Id':'upid','Protein count':'proteome_lenght'})
big_df = big_df.merge(proteome_df, on='upid').drop_duplicates(ignore_index=True)

SAMPLE_TIMES = 10

def get_protein_exist(prot):
    sleep(0.2)

    url = "https://www.ebi.ac.uk/proteins/api/proteins/"
    query = prot.split('_')[0]

    try:
        r = requests.get(url+query, headers={ "Accept" : "application/json"})

    except:
        print(f"Request {prot} not working!")
        return (prot, np.nan)

    else:
        if r.ok:
            return (prot, r.json()['proteinExistence'])


def query_annot_types(prot_list):
    dc = []
    with tqdm(total=len(prot_list), desc="query") as pbar:
        with ThreadPoolExecutor(max_workers=30) as ex:
            
            futures = ex.map(get_protein_exist, prot_list)
            
            for future in futures:
                if future != None:
                    dc.append(future)
                pbar.update(1)

    return {v[0]:v[1] for v in dc}

def plot_plotly_percent(df:pd.DataFrame, nonbig_dict:dict):
    print(df.shape)

    big_values = list(reversed(df.groupby('prot_exist').count().prot_id.to_list() + [0]))
    total_bigps = sum(big_values)
    big_values = [(x/total_bigps)*100 for x in big_values] 
    categories = list(existence_types.values())
    
    fig = go.Figure()
    fig.add_trace(go.Bar(
        name=f'Huge Proteins total:{total_bigps}',
        y=categories, x=big_values,
        opacity=0.75,
        orientation='h'
    ))

    fig.data[0].marker.color = "#5C7EAB"
    fig.data[0].marker.line = dict(width=1, color="#7f7f7f")

    total_nonbigps = sum(nonbig_dict.values())
    nonbig_vals = list(reversed([(x/total_nonbigps)*100 for x in nonbig_dict.values()]))
    print(categories)
    print(big_values)
    print(nonbig_vals)

    fig.add_trace(go.Bar(
        name=f'Non-Huge Proteins total:{total_nonbigps}',
        y=categories, x=nonbig_vals,
        orientation='h',
        opacity=0.75,
    ))

    fig.data[1].marker.color = "#FF7D63"
    fig.data[1].marker.line = dict(width=1, color="#7f7f7f")
    fig.update_layout(
        xaxis_title="% of Dataset",
        barmode='group',
        template='plotly_white',
        title=f"Protein Existence Percent of dataset"
        )
    
    fig.write_html("data/results/protein_existence/protein_existence_level_big_protein_percent.html")


def plot_plotly(df:pd.DataFrame, nondf:pd.DataFrame):
    
    big_values = list(reversed(df.groupby('prot_exist').count().prot_id.to_list() + [0])) 
    categories = list(existence_types.values())
    
    fig = go.Figure()
    fig.add_trace(go.Bar(
        name='Big Proteins',
        y=categories, x=big_values,
        opacity=0.75,
        orientation='h'
    ))
    fig.data[0].marker.color = "#5C7EAB"
    fig.data[0].marker.line = dict(width=1, color="#7f7f7f")

    nonbig_grouped = nondf.groupby(['prot_exist', 'sample']).count()

    for pe, s in nonbig_grouped.index:
        if not ("5: Uncertain",s) in nonbig_grouped.index:
            nonbig_grouped.loc[("5: Uncertain",s),'obs'] = 0
    
    nonbig_avg = []
    nonbig_std = []
    for cat in categories:
        temp = nonbig_grouped.loc[cat]
        nonbig_avg.append(temp.mean().obs)
        nonbig_std.append(temp.std().obs)

    fig.add_trace(go.Bar(
        name='Non-Big Proteins',
        y=categories, x=nonbig_avg,
        orientation='h',
        opacity=0.75,
        error_x=dict(type='data', array=nonbig_std)
    ))

    fig.data[1].marker.color = "#FF7D63"
    fig.data[1].marker.line = dict(width=1, color="#7f7f7f")
    fig.update_layout(barmode='group', template='plotly_white', title=f"Protein Existence (sampled {SAMPLE_TIMES}x)")
        
    fig.write_html("data/results/protein_existence/protein_existence_level_big_protein.html")

def get_nonbigp_existence():
    # Numbers obtained by observing the numbers obtained at the "Protein Existence" section after filtering in the Uniprot web portal
    uniprot_protein_exist_proportion_nonbig = {
        '1: Evidence at protein level': 305_136,
        '2: Evidence at transcript level': 1_333_718,
        '3: Inferred from homology': 74_578_916,
        '4: Predicted': 164_851_138,
        '5: Uncertain': 1_823
        }
    
    if not os.path.exists("data/datasets/uniprotkb-all_exist.tsv"):

        ## Creates a dataset of shape of len equeal to the specified numbers above filled with the keys of the dict
        nonbigp = [pd.DataFrame(np.full((v), k), columns=['prot_exist']) for k, v in uniprot_protein_exist_proportion_nonbig.items()]
        nonbigp = pd.concat(nonbigp, axis=0, ignore_index=True)
        nonbigp['obs'] = np.full_like(nonbigp, '1')
        
    else:
        nonbigp = pd.read_csv("data/datasets/uniprotkb-all_exist.tsv",sep='\t')
    return nonbigp, uniprot_protein_exist_proportion_nonbig

def main():
    df = big_df.copy()
    global existence_types
    existence_types = {'Uncertain':'5: Uncertain','Predicted':'4: Predicted','Inferred from homology':'3: Inferred from homology','Evidence at transcript level':'2: Evidence at transcript level', 'Evidence at protein level' :'1: Evidence at protein level'}

    ## Get annot level for bigP - if NaN assign '4:Predicted'
    print("Doing BigP")
    print(df.shape)
    df['prot_exist'] = df['prot_exist'].map(existence_types)
    print(df.shape)

    ## Get annot level for non-bigP - sample with same size, restrict to the same organisms
    print("Doing Non_BigP")
    nonbigp, nonbigp_dataset  = get_nonbigp_existence()

    plot_plotly_percent(df, nonbigp_dataset)


    ## get the average and std for all
    list_of_idfs = []
    for i in tqdm(range(1,SAMPLE_TIMES+1), desc="sampling"):
        sample = nonbigp.sample(n=len(df), replace=False)
        sample['sample'] = np.full(len(sample), f"sample_{i}")
        list_of_idfs.append(sample)

    non_df = pd.concat(list_of_idfs)
    print(non_df)

    ## plot facet for bigP and NonBigP as before
    print("Plotting")
    plot_plotly(df, non_df)

if __name__ == "__main__":
    main()
