## Imports
import random
import pandas as pd
import numpy as np
from tqdm import tqdm
import plotly.graph_objects as go

## Variables

big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv", sep='\t',
                     usecols=['prot_id','org_id','species','phylum','superkingdom','upid'])
domdf = pd.read_pickle("data/datasets/huge_proteins_pfam.pkl")

N_MAX_PROTEINS_THRESHOLD = 2000
random.seed(1234)

## Helper Functions
def expand_list_values_to_df(pdseries):
    s = set()
    n = len(pdseries)
    for x in pdseries.iloc[:n]:
        [s.add(z) for z in x]
        
    return s

def count_bin_number_of_domains(maindf,coluna='hit_keys'):
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
            z_df.loc[x,hk] += 1
    
    return pd.DataFrame(z_df.sum()).rename(columns={0:'total_found'}).sort_values('total_found',ascending=False).reset_index()

def plot(df,column,taxonomy):
    x_domains = df[column].to_list()
    y_total = df['total_found'].to_list()
    
    fig = go.Bar(name=taxonomy, x=x_domains, y=y_total)

    return fig

def plot_stacked(list_of_go_figures, tax_level):
    fig = go.Figure(data=list_of_go_figures)
    # Change the bar mode
    fig.update_layout(barmode='stack',xaxis_range=[-0.5,25.5],title=f"Most Common Domain in {tax_level.capitalize()} (n_big_proteins<={N_MAX_PROTEINS_THRESHOLD})")
    fig.update_xaxes(categoryorder='total descending')
    fig.write_html(f"data/results/single_domain_superkingdom/domain_single_most_common_stacked_{tax_level}_normbp.html")

def do_all_tax_level(level,subdivision="Bacteria"):
    exclude = ['candidatus thermoplasmatota','candidatus altiarchaeota','candidatus lokiarchaeota','candidate division zixibacteria',np.nan,
                'candidatus absconditabacteria','candidatus omnitrophica','candidatus margulisbacteria','candidatus cloacimonetes','haptista']
    if level == "superkingdom":
        uniques = big_df[level].unique()
    else:
        uniques = big_df[big_df['superkingdom'] == subdivision][level].unique()
    print(uniques)
    go_figures = []
    for tax in uniques:
        if tax not in exclude:
            print(level, tax)
            result_fig = main(domdf,analysis=str(tax),level=level)
            go_figures.append(result_fig)
    plot_stacked(list_of_go_figures=go_figures, tax_level=subdivision)


def main(df:pd.DataFrame, column="hit_keys", analysis="all", level='superkingdom', balance=True):
    if analysis == "all":
        singledomdf = count_bin_number_of_domains(df,coluna=column)
    elif analysis in big_df[level].unique():
        to_use_df = df.sample(n=N_MAX_PROTEINS_THRESHOLD)
        singledomdf = count_bin_number_of_domains(to_use_df,coluna=column)
    
    fig = plot(singledomdf,column=singledomdf.columns[0],taxonomy=analysis)
    return fig

if __name__ == '__main__':
    do_all_tax_level('phylum','Eukaryota')
    do_all_tax_level('phylum','Bacteria')
    do_all_tax_level('phylum','Archaea')
