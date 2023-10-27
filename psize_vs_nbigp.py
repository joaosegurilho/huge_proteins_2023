
import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from ete3 import NCBITaxa

import warnings
warnings.filterwarnings("ignore")

## Variables

big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv", sep='\t',
                     usecols=['prot_id','org_id','species','phylum','superkingdom','upid','assembly_level'])

ncbi = NCBITaxa()

def filter_dataframe(df, clade):
    exclude = [
        "candidatus thermoplasmatota",
        "candidatus altiarchaeota",
        "candidatus lokiarchaeota",
        "candidatus absconditabacteria",
        "candidatus omnitrophica",
        "candidatus margulisbacteria",
        "candidatus cloacimonetes",
        "candidate division zixibacteria"
    ]

    df = df[(df.phylum.notna()) & (~df.phylum.isin(exclude))]

    conditions = (df.superkingdom == clade)
    to_use_df = df[conditions]

    if clade in ("bacteria", "eukaryota"):
        to_use_df = to_use_df[(~to_use_df.phylum.str.startswith("candidat"))]
    
    to_drop = remove_small_proteomes(to_use_df,clade)
    
    to_use_df = to_use_df.set_index('phylum').drop(labels=to_drop).reset_index()

    return to_use_df

def remove_small_proteomes(to_use_df,clade):
    grouped_cnt = to_use_df.groupby('phylum').count()
    to_drop = []
    for tax in to_use_df.phylum.unique():
        if clade == 'archaea':
            break
        else:
            if grouped_cnt.loc[tax, 'species'] <= 5:
                to_drop.append(tax)
            else:
                continue
    return to_drop

def retrive_taxonomy(org):
    org = int(org)
    lineage = ncbi.get_lineage(org)
    trans = {k:v for k,v in ncbi.get_taxid_translator(lineage).items()}
    ranks = {v:k for k,v in ncbi.get_rank(lineage).items()}

    if 'species' not in ranks.keys():
        ranks['species'] = np.nan
        trans[ranks['species']] = np.nan
    if 'phylum' not in ranks.keys():
        ranks['phylum'] = np.nan
        trans[ranks['phylum']] = np.nan
    return trans[ranks['species']], trans[ranks['phylum']], trans[ranks['superkingdom']]

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_bigp_info():
    list_of_subsets = []
    for chunk in os.listdir("data/outputs/uniprot_proteome_size_for_orgs/"):
        if not chunk.split('_')[1] == 'OR':
            continue
        chunk_df = pd.read_csv(os.path.join("data/outputs/uniprot_proteome_size_for_orgs/", chunk), sep='\t')
        list_of_subsets.append(chunk_df)

    bigdf = (
        pd.concat(list_of_subsets, ignore_index=True)
        .rename(columns={"Proteome Id":"upid","Organism Id":"org_id","Protein count":"proteome_lenght"})
        .sort_values('proteome_lenght', ascending=False).drop_duplicates(subset='org_id')
        .drop(columns=['org_id'])
    )

    return bigdf

def get_non_bigp_info():

    nonbigdf = (
        pd.read_csv("data/outputs/uniprot_proteome_size_for_orgs/proteomes_nonbigp_ref_24_03_2023.tsv", sep='\t')
        .rename(columns={"Proteome Id":"upid","Organism Id":"org_id","Protein count":"proteome_lenght"})
        .sort_values('proteome_lenght', ascending=False).drop_duplicates(subset='org_id')
    )

    nonbigdf['n_bigp'] = pd.Series(np.zeros(len(nonbigdf)))

    nonbigdf['species'] = pd.Series(np.full(len(nonbigdf),""))
    nonbigdf['phylum'] = pd.Series(np.full(len(nonbigdf),""))
    nonbigdf['superkingdom'] = pd.Series(np.full(len(nonbigdf),""))
    
    for row, data in nonbigdf.iterrows():
        sp, phy, sk = retrive_taxonomy(data.org_id)
        nonbigdf.at[row,'species'] = sp
        nonbigdf.at[row,'phylum'] = phy
        nonbigdf.at[row,'superkingdom'] = sk

    return nonbigdf

def plot_all(data_f, filterd=True):
    data_f = data_f.sort_values('superkingdom')
    fig = px.scatter(data_f, x='proteome_lenght', y='n_bigp',
                    color="superkingdom",
                    log_x=True,
                    # log_y=True,
                    # trendline='ols',
                    # trendline_color_override='gray',
                    # trendline_options=dict(log_x=True),
                    # trendline_scope = 'overall',
                    hover_name="org_id",
                    hover_data=["proteome_lenght", "n_bigp", "phylum", "species"],
                    labels={
                     "n_bigp": "N. Huge Proteins",
                     "proteome_lenght": "Size of Proteome"}
                    )
    
    fig.update_traces(marker_size=4)

    fig.update_layout(legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            font=dict(
                size=13
            ),
            itemsizing='constant'
        ))
    
    fig.data = fig.data[::-1]
    if not filterd: 
        fig.write_html(f"data/results/psize_vs_bigp_n/proteome_size_vs_number_of_bigp_all_unfilterd.html")
    else:
        fig.write_html(f"data/results/psize_vs_bigp_n/proteome_size_vs_number_of_bigp_all_{filterd}.html")

def plot_all_boxplot(data_f, filterd=True):

    colors = {
        'Archaea':('#87cbff','#5d707f'),
        'Bacteria':('#e5737e','#996b6f'),
        'Eukaryota':('#8bffa3','#6f9978')
    }

    list_of_boxs = []
    data_f = data_f.sort_values('superkingdom')
    for sk in data_f['superkingdom'].unique():
        sk_color = colors[sk]
        grouped_phyla = data_f[data_f.superkingdom == sk].groupby('phylum').count()
        grouped_phyla['median'] = data_f[data_f.superkingdom == sk].groupby('phylum').median()['bigp_probability']
        grouped_phyla['max'] = data_f[data_f.superkingdom == sk].groupby('phylum').max()['bigp_probability']
        grouped_phyla['q3'] = data_f[data_f.superkingdom == sk].groupby('phylum').quantile(q=0.75)['bigp_probability']
        sorted_phyla = grouped_phyla.sort_values(['median','q3','max'], ascending=True)
        
        for phy in sorted_phyla.index:
            if phy.startswith("candidat") and sk != 'archaea':
                continue
            elif len(data_f[data_f.phylum == phy].org_id.unique()) < 5:
                continue
            
            list_of_probs = data_f[data_f.phylum == phy]['bigp_probability'].to_list()
            temp_box = go.Box(
                y=list_of_probs,
                name=phy,
                marker_color=sk_color[0],
                marker_size=3.5,
                boxpoints='all', # no data points
                line_color=sk_color[1],
                line_width=0.8,
                )
            list_of_boxs.append(temp_box)
    fig = go.Figure(list_of_boxs)


    fig.update_layout(template="plotly_white",
                        showlegend=False,
                        xaxis_title="Phylum",
                        yaxis_title="N. Huge Proteins / Size of Proteome (%)",
                        yaxis_tickfont_size=15,
                        yaxis_gridwidth=0.8,
                        yaxis_gridcolor='LightGray'
                        )

    if not filterd:
        fig.write_html(f"data/results/psize_vs_bigp_n/boxplot_proteome_size_vs_number_of_bigp_all_phyla_and_sk_unfiltered.html")
    else:
        fig.write_html(f"data/results/psize_vs_bigp_n/boxplot_proteome_size_vs_number_of_bigp_all_phyla_and_sk_{filterd}.html")

def query_uniprot_other_proteins(orgs, inclusive=True):  
    from .utils import query_uniprot
    logic = 'OR' if inclusive else 'NOT'
    for i, chunk in enumerate(chunks(orgs, 500)):
        query = f"%29%20{logic}%20%28".join([f"organism_id%3A{x}" for x in chunk])
        if not inclusive: query = f"%29%20{logic}%20%28" + query
        url = f'https://rest.uniprot.org/proteomes/search?fields=upid%2Corganism_id%2Cprotein_count&format=tsv&query=%28%28{query}%29%29&size=500'
        # print(url)
        outpath = f"data/outputs/uniprot_proteome_size_for_orgs/proteomes_{logic}_{i}.tsv"
        query_uniprot.run_query(outpath, url)

def remove_outlier_proteomes(df: pd.DataFrame):
    avg_df = df.groupby('phylum').mean()['proteome_lenght']
    std_df = df.groupby('phylum').std()['proteome_lenght']
    labels_to_drop = []
    for phy in df.phylum.dropna().unique():
        phy_df = df[df.phylum == phy]
        mean = avg_df.loc[phy]
        std = std_df.loc[phy]
        labels = phy_df[phy_df.proteome_lenght < mean-std].index # filter based on mean+std proteome len
        labels = list(labels) + list(phy_df.index) if len(phy_df) < 5 else labels # filter based on number of proteomes in the phyla
        for label in labels:
            labels_to_drop.append(label)
    return df.drop(labels=labels_to_drop)

def remove_plasmids_proteomes(df: pd.DataFrame):
    labels_to_drop = df[df.assembly_level.str.startswith('Plasmid',na=False)].index
    return df.drop(labels=labels_to_drop)

def main():
    df = big_df.copy()
    grouped = df.groupby('upid').count()['prot_id']
    df['n_bigp'] = df.upid.apply(lambda x: grouped.loc[x])
    df['protein_type'] = df.prot_id.apply(lambda _: 'big_protein')
    df = df.drop(columns=['prot_id'])

    orgs_toquery = list(df.org_id.unique())
    if not os.path.exists("data/outputs/uniprot_proteome_size_for_orgs/proteomes_OR_0.tsv"):
        query_uniprot_other_proteins(orgs_toquery)
        df = df.merge(get_bigp_info(),on='upid')
    else:
        df = df.merge(get_bigp_info(),on='upid')
    df = remove_outlier_proteomes(df)

    nonbigdf = get_non_bigp_info()
    nonbigdf = nonbigdf[~nonbigdf.org_id.isin(orgs_toquery)]
    nonbigdf['protein_type'] = nonbigdf.org_id.apply(lambda _: 'non_big_protein')
    nonbigdf = remove_outlier_proteomes(nonbigdf)
    
    all_df = pd.concat([df,nonbigdf])
    ## merge similar named phyla
    all_df.phylum.replace('Planctomycetota','Planctomycetes',inplace=True)
    all_df.phylum.replace('Bacteroidota','Bacteroidetes',inplace=True)
    all_df['bigp_probability'] = (all_df['n_bigp'] / all_df['proteome_lenght']) * 100
    all_df = all_df.drop_duplicates(subset='org_id')

    print(all_df['bigp_probability'].describe())
    print(all_df[all_df.phylum == 'Elusimicrobia'].sort_values('bigp_probability', ascending=False))
    print(all_df[all_df['bigp_probability'] > 1.5])

    plot_all(all_df, False)
    plot_all_boxplot(all_df, False)

    plot_all(remove_plasmids_proteomes(all_df), "filter_palsmids")
    plot_all_boxplot(remove_plasmids_proteomes(all_df), "filter_palsmids")
    
    plot_all(remove_outlier_proteomes(all_df),"filter_outliers")
    plot_all_boxplot(remove_outlier_proteomes(all_df), "filter_outliers")

    plot_all(remove_outlier_proteomes(remove_plasmids_proteomes(all_df)),"filter_plasmids_and_outliers")
    plot_all_boxplot(remove_outlier_proteomes(remove_plasmids_proteomes(all_df)), "filter_plasmids_and_outliers")

if __name__ == "__main__":
    # clade_name = sys.argv[1]
    main()

    # desired_proteomes = list(main_df.org_id.astype('str').unique())
    # # print(desired_proteomes)
    # # retrive_info(desired_proteomes[:3])
    # d = get_non_bigp_info(desired_proteomes)
    # print(d)
    # print(d.columns)