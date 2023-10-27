
import warnings
warnings.filterwarnings("ignore")

## Imports
from concurrent.futures import ThreadPoolExecutor
from os import path
import os
import subprocess
import pandas as pd
from tqdm import tqdm
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from sqlalchemy import create_engine

from ete3 import NCBITaxa
ncbi = NCBITaxa()

big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv",sep='\t',
                     usecols=['prot_id', 'prot_len','superkingdom'])

sigdf = pd.read_csv("data/datasets/huge_proteins_sig_tmh_dis.tsv", sep='\t',
                    usecols=['prot_id','signalp','tmh','disorder'])
                    # usecols=["NAME / ID","Signal_Peptide","TMHs","Disorder % [Globular Domains]"],

mergedf = (big_df
    .merge(sigdf, on='prot_id')
    .assign(
        dis_percent=lambda x: x['disorder'].str.split(' % ').apply(lambda y: y[0]).astype('float'),
        tmh_topo=lambda x: x['tmh'].str.split(r' \[').apply(lambda y: y[1].rstrip(']')).astype('object'),
        tmh_cnt=lambda x: x['tmh'].str.split(r' \[').apply(lambda y: y[0]).astype('int'),
        signalp_exist=lambda x: x['signalp'].apply(lambda y: True if y != '-' else False),
    )
)

pio.templates.default = "plotly_white"
SAMPLE_TIMES = 5

db_path="sqlite:///data/datasets/all_data.db"

def query_sql(prots_to_query):
    disk_engine = create_engine(db_path)
    con = disk_engine.connect()

    l = "','".join(prots_to_query)
    l_of_ids = f"('{l}')"
            
    non_big_prots = pd.read_sql_query(f"""
                                        SELECT *
                                        FROM all_data
                                        WHERE prot_id IN {l_of_ids};
                                    """,con)

    return non_big_prots

def run_iupred3(row):
    seqfile = os.path.join("data/outputs/disorder_iupred/outputs",f"{row.prot_id}.faa")
    with open(seqfile,'w') as f:
        f.write(f">{row.prot_id};\n{row.prot_seq}\n")
    
    process = subprocess.run(
        f"iupred3.py {seqfile} long",
        shell=True,
        capture_output=True,
        text=True
    )

    textout = process.stdout.rstrip().split('\n')
    list_dis_per_aa = []
    for line in textout:
        if line.startswith('#'):
            textout.pop(textout.index(line))
            continue
        if line == '' or line == ' ':
            textout.pop(textout.index(line))
            continue
        aa_val = float(line.split('\t')[-1])
        list_dis_per_aa.append(aa_val)
    
    os.remove(seqfile)
    return list_dis_per_aa

def calc_disord_percent(list_dis_per_aa):
    total_aa = len(list_dis_per_aa)
    disorder_aas = len([x for x in list_dis_per_aa if x > 0.5]) #only consider aa disorder if above 0.5
    disorder_percentage = round((disorder_aas /total_aa) * 100, ndigits=2)
    return disorder_percentage
        
def get_dis_nonbigp_offline(df, output_dir):
    print("Doing Non_BigP")

    non_bigdf = pd.read_csv("data/all_data.csv")
    # big_orgs = df.org_id.unique()
    # non_bigdf = non_bigdf[non_bigdf.org_id.isin(big_orgs)] # just for big prot orgs or no?

    ## repeat N times 
    for i in tqdm(range(1,SAMPLE_TIMES+1), desc="sample"):

        #sample from from the all data dataset
        non_bigdf_sample = non_bigdf.sample(n=len(big_df),replace=False)

        non_bigps = non_bigdf_sample.prot_id.iloc[:].to_list()

        seq_df = query_sql(non_bigps)
        
        list_of_dispercent = []
        for _, row in (t := tqdm(seq_df.iterrows(), total=len(seq_df))):
            prot = row.prot_id
            t.set_description(prot)
            list_dis_per_aa = run_iupred3(row)
            if len(list_dis_per_aa) == 0:
                continue
            disorder_percentage = calc_disord_percent(list_dis_per_aa)
            list_of_dispercent.append(disorder_percentage)
        
        seq_df['dis_percent'] = pd.Series(list_of_dispercent)
        seq_df.to_csv(f"data/outputs/disorder_iupred/nonbigp_disorder_{i}.tsv", sep='\t')

def retrive_taxonomy(tax_id):
    tax_id = int(tax_id)
    lineage = ncbi.get_lineage(tax_id)
    ranks = {v:k for k,v in ncbi.get_rank(lineage).items() if v in ['superkingdom','phylum']}
    trans = {k:v.lower() for k,v in ncbi.get_taxid_translator(ranks.values()).items()}
    final = {k:trans[v] for k,v in ranks.items()}
    return final['superkingdom']

def plot_dis_nonbigp(output_dir):
    fig = make_subplots(
        rows=2,cols=1,
        vertical_spacing=0.04
    )
    colors = {
        'bacteria': '#EF553B',
        'eukaryota': '#00cc96'
    }
    samples_output_path = os.path.relpath("data/outputs/disorder_iupred")
    for sample in os.listdir(samples_output_path):
        sample_path = os.path.join(samples_output_path, sample)
        sdf = pd.read_csv(sample_path, sep='\t', usecols=['prot_id','prot_len','org_id','dis_percent'])
        sdf.dropna(subset='dis_percent', inplace=True)
        sdf['superkingdom'] = sdf.org_id.apply(retrive_taxonomy)
        sdf = sdf[sdf.superkingdom.isin(['bacteria', 'eukaryota'])]
        # print(sdf)

        y_b = sdf[sdf.superkingdom == 'bacteria']['dis_percent']
        fig.add_trace(go.Histogram(y=y_b, histfunc="count", marker_color=colors['bacteria']), row=1, col=1)

        y_e = sdf[sdf.superkingdom == 'eukaryota']['dis_percent']
        fig.add_trace(go.Histogram(y=y_e, histfunc="count", marker_color=colors['eukaryota']), row=2, col=1)

    
    fig.update_layout(barmode='overlay', showlegend=False) 
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(nticks=3, showgrid=True)
    fig.update_traces(opacity=0.9)
    fig.add_annotation(x=-0.021,y=0.7, text="Prokaryotes", textangle=-90, xref="paper", yref="paper", font_size=13)
    fig.add_annotation(x=-0.021,y=0.2, text="Eukaryotes", textangle=-90, xref="paper", yref="paper", font_size=13)
    # fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1].capitalize()))
    fig.add_annotation(x=-0.053,y=0.43,
                    text="Disorder %", textangle=-90,
                    xref="paper", yref="paper",
                    font_size=16)
        
    n_sample = sample.split('_')[0] #.rstrip('.tsv').replace('disorder_','')
    output = path.join(output_dir, f"disorder_distribution_{n_sample}.html")
    fig.write_html(output)


def plot_dis(df: pd.DataFrame, output_dir):
    data = df.copy()
    mapping = {
        'Archaea':'Prokaryota',
        'Bacteria':'Prokaryota',
        'Eukaryota':'Eukaryota'
    }

    color_mapping = {
        'Archaea':'#636efa',
        'Bacteria':'#EF553B',
        'Eukaryota':'#00cc96',
    }

    data['domain'] = data.superkingdom.map(mapping)
    print(data)
    fig = px.scatter(data, x='prot_len', y='dis_percent',
                    color="superkingdom", symbol="superkingdom",
                    facet_row='domain', marginal_y="histogram")

    fig.update_layout(
        # legend=dict(
        #     traceorder="normal"
        # ),
        xaxis_title="Protein Length",
        legend_title="Superkingdom",
    )

    fig.layout.annotations[0].x = - 0.1
    fig.layout.annotations[0].textangle = -90
    fig.layout.annotations[0].font=dict(size=12)
    fig.layout.annotations[1].x = - 0.1
    fig.layout.annotations[1].textangle = -90
    fig.layout.annotations[1].font=dict(size=12)

    for plot in fig.data:
        if plot.type == 'histogram':
            plot.marker.color = color_mapping[plot.name]
        if plot.type == 'scattergl':
            plot.marker.symbol = 'diamond'
            plot.marker.size = 4 if plot.name == 'Archaea' else 2.5
            plot.marker.line = dict(width=0.2, color="#7f7f7f")
            plot.marker.color = color_mapping[plot.name]


    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    
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

    fig.for_each_yaxis(lambda y: y.update(title = ""))
    fig.add_annotation(x=-0.1,y=0.43,
                    text="Disorder %", textangle=-90,
                    xref="paper", yref="paper",
                    font_size=15)

    output = path.join(output_dir, f"disorder_distribution.html")
    fig.write_html(output)
    print(fig)


def plot_tmhms(data, output_dir):
    compress_data = data.copy()
    compress_data['tmh_cnt'] = compress_data.tmh_cnt.apply(lambda x: str(x) if x in [0,1,2] else "3+") 

    group_labels = list(sorted(compress_data.tmh_cnt.unique(), reverse=True))
    
    hist_compress_data = [compress_data[compress_data.tmh_cnt == n].prot_len.to_list() for n in group_labels]
    colors = ['#2ca02c', '#1f77b4', '#fe7e0e', '#9467bd' ]
    
    fig = ff.create_distplot(hist_compress_data, group_labels, show_hist=True, show_curve=True, colors=colors)
    fig.update_xaxes(
        type="log",
    )

    fig.update_layout(
        legend=dict(
            traceorder="reversed",
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        ),
        xaxis_title="Protein Length (Log)",
        yaxis_title="Density",
        legend_title="Number of TMHMs",
    )

    output = path.join(output_dir, f"tmhms_distribution_displot.html")
    fig.write_html(output)

def plot_signal(data, output_dir):

    data['signalp_exist'] = data.signalp_exist.apply(lambda x: "Yes" if x == True else "No")
    
    fig = px.histogram(
                data, x='prot_len',
                color='signalp_exist',
                barmode='group'
                )

    fig.update_yaxes(
        type="log"
    )

    fig.update_layout(
        # legend=dict(
        #     traceorder="normal"
        # ),
        xaxis_title="Protein Length",
        yaxis_title="Count",
        legend_title="Signal Exists",
    )

    output = path.join(output_dir, f"signal_distribution.html")
    fig.write_html(output)

def main(output_dir):
    print("Plotting...")
    plot_dis(mergedf, output_dir)
    plot_tmhms(mergedf, output_dir)
    plot_signal(mergedf, output_dir)
    print("Done")

if __name__ == "__main__":
    main("data/results/signal_tmh_dis")