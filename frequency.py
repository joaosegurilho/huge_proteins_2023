from time import sleep
import warnings

import requests
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.io as pio
import os
from tqdm import tqdm
from ete3 import NCBITaxa
import pickle
ncbi = NCBITaxa()

big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv", sep='\t',
                    usecols=['prot_id', 'prot_len','superkingdom','entry_name'])

template = "plotly_white"
pio.templates.default = template

def doBigP():
    print("Doing BigP")
    df = big_df.copy().sort_values('superkingdom')

    fig = px.histogram(
        df,
        x="prot_len",
        color="superkingdom",
        marginal="rug",
        hover_data=['prot_id', 'prot_len','superkingdom','entry_name']
        )
    fig.update_layout(
        xaxis_title="Protein Length",
        yaxis_title="Count",
        xaxis_range=[5000, 46000])

    # Highlight proteins of interest
    proteins_of_interest = [
        ("Q8WZ42",1.2), #TITIN_HUMAN
        ("A0A5A9P0L4",1.2), #Euka largest
        ("A0A410P257",1.1), #Bac largest
        ("A0A842W954",1), #Archaea largest
    ]

    for bp, arr_h in proteins_of_interest:
        data = df[df.prot_id == bp]
        x_coord = data['prot_len'].iloc[0]
        print(x_coord)
        test_h = 0.75
        if arr_h == 1.2:
            test_h = 0.65
        elif bp == "A0A410P257_9BACT":
            test_h = 0.5
        fig.add_annotation(x=x_coord,
                    text=data['entry_name'].iloc[0],
                    showarrow=True,
                    yref="y domain",
                    # The arrow head will be 40% along the y axis, starting from the bottom
                    y=arr_h,
                    ayref='y domain',
                    ax=0.5,
                    ay=test_h,
                    arrowhead=3)

    fig.write_html("data/results/frequency/huge_proteins_frequency.html")

#  ------------------------------------------------------------------------ #

def convert_to_df(d: dict):
    bins = [0] + [x[1] for x in d.keys()]
    counts = [0] + list(d.values())
    return pd.DataFrame({'Count':counts,'bins':bins})

def generate_sizes(step=100):
    def request_size(bounds):
        requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&seqLength={bounds[0]}-{bounds[1]}"
        try:
            r = requests.get(requestURL, headers={ "Accept" : "application/json"})
        except:
            print(f"Failed for range:{bounds}")
        else:
            if r.ok:
                return int(r.headers['x-pagination-totalrecords'])

    if not os.path.exists("data/outputs/protein_size_bounds.pkl"):
        d = {}
        prev_i = 0
        for i in tqdm(range(step, big_df.prot_len.sort_values(ascending=False).iloc[0]+step, step)):
            bounds = (prev_i,i)
            d[bounds] = request_size(bounds)
            if i % 100 == 0:
                sleep(0.5)
            prev_i += step
        print(d)
        with open("data/outputs/protein_size_bounds.pkl",'wb') as f: pickle.dump(d,f,pickle.HIGHEST_PROTOCOL)
    else:
        with open("data/outputs/protein_size_bounds.pkl",'rb') as f: d = pickle.load(f)
    return convert_to_df(d)

def quantile(df: pd.DataFrame, prob:float = 0.75):
    bins = df.bins
    counts = df.Count
    total_count = np.sum(counts)

    # Calculate cumulative counts
    cumulative_counts = np.cumsum(counts)

    # Normalize cumulative counts
    ecdf = cumulative_counts / total_count

    # Find the bin boundary where the ECDF equals prob
    bin_idx = np.argmax(ecdf >= prob)

    # Estimate the quantile using linear interpolation
    if bin_idx == 0:
        quantile = bins[0]
    else:
        bin_width = bins[bin_idx] - bins[bin_idx - 1]
        prob_within_bin = (prob - ecdf[bin_idx - 1]) / (ecdf[bin_idx] - ecdf[bin_idx - 1])
        quantile = bins[bin_idx - 1] + prob_within_bin * bin_width

    return quantile

def doAll():
    print("Querying...")
    all_df = generate_sizes()
    print(all_df)
    print(all_df.describe())

    print("Plotting...")
    fig_all = px.histogram(all_df, x="bins",y="Count", nbins=len(all_df.bins))
    fig_all.update_layout(
        xaxis_title="Protein Length",
        yaxis_title="Count")

    average = (all_df.bins * all_df.Count).sum() / all_df.Count.sum()

    print("average",average)
    print("q75",quantile(all_df))
    print("q80",quantile(all_df, 0.8))
    print("q95",quantile(all_df, 0.95))
    print("q99",quantile(all_df, 0.99))
    print("q999",quantile(all_df, 0.999))
    print("n_prots > than q95",sum(all_df[all_df.bins >= quantile(all_df,0.95)]['Count']))
    print("n_prots > than q99",sum(all_df[all_df.bins >= quantile(all_df,0.99)]['Count']))
    print("n_prots > than 5k",sum(all_df[all_df.bins >= 5000]['Count']))

    huge_p = all_df[all_df.bins >= 5000]
    print(huge_p.Count.sum())

    fig_all.add_vline(
        x=average,
        annotation_text=f"Average Length<br>{round(average)} aa",
        annotation_position="top",
        annotation_font_size=10, 
        line_width=2, 
        line_dash="dot", 
        line_color="DarkGrey"
        )

    fig_all.add_vline(
        x=5000, 
        annotation_text="Huge<br>Proteins", 
        annotation_position="top", 
        annotation_font_size=10, 
        # annotation_yshift=-5, 
        # annotation_textangle=20, 
        line_width=2, 
        line_dash="dot", 
        line_color="DarkGrey"
        )

    print("Saving...")
    # fig_all.write_json("figure.json")
    fig_all.write_html("data/results/frequency/all_proteins_frequency.html")
    print("Done.")


if __name__ == "__main__":
    # doAll()
    doBigP()



