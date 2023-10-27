import os
import pandas as pd
import numpy as np
import requests
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv", sep='\t',
                     usecols=['prot_id','prot_len','org_id','superkingdom','upid','assembly_level'])


def get_proteome_id(prot):

    if prot.startswith('HMU'):
        return (prot, np.nan)

    url = "https://www.ebi.ac.uk/proteins/api/proteins/"
    query = prot.split('_')[0]

    try:
        r = requests.get(url+query, headers={ "Accept" : "application/json"})

    except:
        print(f"Request prot:{prot} not working!")
        return (prot, np.nan)

    else:
        if r.ok:
            proteome = np.nan
            try:
                for evid in r.json()['keywords']:
                    if evid['evidences'][0]['source']['name'] == 'Proteomes':
                        proteome = evid['evidences'][0]['source']['id']
                
            except:

                try:
                    for evid in r.json()['references'][0]['source']['strain'][0]['evidences']:
                        if evid['source']['name'] == 'Proteomes':
                            proteome = evid['source']['id']
                
                except:
                    print(f"Request {prot} no proteome.")
                    return (prot, proteome)

                else:
                    return (prot, proteome)
            
            else:
                return (prot, proteome)

def get_busco_score(proteome):
    url = "https://www.ebi.ac.uk/proteins/api/proteomes/"

    try:
        r = requests.get(url+proteome, headers={ "Accept" : "application/json"})

    except:
        return (proteome, np.nan)

    else:
        if r.ok:
            for score in r.json()['scores']:
                if score['name'] == 'busco':
                    busco = score['property'][1]['value']
            return (proteome, busco)


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def query_busco_scores(prot_list):
    from utils import query_uniprot
    for chunk in chunks(prot_list, 500):
        query = '%29%20OR%20%28'.join([f'upid%3A{x}' for x in chunk])
        url = f'https://rest.uniprot.org/proteomes/search?fields=upid%2Cbusco&format=tsv&query=%28%28{query}%29%29&size=500'
        query_uniprot.run_query("data/datasets/uniprotkb-proteomes-busco.tsv", url)

def clean_busco(df: pd.DataFrame):
    for idx, row in df.iterrows():
        if isinstance(row.busco, str) and row.busco.startswith('C:'):
            df.loc[idx, 'busco'] = row.busco.split(':')[1].rstrip('%[S')
        else:
            df.loc[idx, 'busco'] = row.busco
    return df

def read_busco_scores():
    busco_df = (
        pd.read_csv("data/datasets/uniprotkb-proteomes-busco.tsv", sep='\t')
        .rename(columns={'Proteome Id':'upid','BUSCO':'busco'})
        .set_index('upid')
        .dropna()
        .pipe(clean_busco)
        .astype('float16')
    )
    return busco_df

def clean_completeness(text: str):
    if text.startswith("Plasmid"):
        text = "Plasmid"
    elif text.startswith("Chromosome"):
        text = "Chromosome"
    elif text.__contains__("Linkage group"):
            text = "Linkage Group"
    elif text.startswith("LG"):
        text = "Linkage Group"
    elif text.startswith("Contig"):
        text = "Contig"
    
    return text


def plot_bigp_plen(df):
    fig = make_subplots(
        rows=2,cols=3,
        subplot_titles=tuple(df.superkingdom.unique()),
        vertical_spacing=0.04,
        horizontal_spacing=0.03,
    )

    df_use = df.copy()
    colors = {
        'n_bigp':'#6DA96C',
        'prot_len':'#E27D52'
    }
    for i, sk in enumerate(df.superkingdom.unique(), start=1):
        for j, analysis in enumerate(['n_bigp','prot_len'], start=1):

            x = df_use[df_use.superkingdom == sk].busco.to_list() 
            
            y= df_use[df_use.superkingdom == sk][analysis].to_list() 
            fig.add_trace(go.Scatter(
                    x=x,
                    y=y,
                    mode='markers',
                    marker=dict(
                        color=[colors[analysis] for _ in range(len(y))],
                    )
                ), row=j, col=i)

    fig.update_layout(
        title="N.Big proteins per proteome and Protein length vs Busco Completeness",
        showlegend=False
    )

    fig.update_yaxes(matches=None)
    fig.update_yaxes(title_text="N. Big proteins <br> per Proteome", row=1, col=1)
    fig.update_yaxes(title_text="Protein Length", row=2, col=1)
    fig.update_xaxes(matches=None)

    fig.update_traces(marker_size=2.5)
    fig.for_each_annotation(lambda a: a.update(text=a.text.capitalize(), font_size=18))
    
    fig.add_annotation(x=0.5,y=-0.2,
                    text="Busco Completeness",
                    xref="paper", yref="paper",
                    font_size=15,
                    showarrow=False)

    fig.write_html("data/results/genome_quality/nbigp_and_plen_vs_busco_proteome_quality.html")

def plot(df):
    fig = px.scatter(df,
                    x='n_bigp',
                    y='busco',
                    facet_col='superkingdom',
                    hover_data=['completeness', 'busco', 'n_bigp','proteomes'])

    fig.update_layout(
        title="N.Big proteins per proteome vs Busco Completeness",
        yaxis_title="Busco Completeness Score"
    )

    fig.update_xaxes(matches=None)

    fig.update_traces(marker_size=5)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split('=')[1].capitalize(), font_size=18))
    
    fig.for_each_xaxis(lambda x: x.update(title = ""))
    fig.add_annotation(x=0.5,y=-0.08    ,
                    text="N. Big Proteins per Proteome",
                    xref="paper", yref="paper",
                    font_size=15,
                    showarrow=False)

    fig.write_html("data/results/genome_quality/nbigp_vs_busco_proteome_quality.html")

def plot_bar(df):
    fig = px.histogram(df,
                x='n_bigp',
                y='busco',
                color='busco',
                facet_col='superkingdom',
                histfunc='avg',
                hover_data=['completeness', 'busco', 'n_bigp','proteomes'])

    fig.update_layout(
        title="N.Big proteins per proteome vs Busco Completeness",
        yaxis_title="Busco Completeness Score",
    )

    fig.update_xaxes(matches=None)
    fig.update_yaxes(matches=None)


    fig.write_html("data/results/genome_quality/bar_nbigp_vs_busco_proteome_quality.html")

def plot_by_prot_len(df):
    fig = px.scatter(df,
                    x='prot_lenght',
                    y='busco',
                    color='completeness',
                    facet_col='superkingdom',
                    hover_data=['completeness', 'busco', 'prot_id','proteomes'])
    
    fig.update_layout(
        title="Protein Length vs Busco",
        yaxis_title="Busco Completeness Score"
        )
    fig.update_traces(marker_size=4)

    fig.write_html("data/results/genome_quality/bp_proteins_vs_busco_proteome_quality.html")

def plot_by_genome_quality(df):
    fig = px.box(df,
                x='completeness',
                y='prot_lenght',
                color='completeness',
                facet_col='superkingdom',
                hover_data=['completeness', 'busco', 'prot_id','proteomes'])
    
    fig.update_xaxes(matches=None, title='Genome Component')
    fig.update_layout(
        title="Protein Length vs Assigned Genome Quality",
        template="plotly_white",
        yaxis_title="Protein Length",
        showlegend=False,
        )
    
    fig.update_traces(marker_size=5,marker_line=dict(width=0.5))#,color="#7f7f7f"))
    
    new_categories_euka = list(fig.layout.xaxis3.categoryarray)
    new_categories_euka.pop(new_categories_euka.index('Plasmid'))
    fig.layout.xaxis3.categoryarray = new_categories_euka
    

    fig.write_html("data/results/genome_quality/bp_prot_length_vs_genome_quality.html")

def plot_by_genome_quality_proteomes(df):
    fig = make_subplots(
        rows=1,cols=3,
        subplot_titles=tuple(df.superkingdom.unique()),
        horizontal_spacing=0.04,
        specs=[[{"secondary_y": True}, {"secondary_y": True}, {"secondary_y": True}]]
    )
        
    df_use = df.copy()
    for i, sk in enumerate(df.superkingdom.unique()):

        x = df_use[df_use.superkingdom == sk].assembly_level.to_list() 
        
        y = df_use[df_use.superkingdom == sk].n_bigp.to_list() 
        fig.add_trace(go.Box(x=x,y=y),row=1,col=i+1, secondary_y=False)

        x_values = list(sorted(set(x)))
        y_values = df_use[df_use.superkingdom == sk].groupby(by='assembly_level').count().sort_values(by='assembly_level').prot_id.to_list()
        fig.add_trace(go.Bar(x=x_values,y=y_values,marker=dict(color="#a6c5da",opacity=0.5)), row=1,col=i+1, secondary_y=True)
        
    fig.update_layout(
        title="N.Big proteins per proteome vs Assigned Genome Quality",
        template="plotly_white",
        yaxis_title="N. Big Proteins per Proteome",
        yaxis6=dict(title='N. Proteomes',
                    side='right'),
        showlegend=False
        )
    
    fig.for_each_annotation(lambda a: a.update(text=a.text.capitalize(), font_size=20))

    fig.update_xaxes(matches=None, tickfont=dict(size=16))
    fig.update_yaxes(matches=None, title=dict(font=dict(size=25)),
                     tickfont=dict(size=20), ticklabelposition="inside top")

    fig.add_annotation(x=0.5,y=-0.25,
                    text="Assembly level",
                    xref="paper", yref="paper",
                    font_size=20)
    
    fig.layout.yaxis1.rangemode = 'tozero'
    fig.layout.yaxis1.showline = True
    fig.layout.yaxis1.linewidth = 1

    fig.layout.yaxis2.rangemode = 'tozero'
    fig.layout.yaxis2.showline = True
    fig.layout.yaxis2.linewidth = 1
    fig.layout.yaxis2.linecolor ='LightGray'
    fig.layout.yaxis2.gridcolor ='LightGray'

    fig.layout.yaxis3.rangemode = 'tozero'
    fig.layout.yaxis3.showline = True
    fig.layout.yaxis3.linewidth = 1
    
    fig.layout.yaxis4.rangemode = 'tozero'
    fig.layout.yaxis4.showline = True
    fig.layout.yaxis4.linewidth = 1
    fig.layout.yaxis4.linecolor ='LightGray'
    fig.layout.yaxis4.gridcolor ='LightGray'
    
    fig.layout.yaxis5.rangemode = 'tozero'
    fig.layout.yaxis5.showline = True
    fig.layout.yaxis5.linewidth = 1
    
    fig.layout.yaxis6.rangemode = 'tozero'
    fig.layout.yaxis6.showline = True
    fig.layout.yaxis6.linewidth = 1
    fig.layout.yaxis6.linecolor ='LightGray'
    fig.layout.yaxis6.gridcolor ='LightGray'

  
    fig.write_html("data/results/genome_quality/nbigp_vs_assembly_quality.html")

def main():

    if not os.path.exists("data/datasets/uniprotkb-proteomes-busco.tsv"):
        print("Querying Uniprot...")
        query_busco_scores(big_df.upid.iloc[:])
        busco_df = read_busco_scores()
        df = big_df.merge(busco_df, left_on='upid', right_index=True)
        
        df = big_df.copy()
    else:
        busco_df = read_busco_scores()
        df = big_df.merge(busco_df, left_on='upid', right_index=True)
        print(df)

    df['assembly_level'] = (df['assembly_level']
        .apply(lambda x: x.split(';')[0])
        .apply(clean_completeness)
    )
    
    exclude = [
        'Unplaced LGUn',
        'Unplaced contigs',
        'Genome',
        'Genome assembly',
        'Contig'
    ]
    conditions = (df['assembly_level'].isin(exclude))
    df = df[~conditions]

    print("Plotting...")
    temp = df[['prot_id','upid']].groupby('upid').count().rename(columns={'prot_id':'n_bigp'}).to_dict()['n_bigp']
    df['n_bigp'] = df.upid.map(temp)

    df = df.sort_values('superkingdom')
    plot_bigp_plen(df)
    plot_by_genome_quality_proteomes(df)
    print("Done")

if __name__ == "__main__":
    main()

