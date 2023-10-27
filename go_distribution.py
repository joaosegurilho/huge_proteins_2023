import os
import pandas as pd
import requests
from tqdm import tqdm
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx
import obonet

bigP = pd.read_csv("data/datasets/huge_proteins_dataset.tsv",sep='\t', usecols=['prot_id','gos_id','superkingdom']).dropna(subset='gos_id')
bigP['gos_id'] = bigP['gos_id'].apply(lambda x: x.split('; '))


def correct_gos(go_term):
    try:
        id_to_name[go_term]
    except:
        #print(f"Failed for {go_term}")

        requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
        tool = "/secondaryids"

        r = requests.get(requestURL + go_term + tool, headers={ "Accept" : "application/json"})

        if not r.ok:
            print(r.headers)
            
        response = r.json()['results'][0]
        return response['id']
    else:
        return go_term

def check_obsolete(maindf):
    m_df = maindf.copy()
    
    for go in tqdm(m_df.columns, desc='obsolete:'):
        requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
        tool = "/secondaryids"

        r = requests.get(requestURL + go + tool, headers={ "Accept" : "application/json"})

        if not r.ok:
            print(r.headers)

        response = r.json()['results'][0]

        if response['isObsolete'] == True:
            m_df = m_df.drop(go, axis=1)
    
    return m_df

def get_sublevels(df, gotype, graph):
    #get all the level 2 go terms
    node = name_to_id[mapping_names[gotype]]
    level = [parent for parent, child, key in graph.in_edges(node, keys=True)]
    specify = mapping_specificity[gotype]
    level_updated = []
    for fun in level:
        if fun in specify:
            sub_level =[parent for parent, child, key in graph.in_edges(fun, keys=True)]
            level_updated = level_updated + sub_level
            continue

        level_updated.append(fun)

    #for every protein, check associated go terms    
    function_groups = {}
    for row in tqdm(df.iterrows(), total=len(df), desc='associated'):
        prot_id = row[0] 
        prot_gos = [go for go, v in row[1].items() if v == 1.0]
        
        function_groups[prot_id]=[]
        for go in prot_gos:
            if go in specify:
                #print(id_to_name[go])
                function_groups[prot_id].append(id_to_name[go])
            else:
                try:
                    parents = [superterm for superterm in networkx.descendants(graph, go)]
                except Exception as e:
                    continue
                else:
                    for parent in parents:
                        if parent in level_updated and parent not in function_groups[prot_id]:
                            function_groups[prot_id].append(id_to_name[parent])
    return function_groups

def find_uniques(df) -> pd.DataFrame:
    df = df.groupby(by=['function']).count().sort_values(by='prot_id', ascending=False)
    df_uniques = df.reset_index()
    unique_gos_dict = {}
    for i in tqdm(df_uniques.iterrows(), desc='uniques'):
        s = i[1][0]
        function_list = s.strip('[').strip(']').replace("'","").split(", ")
        counts = int(i[1][1])
        
        ## Get all the unique GOs
        for z in function_list:
            if z in unique_gos_dict:
                unique_gos_dict[z][0] += counts
            else:
                unique_gos_dict[z] = [counts]
    return pd.DataFrame.from_dict(unique_gos_dict).T.sort_values(0,ascending=False)

def plot_together(df,gotype):
    to_use_df = df.copy()
    to_use_df = find_uniques(to_use_df)
    data = to_use_df[to_use_df.index != ''].sort_values(0,ascending=False).head(15).sort_values(0,ascending=True).rename(columns={0:'counts'})

    fig = px.bar(
        data,
        x='counts',
        # facet_col='superkingdom',
        orientation='h',
        labels={
            'counts':"Counts",
            "index":f"Go terms - {gotype[0].upper() + gotype[2].upper()}"
        },
        height=500,
        width=1000
        )
    fig.update_layout(
        template='plotly_white',
        )
    fig.update_traces(marker_color='#f9c8a0',marker_line_color='#A4968A',
                    marker_line_width=1.5)
    fig.write_html(f"data/results/go_distribution/{gotype}_groups_unique.html")

def plot(df, gotype):
    print(df)
    fig = make_subplots(
        rows=3, subplot_titles=tuple(df.superkingdom.unique()))
    
    for i, sk in enumerate(df.superkingdom.unique()):
        to_use_df = df[df.superkingdom == sk]
        to_use_df = find_uniques(to_use_df)
        data = to_use_df[to_use_df.index != ''].sort_values(0,ascending=False).head(15).sort_values(0,ascending=True).rename(columns={0:'counts'})
        print(data)
        fig.add_trace(go.Bar(
            x=data.counts,
            y=data.index,
            orientation='h',
        ), row=i+1, col=1)


    fig.update_layout(template='plotly_white', title=f"Go terms - {gotype[0].upper() + gotype[2].upper()}",showlegend=False)
    fig.update_traces(marker_color='#f9c8a0',marker_line_color='#A4968A',
                    marker_line_width=1.5)
    fig.write_html(f"data/results/go_distribution/{gotype}_groups_unique_by_sk.html")
            
def generate_go_terms_dataset(graph):
    heads = {
        'GO:0003674':'m_function',
        'GO:0008150':'b_process',
        'GO:0005575':'c_component'
    } 

    df = pd.DataFrame(columns=heads.values(), index=bigP.prot_id)
    
    for idx, row in tqdm(bigP.set_index('prot_id').iterrows(), desc='building dataset', total=len(bigP)):
        for go in row.gos_id:
            try:
                descendts = networkx.descendants(graph, go)
            except:
                continue
            gotype = [heads[x] for x in descendts if x in heads.keys()][0]
            if isinstance(df.loc[idx, gotype], list):
                df.loc[idx, gotype].append(go)
            else:
                df.loc[idx, gotype] = go
    
    df.to_csv(f"data/outputs/go_terms/go_terms_dataset.tsv", sep='\t', index=False)

def main():
    url = 'http://purl.obolibrary.org/obo/go.obo'
    graph = obonet.read_obo(url)

    global id_to_name
    global name_to_id
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}

    global mapping_specificity
    mapping_specificity = {
        'b_process':[
            name_to_id['regulation of biological process'],
            name_to_id['multicellular organismal process'],
            name_to_id['developmental process'],
            name_to_id['biological regulation'],
            name_to_id['metabolic process'],
            name_to_id['cellular process'],
            name_to_id['cellular component organization or biogenesis']
            ],
        'm_function':[name_to_id['binding'], name_to_id['catalytic activity']],
        'c_component':[name_to_id['cellular anatomical entity'],name_to_id['protein-containing complex'],name_to_id['intracellular anatomical structure']]
    }

    global mapping_names
    mapping_names = {
        'b_process':'biological_process',
        'm_function':'molecular_function',
        'c_component':'cellular_component'
    }

    for gotype in ['b_process','m_function','c_component']:
        if os.path.exists(f"data/outputs/go_terms/function_groups_{gotype}.pkl"):
            prot_func_df =  pd.read_pickle(f"data/outputs/go_terms/function_groups_{gotype}.pkl")
                
        
        else:
            print(gotype)
            if not os.path.exists(f"data/outputs/go_terms/go_terms_dataset.tsv"):
                generate_go_terms_dataset(graph)
                functiondf = (
                    pd.read_csv(f"data/outputs/go_terms/go_terms_dataset.tsv", sep='\t', usecols=['prot_id',gotype])
                )
                functiondf = pd.crosstab(functiondf['prot_id'],functiondf[gotype])
            else:
                functiondf = (
                    pd.read_csv(f"data/outputs/go_terms/go_terms_dataset.tsv", sep='\t', usecols=['prot_id',gotype])
                )
                functiondf = pd.crosstab(functiondf['prot_id'],functiondf[gotype])
            
            functiondf.columns = [correct_gos(x) for x in tqdm(functiondf.columns, desc='correcting')]
            functiondf = check_obsolete(functiondf)

            function_groups = get_sublevels(functiondf,gotype, graph)

            prot_func_df = bigP.merge(pd.Series(function_groups, name='function'), left_on='prot_id', right_index=True).astype('str')
            prot_func_df.to_pickle(f"data/outputs/go_terms/function_groups_{gotype}.pkl")
                
        plot(prot_func_df, gotype)
        plot_together(prot_func_df, gotype)


if __name__ == "__main__":
    main()