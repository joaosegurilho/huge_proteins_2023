import pandas as pd
from sys import argv, exit
from os import path

def extract_taxonomy(df: pd.DataFrame):
    df = df.copy()
    keep_clades = ['species','genus','family','order','classe','phylum','superkingdom']
    for idx, row in df.iterrows():
        levels = {x.split(' (')[-1].replace(')',''):" ".join(x.split(' (')[:-1]) for x in row.lineage.split(', ') if x.split(' (')[-1].replace(')','') in keep_clades}
        for k, v in levels.items():
            df.loc[idx, k] = v
    return df.drop(columns=['lineage'])

def extract_assembly(df: pd.DataFrame):
    df = df.copy().dropna(subset=['proteome'])
    for idx, row in df.iterrows():
        l = row.proteome.split(': ')
        df.loc[idx, 'upid'] = l[0]
        df.loc[idx, 'assembly_level'] = l[1]
    return df.drop(columns=['proteome'])

if len(argv) > 1:
    dataset_path = path.relpath(argv[1])

# test: "data/datasets/uniprotkb-big5k-2023_03_23-v13_04.tsv"
uniprotkb = (
    pd.read_csv(dataset_path, sep='\t')
    .rename(columns={"Entry":"prot_id","Entry Name":"entry_name","Protein names":"prot_name","Length":"prot_len","Proteomes":"proteome","Taxonomic lineage":"lineage","Organism (ID)":"org_id","Protein existence":"prot_exist","Gene Ontology IDs":"gos_id","Sequence":"prot_seq"})
    .pipe(extract_taxonomy)
    .pipe(extract_assembly)
)

uniprotkb = uniprotkb.drop(labels=uniprotkb.prot_name.str.extract("(isoform)").dropna().index)
uniprotkb.to_csv("data/datasets/huge_proteins_dataset.tsv",sep='\t')