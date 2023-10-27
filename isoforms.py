import os
import re
import pandas as pd
import numpy as np
import requests
from time import sleep
from pprint import pprint
from tqdm import tqdm
from bioservices import UniProt
import pickle

big_df = pd.read_csv("data/datasets/big_proteins_dataset.tsv",sep='\t')
org_meta = pd.read_csv("data/datasets/organism_metadata_dataset.tsv",sep='\t')
merge_df = big_df.merge(org_meta, on='org_id').drop(columns='prot_seq')

u = UniProt(verbose=False, cache=True)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def query_prot_names(prots):
    data = u.get_df(
        prots,
        limit=25,
        columns="Entry Name,Protein names",
        )
    return data.set_index('Entry Name')

def query_isoforms(prots):
    data = u.get_df(
        prots,
        limit=25,
        columns="Entry Name,Alternative products (isoforms)",
        )
    return data.dropna(subset='Alternative products (isoforms)').set_index('Entry Name')

def get_n_of_isoforms(hit):
    if hit != np.nan and not hit.endswith(': '): 
        l_isos = [int(x.split('=')[1]) for x in hit.strip().split(';') if x.strip().startswith('Named')]
        if len(l_isos) > 0:
            return l_isos[0]
        else:
            return 0

def get_isoforms(hit):
    return [x.split('=')[1] for x in hit.strip().split(';') if x.strip().startswith('IsoId')]

def get_iso_len_2(df):
    data = {}
    for isos in tqdm(df.isoforms):
        for isoform_id in isos:
            response = requests.get(f'https://www.uniprot.org/uniprot/{isoform_id}.fasta')
            sequence = ''.join(response.text.split('\n')[1:])
            data[isoform_id] = [sequence]
        sleep(0.5)            
    return pd.DataFrame(data).T.rename(columns={0:'iso_seq'})

def find_in_dataset(df):
    df = df.copy()
    for idx, row in df.iterrows():
        h = big_df[big_df.prot_seq == row.iso_seq]
        df.loc[idx, 'in_dataset'] = not h.empty
        df.loc[idx, 'name_in_dataset'] = '|'.join(h.prot_id.values.tolist()) if not h.empty else np.nan
    return df

def main():
    if not os.path.exists("data/outputs/found_isoforms_of_bigproteins.pkl"):
        found = {}
        for proteins in tqdm(chunks(merge_df.prot_id.to_list(), 25), total=round(len(merge_df)/25)):
            d = query_isoforms(proteins)
            for protein, hit in d.iterrows():
                c = get_n_of_isoforms(hit.iloc[0])
                if c > 0:
                    found[protein] = get_isoforms(hit.iloc[0])
        
        pprint(found)
        with open("data/outputs/found_isoforms_of_bigproteins.pkl",'wb') as f:
            pickle.dump(found, f, pickle.HIGHEST_PROTOCOL)
        pd.DataFrame(found).to_csv("data/datasets/found_isoforms_of_bigproteins.tsv", sep='\t')

    else:
        with open("data/outputs/found_isoforms_of_bigproteins.pkl",'rb') as f:
            found = pickle.load(f)
        founddf = pd.DataFrame(pd.Series(found, name='isoforms'))
        print(founddf)

        orgs_of_inter = merge_df[merge_df.prot_id.isin(founddf.index)]['org_id'].to_list()
        orgs_w_isos_df = merge_df[merge_df.org_id.isin(orgs_of_inter)].set_index('prot_id')
        for prt in orgs_w_isos_df.index:
            # print(prt)
            fnd = found.get(prt)
            has_iso = 1 if isinstance(fnd, list) else 0
            orgs_w_isos_df.loc[prt, 'has_iso'] = has_iso
            orgs_w_isos_df.loc[prt, 'n_iso'] = len(fnd) if has_iso else 0

        orgs_w_isos_df.to_csv("data/outputs/proteins_w_iso.csv", header=True)
        
        if not os.path.exists("data/outputs/found_isoforms_len.pkl"):
            isoslen_df = get_iso_len_2(founddf)
            isoslen_df['iso_len'] = isoslen_df['iso_seq'].apply(lambda x: len(x))
            with open("data/outputs/found_isoforms_len.pkl",'wb') as f:
                pickle.dump(isoslen_df, f, pickle.HIGHEST_PROTOCOL)
        else:
            isoslen_df = pd.read_pickle("data/outputs/found_isoforms_len.pkl")
        
        for idx, isos in founddf.iterrows():
            for iso in isos.values.item():
                isoslen_df.loc[iso, 'original'] = idx
        

        isoslen_df = isoslen_df.reset_index().set_index(['original', 'index'])
        tempdf = isoslen_df[isoslen_df.iso_len > 5000]
        tempdf = find_in_dataset(tempdf).sort_index(level='original').drop(columns='iso_seq')
        print(tempdf)
        tempdf.to_csv("data/outputs/found_isos_w_5k_len.csv", header=True, index=True)
        
if __name__ == "__main__":
    if not os.path.exists("data/outputs/protein_names.pkl"):
        f_names_dfs = []
        for proteins in tqdm(chunks(merge_df.prot_id.iloc[:].to_list(), 25), total=round(len(merge_df)/25)):
            temp_df = query_prot_names(proteins)
            f_names_dfs.append(temp_df)
        protein_names_df = pd.concat(f_names_dfs)
        print(protein_names_df)
        protein_names_df.to_pickle("data/outputs/protein_names.pkl")
    else:
        protein_names_df = pd.read_pickle("data/outputs/protein_names.pkl")
        i_df = protein_names_df['Protein names'].str.extract(r"(isoform)", re.IGNORECASE)
        print(i_df.dropna())