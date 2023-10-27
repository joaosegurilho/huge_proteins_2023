import pandas as pd

big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv", sep='\t').drop(columns='Unnamed: 0')#,usecols=['prot_id','prot_len','org_id','superkingdom','upid','assembly_level'])
clusters_df = pd.read_csv("data/results/clustering/attempt_2/clusters.processed.tsv", sep='\t')
aling_df = (
    pd.read_csv("data/results/clustering/attempt_2/clusters_aligned.tsv", sep='\t',
                        names=['source','target','seq_idt'], usecols=[0,1,3], dtype={'source':'str','target':'str','seq_idt':'float'})
    .assign(source = lambda x: x['source'].str.replace(';',''))
    .assign(target = lambda x: x['target'].str.replace(';',''))
)

doms_df = (
    pd.read_pickle("data/datasets/huge_proteins_pfam.pkl")
    .rename(columns={'hit_keys':'pfam_names', 'hit_description':'pfam_description','coord_on_prot':'coordinates_on_seq'})
)

def explode_df(df):
    _df = df.copy()
    l_of_dfs = []
    for idx, row in _df.iterrows():
        prot_id = [row.prot_id for _ in range(len(row.coordinates_on_seq))]
        pfam_names = []
        pfam_description = []
        pfam_acession = []
        evalue = []
        for i, occ in enumerate(row.ocurrence):
            pfam_names = pfam_names + [row.pfam_names[i] for _ in range(occ)]
            pfam_description = pfam_description + [row.pfam_description[i] for _ in range(occ)]
            pfam_acession = pfam_acession + [row.pfam_acession[i] for _ in range(occ)]
            evalue = evalue + [row.evalue[i] for _ in range(occ)]
        cols = list(_df.columns)
        cols.pop(4)
        d = dict(zip(cols,[prot_id,pfam_names,pfam_description,pfam_acession,evalue,row.coordinates_on_seq]))
        l_of_dfs.append(pd.DataFrame(d))
        if idx > 2:
            break
    # print(l_of_dfs)
    return pd.concat(l_of_dfs)

print(doms_df)
print(doms_df.iloc[0])
print(doms_df.columns)
d_df = doms_df.pipe(explode_df)#, ['pfam_names', 'pfam_description', 'pfam_acession', 'ocurrence', 'evalue', 'coordinates_on_seq'])
d_df['start'],d_df['end'] = d_df['coordinates_on_seq'].str[0], d_df['coordinates_on_seq'].str[1]
d_df = d_df.drop(columns=['coordinates_on_seq'])
d_df.set_index('prot_id').to_csv("data/datasets/final_datasets/huge_proteins_pfam.tsv",sep='\t')


convert_names = {
    'source':'cluster_repr',
    'target':'prot_id',
    'seq_idt':'repr_seq_identity',
    'org_id':'ncbi_org_id',
    'gos_id':"GO IDs"
    }

merged_df = (
    pd.merge(left=clusters_df, right=aling_df, on='target')
    .drop(columns=['source_y']).rename(columns={'source_x':'source'})
    .merge(big_df, left_on='target', right_on='prot_id')
    .drop(columns=['prot_id'])
    .rename(columns=convert_names)
    )

print(merged_df)
merged_df.set_index('prot_id').to_csv("data/datasets/final_datasets/huge_proteins_all_data.tsv",sep='\t')

merged_df = (merged_df
   .merge(doms_df,on='prot_id')
)
merged_df.set_index('prot_id').to_csv("data/datasets/final_datasets/huge_proteins_all_data_w_pfam.tsv",sep='\t')