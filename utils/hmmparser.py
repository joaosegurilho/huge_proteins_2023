import pandas as pd
from tqdm import tqdm
from Bio import SearchIO

def costumHmmerParser(hmmer_out, output):
    t = {}
    for query in tqdm(SearchIO.parse(hmmer_out,"hmmscan3-domtab")):
        t[query.id.replace(';','')] = {
            'hit_keys':query.hit_keys,
            'hit_description':[query[k].description for k in query.hit_keys],
            'pfam_acession':[x.accession for x in query.hits],
            'ocurrence':[len(y) for y in query.hits],
            'evalue':[z.evalue for z in query.hits],
            'coord_on_prot':[(z.query_start,z.query_end) for z in query.fragments]
            }
    t_df = pd.DataFrame(t).transpose().reset_index().rename(columns={'index':'prot_id'}) 
    t_df.to_pickle(output)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='Parse Hmmscan domtblout')
    parser.add_argument('inputpath')
    parser.add_argument('outputpath')
    args = parser.parse_args()
    costumHmmerParser(args.inputpath, args.outputpath)