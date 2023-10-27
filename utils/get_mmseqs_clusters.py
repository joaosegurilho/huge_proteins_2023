import pandas as pd
from os import path
import sys

def main(mmseq_output_df_path):
    mmseq_output_df = pd.read_csv(mmseq_output_df_path, sep='\t', names=['source','target'])
    mmseq_output_df['source'] = mmseq_output_df.source.apply(lambda x: x.rstrip(';'))
    mmseq_output_df['target'] = mmseq_output_df.target.apply(lambda x: x.rstrip(';'))

    mmseq_output_df['cluster'] = mmseq_output_df.groupby('source').ngroup()

    print(mmseq_output_df.sort_values(by="cluster",ascending=True).head(15))
    print(f"Total unique proteins: {len(mmseq_output_df.target.unique())}")
    print(f"Number of clusters: {len(mmseq_output_df.sort_values(by='cluster',ascending=True).cluster.unique())}")
    
    mmseq_output_df_path = ".processed.".join(mmseq_output_df_path.split('.'))
    mmseq_output_df.sort_values(by='cluster',ascending=True).to_csv(mmseq_output_df_path, sep='\t', index=False)
    print(f"Saved in {mmseq_output_df_path}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        mmseq_output_path = path.relpath(sys.argv[1])
        main(mmseq_output_path)
    else:
        print("Need input file.")
        sys.exit()
