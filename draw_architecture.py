
## Imports
import os
import random
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from tqdm import tqdm
from dna_features_viewer import GraphicFeature, GraphicRecord


big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv",sep='\t')
domsP = pd.read_pickle("data/datasets/huge_proteins_pfam.pkl")

## Helper Functions
def pattern_generator(maindf,acession_type="pfam_acession"):
    possible_pfams = set(x for lst in maindf[acession_type].to_list() for x in lst)
    list_of_colors = random.choices(list(mcolors.CSS4_COLORS.keys()),k=len(possible_pfams))
    
    return dict(zip(possible_pfams,list_of_colors))

def deduce_doms(prot_id,acess,occur,coord,colors):
    # input: protien acession, pfam_acessions, occurence for each acession, coordinates of each domain ordered, color for each domain (for consistency)  
    l = []
    cnt = 0
    for o in occur:
        [l.append(acess[cnt]) for _ in range(int(o))]
        cnt +=1

    pat_dict = dict(zip(coord,l)) 
    features = []
    for k,v in pat_dict.items():
        feature = GraphicFeature(start=k[0], end=k[1], strand=0, color=colors[v],label=v) #change strand between +1,-1 or 0 for protein repr
        features.append(feature)
    
    return features

def plot_prot(maindf:pd.DataFrame,save_to,pfam_keys=False):
    codes_df = pd.read_csv("data/outputs/mapping/code2arch.tsv", sep='\t', index_col='arch_id')
    
    if pfam_keys:
        colors = pattern_generator(maindf,acession_type='hit_keys')
    else:
        colors = pattern_generator(maindf)
    for idx, row in tqdm(maindf.iterrows(), total=len(maindf)):
        prot_id = row.prot_id
        prot_len = row.prot_len
        seq = row.prot_seq
        if pfam_keys:
            pfam = row.hit_keys
        else:
            pfam = row.pfam_acession
        coord = row.coord_on_prot
        occur = row.ocurrence
        try:
            arch = codes_df.query(f"prot_id == '{prot_id}'").index[0]
        except:
            print(prot_id)
            print(maindf)
            print(codes_df.query(f"prot_id == '{prot_id}'"))
        
        prot_features = deduce_doms(prot_id,pfam,occur,coord,colors)
        
        record = GraphicRecord(sequence=seq, sequence_length=prot_len, features=prot_features)
        
        ax, _ = record.plot(with_ruler=False)
        ax.set_title(f"protein:{prot_id} | arch:{arch}")
        ax.figure.tight_layout()
        save_local = f"{save_to}/{arch}_{prot_id}_repr.png"
        ax.figure.savefig(save_local)
        plt.close()
    print(f"Saved repr to: {save_to}")

def main_draw_controller(POI,draw_archs=True,save_dir="architectures"):
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    fusedf = big_df.merge(domsP, on='prot_id').drop(columns=['Unnamed: 0'])
    print(fusedf)
    print(fusedf.columns)

    if draw_archs:
        df_to_send = fusedf.query(f"prot_id in {POI}")#.query(f"prot_id in {POI}")
        print(len(df_to_send))
        plot_prot(df_to_send,save_to=save_dir,pfam_keys=True)
    else:
        plot_prot(fusedf.query(f"prot_id in {POI}"),save_to="big_protein_representation")

if __name__ == "__main__":
    # POI_test = ['A0A2Z3H2M1_9BACT',
    #             'M5TU89_9BACT',
    #             'A0A5B9R051_9BACT']

    POI_all = big_df.prot_id.to_list()
    print(POI_all[0])
    print(len(POI_all))
    with open("missing_archs.txt", 'r') as f:
        POI_missing = [x.rstrip() for x in f.readlines()]
        
        print(POI_missing[0])
    
    main_draw_controller(POI_missing,draw_archs=True,save_dir="data/results/architectures")

