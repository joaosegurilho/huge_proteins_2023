
import pandas as pd
from tqdm import tqdm
import os

from goatools import obo_parser, go_enrichment
from goatools.utils import read_geneset
from goatools.anno.idtogos_reader import IdToGosReader

import plotly.express as px

big_df = pd.read_csv("data/datasets/huge_proteins_dataset.tsv", sep='\t',
                    usecols=['prot_id','phylum','superkingdom','gos_id'])


def write_all_genesets(df: pd.DataFrame, path_of_net):
    enr_go_dir = os.path.join(path_of_net,"enrichment")
    if not os.path.exists(enr_go_dir):
            os.mkdir(enr_go_dir)

    #backgorund pop
    study_background_path = os.path.join(enr_go_dir, f"study_background.txt")
    with open(study_background_path,"w") as study_background_f:
        for prot in df.prot_id.to_list():
            study_background_f.write(f"{prot}\n")

    #pop annoted 
    annot_pop_path = os.path.join(enr_go_dir, f"annot_pop.txt")
    with open(annot_pop_path,"w") as annot_pop_f:
        for i, prot, gos in df[['prot_id','gos_id']].fillna(value={'gos_id':''}).itertuples():
            gos = gos.replace(' ','')
            annot_pop_f.write(f"{prot}\t{gos}\n")
    
    #study pop
    for cluster in df.cluster.unique():
        subset = df[df.cluster == cluster]

        temp_path = os.path.join(enr_go_dir, f"study_of_{cluster}")
        if not os.path.exists(temp_path):
            os.mkdir(temp_path)

        study_pop_path = os.path.join(temp_path, f"study_pop_{cluster}.txt")
        with open(study_pop_path,"w") as study_pop_f:
            for prot in subset.prot_id.to_list():
                study_pop_f.write(f"{prot}\n")
        
    return study_background_path, annot_pop_path, enr_go_dir

def perform_enrichment_in_cluster(cluster_num, g_enrich_main, enrich_dir):

    study_path = os.path.join(enrich_dir,f"study_of_{cluster_num}")

    sample_path = os.path.join(study_path, f"study_pop_{cluster_num}.txt")
    prot_of_interest = read_geneset(sample_path)

    study = g_enrich_main.run_study(prot_of_interest)

    def divide_tuples(tupl):
        x, y = tupl
        return (x/y)*100

    map_enr_type = {'e':'enriched','p':'purified'}

    study_df = pd.DataFrame(columns=["GO", "NS", "name", "enriched_purified", "num_in_study" , "ratio_in_study%", "ratio_in_pop%", "p_uncorrected", "bonferroni", "benjamimi_hochberg"])

    for i, g in enumerate(study):
        study_df.loc[i] = [g.GO] + [g.NS, g.name, map_enr_type[g.enrichment], g.ratio_in_study[0], divide_tuples(g.ratio_in_study), divide_tuples(g.ratio_in_pop), g.p_uncorrected, g.p_bonferroni, g.p_fdr_bh]

    enrich_df_path = os.path.join(enrich_dir, f"study_of_{cluster_num}", f"dataframe_{cluster_num}_enriched_purified.csv")
    study_df.to_csv(enrich_df_path, index=False)

    plot_enrichment(study_df, study_path)

def plot_enrichment(df, study_save_to):
    condition_to_filter = (df.num_in_study != 0) & (df.bonferroni <= 0.05)
    n_df = df[condition_to_filter]

    if len(n_df) < 1:
        condition_to_filter = (df.bonferroni <= 0.05)
        n_df = df[condition_to_filter]
        if len(n_df) < 1:
            n_df = df

    fig = px.bar(
        n_df,
        x='num_in_study',
        y='name',
        color='enriched_purified',
        facet_row='NS',
        orientation='h',
        hover_data=['benjamimi_hochberg','bonferroni','p_uncorrected'],
        template='seaborn'
    )

    fig.update_yaxes(matches=None, tickfont_size=9)
    fig.update_yaxes(tickfont_size=9)
    
    figure_path = os.path.join(study_save_to, f"enrichment_of_go_terms_{study_save_to.split('_')[-1]}.html")
    fig.write_html(figure_path)

def main_enr(cluster_file, obo_path):
        
    clusterdf = pd.read_csv(cluster_file, sep='\t')
    merged_df = clusterdf.merge(big_df, left_on='target', right_on='prot_id').drop(columns=['target'])

    background_path, annot_path, enrich_dir = write_all_genesets(merged_df, os.path.dirname(cluster_file))

    obo_dag = obo_parser.GODag(obo_path)
    
    background_pop = read_geneset(background_path)

    annot_obj = IdToGosReader(annot_path, godag=obo_dag)
    go_associassion = annot_obj.get_id2gos()

    g_bonferroni = go_enrichment.GOEnrichmentStudy(
        background_pop,
        go_associassion,
        obo_dag,
        propagate_counts=False,
        alpha=0.05,
        methods=['bonferroni','fdr_bh'],
        pvalcalc='fisher_scipy_stats'
        )

    print("Performing the enrichment study for each cluster...")
    for cluster in tqdm(clusterdf.cluster.unique()):
        perform_enrichment_in_cluster(cluster, g_bonferroni, enrich_dir)
    
    print("Finished with enrichment.")


if __name__ == "__main__":
    main_enr("data/results/clustering/attempt_2/clusters.processed.tsv", "data/outputs/go_terms/go-basic.obo")