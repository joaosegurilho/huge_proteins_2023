import pandas as pd
import os

eggnong_mapp = {
    "A":"RNA processing and modification",
    "B":"Chromatin structure and dynamics",
    "C":"Energy production and conversion",
    "D":"Cell cycle control, cell division, chromosome partitioning",
    "E":"Amino acid transport and metabolism",
    "F":"Nucleotide transport and metabolism",
    "G":"Carbohydrate transport and metabolism",
    "H":"Coenzyme transport and metabolism",
    "I":"Lipid transport and metabolism",
    "J":"Translation, ribosomal structure and biogenesis",
    "K":"Transcription",
    "L":"Replication, recombination and repair",
    "M":"Cell wall/membrane/envelope biogenesis",
    "N":"Cell motility",
    "O":"Posttranslational modification, protein turnover, chaperones",
    "P":"Inorganic ion transport and metabolism",
    "Q":"Secondary metabolites biosynthesis, transport and catabolism",
    "R":"General function prediction only",
    "S":"Function unknown",
    "T":"Signal transduction mechanisms",
    "U":"Intracellular trafficking, secretion, and vesicular transport",
    "V":"Defense mechanisms",
    "W":"Extracellular structures",
    "Y":"Nuclear structure",
    "Z":"Cytoskeleton"
}

def add_cog_desc(df:pd.DataFrame):
    data = df[["cog_cat"]].copy()
    data['cog_desc'] = data.cog_cat.map(eggnong_mapp)
    if len(data[data.isna()]) > 0:
        for idx, row in data[data.cog_desc.isna()].iterrows():
            if row.cog_cat == '-':
                data.loc[idx, 'cog_desc'] = ""
            elif isinstance(row.cog_cat, str):
                val = "|".join([eggnong_mapp.get(char) for char in row.cog_cat])
                data.loc[idx, 'cog_desc'] = val
            # print(row)
    return data["cog_desc"]

def main():
    eggnog_output = "data/outputs/eggnog/outputs"
    to_concat = []
    for f in os.listdir(eggnog_output):
        if f.endswith('.tsv'):
            f = os.path.join(eggnog_output,f)
            _df = pd.read_csv(f, sep='\t', skiprows=4, usecols=['#query','COG_category','max_annot_lvl'])
            _df = _df.drop(_df.tail(n=3).index)            
            to_concat.append(_df)

    df = pd.concat(to_concat, ignore_index=True).rename(columns={"#query":"prot_id","COG_category":"cog_cat","max_annot_lvl":"deepst_cog"})
    df['cog_desc'] = df.pipe(add_cog_desc)
    df.to_csv("data/datasets/eggnog_categories.tsv",sep='\t')
if __name__ == "__main__":
    main()