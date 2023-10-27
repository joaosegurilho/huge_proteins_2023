import sys
import pandas as pd

def main(input_file,output_file, filter=False):
	if not filter:
		to_use_df = pd.read_csv(input_file, sep='\t', usecols=['prot_id','prot_seq'])
	else:
		org_meta = pd.read_csv("/data/datasets/organism_metadata_dataset.tsv", sep='\t')
		big_df = pd.read_csv(input_file, sep='\t')
		if filter in ["eukaryota","bacteria","archaea"]:
			level = "superkingdom"
		else:
			level = "phylum"
		list_of_orgs = org_meta.query(f"{level} == '{filter}'").org_id.to_list() # Filter by taxnomic level desired
		list_of_prots = big_df.query(f"org_id in {list_of_orgs}").prot_id.to_list() # Get the proteins in that taxonomy level
		to_use_df = big_df.query(f"prot_id in {list_of_prots}")

	#print(df.columns)
	with open(output_file, 'w') as fasta_file:
		for idx, row in to_use_df.iterrows():
			fasta_file.write(f">{row.prot_id}\n{row.prot_seq}\n")
			# print(prot_id)
			

if __name__ == '__main__':
	if len(sys.argv) > 1:
		input_file = sys.argv[1]
		output_file = sys.argv[2]
		try:
			filter_by = sys.argv[3]
		except:
			main(input_file,output_file)
		else:
			main(input_file,output_file,filter_by)
			
	else:
		print("Need an output_file")