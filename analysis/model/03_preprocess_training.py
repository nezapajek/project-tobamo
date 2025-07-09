from utils import *

refs_path = sys.argv[1]
contig_fasta_path = sys.argv[2]
orf_fasta_path = sys.argv[3]
pairwise_path = sys.argv[4]
output_dir = sys.argv[5]

# Check if files exist and are in the correct format
check_file_exists(refs_path, "Excel")
check_file_exists(pairwise_path, "CSV")
check_file_exists(contig_fasta_path, "FASTA")
check_file_exists(orf_fasta_path, "FASTA")
os.makedirs(f"results/{output_dir}", exist_ok=True)

# load data
refs = pd.read_excel(refs_path)
pairwise = pd.read_csv(pairwise_path)

# load data, filter out parent amino acid refs
print("Filtering out parent amino acid refs")
pairwise_minus, removed_refs = filter_pairwise(pairwise, refs)

# make mappers from refs
print("Making mappers")
aa_id2type, aa_id2prot, inv_type_acc_dict_path, inv_nt_virusname2id_path = make_mappers(refs)

# aggregate and pivot
print("Aggregating and pivoting")
agg = aggregate_df(pairwise_minus, aa_id2type, aa_id2prot)
pvtd = pivot_df(agg)

# add basic info
print("Adding info")
basic = add_info_basic(pvtd, orf_fasta_path, snakemake=False)

print("Adding info2")
# add info specific to training
input_df = add_info_to_training_input_df(basic, orf_fasta_path, inv_nt_virusname2id_path, inv_type_acc_dict_path)

# save
print("Saving")
input_df.to_csv(f"results/{output_dir}/training_input.csv", index=False)
removed_refs.to_csv(f"results/{output_dir}/removed_refs.csv", index=False)
