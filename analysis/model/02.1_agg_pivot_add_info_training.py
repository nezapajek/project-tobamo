from utils import *

contigs = sys.argv[1]
outdir = sys.argv[2]

df = pd.read_csv(f"results/{outdir}/pairwise_aln.csv")
out = f"results/{outdir}/test_input_df.csv"

aa_id2type = mpu.io.read(f"data/aa_id2type.json")
aa_id2prot = mpu.io.read(f"data/aa_id2prot.json")
inv_nt_species2id_path = "data/inv_nt_dict.json"
inv_type_acc_dict_path = "data/inv_type_acc_dict.json"
orf_fasta_path = f"results/{outdir}/orfs/combined_orfs.fasta"

if os.path.exists(out):
    print(f"Skipping processing as {out} already exists.")
else:
    agg = aggregate_df(df, aa_id2type, aa_id2prot)
    pvtd = pivot_df(agg)
    pvtd_plus = add_info_to_input_df_v2(pvtd, orf_fasta_path, contigs)
    add_info = add_info_to_training_input_df(pvtd_plus, orf_fasta_path, inv_nt_species2id_path, inv_type_acc_dict_path)
    add_info.to_csv(out, index=False)
    print("Processed successfully")
