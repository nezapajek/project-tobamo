from utils_v10 import *

contigs = sys.argv[1]
outdir = sys.argv[2]

df = pd.read_csv(f"results/{outdir}/pairwise_aln.csv")
out = f"results/{outdir}/test_input_df.csv"

aa_id2type = mpu.io.read(f"data/aa_id2type.json")
aa_id2prot = mpu.io.read(f"data/aa_id2prot.json")

if os.path.exists(out):
    print(f"Skipping processing as {out} already exists.")
else:
    agg = aggregate_df(df, aa_id2type, aa_id2prot)
    pvtd = pivot_df(agg)
    add_info = add_info_to_input_df_v2(pvtd, f"results/{outdir}/orfs/combined_orfs.fasta", contigs)
    add_info.to_csv(out, index=False)
    print("Processed successfully")
