python scripts/preprocess.py ../data/tobamo/reference_database.xlsx results/training/orfs/combined_orfs.fasta results/training/pairwise_aln.csv training --train --contigs results/training/sampling/2025-07-11_sampled_contigs_30.fasta

# update path/to/pairwise.csv
python scripts/preprocess.py ./data/tobamo/reference_database.xlsx results/snakemake/orfs/combined_orfs.fasta path/to/pairwise.csv snakemake --test