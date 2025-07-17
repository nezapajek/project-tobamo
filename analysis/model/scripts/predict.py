from utils import *
import argparse
import os


def get_orfs_pairwise(train_or_test):
    # Implement the logic to extract ORFs and perform pairwise alignment
    contigs_path = "../data/contigs/contigs_all_deduplicated.fasta"
    outdir = "results/orfs" if train_or_test == "train" else "results/testing/orfs"

    # Check if the output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Call the existing script or function to get ORFs and pairwise alignments
    os.system(f"python scripts/03_getorfs_pairwise_aln.py {contigs_path} {outdir} forward")


def agg_pivot(train_or_test):
    # Implement the logic to aggregate and pivot the data
    refs_path = "../data/references/reference_database.xlsx"
    pairwise_path = f"results/{'training' if train_or_test == 'train' else 'testing'}/pairwise_aln.csv"
    orf_fasta_path = f"results/orfs/combined_orfs.fasta"
    output_dir = "training" if train_or_test == "train" else "testing"

    os.makedirs(f"results/{output_dir}", exist_ok=True)

    # Call the existing script or function to aggregate and pivot
    os.system(f"python scripts/02_preprocess_training.py {refs_path} {orf_fasta_path} {pairwise_path} {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Preprocessing pipeline for ORF extraction and data aggregation.")
    parser.add_argument("--getorfs-pairwise", action="store_true", help="Extract ORFs and perform pairwise alignment.")
    parser.add_argument("--agg-pivot", action="store_true", help="Aggregate and pivot the data.")
    parser.add_argument("--train", action="store_true", help="Specify if the data is for training.")
    parser.add_argument("--test", action="store_true", help="Specify if the data is for testing.")

    args = parser.parse_args()

    if args.getorfs_pairwise:
        if args.train:
            get_orfs_pairwise("train")
        elif args.test:
            get_orfs_pairwise("test")
        else:
            print("Please specify --train or --test.")

    if args.agg_pivot:
        if args.train:
            agg_pivot("train")
        elif args.test:
            agg_pivot("test")
        else:
            print("Please specify --train or --test.")


if __name__ == "__main__":
    main()
