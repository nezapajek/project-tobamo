import ast
import glob
import os
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from itertools import product

import mpu
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import dump, load
from sklearn.discriminant_analysis import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, f1_score, precision_score, recall_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm


# translate sequence in all three reading frames
def translate_seq(seq_rec: SeqRecord):
    seq = seq_rec.seq
    # Check if the original sequence is divisible by 3, and trim if necessary
    if len(seq) % 3 != 0:
        seq = seq[: -(len(seq) % 3)]
    prot1 = seq.translate()  # First reading frame

    # Check if the sequence starting from the 2nd nucleotide is divisible by 3, and trim if necessary
    seq2 = seq[1:]
    if len(seq2) % 3 != 0:
        seq2 = seq2[: -(len(seq2) % 3)]
    prot2 = seq2.translate()  # Second reading frame

    # Check if the sequence starting from the 3rd nucleotide is divisible by 3, and trim if necessary
    seq3 = seq[2:]
    if len(seq3) % 3 != 0:
        seq3 = seq3[: -(len(seq3) % 3)]
    prot3 = seq3.translate()  # Third reading frame

    protein1 = SeqRecord(prot1, id=f"{seq_rec.id}_aa_frame1", description=seq_rec.description)
    protein2 = SeqRecord(prot2, id=f"{seq_rec.id}_aa_frame2", description=seq_rec.description)
    protein3 = SeqRecord(prot3, id=f"{seq_rec.id}_aa_frame3", description=seq_rec.description)

    return {protein1.id: protein1, protein2.id: protein2, protein3.id: protein3}


def reverse_complement_record(record):
    """
    Returns the reverse complement of a SeqRecord, retaining metadata.

    Parameters:
        record (SeqRecord): The input SeqRecord object.

    Returns:
        SeqRecord: A new SeqRecord object representing the reverse complement.
    """
    rc = record.reverse_complement()
    rc.id = record.id
    rc.name = record.name
    rc.description = record.description
    return rc


def translate_seq_6frames(seq_rec: SeqRecord):
    seq = seq_rec.seq
    # Forward strand translations
    # Check if the original sequence is divisible by 3, and trim if necessary
    if len(seq) % 3 != 0:
        seq = seq[: -(len(seq) % 3)]
    prot1 = seq.translate()  # First reading frame

    # Check if the sequence starting from the 2nd nucleotide is divisible by 3, and trim if necessary
    seq2 = seq[1:]
    if len(seq2) % 3 != 0:
        seq2 = seq2[: -(len(seq2) % 3)]
    prot2 = seq2.translate()  # Second reading frame

    # Check if the sequence starting from the 3rd nucleotide is divisible by 3, and trim if necessary
    seq3 = seq[2:]
    if len(seq3) % 3 != 0:
        seq3 = seq3[: -(len(seq3) % 3)]
    prot3 = seq3.translate()  # Third reading frame

    # Reverse strand translations
    rev_seq = reverse_complement_record(seq_rec).seq
    # Check if the reverse sequence is divisible by 3, and trim if necessary
    if len(rev_seq) % 3 != 0:
        rev_seq = rev_seq[: -(len(rev_seq) % 3)]
    rev_prot1 = rev_seq.translate()  # First reading frame

    # Check if the reverse sequence starting from the 2nd nucleotide is divisible by 3, and trim if necessary
    rev_seq2 = rev_seq[1:]
    if len(rev_seq2) % 3 != 0:
        rev_seq2 = rev_seq2[: -(len(rev_seq2) % 3)]
    rev_prot2 = rev_seq2.translate()  # Second reading frame

    # Check if the reverse sequence starting from the 3rd nucleotide is divisible by 3, and trim if necessary
    rev_seq3 = rev_seq[2:]
    if len(rev_seq3) % 3 != 0:
        rev_seq3 = rev_seq3[: -(len(rev_seq3) % 3)]
    rev_prot3 = rev_seq3.translate()  # Third reading frame

    protein1 = SeqRecord(prot1, id=f"{seq_rec.id}_aa_frame1", description=seq_rec.description)
    protein2 = SeqRecord(prot2, id=f"{seq_rec.id}_aa_frame2", description=seq_rec.description)
    protein3 = SeqRecord(prot3, id=f"{seq_rec.id}_aa_frame3", description=seq_rec.description)
    rev_protein1 = SeqRecord(rev_prot1, id=f"{seq_rec.id}_aa_frame4", description=seq_rec.description)
    rev_protein2 = SeqRecord(rev_prot2, id=f"{seq_rec.id}_aa_frame5", description=seq_rec.description)
    rev_protein3 = SeqRecord(rev_prot3, id=f"{seq_rec.id}_aa_frame6", description=seq_rec.description)

    return {
        protein1.id: protein1,
        protein2.id: protein2,
        protein3.id: protein3,
        rev_protein1.id: rev_protein1,
        rev_protein2.id: rev_protein2,
        rev_protein3.id: rev_protein3,
    }


# run orfipy
def run_orfipy(contig, log_file_dir):
    with tempfile.TemporaryDirectory() as tempdir:
        fasta_path = f"{tempdir}/input.fasta"
        SeqIO.write(contig, fasta_path, "fasta")
        os.makedirs(log_file_dir, exist_ok=True)
        log_file_path = os.path.join(log_file_dir, "log_file")
        os.system(
            f"orfipy {fasta_path} --min 300 --pep PEP --bed BED --between-stops --outdir {tempdir} > {log_file_path} 2>&1"
        )
        orfs = SeqIO.to_dict(SeqIO.parse(f"{tempdir}/PEP", "fasta"))
        return orfs


# pairwise aln process in chunks
def process_in_chunks(args_list, chunk_size: int, max_workers: int):
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for i in range(0, len(args_list), chunk_size):
            chunk = args_list[i : i + chunk_size]
            results.extend(tqdm(executor.map(compute_identity_score_orfxref, chunk), total=len(chunk)))
    return results


# compute identity score
def compute_identity_score_orfxref(args):
    s1, s2 = args

    for aln in align_mafft(s1.seq, s2.seq):
        score, M, N, aln_len, seqA_len, seqB_len, n_aln_len, gaps, gap_ratio, aln_orf_len = score_alignment(aln)
        break
    return {
        "orf_name": s1.id.replace("=", "_"),
        "ref_name": s2.id,
        "identity_score": score,
        "M": M,
        "N": N,
        "aln_len": aln_len,
        "orf_len": seqA_len,
        "ref_len": seqB_len,
        "N/aln_len": n_aln_len,
        "gap_openings": gaps,
        "gap_ratio": gap_ratio,
        "aln_orf_len": aln_orf_len,
    }


# align sequences with MAFFT
def align_mafft(s1, s2):
    with tempfile.TemporaryDirectory() as tempdir:
        fin_name = tempdir + "/input.fasta"
        records = [SeqRecord(Seq(s1), id="seqA"), SeqRecord(Seq(s2), id="seqB")]
        SeqIO.write(records, fin_name, "fasta")
        ec = os.system(f"mafft --quiet {fin_name} > {tempdir}/output.fasta")
        if ec != 0:
            raise RuntimeError("Error in MAFFT")
        aln = AlignIO.read(f"{tempdir}/output.fasta", "fasta")
        aln.seqA = str(aln._records[0].seq.upper())
        aln.seqB = str(aln._records[1].seq.upper())
    return [aln]


# score alignment
def score_alignment(aln):
    M = sum(n1 != n2 for n1, n2 in zip(aln.seqA, aln.seqB) if n1 != "-" and n2 != "-")
    N = sum(1 for n1, n2 in zip(aln.seqA, aln.seqB) if n1 != "-" and n2 != "-")
    seqA_len = sum(1 for n1 in aln.seqA if n1 != "-")
    seqB_len = sum(1 for n2 in aln.seqB if n2 != "-")
    aln_len, aln_orf_len = calculate_alignment_length(aln)
    gaps = gap_openings(aln)
    score = 1 - M / N
    return score, M, N, aln_len, seqA_len, seqB_len, N / aln_len, gaps, gaps / aln_len, aln_orf_len


# calculate alignment length
def calculate_alignment_length(aln):
    match_positions = [i for i, (n1, n2) in enumerate(zip(aln.seqA, aln.seqB)) if n1 != "-" and n2 != "-"]
    if match_positions:
        first_match = match_positions[0]
        last_match = match_positions[-1] + 1
        aln_len = last_match - first_match
    else:
        aln_len = 0
    aln_orf_len = aln_len - sum(1 for n1 in aln.seqA[match_positions[0] : match_positions[-1] + 1] if n1 == "-")
    return aln_len, aln_orf_len


# calculate gap openings
def gap_openings(aln):
    gapsA = sum(1 for i in range(len(aln.seqA) - 1) if aln.seqA[i] == "-" and aln.seqA[i - 1] != "-")
    gapsB = sum(1 for i in range(len(aln.seqB) - 1) if aln.seqB[i] == "-" and aln.seqB[i - 1] != "-")

    if aln.seqA.endswith("-") and gapsA != 0:
        gapsA -= 1
    if aln.seqB.endswith("-") and gapsB != 0:
        gapsB -= 1

    gaps = gapsA + gapsB
    return gaps


def aggregate_df(df, aa_id2type, aa_id2prot):
    df = df.copy()
    # Apply mappings to create new columns
    df["ref_type"] = df["ref_name"].map(aa_id2type)
    df["ref_prot_type"] = df["ref_name"].map(aa_id2prot)

    # Group by 'orf_name', 'ref_type', and 'ref_prot_type' and aggregate min and max values
    grouped = df.groupby(["orf_name", "ref_type", "ref_prot_type"]).agg(
        {
            "identity_score": ["min", "max", "mean"],
            "gap_openings": ["min", "max", "mean"],
            "gap_ratio": ["min", "max", "mean"],
            "N/aln_len": ["min", "max", "mean"],
            "aln_orf_len": ["min", "max", "mean"],
            "M": ["min", "max", "mean"],
            "N": ["min", "max", "mean"],
            "aln_len": ["min", "max", "mean"],
        }
    )

    # Flatten the multi-index column names
    grouped.columns = ["{}_{}".format(col[0], col[1]) for col in grouped.columns]

    # Reset index to make 'orf_name' a column again
    new_grouped = grouped.reset_index()

    assert not new_grouped.isnull().values.any(), "There are NaN values in the dataframe"

    # return new_grouped
    return new_grouped


def pivot_df(
    df,
    orf_name_col="orf_name",
    ref_type_col="ref_type",
    ref_prot_type_col="ref_prot_type",
):
    ref_types = df[ref_type_col].unique()
    ref_prot_types = df[ref_prot_type_col].unique()
    groups = list(product(ref_types, ref_prot_types))

    df_groups = {}
    for g in groups:
        temp_df = df.copy()
        temp_df = temp_df[(temp_df[ref_type_col] == g[0]) & (temp_df[ref_prot_type_col] == g[1])]
        temp_df = temp_df.set_index(orf_name_col)
        temp_df = temp_df.drop(columns=[ref_type_col, ref_prot_type_col])
        temp_df.columns = [f"{g[0]}_{g[1]}_{c}" for c in temp_df.columns]
        df_groups[g] = temp_df.copy()

    df_pivoted = df_groups[groups[0]]
    for g in groups[1:]:
        df_pivoted = df_pivoted.join(df_groups[g], how="outer")

    df_pivoted = df_pivoted.reset_index()

    return df_pivoted


def add_info_to_training_input_df(
    input_df, inv_nt_species2id_path: str, inv_type_acc_dict_path: str, orf_fasta_path: str
):

    if "contig_name" not in input_df.columns:
        input_df["contig_name"] = input_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)")

    # add nt_id info
    input_df["nt_id"] = input_df["orf_name"].str.extract(r"(.*)_start")

    # add species info
    inv_nt_species2id = mpu.io.read(inv_nt_species2id_path)
    input_df["species"] = input_df["nt_id"].map({k: v[0] for k, v in inv_nt_species2id.items()})

    # add orf_type info
    inv_type_acc_dict = mpu.io.read(inv_type_acc_dict_path)
    input_df["orf_type"] = input_df["nt_id"].map({k: v[0] for k, v in inv_type_acc_dict.items()})

    # add strand orientation info
    with open(orf_fasta_path, "r") as file:
        orf_dict = {
            record.description.split()[0]: "REVERSE" if "(-)" in record.description else "FORWARD"
            for record in SeqIO.parse(file, "fasta")
        }
    input_df["strand"] = input_df["orf_name"].map(orf_dict)

    return input_df


def add_info_to_input_df_v2(input_df, orf_fasta_path: str, contig_fasta_path: str):
    input_df["contig_name"] = input_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)", expand=False)

    # ADD 'orf_len'
    orf_lens = {k.replace("=", "_"): len(v.seq) for k, v in SeqIO.to_dict(SeqIO.parse(orf_fasta_path, "fasta")).items()}
    input_df["orf_len"] = input_df["orf_name"].map(orf_lens)
    # ADD 'contig_length'
    contig_lens = {
        k.replace("=", "_"): len(v.seq) for k, v in SeqIO.to_dict(SeqIO.parse(contig_fasta_path, "fasta")).items()
    }
    input_df["contig_length"] = input_df["contig_name"].map(contig_lens)

    # add orf_acc_method info
    input_df["orf_acc_method"] = np.where(input_df["orf_name"].str.contains("frame"), "translate", "orfipy")
    orf_acc_mapping = {"orfipy": 0, "translate": 1}
    input_df["orf_acc_method"] = input_df["orf_acc_method"].map(orf_acc_mapping)

    # add numbers of stop codons info
    all_orfs_stops = {
        k.replace("=", "_"): int(v.seq.count("*"))
        for k, v in SeqIO.to_dict(SeqIO.parse(orf_fasta_path, "fasta")).items()
    }  # NAMING '=' to '_'

    # Check if all dictionary keys are in the 'orf_name' column of input_df
    missing_keys = set(all_orfs_stops.keys()) - set(input_df["orf_name"])

    if missing_keys:
        print(f"Error: The following ORF names are missing in the 'orf_name' column of input_df: {missing_keys}")
        sys.exit(1)

    # Map the number of stop codons to the 'orf_name' column
    input_df["stop_codons"] = input_df["orf_name"].map(all_orfs_stops)
    return input_df


def add_info_to_training_input_df(
    input_df, orf_fasta_path: str, inv_nt_species2id_path: str, inv_type_acc_dict_path: str
):

    if "contig_name" not in input_df.columns:
        input_df["contig_name"] = input_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)")

    input_df["nt_id"] = input_df["orf_name"].str.extract(r"(.*)_start")

    # add species info
    inv_nt_species2id = mpu.io.read(inv_nt_species2id_path)
    input_df["species"] = input_df["nt_id"].map({k: v[0] for k, v in inv_nt_species2id.items()})

    # add orf_type info
    inv_type_acc_dict = mpu.io.read(inv_type_acc_dict_path)
    input_df["orf_type"] = input_df["nt_id"].map({k: v[0] for k, v in inv_type_acc_dict.items()})

    # add strand orientation info
    with open(orf_fasta_path, "r") as file:
        orf_dict = {
            record.description.split()[0]: "REVERSE" if "(-)" in record.description else "FORWARD"
            for record in SeqIO.parse(file, "fasta")
        }
    input_df["strand"] = input_df["orf_name"].map(orf_dict)

    return input_df


def most_extreme_prob(group):
    return group.loc[np.abs(group["probability"] - 0.5).idxmax()]


# contig sampling function
def sampling(seq_record, n, rng: np.random.Generator):
    lens_freq = mpu.io.read("/home/tobamo/analize/model-tobamo/data/lens_freq.json")
    trimmed_seqs = []
    random_lengths = rng.choice(a=list(lens_freq.keys()), p=list(lens_freq.values()), size=n).astype(int)
    seq_length = len(seq_record.seq)
    for rand_length in random_lengths:
        if rand_length < seq_length:
            start = rng.integers(0, seq_length - rand_length)
            trimmed_seq = seq_record.seq[start : start + rand_length]
            trimmed_seqs.append(
                SeqRecord(
                    trimmed_seq,
                    id=f"{seq_record.id}_start-{start}_len-{rand_length}",
                    description=seq_record.description,
                )
            )
        else:
            # If the random length is greater than the sequence length, keep the original sequence
            trimmed_seqs.append(
                SeqRecord(
                    seq_record.seq,
                    id=f"{seq_record.id}_start-0_len-{len(seq_record.seq)}_NT",
                    description=seq_record.description,
                )
            )
    contigs_dict = {trimmed_seq.id: trimmed_seq for trimmed_seq in trimmed_seqs}
    return contigs_dict, len(contigs_dict)


# subsample contigs to the same number as the genome with least contigs
def subsample_contigs_dicts(list_of_dicts: list, selected_n: int, rng: np.random.Generator):
    # Create a combined dictionary to hold the final result
    combined_dict = {}
    os.makedirs("training/sampling_reports", exist_ok=True)
    report_filename = f"training/sampling_reports/sampled_contigs_report.txt"

    # Open a file to write the report
    with open(report_filename, "w") as report_file:
        # Write header for the report
        report_file.write("Contig Subsampling Report\n")
        report_file.write("==================\n")

        # List to store the sizes of dictionaries
        sizes = []

        # Iterate over each dictionary and subsample
        for _, d in enumerate(list_of_dicts):
            dict_size = len(d)  # Get the size of the current dictionary
            sizes.append(dict_size)  # Add the size to the list

            # Calculate the sample size (ensure it doesn't exceed the dictionary size)
            sample_size = min(selected_n, dict_size)

            # Randomly sample `sample_size` keys
            sampled_keys = rng.choice(list(d.keys()), size=sample_size, replace=False)

            # Add the sampled items to the combined dictionary
            for key in sampled_keys:
                combined_dict[key] = d[key]

        # Find the minimum size across all dictionaries
        report_file.write(f"Selected sampling_n: {selected_n}\n")
        report_file.write(f"All contig lenghts: {sizes}\n")
        report_file.write(f"Median dictionary size: {np.median(sizes)}\n")
        report_file.write(f"Mean dictionary size: {np.mean(sizes)}\n")
        report_file.write(f"Minimum dictionary size: {min(sizes)}\n")

    return combined_dict


def process_species(species, input_df):
    print(f"\nProcessing for species: {species}")

    # make sure 'orf_name' is the index
    if input_df.index.name != "orf_name":
        input_df = input_df.set_index("orf_name")

    # remove species from training data
    test_df_all = input_df[input_df["species"] == species]

    # Only forward strand for training
    train_data_fwd = input_df[input_df["strand"] == "FORWARD"]

    # Feature extraction for training
    X_train = train_data_fwd.drop(columns=["orf_type", "strand", "species", "nt_id", "contig_name"])
    y_train = (train_data_fwd["orf_type"] == "tobamo").astype(
        int
    )  # Convert 'orf_type' to binary (1 for 'tobamo', 0 otherwise)

    # Standardize training and test data
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)

    # Initialize the model
    model = RandomForestClassifier(n_estimators=125, max_depth=40, n_jobs=-1)

    # Fit the model
    model.fit(X_train, y_train)

    # Return the model and test data for further evaluation
    return model, test_df_all


def process_df(df, input_df):
    df["orf_name"] = df["orf_name"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    df["ground_truth"] = df["ground_truth"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    df["prediction"] = df["prediction"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    df["probability"] = df["probability"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    columns_to_explode = ["orf_name", "ground_truth", "prediction", "probability"]
    res_df = df.explode(columns_to_explode).reset_index(drop=True)
    res_df["true_positive"] = (res_df["ground_truth"] == 1) & (res_df["prediction"] == 1)
    res_df["false_positive"] = (res_df["ground_truth"] == 0) & (res_df["prediction"] == 1)
    res_df["true_negative"] = (res_df["ground_truth"] == 0) & (res_df["prediction"] == 0)
    res_df["false_negative"] = (res_df["ground_truth"] == 1) & (res_df["prediction"] == 0)

    res_df["prob_0"] = res_df["probability"].apply(lambda x: x[0])
    res_df["prob_1"] = res_df["probability"].apply(lambda x: x[1])
    res_df["max_prob"] = res_df[["prob_0", "prob_1"]].max(axis=1)

    res_df = res_df.copy()
    res_df.loc[:, "match"] = res_df["ground_truth"] == res_df["prediction"]

    if input_df.index.name != "orf_name":
        input_df = input_df.set_index("orf_name")
    orf_orientation_mapper = input_df["strand"].to_dict()
    res_df["strand"] = res_df["orf_name"].map(orf_orientation_mapper)

    contig_length_mapper = input_df["contig_length"].to_dict()
    res_df["contig_length"] = res_df["orf_name"].map(contig_length_mapper)

    confusion_matrix = res_df[["true_positive", "true_negative", "false_positive", "false_negative"]].sum()

    return res_df, confusion_matrix


def select_best_orf(res_df):
    res_df["best_orf"] = res_df.groupby("contig")["max_prob"].transform(lambda x: x == x.max()).astype(int)
    res_df = res_df[res_df["best_orf"] == 1]
    res_df = res_df.copy()
    res_df["rank"] = res_df.groupby("contig")["max_prob"].rank(method="first", ascending=False)
    res_df = res_df[res_df["rank"] == 1]
    return res_df
