import os
import re
import tempfile
import pandas as pd

from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def parse_getorf_fasta_to_df(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((record.id, str(record.seq)))

    df = pd.DataFrame(sequences, columns=["orf_name", "seq"])

    def remove_last_number(s):
        return re.sub(r"_(\d+)$", "", s)

    df["contig"] = df["orf_name"].apply(lambda x: remove_last_number(x))
    df["SRR"] = df["contig"].str.extract(r"_([A-Za-z0-9]+)$")
    df["seq_len"] = df["seq"].apply(len)

    columns = ["SRR", "contig", "orf_name", "seq", "seq_len"]
    df = df[columns]
    return df


def parse_fasta(file_path):
    return SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))


def compute_identity_score_refxref(args):
    s1, s2 = args

    for aln in align_mafft(s1.seq, s2.seq):
        score, M, N, aln_len, seqA_len, seqB_len, n_aln_len, gaps, gap_ratio, aln_orf_len = score_alignment(aln)
        break
    return {
        "identity_score": score,
        "ref_name1": s1.id,
        "ref_name2": s2.id,
        "M": M,
        "N": N,
        "aln_len": aln_len,
        "ref_len1": seqA_len,
        "ref_len2": seqB_len,
        "N/aln_len": n_aln_len,
        "gap_openings": gaps,
        "gap_ratio": gap_ratio,
        "aln_orf_len": aln_orf_len,
    }


def compute_identity_score_orfxref(args):
    s1, s2 = args

    for aln in align_mafft(s1.seq, s2.seq):
        score, M, N, aln_len, seqA_len, seqB_len, n_aln_len, gaps, gap_ratio, aln_orf_len = score_alignment(aln)
        break
    return {
        "identity_score": score,
        "orf_name": s1.id.replace("=", "_"),
        "ref_name": s2.id,
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


def align_mafft(s1, s2):
    with tempfile.TemporaryDirectory() as tempdir:
        fin_name = tempdir + "/input.fasta"
        records = [SeqRecord(Seq(s1), id="seqA"), SeqRecord(Seq(s2), id="seqB")]
        SeqIO.write(records, fin_name, "fasta")
        os.system(f"mafft --quiet {fin_name} > {tempdir}/output.fasta")
        aln = AlignIO.read(f"{tempdir}/output.fasta", "fasta")
        aln.seqA = str(aln._records[0].seq.upper())
        aln.seqB = str(aln._records[1].seq.upper())
    return [aln]


def score_alignment(aln):
    M = sum(n1 != n2 for n1, n2 in zip(aln.seqA, aln.seqB) if n1 != "-" and n2 != "-")
    N = sum(1 for n1, n2 in zip(aln.seqA, aln.seqB) if n1 != "-" and n2 != "-")
    seqA_len = sum(1 for n1 in aln.seqA if n1 != "-")
    seqB_len = sum(1 for n2 in aln.seqB if n2 != "-")
    aln_len, aln_orf_len = calculate_alignment_length(aln)
    gaps = gap_openings(aln)
    score = 1 - M / N
    return score, M, N, aln_len, seqA_len, seqB_len, N / aln_len, gaps, gaps / aln_len, aln_orf_len


def gap_openings(aln):
    gapsA = sum(1 for i in range(len(aln.seqA) - 1) if aln.seqA[i] == "-" and aln.seqA[i - 1] != "-")
    gapsB = sum(1 for i in range(len(aln.seqB) - 1) if aln.seqB[i] == "-" and aln.seqB[i - 1] != "-")

    if aln.seqA.endswith("-") and gapsA != 0:
        gapsA -= 1
    if aln.seqB.endswith("-") and gapsB != 0:
        gapsB -= 1

    gaps = gapsA + gapsB
    return gaps


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


def df_ref_edit(results):
    df = pd.DataFrame(results)
    df = df[
        [
            "ref_name1",
            "ref_name2",
            "identity_score",
            "gap_openings",
            "gap_ratio",
            "N/aln_len",
            "aln_orf_len",
            "M",
            "N",
            "aln_len",
            "ref_len1",
            "ref_len2",
        ]
    ]
    return df


def df_orf_edit(results):
    df = pd.DataFrame(results)
    df["contig_name"] = [remove_last_number(name) for name in df["orf_name"]]
    df["SRR"] = df["contig_name"].str.extract(r"_([A-Za-z0-9]+)$")
    df = df[
        [
            "SRR",
            "contig_name",
            "orf_name",
            "ref_name",
            "identity_score",
            "gap_openings",
            "gap_ratio",
            "N/aln_len",
            "aln_orf_len",
            "M",
            "N",
            "aln_len",
            "orf_len",
            "ref_len",
        ]
    ]
    return df


def remove_last_number(s):
    return re.sub(r"_(\d+)$", "", s)


def compute_identity_score_print_aln(args):
    s1, s2 = args
    retl = []

    aln = align_mafft(s1.seq, s2.seq)[0]
    a = "".join([" ", "|"][n1 == n2] for n1, n2 in zip(aln.seqA, aln.seqB))
    aln_str = "\n".join([aln.seqA, a, aln.seqB])
    score, M, N, aln_len, _, _, n_aln_len, gaps, _ = score_alignment(aln)
    retl = [
        "**********************************************************************************\n",
        s1.id + "\n",
        s2.id + "\n",
        f"score: {score}" + "\n",
        f"M: {M}" + "\n",
        f"N: {N}" + "\n",
        f"aln_len: {aln_len}" + "\n",
        f"N/aln_len: {n_aln_len}" + "\n",
        f"gap_openings: {gaps}",
        "\n\n",
        aln_str,
        "\n\n",
    ]
    return retl
