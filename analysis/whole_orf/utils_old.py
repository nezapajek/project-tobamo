import os
import re
import tempfile
import pandas as pd

from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def parse_fasta(file_path):
    return SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))


def compute_identity_score_refxref(args):
    s1, s2 = args

    for aln in align_mafft(s1.seq, s2.seq):
        score, M, N, aln_len, orf_len, ref_len, n_orflen, gaps = score_alignment(aln)
        break
    return {
        "identity_score": score,
        "ref_name1": s1.id,
        "ref_name2": s2.id,
        "M": M,
        "N": N,
        "aln_len": aln_len,
        "orf_len": orf_len,
        "ref_len": ref_len,
        "N/orf_len": n_orflen,
        "gap_openings": gaps,
    }


def compute_identity_score_orfxref(args):
    s1, s2 = args

    for aln in align_mafft(s1.seq, s2.seq):
        score, M, N, aln_len, orf_len, ref_len, n_orflen, gaps = score_alignment(aln)
        break
    return {
        "identity_score": score,
        "orf_name": s1.id.replace("=", "_"),
        "ref_name": s2.id,
        "M": M,
        "N": N,
        "aln_len": aln_len,
        "orf_len": orf_len,
        "ref_len": ref_len,
        "N/orf_len": n_orflen,
        "gap_openings": gaps,
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
    orf_len = sum(1 for n1 in aln.seqA if n1 != "-")
    aln_len = sum(1 for n1, n2 in zip(aln.seqA, aln.seqB))
    ref_len = sum(1 for n2 in aln.seqB if n2 != "-")
    gaps = gap_openings(aln)
    score = 1 - M / N
    return score, M, N, aln_len, orf_len, ref_len, N / orf_len, gaps


def gap_openings(aln):
    gapsA = sum(1 for i in range(len(aln.seqA) - 1) if aln.seqA[i] == "-" and aln.seqA[i - 1] != "-")
    gapsB = sum(1 for i in range(len(aln.seqB) - 1) if aln.seqB[i] == "-" and aln.seqB[i - 1] != "-")

    if aln.seqA.endswith("-") and gapsA != 0:
        gapsA -= 1
    if aln.seqB.endswith("-") and gapsB != 0:
        gapsB -= 1

    gaps = gapsA + gapsB
    return gaps


def df_edit(results):
    df = pd.DataFrame(results)
    df = df[
        [
            "ref_name1",
            "ref_name2",
            "identity_score",
            "gap_openings",
            "N/orf_len",
            "M",
            "N",
            "aln_len",
            "orf_len",
            "ref_len",
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
            "N/orf_len",
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
