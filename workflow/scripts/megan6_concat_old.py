import sys
import pandas as pd
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser

# python megan6_concat.py {input.csv} {input.fasta} {input.info_out} {output}

csv = sys.argv[1]
fasta_file = sys.argv[2]
diamond_info = sys.argv[3]
results = sys.argv[4]


# define function read_fasta
def read_fasta(file_path, columns):
    with open(file_path) as fasta_file:
        records = []
        for title, sequence in SimpleFastaParser(fasta_file):
            record = []
            record.append(title)
            record.append(sequence)
            records.append(record)
    return pd.DataFrame(records, columns=columns)


# read csvs
raw = pd.read_csv(csv, sep="\t", names=["qseqid", "taxonomy"])
raw_fasta = read_fasta(fasta_file, ["qseqid", "seq"])
df_info = pd.read_csv(
    diamond_info,
    sep="\t",
    names=[
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ],
)

# extract SRR
cols = ["SRR"]
contig_split = pd.DataFrame([re.findall("_([A-Za-z0-9]+)$", s) for s in raw.qseqid], columns=cols)

# mega merge
df = pd.merge(pd.merge(pd.concat([raw, contig_split], axis=1), df_info, on="qseqid"), raw_fasta, on="qseqid")

# save csv
if not df.empty:
    df.to_csv(results, index=False)
else:
    with open(results, "wt") as fout:
        pass
