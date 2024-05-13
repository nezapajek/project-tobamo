### ORFIPY  (version 0.0.4)

conda create -n orfipy -c bioconda orfipy
conda activate orfipy

orfipy data/test_subset.fasta

orfipy data/test_contigs_2024-04-09_non_cellular_.fasta --min 300 --pep PEP --bed BED --between-stops --outdir results/