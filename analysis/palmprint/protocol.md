    cd project-tobamo/analysis/palmprint
    mkdir results/palmscan

GETORF
    getorf -sequence ../data/contigs/contigs_non_cellular_filtered.fasta -outseq results/getorf/getorf_output_find1.fasta -table 0 -find 1 -minsize 270

    getorf -sequence ../data/contigs/contigs_non_cellular_filtered.fasta -outseq results/getorf/getorf_output_find0.fasta -table 0 -find 0 -minsize 270

ORFIPY
    orfipy ../data/contigs/contigs_non_cellular_filtered.fasta --min 270 --pep PEP --bed BED --between-stops --outdir results/orfipy

(getorf -find0 == orfipy --between-stops)

PALMSCAN
docker run -it --rm -v "$PWD":/usr/src/palmscan/data -w /usr/src/palmscan/data palmscan -search_pssms results/getorf/getorf_output_find1.fasta -tsv results/palmscan/palmscan_hits_find1.tsv -report_pssms results/palmscan/palmscan_report_find1 -fasta results/palmscan/palmscan_pp_find1.fa

docker run -it --rm -v "$PWD":/usr/src/palmscan/data -w /usr/src/palmscan/data palmscan -search_pssms results/getorf/getorf_output_find0.fasta -tsv results/palmscan/palmscan_hits_find0.tsv -report_pssms results/palmscan/palmscan_report_find0 -fasta results/palmscan/palmscan_pp_find0.fa