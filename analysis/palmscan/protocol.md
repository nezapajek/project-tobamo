0. snakemake output: 'project-tobamo/results/megan6_results_combined.csv'
1. make contig fasta from snakemake output
    01_make_contig_fasta.ipynb
2.1 run getorf (see requirements.txt)
    bash 02_getorf.sh
2.2 getorf stats
    02_getorf_stats.ipynb
3.1 run palmscan (see requirements.txt)
    bash 03_palmscan.sh
3.2 make palmscan positive contig fasta
    03_make_palmscan_pos_contig.ipynb
4. make palmscan positive contig fasta
    04_make_palmscan_pos_contig.ipynb