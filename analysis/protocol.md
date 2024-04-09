0. snakemake output: 'project-tobamo/results/megan6_results_combined.csv'
1. make contig fasta from snakemake output
    01_make_contig_fasta.ipynb
2. run getorf (see requirements.txt)
    bash 02_getorf.sh
3. run palmscan (see requirements.txt)
    bash 03_palmscan.sh
4. make palmscan positive contig fasta
    04_make_palmscan_pos_contig.ipynb