#!/bin/bash
docker run -it --rm -v "$PWD":/usr/src/palmscan/analysis_data -w /usr/src/palmscan/analysis_data palmscan -search_pssms analysis_data/Virgaviridae_aln.txt -tsv analysis_results/virga_aln/palmscan_hits.tsv -report_pssms analysis_results/virga_aln/palmscan_report -fasta analysis_results/virga_aln/palmscan_pp.fa