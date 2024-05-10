#!/bin/bash
blastn -query results/blast/other_contigs_139.fasta -db nt -out results/blast/blast_out_non_mv_only.tsv -outfmt 6 -num_threads 40 -max_target_seqs 10 -evalue 0.0001