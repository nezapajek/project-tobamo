1. fitted a curve 

notebook: notebooks/01_fit_distribution_curve.ipynb
output: results/training/sampling/fitted_curve_lens_freq.json

2. sample references

python scripts/02_sample_refs.py <path/to/reference.fasta> <out_dir_name> <sampling_num> <subsampling_num> <path/to/lens_freq.json>
python scripts/02_sample_refs.py ../data/tobamo/reference_nukleotidne.fasta training 300 30 results/training/sampling/fitted_curve_lens_freq.json

3. pairwise alignment

time python 03_getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation>
time python scripts/03_getorfs_pairwise_aln.py results/training/sampling/2025-07-11_sampled_contigs_30.fasta training unknown

4. preprocess training data 

python 04_preprocess_training.py <path/to/reference_database.xlsx> <path/to/sampled_contigs.fasta> <path/to/orf.fasta> <path/to/pairwise_aln.csv> <output_dir>
python scripts/04_preprocess_training.py ../data/tobamo/reference_database.xlsx results/training/sampling/2025-07-11_sampled_contigs_30.fasta results/training/orfs/combined_orfs.fasta results/training/pairwise_aln.csv training

5. training TEST