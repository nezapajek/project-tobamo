# using Conda
conda create --name <env_name> python=3.10 \
conda activate <env_name> \
pip install -r requirements.txt \
conda install bioconda::orfipy \
conda install -c bioconda mafft=7.520 \


# training the model
python 00_sample_refs.py <path/to/reference.fasta> <out_dir_name> <sampling_num> <subsampling_num> \
python 01_getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation> \
python 02.1_agg_pivot_add_info_training.py <path/to/contig.fasta> <out_dir_name> \
python 03.1_train_and_evaluate.py <out_dir_name> \
python 03.2_train.py <out_dir_name> <subsampling_num> \

# preprocessing and predicting the contigs
python 01_getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation> 
python 02.2_agg_pivot_add_info.py <path/to/contig.fasta> <out_dir_name> 
python 03.3_predict_and_report.py <path/to/input_df> <out_dir_name> #change dataset specific details in script

# example
python 00_sample_refs.py data/virga_nt.fasta training 100 25 \
python 01_getorfs_pairwise_aln.py results/training/sampled_contigs/sampled_contigs_25.fasta training unknown \
python 02.1_agg_pivot_add_info_training.py results/training/sampled_contigs/sampled_contigs_25.fasta training \
python 03.1_train_and_evaluate.py training \
python 03.2_train.py training 25