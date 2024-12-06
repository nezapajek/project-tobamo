# Tobamovirus Contig Classification with Random Forest

This repository provides Python code for training and using a Random Forest Classifier to predict Tobamoviral sequences in contigs. The workflow involves training the model on a reference dataset and then using it to classify contigs from a separate query set.

**Key Features:**

* Leverages the power of Random Forests for robust classification.
* Preprocesses contigs for optimal feature extraction.
* Supports contig orientation handling (unknown, forward, or reverse).
* Provides training and prediction scripts.

**Installation**

1.  **Create a conda environment:**

```bash
# using Conda
conda create --name <env_name> python=3.10
conda activate <env_name> 
pip install -r requirements.txt 
conda install bioconda::orfipy 
conda install -c bioconda mafft=7.520 
```
2. **TRAINING THE MODEL**

The code follows a modular approach with several scripts for different stages. Here's a breakdown of the steps:

2.1. Simulating Sequencing and Assembly for Training Data

We aimed to mimic the distribution of contig lengths typically observed in real-world sequencing and de novo genome assembly pipelines, such as our Snakemake pipeline. To achieve this, we fragmented reference genomes into contigs, ensuring that the resulting contig lengths closely matched the empirical distribution. This approach allowed us to generate a training dataset that more accurately reflects the characteristics of actual sequencing data.

```bash
python 00_sample_refs.py <path/to/reference.fasta> <out_dir_name> <sampling_num> <subsampling_num> 
```

2.2. Find ORFs and pairwise align them to our reference database using MAFFT

```bash
python 01_getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation>
```

2.3 Aggregating (min, max, mean) per featre and pivoting the table into desired format for further analysis
```bash
python 02.1_agg_pivot_add_info_training.py <path/to/contig.fasta> <out_dir_name> 
```

2.4. Training and evaluating the model
```bash
python 03.1_train_and_evaluate.py <out_dir_name> 
python 03.2_train.py <out_dir_name> <subsampling_num> 
```

3. **USING THE MODEL**

3.1.  preprocessing the contigs 
```bash
python 01_getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation> \
python 02.2_agg_pivot_add_info.py <path/to/contig.fasta> <out_dir_name> \
```

```bash
python 03.3_predict_and_report.py <path/to/input_df> <out_dir_name> #change dataset specific details in script
```

4. **example usage**
```bash
python 00_sample_refs.py data/virga_nt.fasta training 100 25 \
python 01_getorfs_pairwise_aln.py results/training/sampled_contigs/sampled_contigs_25.fasta training unknown \
python 02.1_agg_pivot_add_info_training.py results/training/sampled_contigs/sampled_contigs_25.fasta training \
python 03.1_train_and_evaluate.py training \
python 03.2_train.py training 25
```