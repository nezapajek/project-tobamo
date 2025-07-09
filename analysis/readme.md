# Analysis

This folder contains the analysis pipeline and related code for preprocessing, training, and evaluating models, as well as for exploring and summarizing contig data. It is organized into X subfolders, each focused on a distinct part of the workflow.

## Folder Structure

### `contigs_report/`
This folder contains scripts and notebooks used to inspect and process the Snakemake output contigs. It includes steps for:

- Filtering contigs (e.g., removing those associated with *cellular organisms*)
- Excluding duplicated or low-quality runs (e.g., SRR6846476)
- Generating FASTA files for downstream use
- Creating supplementary data tables for reporting

### `model/`
This folder contains code for the machine learning pipeline, including:

- Data preprocessing and feature extraction
- Model training
- Evaluation and performance testing

The model is trained using contig features and tested on curated datasets.

### `palmprint/`
Contains scripts and tools for analyzing contig *palmprints*. These are specific sequence motifs or patterns used for viral classification or similarity-based filtering.

### `data/`
This folder stores .fasta files of snakemake contigs