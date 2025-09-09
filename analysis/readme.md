# Analysis Pipeline

This folder contains the post-processing analysis pipeline and related code for preprocessing, training, and evaluating models, as well as for exploring and summarizing contig data from the main Snakemake workflow.

## Overview

The analysis pipeline is organized into 4 main components, each focused on a distinct part of the post-workflow analysis:

```
analysis/
├── contigs_report/    # Contig inspection and filtering
├── model/            # Machine learning pipeline 
```

## Folder Structure

### `contigs_report/`

**Purpose:** Inspect and process Snakemake output contigs for quality control and filtering.

**Key Functions:**
- Filter contigs to remove those associated with *cellular organisms*
- Exclude duplicated or low-quality runs (e.g., problematic samples like SRR6846476)
- Generate filtered FASTA files for downstream analysis
- Create supplementary data tables and summary reports


### `model/`

**Purpose:** Machine learning pipeline for automated viral sequence classification and prediction.

**Key Components:**
- **Data preprocessing:** Feature extraction from contig sequences
- **Model training:** Training classifiers on curated viral/non-viral datasets  
- **Evaluation:** Performance testing using cross-validation and held-out test sets
- **Prediction:** Classification of new contigs from the Snakemake pipeline

**Input:** Processed contigs from `contigs_report/` and `data/`
**Output:** Trained models, performance metrics, and classification predictions

See [model/readme.md](model/readme.md) for detailed information.


### Prerequisites

Install analysis dependencies:

```bash
# Option 1: Using pip
pip install -r analysis_requirements.txt

# Option 2: Using conda  
conda env create -f analysis_conda-requirements.txt
conda activate tobamo-analysis
```