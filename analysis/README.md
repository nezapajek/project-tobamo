# Analysis Pipeline

This folder contains the post-processing analysis pipeline and related code for preprocessing, training, and evaluating models, as well as for exploring and summarizing contig data from the main Snakemake workflow.

## Overview

The analysis pipeline is organized into 8 top-level components covering core analysis workflows and supporting data/reference resources:

```
analysis/
├── contigs_report/         # Contig inspection and filtering
├── model/                  # Machine learning pipeline 
├── palmprint/              # Viral palmprint domain identification
├── clustering/             # Contig clustering for phylogenetic analysis
├── phylogenetic_placement/ # Phylogenetic placement using EPA-ng
├── data/                   # Shared datasets and curated inputs
├── references/             # Reference sequence preparation workflows
└── supplementary_data/     # Supplementary tables and metadata assets
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

See [model/README.md](model/README.md) for detailed information.


### `palmprint/`

**Purpose:** Identify viral palmprint domains (conserved RNA-dependent RNA polymerase regions) in contig sequences.

**Key Components:**
- **ORF prediction:** Extract open reading frames using GETORF (EMBOSS) or ORFIPY
- **Domain scanning:** Search for viral palmprints using PalmScan and position-specific scoring matrices (PSSMs)
- **Sequence analysis:** Generate detailed reports and extract palmprint-containing sequences
- **Docker integration:** Containerized PalmScan workflow for reproducible analysis

**Input:** Nucleotide sequences from contigs (converted to amino acids via ORF prediction)
**Output:** TSV hit tables, detailed reports, and FASTA files with identified palmprint sequences

See [palmprint/README.md](palmprint/README.md) and [palmprint/protocol.md](palmprint/protocol.md) for detailed workflows.


### `clustering/`

**Purpose:** Cluster contigs obtained from the Snakemake workflow for downstream phylogenetic analysis.

**Key Functions:**
- Group similar contigs based on sequence similarity
- Prepare clustered sequences for phylogenetic tree construction
- Generate baseline for selecting representative sequences from clusters for further analysis
- Support downstream evolutionary and comparative genomics studies

**Input:** Contigs from Snakemake workflow output
**Output:** Clustered contig groups ready for phylogenetic analysis


### `phylogenetic_placement/`

**Purpose:** Phylogenetic placement of viral contigs onto reference trees using EPA-ng.

**Key Components:**
- **Reference preparation:** Align reference sequences and build phylogenetic tree with MAFFT and IQ-TREE
- **Phylogenetic placement:** Place query contigs onto reference tree using EPA-ng
- **Visualization:** Convert placement results to annotated trees with gappa
- **Orientation handling:** Tools to test both forward and reverse complement orientations

**Input:** Viral contigs from Snakemake workflow and curated reference sequences
**Output:** Phylogenetic placement results (.jplace) and annotated tree visualizations

See [phylogenetic_placement/README.md](phylogenetic_placement/README.md) and [phylogenetic_placement/protocol.md](phylogenetic_placement/protocol.md) for detailed workflows.


### `data/`

**Purpose:** Central storage for analysis-ready datasets, curated sequence collections, and domain-scientist-provided inputs shared across workflows.

**Key Contents:**
- `contigs/` and `serratus/` data resources
- `domain_sci_input/` curated input sets
- `non-virga_representatives/` model control datasets
- `tobamo/` project-specific sequence and metadata resources


### `references/`

**Purpose:** Workflows and scripts for generating and curating reference representatives used in downstream model and phylogenetic analyses.

**Key Functions:**
- Build curated Virgaviridae representative sets
- Generate non-virga control representatives
- Produce reference outputs for downstream placement and evaluation

See [references/README.md](references/README.md) for workflow details.


### `supplementary_data/`

**Purpose:** Supplementary deliverables and metadata support for reporting and manuscript-ready outputs.

**Key Contents:**
- `raw_files/` source supplementary files
- `files/` generated/organized supplementary artifacts
- notebook utilities for supplementary report assembly and variable dictionary generation
