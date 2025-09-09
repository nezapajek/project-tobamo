# Snakemake workflow: `project-tobamo`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/nezapajek/project-tobamo/workflows/Tests/badge.svg?branch=main)](https://github.com/nezapajek/project-tobamo/actions?query=branch%3Amain+workflow%3ATests)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

## Overview

A comprehensive Snakemake workflow for the preparation of a curated catalogue of sequences of possible new tobamoviruses by scanning large accumulated datasets from different metagenomics data repositories.

**Project Website:** http://projects.nib.si/tobamo/

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Workflow](#workflow)
- [Output](#output)
- [Analysis](#analysis)
- [Documentation](#documentation)
- [Citation](#citation)

## Documentation

### 📚 Complete Documentation

- **🚀 [Quick Start Guide](QUICKSTART.md)** - Get running in 15 minutes
- **💾 [Installation Guide](INSTALLATION.md)** - Detailed setup instructions
- **⚙️ [Configuration Guide](config/README.md)** - Sample and parameter configuration
- **🔧 [Workflow Rules](workflow/RULES.md)** - Detailed workflow documentation
- **❓ [FAQ](FAQ.md)** - Frequently asked questions and troubleshooting
- **📊 [Analysis Pipeline](analysis/readme.md)** - Post-processing and machine learning

## Quick Start

```bash
# 1. Clone repository
git clone https://github.com/nezapajek/project-tobamo.git
cd project-tobamo

# 2. Install dependencies
conda env create -f workflow/envs/snakemake.yaml
conda activate tobamo-snakemake

# 3. Configure samples
# Edit config/samples_test.tsv with your SRA accessions

# 4. Test run
snakemake -n --configfile config/config.yaml

# 5. Full run
snakemake --use-conda -c32 -p -k
```

## Installation

### Prerequisites

- [Conda/Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [Snakemake ≥6.3.0](https://snakemake.readthedocs.io/)
- Minimum 32GB RAM recommended
- ~500GB free disk space for full analysis

### Environment Setup

```bash
# Create and activate Snakemake environment
conda env create -f workflow/envs/snakemake.yaml
conda activate tobamo-snakemake

# Optional: Install analysis dependencies for post-processing
pip install -r analysis_requirements.txt
```

## Configuration

### Sample Configuration

Edit `config/samples_*.tsv` files to specify your SRA accessions:

- `samples_debug.tsv` - Small debug dataset (2 samples, 1 single-end (SE) + 1 paired-end (PE))
- `samples_test.tsv` - Small test dataset (253 test samples)
- `samples_all.tsv` - Complete dataset (278 samples - 253 test + 25 control samples)

### Database Setup TK

Required databases are automatically downloaded, but you can pre-download:

```bash
# Download nr.dmnd database (~100GB)
wget -P resources/ https://example.com/nr.dmnd

# Download tpdb2.dmnd database
wget -P resources/ https://example.com/tpdb2.dmnd
```

## Usage

### Basic Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=nezapajek%2Fproject-tobamo).

### Command Examples

```bash
# Dry run to check workflow
snakemake -n

# Run with test samples
snakemake --use-conda -c8 --configfile config/config.yaml

# Full production run
time snakemake --use-conda -c32 -p -k > output_$(date +%Y-%m-%d).txt 2>&1

# Generate workflow report
snakemake --report report.html
```

## Workflow

The workflow consists of 5 main steps:

1. **Quality Control and Trimming** 
   - Tool: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
   - Removes adapters and low-quality sequences
   - Handles both paired-end and single-end reads

2. **De Novo Assembly**
   - Tools: [MEGAHIT](https://www.metagenomics.wiki/tools/assembly/megahit) and [SPAdes](https://github.com/ablab/spades)
   - Assembles trimmed reads into contigs
   - Optimized for metagenomic data

3. **Viral Protein Search**
   - Tool: [Diamond](https://bio.tools/diamond) BLASTx
   - Searches against tobamoviral protein sequences (tpdb2.dmnd)
   - Identifies potential viral contigs

4. **Taxonomic Search**
   - Tool: [Diamond](https://bio.tools/diamond) BLASTx
   - Searches against NCBI non-redundant (NR) protein database
   - Provides broader taxonomic context

5. **Taxonomy Assignment**
   - Tool: [MEGAN6](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/)
   - Assigns taxonomic classifications
   - Filters and curates results

### Workflow Diagram

```
SRA Data → Trimming → De novo assembly → Filtering → Diamond (tpdb2) → Diamond (NR) → MEGAN6 
```

## Output

### Main Output Files

- **`results/megan6_results_combined.csv`** - Combined results from all samples with taxonomic classifications
- **`results/{accession}/09_{accession}_megan6_results.csv`** - Individual sample results with:
  - Contig sequences and metadata
  - Viral protein search hits (tpdb2 database)
  - NR database taxonomic assignments
  - Quality scores and e-values
  - Contig length and coverage information
- **Intermediate files** - Complete processing chain preserved for reproducibility and debugging

### Output Structure

```
results/
├── megan6_results_combined.csv          # Main combined results from all samples
├── megan6_results_combined_add_nr_taxa.csv  # Combined results with NR taxonomy
└── {SAMPLE_ID}/                          # Individual sample directories (SRR*/ERR*/DRR*)
    │
    ├── 01_*_trim_*.done                  # Trimming completion flags
    ├── 01_*_trim_single.fq.gz            # Single-end trimmed reads (SE samples)
    ├── 01_*_trim_*_paired.fq.gz          # Paired trimmed reads (PE samples: R1, R2)
    ├── 01_*_trim_*_unpaired.fq.gz        # Unpaired trimmed reads (PE samples: R1, R2)
    │
    ├── 02_*_benchmark_*.txt               # Performance benchmarks
    ├── 02_*_spades_isolate_contigs.fasta # SPAdes assembly contigs
    ├── 02_*_megahit_contigs.fasta        # MEGAHIT assembly contigs  
    ├── 02_*_*.no_isolate                 # SPAdes failure flags
    ├── 02_*_*.no_megahit                 # MEGAHIT failure flags
    │
    ├── 03_*_contigs_combined.fasta       # Combined SPAdes + MEGAHIT contigs
    ├── 03_*_contigs_unique.fasta         # Deduplicated combined contigs
    │
    ├── 04_*_contigs_add_accession.fasta  # Contigs with sample ID prefixes
    ├── 04_*_contigs_filtered.fasta       # Length-filtered contigs (600-8000 bp)
    │
    ├── 05_*_benchmark_diamond_tpdb2.txt  # Diamond tpdb2 benchmark
    ├── 05_*_diamond_tpdb2.daa            # Tobamovirus protein search results
    │
    ├── 06_*_benchmark_diamond_nr.txt     # Diamond NR benchmark
    ├── 06_*_diamond_info.tsv             # tpdb2 search summary (tabular)
    ├── 06_*_diamond_tpdb2_selected.fasta # Contigs with viral hits (for NR search)
    ├── 06_*_diamond_nr.daa               # NR protein database search results
    ├── 06_*_diamond_nr_info.tsv          # NR search summary (tabular)
    │
    ├── 07_*_benchmark_meganizer_tpdb2.txt # MEGAN processing benchmark
    ├── 07_*_meganizer_tpdb2.daa          # MEGAN-processed search results
    │
    ├── 08_*_meganizer_tpdb2_read_classification.tsv  # Per-contig taxonomic classification
    ├── 08_*_meganizer_tpdb2_class_count.tsv          # Taxonomic summary counts
    │
    └── 09_*_megan6_results.csv           # Final integrated results (main output)
```

## Analysis

The `analysis/` folder contains post-processing scripts and notebooks:

- **`contigs_report/`** - Contig filtering and metadata analysis
- **`model/`** - Machine learning pipeline for viral classification  
- **`palmprint/`** - Sequence motif analysis tools
- **`data/`** - Processed FASTA files from Snakemake output and input by domain scientists
See [analysis/readme.md](analysis/readme.md) for detailed information.

### Log Files

Check logs in the `logs/` directory for detailed error information:
- `logs/trim_pe/` - Trimming logs
- `logs/megahit_*/` - Assembly logs  
- `logs/diamond_*/` - Database search logs

## Citation

If you use this workflow in a paper, please cite:

- This repository: `https://github.com/nezapajek/project-tobamo`
- The workflow DOI: [Add DOI when available]
- Related publication: [Add publication when available]

## License

This project is licensed under [LICENSE](LICENSE).

## Contact

For questions and support, please open an issue on GitHub or contact [contact information].