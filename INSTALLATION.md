# Installation Guide

This guide provides detailed installation instructions for the project-tobamo workflow.

## System Requirements

### Hardware Requirements
- **RAM:** Minimum 32GB, recommended 64GB+
- **Storage:** ~500GB free space for full analysis
- **CPU:** Multi-core processor (8+ cores recommended)
- **Network:** High-speed internet for database downloads

### Software Requirements
- **Operating System:** Linux (tested on Ubuntu 18.04+, CentOS 7+)
- **Python:** 3.8-3.10
- **Conda/Miniconda:** Latest version
- **Git:** For repository cloning

## Step-by-Step Installation

### 1. Install Miniconda (if not already installed)

```bash
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install Miniconda
bash Miniconda3-latest-Linux-x86_64.sh

# Restart shell or source bashrc
source ~/.bashrc
```

### 2. Clone Repository

```bash
git clone https://github.com/nezapajek/project-tobamo.git
cd project-tobamo
```

### 3. Create Conda Environment

```bash
# Create environment with Snakemake
conda create -n tobamo-snakemake snakemake=7.32.4 python=3.10 -c conda-forge -c bioconda

# Activate environment
conda activate tobamo-snakemake
```

### 4. Install Additional Dependencies

**Note:** Snakemake will automatically handle all tool dependencies. Each workflow rule creates its own conda environment as defined in `workflow/envs/`. You don't need to manually install individual tools.

Dependencies are automatically managed by Snakemake when using the `--use-conda` flag

Optional: Install analysis dependencies for post-processing only.
See [analysis/readme.md](analysis/readme.md) for detailed information.
```bash
pip install -r analysis_requirements.txt
```

The workflow automatically installs and manages:
- Trimmomatic (quality control)
- MEGAHIT and SPAdes (assembly)
- Diamond (sequence search)
- MEGAN tools (taxonomic classification)
- BioPython and other utilities

### 5. Database Setup

#### Required Databases

**Note:** Databases need to be downloaded manually

1. **NCBI BLAST Database**
```bash
# Create database directory
mkdir -p blast_db

# Download all BLAST databases
wget --directory-prefix=blast_db --cut-dirs=2 -Anr* ftp://ftp.ncbi.nlm.nih.gov/blast/db/*
```

2. **Tobamovirus Protein Database (tpdb2)**
```bash
# Build the tobamovirus database tpdb2
diamond makedb --in <path/to/tobamo_proteins.fasta> -d <path/to/output/tpdb2.dmnd>

# Copy to resources directory
cp <path/to/tpdb2.dmnd> resources/tpdb2.dmnd
```

3. **MEGAN Mapping Database**
```bash
# Download MEGAN mapping files
wget -P resources/ https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-map-Feb2022.db
```

### 6. Test Installation

```bash
# Test with dry run
snakemake -n --configfile config/config.yaml

# Test with small dataset
cp config/samples_debug.tsv config/samples.tsv
snakemake --use-conda -c4 -p --configfile config/config.yaml
```

## Environment Management

### Creating Analysis Environment

If you plan to use the analysis pipeline:

```bash
# Option 1: From conda requirements
conda env create -f analysis_conda-requirements.txt -n tobamo-model

# Option 2: From pip requirements
conda create -n tobamo-model python=3.10
conda activate tobamo-model
pip install -r analysis_requirements.txt
```

### Switching Environments

```bash
# For running Snakemake workflow
conda activate tobamo-snakemake

# For running analysis scripts
conda activate tobamo-model
```

### Verification Steps

After installation, verify everything works:

```bash
# 1. Check Snakemake installation
snakemake --version

# 2. Test workflow syntax (this will also validate conda environments)
snakemake -n --configfile config/config.yaml

# 3. Run test dataset (tools will be installed automatically)
snakemake --use-conda -c2 -np --configfile config/config.yaml
```

**Note:** Individual tools (trimmomatic, megahit, spades, diamond, etc.) don't need to be manually verified since Snakemake manages them automatically in isolated conda environments.

## Performance Optimization

### For High-Performance Computing (HPC)

```bash
# Use cluster execution
snakemake --cluster "sbatch -c {threads} --mem={resources.mem_mb}" \
          --jobs 100 \
          --use-conda

# Or with profiles
snakemake --profile slurm
```

### For Local Workstations

```bash
# Optimize thread usage
snakemake --use-conda -c$(nproc) --resources mem_mb=32000

# Monitor resource usage
htop  # or similar system monitor
```

## Next Steps

After successful installation:

1. Review [Configuration Guide](config/README.md)
2. Run test workflow with small dataset
3. Configure your own sample list
4. Execute full analysis pipeline
