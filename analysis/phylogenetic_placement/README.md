# Phylogenetic Placement Environment

This directory contains the setup for a dedicated conda environment for phylogenetic placement analysis of tobamovirus sequences.

## Environment Specifications

- **Python**: 3.10
- **Bioinformatics Tools**:
  - INFERNAL 1.1.5 (RNA alignment inference)
  - pplacer 1.1.alpha20 (phylogenetic placement)
  - EPA-ng 0.3.8 (evolutionary placement algorithm)
  - IQ-TREE (phylogenetic inference)
  - MAFFT (multiple sequence alignment)
  - FastTree (approximate phylogenetic trees)
- **Python Packages**:
  - BioPython (sequence analysis)
  - pandas, numpy, matplotlib, seaborn (data analysis)
  - applese (phylogenetic placement analysis)
  - Jupyter notebook support

## Installation Options

### Option 1: Using the setup script (recommended)
```bash
cd /home/tobamo/analize/project-tobamo/analysis/phylogenetic_placement
./setup_environment.sh
```

### Option 2: Using conda environment file
```bash
cd /home/tobamo/analize/project-tobamo/analysis/phylogenetic_placement
conda env create -f environment.yml
```

### Option 3: Manual installation
```bash
# Create environment
conda create -n phylo_placement python=3.10 -y
conda activate phylo_placement

# Add channels
conda config --add channels bioconda
conda config --add channels conda-forge

# Install packages
conda install -c bioconda infernal=1.1.5 pplacer=1.1.alpha20 epa-ng=0.3.8 -y
conda install biopython pandas numpy matplotlib seaborn jupyter -y
pip install applese
```

## Verification

After installation, verify everything is working:
```bash
./verify_environment.sh
```

## Usage

1. Activate the environment:
```bash
conda activate phylo_placement
```

2. Start Jupyter notebook:
```bash
jupyter notebook
```

3. Select the "Phylogenetic Placement (Python 3.10)" kernel in your notebooks

## Tools Overview

- **INFERNAL**: For creating and using covariance models for RNA alignment
- **pplacer**: Maximum likelihood phylogenetic placement of sequences
- **EPA-ng**: Next-generation evolutionary placement algorithm
- **applese**: Python package for phylogenetic placement analysis and visualization

## Project Structure

```
phylogenetic_placement/
├── requirements.txt          # Python package requirements
├── environment.yml          # Conda environment specification
├── setup_environment.sh     # Automated setup script
├── verify_environment.sh    # Installation verification
├── README.md               # This file
├── data/                   # Input data files
└── pplacer/                # pplacer-specific analysis
    ├── pplacer.ipynb       # Main analysis notebook
    └── results/            # Output files
```

## Troubleshooting

If installation fails:

1. **Channel issues**: Ensure bioconda and conda-forge channels are added
2. **Version conflicts**: Try installing tools individually
3. **applese installation**: This requires pip and may need C++ compiler
4. **Permissions**: Ensure you have write access to conda environments

For specific tool issues, consult:
- [INFERNAL documentation](http://eddylab.org/infernal/)
- [pplacer documentation](https://matsen.fhcrc.org/pplacer/)
- [EPA-ng documentation](https://github.com/Pbdas/epa-ng)