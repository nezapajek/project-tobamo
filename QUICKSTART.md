# Quick Start Guide

Get up and running with project-tobamo in 15 minutes!

## Prerequisites

- Linux system with 32GB+ RAM
- Conda/Miniconda installed
- ~50GB free disk space (for test run)

## 1. Installation (5 minutes)

```bash
# Clone repository
git clone https://github.com/nezapajek/project-tobamo.git
cd project-tobamo

# Create conda environment
conda create -n tobamo-snakemake snakemake=7.32.4 python=3.10 -c conda-forge -c bioconda
conda activate tobamo-snakemake
```

## 2. Quick Test (10 minutes)

```bash
# Run with test dataset (few samples)
cp config/samples_test.tsv config/current_samples.tsv

# Edit config to use test samples
echo "samples: config/current_samples.tsv" > config/config.yaml

# Dry run to check workflow
snakemake -n --configfile config/config.yaml

# Run test (should complete in ~10 minutes)
snakemake --use-conda -c4 -p --configfile config/config.yaml
```

## 3. View Results

```bash
# Check main results file
head results/megan6_results_combined.csv

# Check individual sample results
ls results/*/09_*_megan6_results.csv
```

## 4. Full Production Run

```bash
# Use complete dataset (279 samples)
echo "samples: config/samples_all.tsv" > config/config.yaml

# Run full analysis (will take days)
nohup snakemake --use-conda -c32 -p --configfile config/config.yaml > output.log 2>&1 &

# Monitor progress
tail -f output.log
```

## Next Steps

1. **Explore Results:** Check the analysis pipeline in `analysis/`
2. **Customize:** Modify sample lists in `config/`
3. **Scale Up:** Use HPC cluster for large datasets
4. **Analysis:** Run post-processing scripts for publication

## Troubleshooting

### Quick Fixes

```bash
# Memory issues - reduce threads
snakemake --use-conda -c8 -p --configfile config/config.yaml

# Restart failed jobs
snakemake --use-conda -c32 -p --configfile config/config.yaml --rerun-incomplete

# Check specific sample
snakemake --use-conda -c4 -p results/SRR1234567/09_SRR1234567_megan6_results.csv
```

### Common Issues

1. **"Command not found"** → Activate conda environment
2. **"Out of memory"** → Reduce thread count or increase system RAM
3. **"Database not found"** → Databases download automatically on first run

## Getting Help

- 📖 **Full Documentation:** [README.md](README.md)
- 🛠️ **Installation Issues:** [INSTALLATION.md](INSTALLATION.md)  
- ⚙️ **Configuration:** [config/README.md](config/README.md)
- 🐛 **Bug Reports:** [GitHub Issues](https://github.com/nezapajek/project-tobamo/issues)