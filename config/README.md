# Configuration Guide

This directory contains configuration files for the tobamo virus detection workflow.

## Configuration Files

### `config.yaml`

Main configuration file that specifies:
- Sample file location (`samples: config/samples_all.tsv`)
- Additional workflow parameters

```yaml
samples: config/samples_all.tsv  # Path to sample list
# Add other configuration parameters as needed
```

### Sample Files

Different sample files are provided for various use cases:

| File | Description | Samples | Use Case |
|------|-------------|---------|----------|
| `samples_debug.tsv` | Debug dataset | 2 | Troubleshooting |
| `samples_12.tsv` | Small dataset | 12 | Development and debugging |
| `samples_test.tsv` | Test dataset | 253 | Test samples |
| `samples_all.tsv` | Complete dataset | 279 | Test and control samples |

### Sample File Format

Sample files are tab-separated with a single column header:

```tsv
samples
SRR1234567
ERR2345678
DRR3456789
```

**Requirements:**
- First line must be `samples` (header)
- One SRA accession per line
- Supported prefixes: SRR, ERR, DRR
- No empty lines or comments

## Usage Examples

### Basic Configuration

1. **Choose appropriate sample file:**
   ```bash
   # For testing
   cp config/samples_test.tsv config/my_samples.tsv
   
   # For production
   cp config/samples_all.tsv config/my_samples.tsv
   ```

2. **Edit config.yaml:**
   ```yaml
   samples: config/my_samples.tsv
   ```

### Custom Sample List

1. **Create custom sample file:**
   ```bash
   echo "samples" > config/custom_samples.tsv
   echo "SRR1234567" >> config/custom_samples.tsv
   echo "ERR2345678" >> config/custom_samples.tsv
   ```

2. **Update configuration:**
   ```yaml
   samples: config/custom_samples.tsv
   ```

## Validation

Before running the workflow, validate your configuration:

```bash
# Check sample file format
snakemake -n --configfile config/config.yaml

# Validate specific samples exist in SRA
snakemake --use-conda -n -R download_sra
```

## Advanced Configuration

For advanced users, additional parameters can be added to `config.yaml`:

```yaml
samples: config/samples_all.tsv

# Example additional parameters
assembly:
  megahit_memory: 0.9  # Memory fraction for MEGAHIT
  spades_memory: 500   # Memory limit in GB for SPAdes

diamond:
  sensitivity: "ultra-sensitive"  # Diamond sensitivity
  evalue: 1e-5                   # E-value threshold

megan:
  min_score: 50        # Minimum bit score
  max_expected: 0.01   # Maximum expected value
```

## Troubleshooting

### Common Configuration Issues

1. **Invalid sample format:**
   - Ensure header is exactly `samples`
   - Check for extra spaces or tabs
   - Verify SRA accession format

2. **File path issues:**
   - Use relative paths from project root
   - Ensure sample files exist before running

3. **Memory configuration:**
   - Adjust memory settings for your system
   - Monitor resource usage during runs
