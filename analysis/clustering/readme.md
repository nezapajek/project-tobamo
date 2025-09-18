# Clustering

This folder contains the scripts for clustering contigs obtained from the Snakemake workflow for downstream phylogenetic analysis.

## Installation

### Option 1: Manual Installation

Install the required software and python packages:
- BLAST+ (v.2.15.0 or later)
- Python 3.11

Python packages:
- pandas 
- networkx  
- matplotlib
- biopython 

### Option 2: Installation with conda (recommended)

```bash
conda env create -f clustering_env.yaml
conda activate clustering_env
```


## Usage / Workflow

### 2.1 Run all-vs-all BLASTN

Define input/output files in the script. 

**Inputs:**
- fasta file with contigs to cluster (FASTA_FILE), .fasta

**Outputs:**
- blast database created from contigs (DB_NAME), multiple blastdb files
- file with complete BLASTN results (OUTPUT_FILE), .tsv

Run the blastn.sh script with

```bash
./blast.sh
```

### 2.2 Remove BLASTN self-hits and duplicated hits 

Define input/output files in the script. 

**Inputs:**
- file with complete BLASTN results (INPUT_FILE), .tsv

**Outputs:**
- file with cleaned BLASTN results (OUTPUT_FILE), .tsv

Run the filter_self_hits.py script with:

```bash
python filter_self_hits.py
```

### 2.3 Filter BLASTN results for specified parameters

Define input/output files in the script.

**Inputs:**
- file with cleaned BLASTN results (INPUT_FILE), .tsv

**Outputs:**
- file with edges for the network graph (OUTPUT_EDGES), .txt
- file with filtered BLASTN results (OUTPUT_TABLE), .txt

Define parameters for filtering:
- MIN_IDENTITY = minimal % identity to retain a BLASTN hit
- MIN_LENGTH = minimal aligned fraction of the contigs needed to retain a BLASTN hit

Both criteria must be met for the hits to be retained in the filtered BLASTN results. 

Run the filter_blast_results.py script with:

```bash
python filter_blast_results.py
```

### 2.4 Create a network graph and extract clusters from the graph

Define input/output files in the script.

**Inputs:**
- fasta file with contigs to cluster (FASTA_FILE), .fasta
- file with edges for the network graph (OUTPUT_EDGES), .txt

**Outputs:**
- file with a list of contigs and corresponding cluster ids (CLUSTER_FILE), .tsv

Run the clustering.py script with:

```bash
python clustering.py
```


## Contact

For questions and support, contact lana.vogrinec@nib.si