# **Tobamovirus Contig Classification with Random Forest**

This repository provides Python scripts for training and using a **Random Forest Classifier** and **Logistic Regression** to predict **Tobamoviral sequences**. The workflow includes simulating training data, training the model, and using it to classify query contigs. The process ensures a realistic approach to handling sequencing data and contig assembly.

---

## **Table of Contents**
1. [**Installation**](#installation)
2. [**Workflow Overview**](#workflow-overview)
3. [**Training the Model**](#training-the-model)
4. [**Using the Model for Classification**](#using-the-model-for-classification)
5. [**Example Usage**](#example-usage)
6. [**Notes**](#notes)
7. [**Contact**](#contact)

---

## **1. Installation**
Before running the scripts, you need to set up the environment and install the required tools. We recommend using **Conda** to manage dependencies.  

Run the following commands to set up the environment:  

```bash
# Create and activate the environment
conda env create -f analysis_environment.yml
conda activate tobamo-model
```
Note: You can rename environment in analysis_environment.yml file.

---

## **2. Workflow Steps**

The workflow is divided into two main parts:
  - **A. Training the Model**
  - **B. Using the Model for Classification**

---

### **A. Training the Model**

This process involves several steps. Each step corresponds to a specific script that performs part of the pipeline.

#### **Step 1: Fitting the curve on Snakemake output data for weighted random sampling**

First we take a look inside the Snakemake pipeline output. This part demands some manual checkup and needs to be tailored for each study. In our case, we removed contigs that had hits on *cellular organisms* and contigs from *SRR6846476*, after consulting domain scientists. We then fitted a curve on contig length distribution of the selected contigs, which we'll later use for random weighted sampling of reference genomes, to generate training data.

Example implementation is available in [`notebooks/01_fit_distribution_curve.ipynb`](notebooks/00_fit_distribution_curve.ipynb).

#### **Step 2: Simulating Sequencing and Assembly for Training Data**

This step fragments reference genomes to generate contigs with realistic lengths, ensuring that the training data resembles actual sequencing data.

**Command**
```bash
python scripts/02_sample_refs.py <path/to/reference.fasta> <out_dir_name> <sampling_num> <subsampling_num> <path/to/lens_freq.json>
```
**Arguments**
| Argument                | Description                                              |
|-------------------------|----------------------------------------------------------|
| `<path/to/reference.fasta>` | Path to the reference FASTA file.                       |
| `<out_dir_name>`        | Directory where the output files will be saved.           |
| `<sampling_num>`        | Total number of samples to generate.                      |
| `<subsampling_num>`     | Number of subsamples per reference sequence.              |
| `<path/to/lens_freq.json>`| Path to the lens_freq.json (Step 1 output)              |

#### **Step 3:  Finding ORFs and Pairwise Alignment**

This step identifies **Open Reading Frames (ORFs)** in contigs and performs **pairwise alignment** with reference proteins. It uses **Orfipy**  and **biopython Bio.Seq.Seq.translate** method to detect ORFs and **MAFFT** to perform pairwise alignments against known reference sequences, such as RdRp ORF1, RdRp ORF2, and Coat Protein from species within the family *Virgaviridae*. 

**Command**
```bash
python 03_getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation>
```
**Arguments**
| Argument                  | Description                                                     |
|---------------------------|-----------------------------------------------------------------|
| `<path/to/contig.fasta>`   | Path to the contig FASTA file to be processed.                  |
| `<out_dir_name>`           | Directory where the output files will be saved.                 |
| `<contig_orientation>`     | Orientation of the contigs (e.g., forward, reverse, or unknown).|

note: if you want to use different reference proteins, change "data/all_proteins.fasta" or change path aa_refs in script.


#### **Step 4: Data Processing and Training Input Generation**

This step processes reference data and pairwise alignment results to create a training input dataset. It performs data filtering, aggregation, and enrichment with additional sequence information to prepare the final training dataset.

**Command**
```bash
python 04_preprocess_training.py <path/to/reference_database.xlsx> <path/to/sampled_contigs.fasta> <path/to/orf.fasta> <path/to/pairwise_aln.csv> <output_dir>
```

**Arguments**
| Argument                           | Description                                                     |
|------------------------------------|-----------------------------------------------------------------|
| `<path/to/reference_database.xlsx>`| Path to the Excel file containing reference protein data.       |
| `<path/to/sampled_contigs.fasta>`  | Path to the sampled contigs FASTA file.                        |
| `<path/to/orf.fasta>`              | Path to the ORF FASTA file from previous step.                 |
| `<path/to/pairwise_aln.csv>`       | Path to the CSV file containing pairwise alignment results.     |
| `<output_dir>`                     | Directory name where output files will be saved in results/.    |

**Processing Steps**
1. **File Validation**: Checks that all input files exist and are in correct formats (Excel, CSV, FASTA)
2. **Data Loading**: Loads reference data from Excel and pairwise alignment results from CSV
3. **Filtering**: Removes parent amino acid references from pairwise data
4. **Mapper Creation**: Creates mappings between amino acid IDs, protein types, and virus names
5. **Data Aggregation**: Aggregates pairwise alignment data by protein type and ID
6. **Data Pivoting**: Restructures data for training format
7. **Information Enrichment**: Adds basic sequence information and training-specific metadata
8. **Output Generation**: Saves processed training input and removed references

**Output Files**
- `results/<output_dir>/training_input.csv` - Final training dataset


#### **Step 5A :  Model selection and hyperparameter tuning

This step evaluates multiple classification models to identify the best-performing one for distinguishing "tobamo" ORFs based on genomic features. It uses stratified k-fold cross-validation, subsampled contigs, and standard classification metrics.

**Command**
```bash
python scripts/04A_model_selection.py <path/to/input_data.csv> <path/to/references.xlsx> <path/to/contigs.fasta>
```

**Arguments**
| Argument                  | Description                                                     |
|---------------------------|-----------------------------------------------------------------|
| `<path/to/input_data.csv>`   | Path to Step 3. output (training_input.csv)               |
| `<path/to/references.xlsx>`           | Path to reference excel file                |
| `<path/to/contigs.fasta>`     | Orientation of the contigs (e.g., forward, reverse, or unknown).|

note: if you're using custom reference file, format in the same manner or change scripts accordingly.
note2: best performing model for our data: RandomForestClassifier

#### **Step 5B : Model  


