# **Tobamovirus Contig Classification with Random Forest**

This repository provides Python scripts for training and using a **Random Forest Classifier** to predict **Tobamoviral sequences** in contigs. The workflow includes simulating training data, training the model, and using it to classify query contigs. The process ensures a realistic approach to handling sequencing data and contig assembly.

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

Example implementation is available in [`00_fit_distribution_curve.ipynb`](./00_fit_distribution_curve.ipynb).

#### **Step 2: Simulating Sequencing and Assembly for Training Data**

This step fragments reference genomes to generate contigs with realistic lengths, ensuring that the training data resembles actual sequencing data.

**Command**
```bash
python 00_sample_refs.py <path/to/reference.fasta> <out_dir_name> <sampling_num> <subsampling_num> 
```
**Arguments**
| Argument                | Description                                              |
|-------------------------|----------------------------------------------------------|
| `<path/to/reference.fasta>` | Path to the reference FASTA file.                       |
| `<out_dir_name>`        | Directory where the output files will be saved.           |
| `<sampling_num>`        | Total number of samples to generate.                      |
| `<subsampling_num>`     | Number of subsamples per reference sequence.              |

#### **Step 3:  Finding ORFs and Pairwise Alignment**

This step identifies **Open Reading Frames (ORFs)** in contigs and performs **pairwise alignment** with reference proteins. It uses **Orfipy**  and **biopython Bio.Seq.Seq.translate** method to detect ORFs and **MAFFT** to perform pairwise alignments against known reference sequences, such as RdRp ORF1, RdRp ORF2, and Coat Protein from species within the family *Virgaviridae*.

**Command**
```bash
python 01_getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation>
```
**Arguments**
| Argument                  | Description                                                     |
|---------------------------|-----------------------------------------------------------------|
| `<path/to/contig.fasta>`   | Path to the contig FASTA file to be processed.                  |
| `<out_dir_name>`           | Directory where the output files will be saved.                 |
| `<contig_orientation>`     | Orientation of the contigs (e.g., forward, reverse, or unknown).|