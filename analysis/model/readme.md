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
conda create --name <env_name> python=3.10
conda activate <env_name>

# Install Python dependencies
pip install -r requirements.txt

# Install bioinformatics tools
conda install bioconda::orfipy 
conda install -c bioconda mafft=7.520
```
Note: Replace <env_name> with the name of your Conda environment.

---

## **2. Workflow Steps**

The workflow is divided into two main parts:
  - **A. Training the Model**
  - **B. Using the Model for Classification**

---

### **A. Training the Model**

This process involves several steps. Each step corresponds to a specific script that performs part of the pipeline.

#### **Step 1: Simulating Sequencing and Assembly for Training Data****

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

