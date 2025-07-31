# *Tobamovirus Contig Classification with Random Forest and Logistic Regression*

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

# **1. Installation**
Before running the scripts, you need to set up the environment and install the required tools. We recommend using **Conda** to manage dependencies.  

Run the following commands to set up the environment:  

```bash
# Create and activate the environment
conda env create -f analysis_environment.yml
conda activate tobamo-model
```
Note: You can rename environment in analysis_environment.yml file.

---

# **2. Workflow Steps**

The workflow is divided into two main parts:
  - **A. Training the Model**
  - **B. Using the Model for Classification**

---

## **A. Training the Model**

This process involves several steps. Each step corresponds to a specific script that performs part of the pipeline.

### **Step 1: Fitting the curve on Snakemake output data for weighted random sampling**

First we take a look inside the Snakemake pipeline output. This part demands some manual checkup and needs to be tailored for each study. In our case, we removed contigs that had hits on *cellular organisms* and contigs from *SRR6846476*, after consulting domain scientists. We then fitted a curve on contig length distribution of the selected contigs, which we'll later use for random weighted sampling of reference genomes, to generate training data.

Example implementation is available in [`notebooks/01_fit_distribution_curve.ipynb`](notebooks/00_fit_distribution_curve.ipynb).

### **Step 2: Simulating Sequencing and Assembly for Training Data**

This step fragments reference genomes to generate contigs with realistic lengths, ensuring that the training data resembles actual sequencing data.

**Command**
```bash
python scripts/sample_refs.py <path/to/reference.fasta> <out_dir_name> <sampling_num> <subsampling_num> <path/to/lens_freq.json>
```
**Arguments**
| Argument                | Description                                              |
|-------------------------|----------------------------------------------------------|
| `<path/to/reference.fasta>` | Path to the reference FASTA file.                       |
| `<out_dir_name>`        | Directory where the output files will be saved.           |
| `<sampling_num>`        | Total number of samples to generate.                      |
| `<subsampling_num>`     | Number of subsamples per reference sequence.              |
| `<path/to/lens_freq.json>`| Path to the lens_freq.json (Step 1 output)              |

### **Step 3:  Finding ORFs and Pairwise Alignment**

This step identifies **Open Reading Frames (ORFs)** in contigs and performs **pairwise alignment** with reference proteins. It uses **Orfipy**  and **biopython Bio.Seq.Seq.translate** method to detect ORFs and **MAFFT** to perform pairwise alignments against known reference sequences, such as RdRp ORF1, RdRp ORF2, and Coat Protein from species within the family *Virgaviridae*. 

**Command**
```bash
python getorfs_pairwise_aln.py <path/to/contig.fasta> <out_dir_name> <contig_orientation>
```
**Arguments**
| Argument                  | Description                                                     |
|---------------------------|-----------------------------------------------------------------|
| `<path/to/contig.fasta>`   | Path to the contig FASTA file to be processed.                  |
| `<out_dir_name>`           | Directory where the output files will be saved.                 |
| `<contig_orientation>`     | Orientation of the contigs (e.g., forward or unknown).|

note: if you want to use different reference proteins, change "data/all_proteins.fasta" or change path aa_refs in script.


### **Step 4: Data Processing and Training Input Generation**

This step processes reference data and pairwise alignment results to create a training input dataset. It performs data filtering, aggregation, and enrichment with additional sequence information to prepare the final training dataset.

**Command**
```bash
python scripts/preprocess.py <path/to/reference_database.xlsx> <path/to/orf.fasta> <path/to/pairwise_aln.csv> <output_dir> --train --contigs <path/to/sampled_contigs.fasta>
```

**Arguments**
| Argument                           | Description                                                     |
|------------------------------------|-----------------------------------------------------------------|
| `<path/to/reference_database.xlsx>`| Path to the Excel file containing reference protein data.       |
| `<path/to/orf.fasta>`              | Path to the ORF FASTA file from previous step.                 |
| `<path/to/pairwise_aln.csv>`       | Path to the CSV file containing pairwise alignment results.     |
| `<output_dir>`                     | Directory name where output files will be saved in results/.    |
| `--train`                          | Flag to specify processing for training data.                  |
| `--contigs`                        | Path to the contigs FASTA file (required for training).        |

**Processing Steps**
1. **File Validation**: Checks that all input files exist and are in correct formats (Excel, CSV, FASTA)
2. **Data Loading**: Loads reference data from Excel and pairwise alignment results from CSV
3. **Filtering**: Removes parent amino acid references from pairwise data (training mode only)
4. **Mapper Creation**: Creates mappings between amino acid IDs, protein types, and virus names
5. **Data Aggregation**: Aggregates pairwise alignment data by protein type and ID
6. **Data Pivoting**: Restructures data for model input format
7. **Information Enrichment**: Adds basic sequence information and metadata
8. **Output Generation**: Saves processed input data

**Output Files**
- `results/<output_dir>/training_input.csv` - Final training dataset (when using --train)


### **Step 5: Model Training and Evaluation**

Our machine learning pipeline employs a three-stage approach to ensure optimal classification performance.

#### Model Selection

First, we conduct comprehensive algorithm selection through systematic evaluation of multiple classification models:

- **Algorithms Tested**: Logistic Regression, Random Forest, SVM, Decision Tree, Naive Bayes, and KNN
- **Hyperparameter Optimization**: Grid search across all relevant parameters for each algorithm
- **Evaluation Method**: 5-fold cross-validation with stratified sampling
- **Results**: Random Forest demonstrated superior performance across evaluation metrics

#### Cross-Validation Evaluation

We then compare two different strategies for contig-level prediction:

- **ORF Prediction**: All models use Random Forest for Open Reading Frame (ORF) classification (best performing model from previous step)
- **Contig Prediction Methods**:
  - **Extreme Method**: Uses the most confident ORF prediction score
  - **Histogram Method**: Bins ORF predictions and uses Logistic Regression
- **Validation Process**: Leave-One-Out Cross-Validation (LOOCV) repeated 30 times
- **Performance Assessment**: Comprehensive metrics including accuracy, F1 score, precision, and recall
- **Winner**: Histogram-based approach (bins=10) achieved superior performance

#### Final Model Training

The production model combines the best components from our evaluation:

- **ORF Classifier (Morf)**: Random Forest trained on all available training data
- **Contig Classifier (Mc)**: Logistic Regression using binned ORF prediction probabilities
- **Feature Importance**: Analysis reveals most informative sequence characteristics
- **Serialization**: Models saved as joblib files for deployment in production pipeline

This multi-stage approach ensures robust performance across diverse viral sequence data while maintaining interpretability

**Command**
```bash
python scripts/train_model_pipeline.py <path/to/training_input.csv> <path/to/references.xlsx> <path/to/contigs.fasta> --stage <stage> [options]
```

**Arguments**
| Argument                           | Description                                                     |
|------------------------------------|-----------------------------------------------------------------|
| `<path/to/training_input.csv>`     | Path to the training input CSV from Step 4.                     |
| `<path/to/references.xlsx>`        | Path to the Excel file containing reference protein data.       |
| `<path/to/contigs.fasta>`          | Path to the contigs FASTA file.                                |
| `--stage`                          | Pipeline stage to run (`select`, `evaluate`, or `final`).       |
| `--outdir`                         | Output directory name (default: "default").                     |
| `--iterations`                     | Number of iterations for cross-validation (default: 30).        |
| `--sample_depth`                   | Number of contigs to sample per species (default: 30).          |
| `--seed`                           | Random seed for reproducibility (default: 42).                  |

**Pipeline Stages**

1. **Model Selection** (`--stage select`)
   - Performs grid search to find the best-performing model type
   - Tests multiple models (RandomForest, LogisticRegression, SVM, etc.)
   - Evaluates models using 5-fold cross-validation
   - Saves performance metrics for each model
   
   ```bash
   python scripts/train_model_pipeline.py training_input.csv references.xlsx contigs.fasta --stage select --outdir model_selection
   ```

2. **Cross-Validation Evaluation** (`--stage evaluate`)
   - Performs extensive evaluation using multiple iterations of 5-fold cross-validation
   - Compares two prediction methods:
     - Most extreme probability approach
     - Histogram-based approach (using logistic regression with bins=10)
   - Generates comprehensive performance metrics
   
   ```bash
   python scripts/train_model_pipeline.py training_input.csv references.xlsx contigs.fasta --stage evaluate --iterations 30 --sample_depth 30 --outdir evaluation_results
    ```

3. **Final Model Training** (`--stage final`)
   - Trains the final production model on all training data
   - Uses the best-performing approach (Random Forest + Logistic Regression histogram)
   - Saves all necessary model files for deployment
   
   ```bash
   python scripts/train_model_pipeline.py training_input.csv references.xlsx contigs.fasta --stage final --outdir final_model
   ```

**Output Files**
- **Model Selection**:
  - `results/<outdir>/performance_metrics.csv` - Performance metrics for all tested models
  - `results/<outdir>/best_model.txt` - Information about the best-performing model

- **Cross-Validation Evaluation**:
  - `results/<outdir>/extreme_predictions_results.csv` - Results using most extreme probability method
  - `results/<outdir>/histogram_predictions_results.csv` - Results using histogram-based approach
  - `results/<outdir>/method_comparison.csv` - Performance comparison between methods
  - `results/<outdir>/best_method.txt` - Information about the best-performing method

- **Final Model**:
  - `results/<outdir>/rf_model.joblib` - Trained Random Forest model
  - `results/<outdir>/rf_scaler.joblib` - StandardScaler for feature normalization
  - `results/<outdir>/rf_feature_names.csv` - Feature names used by the model
  - `results/<outdir>/lr_histogram_10_model.joblib` - Trained Logistic Regression histogram model
  - `results/<outdir>/feature_importances.csv` - All feature importances ranked
  - `results/<outdir>/top_20_features.csv` - Top 20

  ## **B. Using the Model for Classification**

Once the model has been trained and finalized, you can use it to classify new query contigs.

Step 1: Preprocessing Query Contigs

You must process query contigs through ORF detection, pairwise alignment, and training input preparation. These steps replicate parts of the training pipeline.

1.1: ORF Detection and Pairwise Alignment

```bash
# 1. ORF detection and alignment
python scripts/getorfs_pairwise_aln.py ../data/contigs/contigs_all_deduplicated.fasta <out_dir_name> <contig_orientation>
```
**Output File**
- `results/<output_dir>/pairwise_aln.csv` - pairwise alignment metrics table of test data

1.2. Preprocessing of test data (same as training, but with different name parsing and without ground truth)

```bash
# 2. Generate processed input for prediction
python scripts/preprocess.py <path/to/reference_database.xlsx> <path/to/orf.fasta> <path/to/pairwise_aln.csv> <output_dir> --test
```

**Arguments for Test Processing**
| Argument                           | Description                                                     |
|------------------------------------|-----------------------------------------------------------------|
| `<path/to/reference_database.xlsx>`| Path to the Excel file containing reference protein data.       |
| `<path/to/orf.fasta>`              | Path to the ORF FASTA file from previous step.                 |
| `<path/to/pairwise_aln.csv>`       | Path to the CSV file containing pairwise alignment results.     |
| `<output_dir>`                     | Directory name where output files will be saved in results/.    |
| `--test`                           | Flag to specify processing for test/query data.                |

**Output File**
- `results/<output_dir>/testing_input.csv` - Processed data for prediction

Step 2: Predicting Tobamovirus Contigs

Once you have the processed input features and the trained models, run the prediction script:

```bash
python scripts/predict_contigs.py results/<output_dir>/testing_input.csv results/final_model --outdir predictions --bin-num 10
```

**Arguments**
| Argument                     | Description                                                                 |
|------------------------------|-----------------------------------------------------------------------------|
| `<testing_input_df.csv>`     | Path to the processed input CSV file from preprocessing step.              |
| `<model_dir>`                | Directory containing all model files (RF model, scaler, LR model, etc.).    |
| `--outdir`                   | Name of output directory for prediction results (default: "predictions").   |
| `--bin-num`                  | Number of bins for histogram approach (default: 10).                        |

The script expects these files in the model directory with standard names:
- `rf_model.joblib` - Trained Random Forest model
- `rf_scaler.joblib` - StandardScaler for feature normalization
- `rf_feature_names.csv` - Feature names used by the model
- `lr_histogram_<bin-num>_model.joblib` - Trained Logistic Regression histogram model

**Output Files**
After running the script, the following files will be saved in `results/<outdir>/`:

| Filename                   | Description                                                  |
|----------------------------|--------------------------------------------------------------|
| `orf_predictions.csv`      | ORF-level predictions with probability scores.               |
| `contig_predictions.csv`   | Final contig-level predictions with class and probability.   |

---