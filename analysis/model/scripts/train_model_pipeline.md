# 1. Model selection - to find the best base model via grid search:
python scripts/train_model_pipeline.py input_data.csv refs.xlsx contigs.fasta --stage select --outdir model_selection

# 2. Cross-validation evaluation - 30 iterations of 5-fold CV with both methods:
python train_model_pipeline.py input_data.csv refs.xlsx contigs.fasta --stage evaluate --iterations 30 --sample_depth 30 --outdir evaluation_results

# 3. Train final model - using RF + LR histogram with bins=10:
python train_model_pipeline.py input_data.csv refs.xlsx contigs.fasta --stage final --outdir final_model