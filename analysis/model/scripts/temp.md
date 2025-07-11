# For model selection and hyperparameter tuning:
python train_model_pipeline.py input_data.csv refs.xlsx contigs.fasta --stage select --outdir tuning_results

# For model evaluation with cross-validation:
python train_model_pipeline.py input_data.csv refs.xlsx contigs.fasta --stage evaluate --outdir cv_results

# For histogram-based evaluation:
python train_model_pipeline.py input_data.csv refs.xlsx contigs.fasta --stage hist_test --train_depth 30 --test_depth 30 --iterations 5 --outdir hist_eval

# For training the final model:
python train_model_pipeline.py input_data.csv refs.xlsx contigs.fasta --stage final --outdir final_model