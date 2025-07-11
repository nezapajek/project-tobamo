#!/usr/bin/env python3
# filepath: /home/tobamo/analize/project-tobamo/analysis/model/scripts/train_model_pipeline.py

from utils import *
import argparse
import os
import sys
import pandas as pd
from joblib import dump
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, f1_score, roc_auc_score, recall_score, precision_score
)

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Train and evaluate tobamovirus classification models")
    parser.add_argument("input_df", help="Path to input DataFrame CSV")
    parser.add_argument("refs", help="Path to references Excel file")
    parser.add_argument("contigs", help="Path to contigs FASTA file")
    parser.add_argument("--outdir", default="default", help="Output directory name")
    parser.add_argument("--stage", choices=["select", "evaluate", "hist_test", "final"], 
                        default="final", help="Pipeline stage to run")
    parser.add_argument("--train_depth", type=int, default=30, 
                        help="Number of contigs to sample per species for training")
    parser.add_argument("--test_depth", type=int, default=30, 
                        help="Number of contigs to sample per species for testing")
    parser.add_argument("--iterations", type=int, default=1, 
                        help="Number of iterations for histogram testing")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    return parser.parse_args()

def check_inputs(input_df_path, refs_path, contig_fasta_path):
    """Check if input files exist and are in correct format"""
    check_file_exists(input_df_path, "CSV")
    check_file_exists(refs_path, "Excel")
    check_file_exists(contig_fasta_path, "FASTA")
    
    return pd.read_csv(input_df_path, index_col=0), pd.read_excel(refs_path), \
           SeqIO.to_dict(SeqIO.parse(contig_fasta_path, "fasta"))

def model_selection(input_df, refs, contigs, outdir="training/tuning", random_seed=42):
    """Perform model selection and hyperparameter tuning"""
    print("Starting model selection and hyperparameter tuning...")
    
    # Define the models and their parameter grids
    models = {
        "LogisticRegression": {
            "model": LogisticRegression(max_iter=1000),
            "params": {"C": [0.1, 1, 10], "solver": ["liblinear", "lbfgs", "saga"]},
        },
        "RandomForest": {
            "model": RandomForestClassifier(),
            "params": {"n_estimators": [10, 50, 100, 150, 200], "max_depth": [5, 10, 20, 40, None]},
        },
        "SVM": {"model": SVC(), "params": {"C": [10, 50, 100], "kernel": ["linear", "rbf", "poly"]}},
        "DecisionTree": {"model": DecisionTreeClassifier(), "params": {"max_depth": [None, 5, 10]}},
        "NaiveBayes": {"model": GaussianNB(), "params": {}},
        "KNN": {"model": KNeighborsClassifier(), "params": {"n_neighbors": [3, 5, 7, 15]}},
    }
    
    # Create 5-fold CV split
    folds = stratified_kfold_split(refs, n_splits=5, random_state=random_seed)
    
    performance_metrics_list = []
    
    for idx, (train_refs, test_refs) in enumerate(folds):
        print(f"Processing fold {idx+1}/5")
        
        # Prepare TRAIN and TEST data
        train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)
        
        # Subsample contigs
        selected_training_contigs = subsample_contigs(
            train_refs, contigs, num=30, 
            output_dir=f"results/{outdir}/{idx}_train_report", 
            random_seed=random_seed
        )
        selected_test_contigs = subsample_contigs(
            test_refs, contigs, num=30, 
            output_dir=f"results/{outdir}/{idx}_test_report", 
            random_seed=random_seed
        )
        
        # Filter datasets to selected contigs
        train = train_all[train_all["contig_name"].isin(selected_training_contigs)]
        test = test_all[test_all["contig_name"].isin(selected_test_contigs)]
        
        # Use only forward ORFs for training
        train = train[train["strand"] == "FORWARD"]
        
        # Prepare features and target
        X_train = train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
        y_train = (train["orf_type"] == "tobamo").astype(int)
        X_test = test.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
        y_test = (test["orf_type"] == "tobamo").astype(int)
        
        # Standardize data
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
        
        # Find best model
        best_model, best_score, best_params, best_model_name = None, -float("inf"), None, None
        
        for model_name, model_info in models.items():
            print(f"  Testing {model_name}...")
            model, params = model_info["model"], model_info["params"]
            
            grid_search = GridSearchCV(model, param_grid=params, n_jobs=32, cv=5, scoring="accuracy")
            grid_search.fit(X_train, y_train)
            
            if grid_search.best_score_ > best_score:
                best_model = grid_search.best_estimator_
                best_score = grid_search.best_score_
                best_params = grid_search.best_params_
                best_model_name = model_name
        
        # Evaluate best model
        y_pred = best_model.predict(X_test)
        
        metrics = {
            "fold": idx,
            "model": best_model_name,
            "accuracy": accuracy_score(y_test, y_pred),
            "auc_roc": roc_auc_score(y_test, y_pred),
            "f1_score": f1_score(y_test, y_pred),
            "recall": recall_score(y_test, y_pred),
            "specificity": recall_score(y_test, y_pred, pos_label=0),
            "precision": precision_score(y_test, y_pred),
            "negative_predictive_value": precision_score(y_test, y_pred, pos_label=0),
            "best_params": best_params,
        }
        
        performance_metrics_list.append(metrics)
    
    # Save results
    os.makedirs(f"results/{outdir}", exist_ok=True)
    pd.DataFrame(performance_metrics_list).to_csv(f"results/{outdir}/performance_metrics.csv")
    print("Model selection completed!")

def evaluate_model(input_df, refs, contigs, outdir="training/cross_validation", random_seed=42):
    """Evaluate model performance with 5-fold cross-validation"""
    print("Starting model evaluation...")
    
    os.makedirs(f"results/{outdir}", exist_ok=True)
    
    # Create 5-fold CV split
    folds = stratified_kfold_split(refs, n_splits=5, random_state=random_seed)
    
    dfs_test_morf = []
    fold_final_predictions = []
    
    for idx, (train_refs, test_refs) in enumerate(folds):
        print(f"Processing fold {idx+1}/5")
        
        # Prepare TRAIN and TEST data
        train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)
        
        # Subsample contigs
        selected_training_contigs = subsample_contigs(
            train_refs, contigs, num=30, 
            output_dir=f"results/{outdir}/{idx}_train_report", 
            random_seed=random_seed
        )
        selected_test_contigs = subsample_contigs(
            test_refs, contigs, num=30, 
            output_dir=f"results/{outdir}/{idx}_test_report", 
            random_seed=random_seed
        )
        
        # Filter datasets to selected contigs
        train = train_all[train_all["contig_name"].isin(selected_training_contigs)]
        test = test_all[test_all["contig_name"].isin(selected_test_contigs)]
        
        # Use only forward ORFs for training
        train = train[train["strand"] == "FORWARD"]
        
        # Train RandomForest model
        morf, sorf, _ = train_rf_on_all_data(train)
        
        # Make predictions on test data
        test_orf_predictions = predict_orfs(test, morf, sorf, refs=True)
        dfs_test_morf.append(test_orf_predictions)
        
        # Add most extreme prediction
        most_extreme = filter_extreme_probability(test_orf_predictions, idx, refs_=True)
        fold_final_predictions.append(most_extreme)
    
    # Save fold results
    final_predictions = pd.concat(fold_final_predictions)
    final_predictions.to_csv(f"results/{outdir}/final_predictions.csv")
    
    test_morf = pd.concat(dfs_test_morf)
    test_morf.to_csv(f"results/{outdir}/test_morf.csv")
    
    # Train model on all data and export
    print("Training final model on all data...")
    subsampled = subsample_contigs(
        refs, contigs, num=30, 
        output_dir=f"results/{outdir}/all_report", 
        random_seed=random_seed
    )
    input_df_subsampled = input_df[input_df["contig_name"].isin(subsampled)]
    all_train = input_df_subsampled[input_df_subsampled["strand"] == "FORWARD"]
    morf_final, sorf_final, _ = train_rf_on_all_data(all_train)
    
    # Save the model
    os.makedirs(f"results/{outdir}/model", exist_ok=True)
    dump(morf_final, f"results/{outdir}/model/rf_model.joblib")
    dump(sorf_final, f"results/{outdir}/model/rf_scaler.joblib")
    
    # Save feature importances
    feature_importances = morf_final.feature_importances_
    feature_names = all_train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"]).columns
    feature_importances_df = pd.DataFrame({"Feature": feature_names, "Importance": feature_importances})
    feature_importances_df.to_csv(f"results/{outdir}/feature_importances.csv")
    
    print("Model evaluation completed!")

def histogram_evaluation(input_df, refs, contigs, outdir="training/hist_test", 
                         train_depth=30, test_depth=30, iterations=30, random_seed=42):
    """Perform extensive evaluation with histogram approach"""
    print(f"Starting histogram evaluation with {iterations} iterations...")
    
    os.makedirs(f"results/{outdir}/sampling_{train_depth}_{test_depth}", exist_ok=True)
    
    # Empty list to store final predictions
    fold_final_predictions = []
    
    # Iterate multiple times
    for i in range(iterations):
        print(f"Starting iteration {i+1}/{iterations}")
        iter_seed = random_seed + i
        
        # Create 5 folds
        folds = stratified_kfold_split(refs, n_splits=5, random_state=iter_seed)
        
        for idx, (train_refs, test_refs) in enumerate(folds):
            # Prepare TRAIN and TEST data
            train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)
            
            # Subsample contigs
            selected_training_contigs = subsample_contigs(
                train_refs, contigs, num=train_depth, 
                output_dir=f"results/{outdir}/sampling_{train_depth}_{test_depth}/{i}_{idx}_train_report", 
                random_seed=iter_seed
            )
            selected_test_contigs = subsample_contigs(
                test_refs, contigs, num=test_depth, 
                output_dir=f"results/{outdir}/sampling_{train_depth}_{test_depth}/{i}_{idx}_test_report", 
                random_seed=iter_seed
            )
            
            # Filter datasets
            train = train_all[train_all["contig_name"].isin(selected_training_contigs)]
            test = test_all[test_all["contig_name"].isin(selected_test_contigs)]
            
            # Use only forward ORFs for training
            train = train[train["strand"] == "FORWARD"]
            
            # Train RF model and make predictions
            morf_predictions = train_rf_and_predict(train)
            
            # Train LR model with histogram approach
            bins = [10]  # Can expand with more bin sizes: [10, 15, 20, 25, 30]
            mc_dict = train_lr_and_predict_hist_test(morf_predictions, bins)
            
            # Train RF on all data
            morf, sorf, _ = train_rf_on_all_data(train)
            
            # Test phase
            test_orf_predictions = predict_orfs(test, morf, sorf, refs=True)
            
            # Predict using LR models
            for mc_name_n, mc in mc_dict.items():
                mc_name = mc_name_n.split("_")[0]
                num = int(mc_name_n.split("_")[1])
                final_predictions = predict_contigs_hist_test(test_orf_predictions, mc, idx, mc_name, num, refs_=True)
                final_predictions["n"] = num
                final_predictions["treshold"] = 0.5
                final_predictions["random_seed"] = iter_seed
                fold_final_predictions.append(final_predictions)
    
    # Save results
    print("Saving results...")
    fold_final_predictions_df = pd.concat(fold_final_predictions)
    fold_final_predictions_df["contig_length"] = (
        fold_final_predictions_df["contig_name"].str.extract(r"len-(\d+)").astype(int)
    )
    fold_final_predictions_df["predicted_class"] = np.where(
        fold_final_predictions_df["prob_1"] >= fold_final_predictions_df["treshold"], 1, 0
    )
    fold_final_predictions_df.to_csv(
        f"results/{outdir}/hist_test_fold_final_predictions_{train_depth}_{test_depth}.csv", 
        index=False
    )
    print("Histogram evaluation completed!")

def train_final_model(input_df, refs, contigs, outdir="final_model", random_seed=42):
    """Train the final model on all data"""
    print("Training final models...")
    
    os.makedirs(f"results/{outdir}/sampling_report", exist_ok=True)
    
    # Keep only training data
    train_refs = refs.loc[refs["training"] == 1]
    train_refs.loc[:, "sampling_prob"] = train_refs["length"] / train_refs.groupby("virus_name")["length"].transform("sum")
    
    # Subsample contigs
    selected_contigs = subsample_contigs(
        train_refs, contigs, num=30, 
        output_dir=f"results/{outdir}/sampling_report", 
        random_seed=random_seed
    )
    
    # Filter to selected contigs
    train = input_df[input_df["contig_name"].isin(selected_contigs)]
    
    # Use only forward ORFs for training
    train_fwd = train[train["strand"] == "FORWARD"]
    
    # Train RandomForest model
    morf, sorf, features = train_rf_on_all_data(train_fwd)
    
    # Save models
    dump(morf, f"results/{outdir}/random_forest_model.joblib")
    dump(sorf, f"results/{outdir}/random_forest_scaler.joblib")
    
    # Save feature names
    pd.Series(features).to_csv(f"results/{outdir}/random_forest_feature_names.csv", index=False, header=False)
    
    # Save feature importances
    feature_importances = morf.feature_importances_
    feature_names = train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"]).columns
    feature_importances_df = pd.DataFrame({"Feature": feature_names, "Importance": feature_importances})
    feature_importances_df.to_csv(f"results/{outdir}/final_model_feature_importance.csv", index=False)
    
    # Train RF on subset and get predictions with LOOCV for LR model training
    morf_predictions = train_rf_and_predict(train)
    
    # Train Logistic Regression model on LOOCV predictions
    mc_dict = train_lr_and_predict(morf_predictions, methods=["histogram"])
    
    # Save LR models
    for mc_name, mc in mc_dict.items():
        dump(mc, f"results/{outdir}/logistic_regression_model_{mc_name}.joblib")
    
    print("Final models trained and saved!")

def main():
    args = parse_args()
    print(f"Starting pipeline at stage: {args.stage}")
    
    # Check inputs and load data
    input_df, refs, contigs = check_inputs(args.input_df, args.refs, args.contigs)
    print("Input data loaded successfully")
    
    # Run the requested pipeline stage
    if args.stage == "select":
        model_selection(input_df, refs, contigs, 
                        outdir=f"{args.outdir}/tuning", random_seed=args.seed)
    
    elif args.stage == "evaluate":
        evaluate_model(input_df, refs, contigs, 
                       outdir=f"{args.outdir}/cross_validation", random_seed=args.seed)
    
    elif args.stage == "hist_test":
        histogram_evaluation(input_df, refs, contigs, 
                            outdir=f"{args.outdir}/hist_test",
                            train_depth=args.train_depth, 
                            test_depth=args.test_depth,
                            iterations=args.iterations, 
                            random_seed=args.seed)
    
    elif args.stage == "final":
        train_final_model(input_df, refs, contigs, 
                         outdir=args.outdir, random_seed=args.seed)
    
    print("Pipeline completed successfully!")

if __name__ == "__main__":
    main()