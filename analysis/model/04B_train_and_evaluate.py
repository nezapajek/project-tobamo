from utils import *

input_df_path = sys.argv[1]
refs_path = sys.argv[2]
contig_fasta_path = sys.argv[3]
random_seed = 42

# Check if files exist and are in the correct format
check_file_exists(input_df_path, "CSV")
check_file_exists(refs_path, "Excel")
check_file_exists(contig_fasta_path, "FASTA")

os.makedirs("results/training/cross_validation", exist_ok=True)

# load data
input_df = pd.read_csv(input_df_path, index_col=0)
refs = pd.read_excel(refs_path)
contigs = SeqIO.to_dict(SeqIO.parse(contig_fasta_path, "fasta"))

# Create 5 FOLDS of TRAIN(80%) and TEST(20%)
folds = stratified_kfold_split(refs, n_splits=5, random_state=random_seed)

dfs_test_morf = []
fold_final_predictions = []
performance_metrics_list = []

for idx, (train_refs, test_refs) in enumerate(folds):
    # Prepare TRAIN and TEST data
    train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)

    # subsample contigs
    selected_training_contigs = subsample_contigs(
        train_refs,
        contigs,
        num=30,
        output_dir=f"results/training/cross_validation/{idx}_train_report",
        random_seed=random_seed,
    )
    selected_test_contigs = subsample_contigs(
        test_refs,
        contigs,
        num=30,
        output_dir=f"results/training/cross_validation/{idx}_test_report",
        random_seed=random_seed,
    )

    # Filter train and test to keep only contig_name values that are in selected_training_contigs and selected_test_contigs respectively
    train = train_all[train_all["contig_name"].isin(selected_training_contigs)]
    test = test_all[test_all["contig_name"].isin(selected_test_contigs)]

    # check
    assert set(train["contig_name"]).issubset(set(selected_training_contigs))
    assert set(test["contig_name"]).issubset(set(selected_test_contigs))

    # Use only forward ORFs for training
    train = train[train["strand"] == "FORWARD"]

    ################################################ TRAIN PHASE

    # Train RandomForest model on all data
    morf, sorf = train_rf_on_all_data(train)

    ##################################### TEST phase

    # Make orf predictions on test data
    test_orf_predictions = predict_orfs(test, morf, sorf, refs=True)
    dfs_test_morf.append(test_orf_predictions)

    # Add most extreme prediction
    most_extreme = filter_extreme_probability(test_orf_predictions, idx, refs_=True)
    fold_final_predictions.append(most_extreme)


# Save fold results
final_predictions = pd.concat(fold_final_predictions)
final_predictions.to_csv("results/training/final_predictions.csv")

test_morf = pd.concat(dfs_test_morf)
test_morf.to_csv("results/training/test_morf.csv")

performance_metrics = pd.DataFrame(performance_metrics_list)
performance_metrics.to_csv("results/training/performance_metrics.csv")

# Train model on all data and export
subsampled = subsample_contigs(
    refs, contigs, num=30, output_dir="results/training/cross_validation/all_report", random_seed=random_seed
)
input_df_subsampled = input_df[input_df["contig_name"].isin(subsampled)]
all_train = input_df_subsampled[input_df_subsampled["strand"] == "FORWARD"]
morf_final, sorf_final = train_rf_on_all_data(all_train)

# Save the model and scaler
os.makedirs("results/training/model", exist_ok=True)
dump(morf_final, "results/training/model/rf_model.joblib")
dump(sorf_final, "results/training/model/rf_scaler.joblib")

# Calculate and save feature importances
feature_importances = morf_final.feature_importances_
feature_names = all_train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"]).columns
feature_importances_df = pd.DataFrame({"Feature": feature_names, "Importance": feature_importances})

feature_importances_df.to_csv("results/training/feature_importances.csv")
