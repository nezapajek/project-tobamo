from utils import *

input_df_path = sys.argv[1]
refs_path = sys.argv[2]
contig_fasta_path = sys.argv[3]
train_depth = int(sys.argv[4])  # number of contigs to sample per species
test_depth = int(sys.argv[5])  # number of contigs to sample per species in test

print("Checking input data")
# Check if files exist and are in the correct format
check_file_exists(input_df_path, "CSV")
check_file_exists(refs_path, "Excel")
check_file_exists(contig_fasta_path, "FASTA")

# load data
input_df = pd.read_csv(input_df_path, index_col=0)
refs = pd.read_excel(refs_path)
contigs = SeqIO.to_dict(SeqIO.parse(contig_fasta_path, "fasta"))

os.makedirs(f"results/training/hist_test/sampling_{train_depth}_{test_depth}", exist_ok=True)

# empty list to store final predictions
fold_final_predictions = []

# iterate 30 times
for i in range(30):
    print(f"Starting iteration {i}")
    random_seed = 42 + i

    # Create 5 FOLDS of TRAIN(80%) and TEST(20%)
    folds = stratified_kfold_split(refs, n_splits=5, random_state=random_seed)

    print("Starting 5-fold cross-validation")
    for idx, (train_refs, test_refs) in enumerate(folds):
        # Prepare TRAIN and TEST data
        train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)

        # subsample contigs
        selected_training_contigs = subsample_contigs(
            train_refs,
            contigs,
            num=train_depth,
            output_dir=f"results/training/hist_test/sampling_{train_depth}_{test_depth}/{i}_{idx}_train_report",
            random_seed=random_seed,
        )
        selected_test_contigs = subsample_contigs(
            test_refs,
            contigs,
            num=test_depth,
            output_dir=f"results/training/hist_test/sampling_{train_depth}_{test_depth}/{i}_{idx}_test_report",
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

        # Train RandomForest model on subtrain data and make predictions for subtest data
        morf_predictions = train_rf_and_predict(train)

        # Train Logistic Regression model on LOOCV predictions - hist test
        bins = [10]  # , 15, 20, 25, 30]
        mc_dict = train_lr_and_predict_hist_test(morf_predictions, bins)

        print("Training RandomForest model")
        # Train RandomForest model on all data
        morf, sorf, _ = train_rf_on_all_data(train)

        ##################################### TEST phase

        # Make orf predictions on test data
        test_orf_predictions = predict_orfs(test, morf, sorf, refs=True)

        print("Predicting contigs using logistic regression models")
        # Predict contigs using logistic regression models
        for mc_name_n, mc in mc_dict.items():
            mc_name = mc_name_n.split("_")[0]
            num = int(mc_name_n.split("_")[1])
            final_predictions = predict_contigs_hist_test(test_orf_predictions, mc, idx, mc_name, num, refs_=True)
            final_predictions["n"] = num
            final_predictions["treshold"] = 0.5
            final_predictions["random_seed"] = random_seed
            fold_final_predictions.append(final_predictions)

# Save results
print("Saving results")
fold_final_predictions_df = pd.concat(fold_final_predictions)
fold_final_predictions_df["contig_length"] = (
    fold_final_predictions_df["contig_name"].str.extract(r"len-(\d+)").astype(int)
)
fold_final_predictions_df["predicted_class"] = np.where(
    fold_final_predictions_df["prob_1"] >= fold_final_predictions_df["treshold"], 1, 0
)
fold_final_predictions_df.to_csv(
    f"results/training/hist_test/hist_test_fold_final_predictions_{train_depth}_{test_depth}.csv", index=False
)
print("Done")
