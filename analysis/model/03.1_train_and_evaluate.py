from utils_v10 import *


outdir = sys.argv[1]

input_df = pd.read_csv(f"results/{outdir}/test_input_df.csv")
orf_fasta_path = f"results/{outdir}/orfs/combined_orfs.fasta"

nt_species2id = mpu.io.read("data/nt_dict.json")

out_model_stats = f"results/{outdir}/model/model_stats.csv"
out_model_results = f"results/{outdir}/model/model_results.csv"

if os.path.exists(out_model_stats) and os.path.exists(out_model_results):
    print("Model already trained and evaluated. Checking results...")
    results_df = pd.read_csv(out_model_results)

else:
    print("Training and evaluating model...")

    # Placeholders to store results
    os.makedirs(f"results/{outdir}/model", exist_ok=True)
    results = []
    model_stats = []

    # LOO
    for species in nt_species2id.keys():
        print(f"\nProcessing species: {species}")

        # make sure 'orf_name' is the index
        if input_df.index.name != "orf_name":
            input_df = input_df.set_index("orf_name")

        # remove species from training data
        test_df_all = input_df[input_df["species"] == species]

        # Only forward strand for training
        train_data_fwd = input_df[input_df["strand"] == "FORWARD"]

        # Feature extraction for training
        X_train = train_data_fwd.drop(columns=["orf_type", "strand", "species", "nt_id", "contig_name"])
        y_train = (train_data_fwd["orf_type"] == "tobamo").astype(
            int
        )  # Convert 'orf_type' to binary (1 for 'tobamo', 0 otherwise)

        # Standardize training and test data
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)

        # Initialize the model
        model = RandomForestClassifier(n_estimators=125, max_depth=40, n_jobs=-1)

        # Train the model on the entire dataset
        model.fit(X_train, y_train)

        # Evaluate the model on each contig in the test species
        for contig in test_df_all["contig_name"].unique():
            test_df = test_df_all[test_df_all["contig_name"] == contig]

            # Feature extraction for testing
            X_test = test_df.drop(columns=["orf_type", "strand", "species", "nt_id", "contig_name"])
            X_test = scaler.transform(X_test)  # Standardize test data
            y_test = (test_df["orf_type"] == "tobamo").astype(int)  # Convert class to binary

            # Make predictions
            y_pred = model.predict(X_test)
            y_prob = model.predict_proba(X_test)

            # Store predictions
            results.append(
                {
                    "species": species,
                    "contig": contig,
                    "orf_name": test_df.index.tolist(),
                    "ground_truth": y_test.tolist(),
                    "prediction": y_pred.tolist(),
                    "probability": y_prob.tolist(),
                }
            )

            # Compute performance metrics
            model_stats.append(
                {
                    "species": species,
                    "contig": contig,
                    "accuracy": accuracy_score(y_test, y_pred),
                    "precision": precision_score(y_test, y_pred, average="weighted", zero_division=0),
                    "recall": recall_score(y_test, y_pred, average="weighted", zero_division=0),
                    "f1_score": f1_score(y_test, y_pred, average="weighted", zero_division=0),
                }
            )

    model_stats_df = pd.DataFrame(model_stats)
    model_stats_df.to_csv(out_model_stats, index=False)
    results_df = pd.DataFrame(results)
    results_df.to_csv(out_model_results, index=False)

# Process results

res_df_out = f"results/{outdir}/model/model_results_processed.csv"
best_out = f"results/{outdir}/model/best_orfs.csv"

if os.path.exists(res_df_out) and os.path.exists(best_out):
    print("Results already processed.")
else:
    print("Processing results...")
    res_df, confusion_matrix = process_df(results_df, input_df)
    res_df.to_csv(res_df_out, index=False)
    confusion_matrix_df = pd.DataFrame(confusion_matrix, columns=["count"])
    confusion_matrix_df.to_csv(f"results/{outdir}/model/confusion_matrix.csv")

    # Select best evidence
    best = select_best_orf(res_df)
    best.to_csv(best_out, index=False)
    cm = best[["true_positive", "true_negative", "false_positive", "false_negative"]].sum()
    cm = pd.DataFrame(cm, columns=["count"])
    cm.to_csv(f"results/{outdir}/model/best_confusion_matrix.csv")
