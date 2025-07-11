from utils import *

input_df_path = sys.argv[1]
refs_path = sys.argv[2]
contig_fasta_path = sys.argv[3]
random_seed = 42
outdir = sys.argv[4]

# Check if files exist and are in the correct format
check_file_exists(input_df_path, "CSV")
check_file_exists(refs_path, "Excel")
check_file_exists(contig_fasta_path, "FASTA")

os.makedirs(f"results/{outdir}/sampling_report", exist_ok=True)

# load data
input_df = pd.read_csv(input_df_path, index_col=0)
refs = pd.read_excel(refs_path)
contigs = SeqIO.to_dict(SeqIO.parse(contig_fasta_path, "fasta"))

# keep only training data
train_refs = refs.loc[refs["training"] == 1]
train_refs.loc[:, "sampling_prob"] = train_refs["length"] / train_refs.groupby("virus_name")["length"].transform("sum")

# subsample 30 contigs per species
selected_contigs = subsample_contigs(
    train_refs, contigs, num=30, output_dir=f"results/{outdir}/sampling_report", random_seed=random_seed
)

# Filter train to keep only contig_name values that are in selected_training_contigs
train = input_df[input_df["contig_name"].isin(selected_contigs)]

# Use only forward ORFs for training
train_fwd = train[train["strand"] == "FORWARD"]
# Train RandomForest model on all data
morf, sorf, features = train_rf_on_all_data(train_fwd)

# Save models
dump(morf, f"results/{outdir}/random_forest_model.joblib")
dump(sorf, f"results/{outdir}/random_forest_scaler.joblib")

# Save feature names to a file for later use
pd.Series(features).to_csv(f"results/{outdir}/random_forest_feature_names.csv", index=False, header=False)

# Calculate feature importance
feature_importances = morf.feature_importances_
feature_names = train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"]).columns
feature_importances_df = pd.DataFrame({"Feature": feature_names, "Importance": feature_importances})
feature_importances_df.to_csv(f"results/{outdir}/final_model_feature_importance.csv", index=False)

# train rf on subset and get predictions with LOOCV for LR model training
morf_predictions = train_rf_and_predict(train)

# Train Logistic Regression model on LOOCV predictions
mc_dict = train_lr_and_predict(morf_predictions, methods=["histogram"])

# Predict contigs using logistic regression models
for mc_name, mc in mc_dict.items():
    dump(mc, f"results/{outdir}/logistic_regression_model_{mc_name}.joblib")
