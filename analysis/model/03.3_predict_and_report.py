from utils import *

input = sys.argv[1]
outdir = sys.argv[2]
scaler = load("data/scaler_n25.joblib")
model = load("data/model_n25.joblib")

input_df = pd.read_csv(input)

if input_df.index.name != "orf_name":
    input_df = input_df.set_index("orf_name")

input_df = input_df.drop(columns=["contig_name"])
X_test = scaler.transform(input_df)
y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:, 1]
results_df = pd.DataFrame({"orf_name": input_df.index, "prediction": y_pred, "probability": y_prob})

# add contig name
results_df["contig_name"] = results_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)", expand=False)

# select best evidence
extreme_results_df = results_df.groupby("contig_name").apply(most_extreme_prob).reset_index(drop=True)

# sort by probability
df = extreme_results_df.sort_values(by="probability", ascending=False)

# add info ###SPECIFIC TO THIS DATASET
# df["SRR"] = df["contig_name"].str.split("_").str[-1]
orf_len_mapper = input_df["orf_len"].to_dict()
contig_len_mapper = input_df["contig_length"].to_dict()
df["orf_len"] = df["orf_name"].map(orf_len_mapper)
df["contig_len"] = df["orf_name"].map(contig_len_mapper)

# Define the desired column order
desired_order = [
    # "SRR",
    "contig_name",
    "contig_len",
    "orf_name",
    "orf_len",
    "prediction",
    "probability",
    # "rdrp_tobamo",
    # "rdrp_og",
    # "cp",
    "orf_len",
]

# Reorder the columns
df = df[desired_order]

df.to_csv(f"results/{outdir}/ranked_contig_predictions.csv", index=False)
print("Done!")
