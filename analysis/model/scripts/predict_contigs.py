from utils import *

model_input_df_path = sys.argv[1]
morf_path = sys.argv[2]
sorf_path = sys.argv[3]
mc_path = sys.argv[4]
feature_names_path = sys.argv[5]
outdir = sys.argv[6]

morf = load(morf_path)
sorf = load(sorf_path)
mc = load(mc_path)
input_df = pd.read_csv(model_input_df_path)
os.makedirs(f"results/{outdir}", exist_ok=True)

cols = pd.read_csv(feature_names_path, header=None)[0].tolist()

# add orf_name in index, and remove unnecessary columns
if input_df.index.name != "orf_name":
    input_df = input_df.set_index("orf_name")
input_df = input_df[cols]

# scale data
X_test = sorf.transform(input_df)
# make ORF predictions
y_pred = morf.predict(X_test)
y_prob = morf.predict_proba(X_test)[:, 1]
results_df = pd.DataFrame({"orf_name": input_df.index, "prediction": y_pred, "prob_1": y_prob})
# add contig name
results_df["contig_name"] = results_df["orf_name"].str.extract(r"(.*)(?:_ORF\.\d+|_aa_frame\d+)", expand=False)
# save ORF predictions
results_df.to_csv(f"results/{outdir}/contig_predictions_all.csv", index=False)

# prepare predictions for contig prediction
bin_df = prepare_bin_df(results_df, "histogram", num_bins=10)
X_test = bin_df
# predict cotnigs
final_predictions = mc.predict(X_test)
final_predictions_prob = mc.predict_proba(X_test)[:, 1]
final_predictions_df = pd.DataFrame(
    {"contig_name": bin_df.index, "predicted_class": final_predictions, "prob_1": final_predictions_prob}
)
# save contig predictions
final_predictions_df.to_csv(f"results/{outdir}/contig_predictions.csv", index=False)
