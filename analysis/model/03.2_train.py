from utils import *

outdir = sys.argv[1]
subsampling_n = int(sys.argv[2])

input_df = pd.read_csv(f"results/{outdir}/test_input_df.csv")

out_model = f"results/{outdir}/model/model_n{subsampling_n}.joblib"
out_scaler = f"results/{outdir}/model/scaler_n{subsampling_n}.joblib"

if os.path.exists(out_model):
    raise FileExistsError(f"{out_model} already exists")
else:
    print(f"Training model with subsampling_n={subsampling_n}")
    if input_df.index.name != "orf_name":
        input_df = input_df.set_index("orf_name")

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

    # Save the scaler and the model
    dump(model, out_model)
    dump(scaler, out_scaler)
    print(f"Model saved to {out_model}")
