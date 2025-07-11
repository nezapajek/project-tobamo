from utils import *

input_df_path = sys.argv[1]
refs_path = sys.argv[2]
contig_fasta_path = sys.argv[3]
random_seed = 42

# Check if files exist and are in the correct format
check_file_exists(input_df_path, "CSV")
check_file_exists(refs_path, "Excel")
check_file_exists(contig_fasta_path, "FASTA")

# load data
input_df = pd.read_csv(input_df_path, index_col=0)
refs = pd.read_excel(refs_path)
contigs = SeqIO.to_dict(SeqIO.parse(contig_fasta_path, "fasta"))

# Define the models and their parameter grids
models = {
    "LogisticRegression": {
        "model": LogisticRegression(max_iter=1000),
        "params": {"C": [0.1, 1, 10], "solver": ["liblinear", "lbfgs", "saga"]},
    },
    "RandomForest": {
        "model": RandomForestClassifier(),
        "params": {"n_estimators": [10, 25, 50, 75, 100, 125, 150, 175, 200], "max_depth": [5, 10, 20, 30, 40, None]},
    },
    "SVM": {"model": SVC(), "params": {"C": [10, 50, 100], "kernel": ["linear", "rbf", "poly"]}},
    "DecisionTree": {"model": DecisionTreeClassifier(), "params": {"max_depth": [None, 5, 10]}},
    "NaiveBayes": {"model": GaussianNB(), "params": {}},
    "KNN": {"model": KNeighborsClassifier(), "params": {"n_neighbors": [3, 5, 7, 15]}},
}


# Create 5 FOLDS of TRAIN(80%) and TEST(20%)
folds = stratified_kfold_split(refs, n_splits=5, random_state=42)

dfs_test_morf = []
fold_final_predictions = []
future_importances_list = []
performance_metrics_list = []

for idx, (train_refs, test_refs) in enumerate(folds):
    # Prepare TRAIN and TEST data
    train_all, test_all = prepare_train_test(input_df, train_refs, test_refs)

    # subsample contigs
    selected_training_contigs = subsample_contigs(
        train_refs, contigs, num=30, output_dir=f"results/training/sampling/{idx}_train_report", random_seed=random_seed
    )
    selected_test_contigs = subsample_contigs(
        test_refs, contigs, num=30, output_dir=f"results/training/sampling/{idx}_test_report", random_seed=random_seed
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

    # Separate features and target
    X_train = train.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
    y_train = (train["orf_type"] == "tobamo").astype(int)
    X_test = test.drop(columns=["orf_type", "strand", "virus_name", "accession", "contig_name"])
    y_test = (test["orf_type"] == "tobamo").astype(int)

    # Standardize training and test data
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Initialize variables to store the best model for the current species
    best_model = None
    best_score = -float("inf")
    best_params = None
    best_model_name = None

    # Loop through each model and its parameters for GridSearchCV
    for model_name, model_info in models.items():
        print(f"Running GridSearchCV for {model_name} on training data...")

        # Get the model and parameters
        model = model_info["model"]
        params = model_info["params"]

        # Initialize GridSearchCV
        grid_search = GridSearchCV(model, param_grid=params, n_jobs=32, cv=5, scoring="accuracy")

        # Fit GridSearchCV to the training data
        grid_search.fit(X_train, y_train)

        # Check if this model performs better than the current best
        if grid_search.best_score_ > best_score:
            best_model = grid_search.best_estimator_
            best_score = grid_search.best_score_
            best_params = grid_search.best_params_
            best_model_name = model_name

    # Evaluate the best model on the test data
    y_pred = best_model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    auc_roc = roc_auc_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    specificity = recall_score(y_test, y_pred, pos_label=0)
    precision = precision_score(y_test, y_pred)
    negative_predictive_value = precision_score(y_test, y_pred, pos_label=0)

    performance_metrics = {
        "fold": idx,
        "model": best_model_name,
        "accuracy": accuracy,
        "auc_roc": auc_roc,
        "f1_score": f1,
        "recall": recall,
        "specificity": specificity,
        "precision": precision,
        "negative_predictive_value": negative_predictive_value,
        "best_params": best_params,
    }

    performance_metrics_list.append(performance_metrics)

# Save results
performance_metrics = pd.DataFrame(performance_metrics_list)
os.makedirs("results/training/tuning", exist_ok=True)
performance_metrics.to_csv("results/training/tuning/performance_metrics.csv")
