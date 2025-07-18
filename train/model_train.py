import os  # or use pathlib.Path for more advanced handling
import pandas as pd
from autogluon.tabular import TabularDataset, TabularPredictor
import zarr
from imblearn.over_sampling import RandomOverSampler
import random
import numpy as np
import torch

# Function to set the seed for reproducibility
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True  # Ensure deterministic behavior
        torch.backends.cudnn.benchmark = False

set_seed(42)  # 42, 123, 2025

def balance_data(X_train, y_train, method='oversampling', random_state=42):
    """
    Balances the training data using oversampling or undersampling.

    Parameters:
    - X_train: Feature DataFrame (without target column).
    - y_train: Target Series or DataFrame.
    - method: 'oversampling' or 'undersampling'.
    - random_state: Random seed for reproducibility.

    Returns:
    - X_train_balanced: Balanced feature DataFrame.
    - y_train_balanced: Balanced target Series.
    """
    if method == 'oversampling':
        ros = RandomOverSampler(random_state=random_state)
        X_train_balanced, y_train_balanced = ros.fit_resample(X_train, y_train)
        print("Applied oversampling to the training set.")
    elif method == 'undersampling':
        rus = RandomUnderSampler(random_state=random_state)
        X_train_balanced, y_train_balanced = rus.fit_resample(X_train, y_train)
        print("Applied undersampling to the training set.")
    else:
        raise ValueError("Invalid method. Choose 'oversampling' or 'undersampling'.")
    
    return X_train_balanced, y_train_balanced

# Function to set the output directory
def set_output_directory(path):
    os.makedirs(path, exist_ok=True)
    return path

OUT_DIR = set_output_directory('/data1/deepLN/train/train_sagemaker_scripts/autoB_cluster1_fullfeature')

input_path = '/data1/deepLN/train/data/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v30.modelB_train_1.zarr'

print("Load zarr directories")
root = zarr.open(input_path, mode='r')
column_names = root['metadata'].attrs['columns']
sample_ids = root['metadata'].attrs['sample_ids']
data = root['data'][:]
train_df = pd.DataFrame(data, columns=column_names, index = sample_ids)
train_df.index.name = 'SAMPLE'
train_df.reset_index(inplace=True)

dt = pd.read_table('/data1/deepLN/train/train_sagemaker_scripts/output/new_agg_v27_modelA_train/kor_sfari_mssng.feature_selection_7L.shap.mean.tsv.gz')

filtered_train_df = train_df[['label', 'SAMPLE'] + dt['Feature'].tolist()]

print("filtered_train_df.shape:", filtered_train_df.shape, flush=True)

# Balance the training data
X = filtered_train_df.drop(['SAMPLE', 'label'], axis=1)  # Features only
y = train_df['label']  # Target labels

X_train_balanced, y_train_balanced = balance_data(X, y, method='oversampling')

# Convert balanced data to a DataFrame
X_train_balanced = pd.DataFrame(X_train_balanced, columns=X.columns)
# Concatenate the balanced features and target as a single DataFrame
X_train_balanced = pd.concat([X_train_balanced, pd.Series(y_train_balanced, name='label')], axis=1)

train = TabularDataset(X_train_balanced)

print("Balanced training dataset prepared for AutoGluon.", flush=True)

predictor = TabularPredictor(
    label='label',
    problem_type='binary',
    eval_metric='roc_auc',
    path=os.path.join(OUT_DIR, "model")
).fit(train, 
      presets='best_quality', 
      num_stack_levels=1,
    #   num_bag_folds=0,
      holdout_frac=(3/17), # 0.3
    #   fit_weighted_ensemble = False,
    #   hyperparameters={'NN_TORCH': {}},
      excluded_model_types=['RF', 'XT', 'GBM', 'KNN', 'LR', 'CAT', 'XGB', 'FASTAI'],
      time_limit=3600*2,
      num_gpus=2)

# Save leaderboard
ld_board = predictor.leaderboard(train, silent=True)
ld_board.to_csv(os.path.join(OUT_DIR, 'leaderboard.tsv'), sep='\t', index=False)
print("Saved leader board.", flush=True)

input_path = '/data1/deepLN/train/data/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v30.modelB_test_1.zarr'
print("Load zarr directories")
root = zarr.open(input_path, mode='r')
column_names = root['metadata'].attrs['columns']
sample_ids = root['metadata'].attrs['sample_ids']
data = root['data'][:]
test_df = pd.DataFrame(data, columns=column_names, index = sample_ids)
test_df.index.name = 'SAMPLE'
test_df.reset_index(inplace=True)

filtered_test_df = test_df[['label', 'SAMPLE'] + dt['Feature'].tolist()]

print("filtered_test_df.shape:", filtered_test_df.shape, flush=True)

model_to_use = predictor.model_best
test = TabularDataset(filtered_test_df.drop(['SAMPLE', 'label'], axis=1))
pred_y = predictor.predict(test, model=model_to_use)

# Save predictions
predictions_df = pd.DataFrame({'SAMPLE': test_df['SAMPLE'], 'Prediction': pred_y})
predictions_df.to_csv(os.path.join(OUT_DIR, 'predictions.tsv'), sep='\t', index=False)

print("Saved predictions.", flush=True)

# Use predict_proba to get prediction probabilities
pred_proba = predictor.predict_proba(test, model=model_to_use)

# Save prediction probabilities
# pred_proba is a DataFrame where each column represents a class probability
pred_proba['SAMPLE'] = sample_ids  # Add SAMPLE column for reference
pred_proba.to_csv(os.path.join(OUT_DIR, 'prediction_probabilities.tsv'), sep='\t', index=False)

print("Saved prediction probabilities.", flush=True)

# Save the best model
best_model = predictor.model_best
print(f"Best model: {best_model}")

results = predictor.evaluate(data=test_df.drop(['SAMPLE'], errors='ignore'), model = best_model, decision_threshold=0.5, display=True)
results_df = pd.DataFrame([results])
results_df.to_csv(os.path.join(OUT_DIR, 'evaluation_results.tsv'), sep="\t", index=False)

print("Saved evaluation.", flush=True)
