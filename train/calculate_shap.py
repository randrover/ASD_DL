import os
import subprocess
import sys
import tarfile
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, accuracy_score, f1_score
import argparse
import random
import datetime
import torch
import shap
import optuna
import torch.nn as nn
from tqdm import tqdm
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import zarr
import os
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, precision_score, recall_score, precision_recall_curve, auc
import wandb
import time
import multiprocessing
from multiprocessing import Pool, cpu_count
from functools import partial
import csv
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
import torch
import zarr
import pandas as pd
from torch.nn.functional import sigmoid
from train_7L import SimpleNN

from sklearn.cluster import KMeans

def define_kmeans_background(X_train, n_clusters=10):
    """
    Use KMeans to define cluster centers as background data for SHAP.
    """
    kmeans = KMeans(n_clusters=n_clusters, random_state=42).fit(X_train)
    return torch.tensor(kmeans.cluster_centers_, dtype=torch.float32)

input_path = 'data/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v27.modelA_train.zarr'

print("Load zarr directories")
root = zarr.open(input_path, mode='r')
column_names = root['metadata'].attrs['columns']
sample_ids = root['metadata'].attrs['sample_ids']
data = root['data'][:]
train_df = pd.DataFrame(data, columns=column_names, index = sample_ids)
train_df.index.name = 'SAMPLE'
train_df.reset_index(inplace=True)

# high_f = pd.read_table('data/kor_sfari_mssng.feature_importance_10L.shap.mean.tsv.gz')
sorted_high_f = pd.read_table('data/kor_sfari_mssng.feature_selection_7L.shap.mean.tsv.gz')

filtered_columns = sorted_high_f['Feature'].tolist() + ['label']
train_df_filtered = train_df[filtered_columns]

print("train_df_filtered.shape:", train_df_filtered.shape, flush=True)

# Balance the training data
device = 'cuda'
X_train = train_df_filtered.drop(['label'], axis=1)  # Features only
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_train)
X_train_tensor=torch.from_numpy(X_scaled.astype('float32')).float().to(device)

zarr_path = 'data/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v27.modelA_test.zarr'

root = zarr.open(zarr_path, mode='r')
header = [name.replace(',', '').replace(' ', '_') for name in list(root['metadata'].attrs['columns'])]
data = root['data'][:]
sample_ids = root['metadata'].attrs['sample_ids']    
header_index = {col: idx for idx, col in enumerate(header)}
df = pd.DataFrame(data, columns=header)

filtered_columns = sorted_high_f['Feature'].tolist() + ['label']
df_filtered = df[filtered_columns]

X_test = df_filtered.drop(columns=['label'])
y_test = df_filtered['label']

# Scale data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_test)

device= 'cuda'
X_test_tensor = torch.from_numpy(X_scaled.astype('float32')).float().to(device)
y_test_tensor = torch.from_numpy(y_test.values.astype('float32')).float().to(device)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model_path = "model/20250430_170444_9311_best_model.pth" # model path
model = SimpleNN(X_test_tensor.shape[1]).to(device)
model.load_state_dict(torch.load(model_path, map_location=device, weights_only=False)) 
model.eval()

output_dir = 'new_agg_v27_modelA_train'

predict_prob_file = os.path.join(output_dir, f'kor_sfari_mssng_predict_prob.test_sample.oversample.tsv.gz')

if predict_prob_file:
    with torch.no_grad():
        outputs = model(X_test_tensor).squeeze()
        y_pred_prob = sigmoid(outputs).cpu().numpy() 
        y_pred_test = (y_pred_prob >= 0.5).astype(int) 

        # 결과를 DataFrame으로 저장
        results_df = pd.DataFrame({
            'SAMPLE': sample_ids, 
            'TrueLabel': y_test_tensor.cpu().numpy(),
            'PredictedLabel': y_pred_test,
            'PredictedProbability': y_pred_prob
        })

        # Calculate metrics
        true_labels = y_test_tensor.cpu().numpy()

        # AUROC
        auroc = roc_auc_score(true_labels, y_pred_prob)

        # AUPRC
        precision, recall, _ = precision_recall_curve(true_labels, y_pred_prob)
        auprc = auc(recall, precision)

        # Accuracy
        accuracy = accuracy_score(true_labels, y_pred_test)

        # F1 Score
        f1 = f1_score(true_labels, y_pred_test)

        # Print metrics
        print(f"AUROC: {auroc:.4f}")
        print(f"AUPRC: {auprc:.4f}")
        print(f"Accuracy: {accuracy:.4f}")
        print(f"F1 Score: {f1:.4f}")

        
        results_df.to_csv(predict_prob_file, sep='\t', index=False)
        print(f"Predict probability saved to {predict_prob_file}", flush=True)
        
feature_names = [col for col in header if col != 'label']



def calculate_importance(model, X_tensor, y_tensor, importance_method='permutation', n_repeats=10,
                         random_state=42, num_processes=None, permutation_split=1, permutation_block=0,
                         X_train_tensor=None):
    if importance_method == 'shap':
        return calculate_shap_values(model=model, X_tensor=X_tensor, num_processes=num_processes,
                                     X_train_tensor=X_train_tensor)
    elif importance_method == 'permutation':
        return custom_permutation_importance(model, X_tensor, y_tensor, n_repeats=n_repeats, random_state=random_state,
                                             num_processes=num_processes, permutation_split=permutation_split, 
                                             permutation_block=permutation_block)
    else:
        raise ValueError("Invalid importance method. Choose 'shap' or 'permutation'.")

def calculate_shap_for_batch(batch_idx, X_batch, model, device, background):
    # Initialize DeepExplainer for the subset
    explainer = shap.DeepExplainer(model, background)
    
    # Calculate SHAP values for the subset
    shap_values = explainer.shap_values(X_batch, check_additivity=False)
    
    # Reshape the SHAP values before returning (shape: [n_samples, n_features])
    reshaped_shap_values = np.squeeze(shap_values)
    
    # Debugging output to check the number of samples in each batch
    print(f"Batch {batch_idx}: SHAP values shape {reshaped_shap_values.shape}")
    
    # Return reshaped SHAP values for the batch
    return reshaped_shap_values

def calculate_shap_values(X_tensor, model, num_processes=1, X_train_tensor=None):
    seed = 42
    torch.manual_seed(seed)
    np.random.seed(seed)

    # Ensure the model is in evaluation mode
    model.eval()

    # Move the model and data to the appropriate device
    device = next(model.parameters()).device
    model = model.to(device)
    X_tensor = X_tensor.to(device)
    
    n_background=10
    background = define_kmeans_background(X_train_tensor.cpu().numpy(), n_clusters=n_background).to(device)
    
    # Split the dataset into smaller batches for multiprocessing
    batches = np.array_split(X_tensor, num_processes)
    
    # Debugging: Print out the number of batches and samples per batch
    print(f"Total batches: {len(batches)}")
    for i, batch in enumerate(batches):
        print(f"Batch {i} size: {batch.shape[0]}")
    
    # Initialize an empty list to store the results
    results = []

    for i, batch in tqdm(enumerate(batches), total=len(batches), desc="Processing SHAP batches"):
        result = calculate_shap_for_batch(i, batch, model, device, background)
        results.append(result)
    
    # Concatenate the reshaped SHAP values from each batch row-wise (vertically)
    shap_values_combined = np.concatenate(results, axis=0)

    squeezed_array = np.squeeze(shap_values_combined)
    df_shap_values = pd.DataFrame(squeezed_array)
    
    return df_shap_values

use_feature_selection=True
n_repeats=10
num_p=30
permutation_split=None
permutation_block=None
importance_method='shap'
output_dir='new_agg_v27_modelA_train'
feature_importance_file=os.path.join(output_dir, f'kor_sfari_mssng.feature_importance.test_sample.shap.tsv.gz')

feature_importance = calculate_importance(
    model, X_test_tensor, y_test_tensor, importance_method=importance_method,
    n_repeats=n_repeats, num_processes=5, 
    permutation_split=None, permutation_block=None,
    X_train_tensor=X_train_tensor
)
# Convert to DataFrame and save as CSV
feature_importance.columns = X_test.columns.tolist()
feature_importance['SAMPLE'] = sample_ids
feature_importance.to_csv(feature_importance_file, sep='\t', index=False, compression='gzip') 
mean_shap_values = feature_importance.drop(columns=['SAMPLE']).abs().mean().reset_index()
mean_shap_df = mean_shap_values.rename(columns={'index': 'Feature', 0: 'MeanSHAP'})
mean_feature_importance_file = feature_importance_file.replace(".tsv.gz", ".mean.tsv.gz")
mean_shap_df.to_csv(mean_feature_importance_file, sep='\t', index=False, compression='gzip')









zarr_path = 'data/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v27.modelB_train.zarr'

root = zarr.open(zarr_path, mode='r')
header = [name.replace(',', '').replace(' ', '_') for name in list(root['metadata'].attrs['columns'])]
data = root['data'][:]
sample_ids1 = root['metadata'].attrs['sample_ids']    
header_index = {col: idx for idx, col in enumerate(header)}

df1 = pd.DataFrame(data, columns=header)

zarr_path = 'data/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v27.modelB_test.zarr'

root = zarr.open(zarr_path, mode='r')
header = [name.replace(',', '').replace(' ', '_') for name in list(root['metadata'].attrs['columns'])]
data = root['data'][:]
sample_ids2 = root['metadata'].attrs['sample_ids']    
header_index = {col: idx for idx, col in enumerate(header)}

df2 = pd.DataFrame(data, columns=header)

df_combined = pd.concat([df1, df2], axis=0).reset_index(drop=True)

sample_ids = sample_ids1 + sample_ids2

filtered_columns = sorted_high_f['Feature'].tolist() + ['label']
df_filtered = df_combined[filtered_columns]

X_test = df_filtered.drop(columns=['label'])
y_test = df_filtered['label']

# Scale data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_test)

X_test_tensor = torch.from_numpy(X_scaled.astype('float32')).float().to(device)
y_test_tensor = torch.from_numpy(y_test.values.astype('float32')).float().to(device)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = SimpleNN(X_test_tensor.shape[1]).to(device)
model.load_state_dict(torch.load(model_path, map_location=device, weights_only=False)) 
model.eval()

use_feature_selection=True
n_repeats=10
num_p=30
permutation_split=None
permutation_block=None
importance_method='shap'
output_dir='new_agg_v27_modelA_train'
feature_importance_file=os.path.join(output_dir, f'kor_sfari_mssng.feature_importance.modelB_sample.shap.tsv.gz')

feature_importance = calculate_importance(
    model, X_test_tensor, y_test_tensor, importance_method=importance_method,
    n_repeats=n_repeats, num_processes=5, 
    permutation_split=None, permutation_block=None,
    X_train_tensor=X_train_tensor
)
# Convert to DataFrame and save as CSV
feature_importance.columns = X_test.columns.tolist()
feature_importance['SAMPLE'] = sample_ids
feature_importance.to_csv(feature_importance_file, sep='\t', index=False, compression='gzip') 
mean_shap_values = feature_importance.drop(columns=['SAMPLE']).abs().mean().reset_index()
mean_shap_df = mean_shap_values.rename(columns={'index': 'Feature', 0: 'MeanSHAP'})
mean_feature_importance_file = feature_importance_file.replace(".tsv.gz", ".mean.tsv.gz")
mean_shap_df.to_csv(mean_feature_importance_file, sep='\t', index=False, compression='gzip')


predict_prob_file = os.path.join(output_dir, f'kor_sfari_mssng_predict_prob.modelB_sample.oversample.tsv.gz')

if predict_prob_file: 
    with torch.no_grad():
        outputs = model(X_test_tensor).squeeze()
        y_pred_prob = sigmoid(outputs).cpu().numpy() 
        y_pred_test = (y_pred_prob >= 0.5).astype(int)  

       
        results_df = pd.DataFrame({
            'SAMPLE': sample_ids,  # 샘플 ID
            'TrueLabel': y_test_tensor.cpu().numpy(),
            'PredictedLabel': y_pred_test,
            'PredictedProbability': y_pred_prob
        })

        # Calculate metrics
        true_labels = y_test_tensor.cpu().numpy()

        # AUROC
        auroc = roc_auc_score(true_labels, y_pred_prob)

        # AUPRC
        precision, recall, _ = precision_recall_curve(true_labels, y_pred_prob)
        auprc = auc(recall, precision)

        # Accuracy
        accuracy = accuracy_score(true_labels, y_pred_test)

        # F1 Score
        f1 = f1_score(true_labels, y_pred_test)

        # Print metrics
        print(f"AUROC: {auroc:.4f}")
        print(f"AUPRC: {auprc:.4f}")
        print(f"Accuracy: {accuracy:.4f}")
        print(f"F1 Score: {f1:.4f}")

        
        results_df.to_csv(predict_prob_file, sep='\t', index=False)
        print(f"Predict probability saved to {predict_prob_file}", flush=True)
