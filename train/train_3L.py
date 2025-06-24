import os
import subprocess
import sys
import tarfile

# Disable SageMaker Debugger
# os.environ['SAGEMAKER_DEBUGGER_ENABLED'] = 'false'
# subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "-q", "wandb", "zarr", "imblearn"])

from sklearn.utils.class_weight import compute_class_weight
import argparse
import random
import datetime
import torch
import shap
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

# 1. Define the SimpleNN class
class SimpleNN(nn.Module):
    def __init__(self, input_size):
        super(SimpleNN, self).__init__()
        self.fc1 = nn.Linear(input_size, 256)
        self.fc2 = nn.Linear(256, 64)
        self.fc3 = nn.Linear(64, 1)  # Output layer without sigmoid

    def forward(self, x):
        x = torch.relu(self.fc1(x)) 
        x = torch.relu(self.fc2(x)) 
        x = self.fc3(x)  # Raw output (logits), no sigmoid
        return x
    
# Function to extract .tar.gz file
def extract_zarr_file(tar_path, extract_path):
    print(f"Extracting {tar_path} to {extract_path}...", flush=True)
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=extract_path)
    print(f"Extraction complete.", flush=True)

def load_sample_meta(file_path):
    """
    Load sample meta information from a TXT file.
    The file should have at least two columns: SAMPLE and the feature (e.g., Sex).
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Sample meta information file not found: {file_path}")

    # Load meta information
    meta_df = pd.read_csv(file_path, sep='\t')

    # Ensure required columns exist
    if 'SAMPLE' not in meta_df.columns:
        raise ValueError("The sample meta file must contain a 'SAMPLE' column.")

    # Check for 'Sex' column and handle missing values
    if 'Sex' in meta_df.columns:
        missing_sex_count = meta_df['Sex'].isna().sum()
        if missing_sex_count > 0:
            print(f"Warning: {missing_sex_count} rows have NaN in the 'Sex' column and will be removed.", flush=True)
            meta_df = meta_df.dropna(subset=['Sex'])

    return meta_df

# Function to set the seed for reproducibility
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True  # Ensure deterministic behavior
        torch.backends.cudnn.benchmark = False

def get_output_path(input_path, suffix="_feature_importance.tsv.gz"):
    base_name = os.path.splitext(os.path.basename(input_path))[0]
    return f"{base_name}{suffix}"

# Define a function to save logs to a CSV file
def save_logs_to_csv(epoch, train_loss, eval_loss, test_acc, test_auroc, test_f1, test_precision, test_recall, log_file, 
                     best_f1, best_threshold_for_f1, best_accuracy, best_threshold_for_acc, AUPRC, overwrite=False):
    mode = 'w' if overwrite else 'a'
    
    with open(log_file, mode=mode, newline='') as file:
        writer = csv.writer(file)
        
        # Write the header if overwriting
        if overwrite:
            writer.writerow(['Epoch', 'Train Loss', 'Eval Loss', 'Test Accuracy', 'Test AUROC', 'Test F1', 'Test Precision', 
                             'Test Recall', 'Best F1', 'Best Threshold for F1', 'Best Accuracy', 'Best Threshold for Accuracy', 'AUPRC'])

        # Write the log data, including the best F1 score, best threshold for F1, best accuracy, best threshold for accuracy, and AUPRC
        writer.writerow([epoch, train_loss, eval_loss, test_acc, test_auroc, test_f1, test_precision, test_recall, 
                         best_f1, best_threshold_for_f1, best_accuracy, best_threshold_for_acc, AUPRC])

def balance_data(X_train, y_train, method, random_state=42):
    if method == 'oversampling':
        ros = RandomOverSampler(random_state=random_state)
        X_train_balanced, y_train_balanced = ros.fit_resample(X_train, y_train)
        print("Applied oversampling to the training set.")
    elif method == 'undersampling':
        rus = RandomUnderSampler(random_state=random_state)
        X_train_balanced, y_train_balanced = rus.fit_resample(X_train, y_train)
        print("Applied undersampling to the training set.")
    return X_train_balanced, y_train_balanced

def train_model(zarr_path, learning_rate=0.001, epochs=50, batch_size=32, threshold=1.0, group_list=None, feature_info_df=None,
                log_file=None, patience=10, device='cpu', male=False, female=False, sex_fixed=False, oversampling=False, undersampling=False,
                stratify_col=None, model_path=None, use_weight=False, sample_meta_df=None):  # Add stratify_col as an argument
    # Start time for training
    start_time = time.time()

    # Initialize wandb
    wandb.init(project="simple_nn_project", config={
        "learning_rate": learning_rate,
        "epochs": epochs,
        "batch_size": batch_size,
        "threshold": threshold,
        "group_list": group_list,
        "male_filter": male,
        "female_filter": female,
        "sex_fixed": sex_fixed
    })

    root = zarr.open(zarr_path, mode='r')
    header = [name.replace(',', '').replace(' ', '_') for name in list(root['metadata'].attrs['columns'])]
    
    # Extract metadata columns
    metadata_columns = root['metadata'].attrs['columns']

    # Save as TSV
    output_path = "metadata_columns.tsv"  # Replace with your desired output file path
    pd.DataFrame(metadata_columns, columns=["Column_Name"]).to_csv(output_path, sep="\t", index=False)
    
    data = root['data'][:]
    sample_ids = root['metadata'].attrs['sample_ids']
    
    header_index = {col: idx for idx, col in enumerate(header)}
    
    del root

    if sample_meta_df is not None:
        # Merge sample meta information with data
        meta_sample_ids = sample_meta_df['SAMPLE'].values

        missing_samples = set(sample_ids) - set(sample_meta_df['SAMPLE'])
        if missing_samples:
            print(f"Warning: {len(missing_samples)} samples in the dataset are missing in the sample_meta_df.", flush=True)

        # Filter data to include only samples present in the meta file
        mask = np.isin(sample_ids, meta_sample_ids)
        data = data[mask]
        sample_ids = np.array(sample_ids)[mask].tolist()

        # Merge the meta information with the data for indexing
        sample_meta_df = sample_meta_df[sample_meta_df['SAMPLE'].isin(sample_ids)]

    # Filter based on sex
    if sample_meta_df is not None and 'Sex' in sample_meta_df.columns:
        # Validate 'Sex' column
        invalid_sex_values = sample_meta_df[~sample_meta_df['Sex'].isin(['Male', 'Female'])]
        if not invalid_sex_values.empty:
            raise ValueError(
                f"Invalid values found in 'Sex' column: {invalid_sex_values['Sex'].unique()}. "
                "Allowed values are 'Male' and 'Female'."
            )
        if args.male:
            print("Filtering for male samples only.", flush=True)
            male_samples = sample_meta_df[sample_meta_df['Sex'] == 'Male']['SAMPLE'].values
            mask = np.isin(sample_ids, male_samples)
            data = data[mask]
            sample_ids = np.array(sample_ids)[mask].tolist()
        if args.female:
            print("Filtering for female samples only.", flush=True)
            female_samples = sample_meta_df[sample_meta_df['Sex'] == 'Female']['SAMPLE'].values
            mask = np.isin(sample_ids, female_samples)
            data = data[mask]
            sample_ids = np.array(sample_ids)[mask].tolist()

    # Optional filtering based on group_list and feature_info_df
    if group_list and feature_info_df is not None:
        if 'Feature' not in feature_info_df.columns or 'Group' not in feature_info_df.columns:
            raise ValueError("The feature info file must contain 'Feature' and 'Group' columns.")

        feature_info_df['Feature'] = feature_info_df['Feature'].str.replace(',', '').str.replace(' ', '_')
        filtered_feature_info = feature_info_df[feature_info_df['Group'].isin(group_list)]
        features_to_keep = filtered_feature_info['Feature'].values

        # Ensure the 'label' column is included
        if 'label' not in features_to_keep:
            features_to_keep = list(features_to_keep) + ['label']  # Add 'label' if missing

        feature_indices = [header_index[feature] for feature in features_to_keep if feature in header_index]
        
        # Update header and header_index based on filtered features
        updated_header = [header[i] for i in feature_indices]
        header_index = {col: idx for idx, col in enumerate(updated_header)}

        # Update label_idx to reflect new header_index
        label_idx = header_index.get('label')
        if label_idx is None:
            raise ValueError("The 'label' column is missing after feature filtering.")
        
        sex_in_features = 'Sex' in features_to_keep

        if 'Sex' in features_to_keep:
            print("Sex is included as a feature for training.", flush=True)
        else:
            print("Sex is excluded as a feature. It will only be used for stratification.", flush=True)

        data = data[:, feature_indices]
        feature_indices = list(range(data.shape[1]))
    else:
        # If no filtering is applied, use the full header
        updated_header = header
        feature_indices = list(range(len(header)))
        sex_in_features = False

    # Apply threshold for dropping rows with NaN values
    threshold_value = int(data.shape[1] * threshold)
    non_nan_counts = np.sum(~np.isnan(data), axis=1)
    mask = non_nan_counts >= threshold_value    
    # Filter rows based on the threshold
    data = data[mask]

    stratify_data = None  # Ensure stratify_data is always initialized
    label_idx = header_index.get('label')

    # Handle stratification column
    if stratify_col and stratify_col in header_index:
        # Get the stratify column index
        stratify_idx = header_index[stratify_col]
        stratify_data = data[:, stratify_idx]  # Extract stratify data
        # Exclude stratify column and label from feature indices
        feature_indices = [i for i in feature_indices if i not in {label_idx, stratify_idx}]
    else:
        # Handle cases when stratify column is not specified or not in header_index
        if sample_meta_df is not None and 'Sex' in sample_meta_df.columns:
            sex_map = dict(zip(sample_meta_df['SAMPLE'], sample_meta_df['Sex']))
            sex_labels = np.array([sex_map.get(sample, np.nan) for sample in sample_ids])
            if pd.isna(sex_labels).any():
                raise ValueError("Some samples in the dataset are missing 'Sex' information in the sample_meta_df.")
            
            feature_indices = [i for i in feature_indices if i not in {label_idx}]
            # Combine 'Sex' and label for stratification if sex_fixed is enabled
            if sex_fixed:
                stratify_data = np.char.add(sex_labels.astype(str), '_' + data[:, label_idx].astype(str))
            else:
                stratify_data = data[:, label_idx]  # Default to label-only stratification
        else:
            # Handle cases when sample_meta_df is not provided or does not include 'Sex'
            if not sex_in_features and 'Sex' in header_index:
                # Exclude 'Sex' column if not used as a feature
                sex_idx = header_index['Sex']
                feature_indices = [i for i in feature_indices if i not in {label_idx, sex_idx}]
            else:
                feature_indices = [i for i in feature_indices if i not in {label_idx}]
        
            # Determine stratify_data based on sex_fixed flag
            if not sex_fixed:
                stratify_data = data[:, label_idx]  # Use label directly for stratification
            else:
                # Combine 'Sex' and label if 'Sex' exists, otherwise use label
                sex_idx = header_index.get('Sex', None)
                stratify_data = (
                    (data[:, sex_idx].astype(str) + '_' + data[:, label_idx].astype(str))
                    if sex_idx is not None else data[:, label_idx]
                )

    if 'label' in updated_header:
        feature_names = [col for col in updated_header if col != 'label']
    else:
        feature_names = updated_header
    
    num_features = len(feature_indices)
    print(f"Number of features being used: {num_features}", flush=True)

    X = data[:, feature_indices]
    y = data[:, label_idx]

    # 3. Data preprocessing
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    if sample_meta_df is not None and 'Set' in sample_meta_df.columns:
        # Ensure only samples in the meta file are used
        sample_meta_df = sample_meta_df[sample_meta_df['SAMPLE'].isin(sample_ids)]
        sample_ids = sample_meta_df['SAMPLE'].values
        mask = np.isin(sample_ids, sample_meta_df['SAMPLE'])
        X = X[mask]
        y = y[mask]

        # Filter by 'Set' column
        train_samples = sample_meta_df[sample_meta_df['Set'].str.lower() == 'train']['SAMPLE'].values
        test_samples = sample_meta_df[sample_meta_df['Set'].str.lower() == 'test']['SAMPLE'].values
        
        print(f"Unique Set values: {sample_meta_df['Set'].unique()}", flush=True)
        print(f"Train samples: {train_samples}", flush=True)
        print(f"Test samples: {test_samples}", flush=True)

        # Create masks for train and test sets
        train_mask = np.isin(sample_ids, train_samples)
        test_mask = np.isin(sample_ids, test_samples)

        # Assign datasets
        X_train = X[train_mask, :]
        y_train = y[train_mask]
        X_test = X[test_mask, :]
        y_test = y[test_mask]
        
        # Assign matching sample IDs for train and test
        s_train = sample_ids[train_mask]
        s_test = sample_ids[test_mask]
    else:
        # Now split the data using this combined stratification column
        # 4. Split the data into training and testing sets
        X_train, X_test, y_train, y_test, s_train, s_test = train_test_split(X, y, sample_ids, test_size=(3/17), random_state=42, stratify=stratify_data)

    # 3. 훈련 세트에서 선택한 방법으로 클래스 비율 맞추기
    if oversampling:
        X_train_balanced, y_train_balanced = balance_data(X_train, y_train, method='oversampling', random_state=42)
    elif undersampling:
        X_train_balanced, y_train_balanced = balance_data(X_train, y_train, method='undersampling', random_state=42)
    else:
        X_train_balanced, y_train_balanced = X_train, y_train


    # 5. Convert the training and testing sets to PyTorch tensors
    X_train_tensor = torch.from_numpy(X_train_balanced).float().to(device)
    y_train_tensor = torch.from_numpy(y_train_balanced).float().to(device)
    X_test_tensor = torch.from_numpy(X_test).float().to(device)
    y_test_tensor = torch.from_numpy(y_test).float().to(device)

    model = SimpleNN(X_train_tensor.shape[1]).to(device)

    if use_weight:
        class_weights = compute_class_weight(
            class_weight='balanced', 
            classes=np.unique(y_train_balanced), 
            y=y_train_balanced
        )
        class_weights = torch.tensor([class_weights[0]/class_weights[1]], dtype=torch.float).to(device)
        
        # For binary classification, class_weights should have size 2 (for class 0 and class 1)
        # BCEWithLogitsLoss expects weights in the shape of (num_classes, ), which is [2] in binary classification
        criterion = nn.BCEWithLogitsLoss(pos_weight=class_weights)
    else:
        criterion = nn.BCEWithLogitsLoss()

    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

    # Initialize early stopping and best model tracking
    best_eval_loss = float('inf')
    best_model_state = None
    best_epoch = 0
    no_improvement_epochs = 0

    # 9. Train the model
    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        for inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs).squeeze()
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()

        avg_train_loss = running_loss / len(train_loader)

        # Evaluate the model on the test set
        model.eval()
        with torch.no_grad():
            outputs = model(X_test_tensor).squeeze()
            eval_loss = criterion(outputs, y_test_tensor).item()
            y_pred_prob = torch.sigmoid(outputs).cpu().numpy()
            y_pred_test = (y_pred_prob >= 0.5).astype(int)
            
            y_pred_test = torch.tensor(y_pred_test)

            # Calculate metrics
            test_acc = accuracy_score(y_test_tensor.cpu(), y_pred_test.cpu())
            test_auroc = roc_auc_score(y_test_tensor.cpu(), y_pred_prob)
            test_f1 = f1_score(y_test_tensor.cpu(), y_pred_test.cpu())
            test_precision = precision_score(y_test_tensor.cpu(), y_pred_test.cpu())
            test_recall = recall_score(y_test_tensor.cpu(), y_pred_test.cpu())

            # Now calculate metrics across thresholds to find the optimal F1 and Accuracy thresholds
            thresholds = np.linspace(0, 1, 101)  # 101 values from 0.0 to 1.0 in steps of 0.01
            f1_scores = []
            acc_scores = []
            
            for threshold in thresholds:
                y_pred_threshold = (y_pred_prob >= threshold).astype(int)
                f1 = f1_score(y_test_tensor.cpu(), y_pred_threshold)
                acc = accuracy_score(y_test_tensor.cpu(), y_pred_threshold)
                f1_scores.append(f1)
                acc_scores.append(acc)

            # Compute Precision-Recall Curve and AUPRC
            precision, recall, _ = precision_recall_curve(y_test, y_pred_prob)
            auprc = auc(recall, precision)
            
            # Find the best F1 score and corresponding threshold
            best_f1 = max(f1_scores)
            best_threshold_for_f1 = thresholds[f1_scores.index(best_f1)]

            # Find the best accuracy and corresponding threshold
            best_acc = max(acc_scores)
            best_threshold_for_acc = thresholds[acc_scores.index(best_acc)]
            
            # Log the best F1 score, accuracy, AUPRC, and thresholds to WandB
            wandb.log({
                "epoch": epoch + 1,
                "train_loss": avg_train_loss,
                "eval_loss": eval_loss,
                "test_accuracy": test_acc,
                "test_auroc": test_auroc,
                "test_f1": test_f1,
                "test_precision": test_precision,
                "test_recall": test_recall,
                "best_f1": best_f1,
                "best_threshold_for_f1": best_threshold_for_f1,
                "best_accuracy": best_acc,
                "best_threshold_for_acc": best_threshold_for_acc,
                "AUPRC": auprc
            })

            # Save logs to CSV file
            # For the first epoch, use overwrite=True to overwrite the existing file
            save_logs_to_csv(epoch + 1, avg_train_loss, eval_loss, test_acc, test_auroc, test_f1, test_precision, test_recall, log_file,
                            best_f1=best_f1, best_threshold_for_f1=best_threshold_for_f1,
                            best_accuracy=best_acc, best_threshold_for_acc=best_threshold_for_acc,
                            AUPRC=auprc, overwrite=(epoch == 0))


        print(f'Epoch {epoch+1}/{epochs}, Train Loss: {avg_train_loss:.4f}, Eval Loss: {eval_loss:.4f}, Test Accuracy: {test_acc:.4f}, Test AUROC: {test_auroc:.4f}, Test F1: {test_f1:.4f}',
              flush=True)

        # Save the best model if improvement is greater than or equal to 0.0001
        if best_eval_loss - eval_loss >= 0.0001:
            best_eval_loss = eval_loss
            best_model_state = model.state_dict()
            best_epoch = epoch + 1  # Update the best epoch number
            no_improvement_epochs = 0  # Reset the counter if improvement happens
        else:
            no_improvement_epochs += 1

        if no_improvement_epochs >= patience:
            print(f"Early stopping at epoch {epoch+1}", flush=True)
            break


    # End time for training
    end_time = time.time()
    
    # Calculate and print the total training time
    training_time = end_time - start_time
    print(f"Total training time: {training_time:.2f} seconds", flush=True)

    # 10. Save the best trained model
    if best_model_state is not None:
        torch.save(best_model_state, model_path)
        print(f"Best model saved as {model_path} with eval loss: {best_eval_loss:.4f} at epoch {best_epoch}", flush=True)

    # 11. Evaluate final model accuracy using the best model
    if best_model_state is not None:
        model.load_state_dict(best_model_state)

    # 11. Evaluate final model accuracy
    with torch.no_grad():
        y_pred_train = torch.sigmoid(model(X_train_tensor).squeeze()).round()
        y_pred_test = torch.sigmoid(model(X_test_tensor).squeeze()).round()

        # Move tensors to CPU before passing to accuracy_score
        train_acc = accuracy_score(y_train_tensor.cpu(), y_pred_train.cpu())
        test_acc = accuracy_score(y_test_tensor.cpu(), y_pred_test.cpu())

    wandb.log({"train_accuracy": train_acc, "test_accuracy": test_acc})

    print(f'Train Accuracy: {train_acc:.4f}', flush=True)
    print(f'Test Accuracy: {test_acc:.4f}', flush=True)

    wandb.finish()

    # Return the test tensors, model, and feature names for further evaluation
    return X_train_tensor, y_train_tensor, X_test_tensor, y_test_tensor, feature_names, s_train, s_test, X_train, y_train

# Ensure the "spawn" start method is used for multiprocessing
multiprocessing.set_start_method('spawn', force=True)

from sklearn.cluster import KMeans

def define_kmeans_background(X_train, n_clusters=10):
    """
    Use KMeans to define cluster centers as background data for SHAP.
    """
    kmeans = KMeans(n_clusters=n_clusters, random_state=42).fit(X_train)
    return torch.tensor(kmeans.cluster_centers_, dtype=torch.float32)

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
    
    # background_data = X_train_tensor.cpu().numpy()[np.random.choice(X_train_tensor.cpu().numpy().shape[0], 1000, replace=False)]
    # background = torch.tensor(background_data, dtype=torch.float32).to(device)
    
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

def calculate_importance_for_column(model, X_tensor, y_tensor, baseline_loss, col, n_repeats=10, random_state=42):
    rng = np.random.RandomState(random_state)
    
    original_col = X_tensor[:, col].clone()
    scores = []
    
    for _ in range(n_repeats):
        # Permute the column
        X_tensor[:, col] = X_tensor[:, col][torch.randperm(X_tensor.size(0))]
        
        # Predict with the permuted data
        new_preds = model(X_tensor).squeeze().detach().cpu().numpy()
        new_loss = np.mean((new_preds - y_tensor.cpu().numpy()) ** 2)
        
        # Calculate the score (importance)
        scores.append(new_loss - baseline_loss)
        
        # Restore the original column
        X_tensor[:, col] = original_col

    return torch.tensor(np.mean(scores), dtype=torch.float32)

def custom_permutation_importance(model, X_tensor, y_tensor, n_repeats=10, random_state=42,
                                  num_processes=None, permutation_split=1, permutation_block=0):
    rng = np.random.RandomState(random_state)
    
    # Calculate baseline loss
    baseline_preds = model(X_tensor).squeeze().detach().cpu().numpy()
    baseline_loss = np.mean((baseline_preds - y_tensor.cpu().numpy()) ** 2)
    
    # Split features into blocks for permutation
    num_features = X_tensor.shape[1]
    block_size = num_features // permutation_split
    start_idx = block_size * permutation_block
    end_idx = start_idx + block_size if permutation_block < permutation_split - 1 else num_features
    
    # Select the column indices for the specified block
    col_indices = list(range(start_idx, end_idx))
    print(f"Calculating permutation importance for features in block {permutation_block + 1}/{permutation_split} (features {start_idx} to {end_idx - 1})")
    
    # Define a partial function for multiprocessing, leaving out 'col'
    partial_func = partial(calculate_importance_for_column, model, X_tensor, y_tensor, baseline_loss, n_repeats=n_repeats, random_state=random_state)
    
    # Use multiprocessing to calculate importances in parallel
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Pass col_indices directly as individual items to partial_func
        results = list(tqdm(pool.imap(partial_func, col_indices), total=len(col_indices), desc="Calculating Permutation Importance"))
    
    # Store the results in the importances tensor
    importances = torch.zeros(num_features, dtype=torch.float32)
    for idx, importance in zip(col_indices, results):
        importances[idx] = importance
    
    return importances.numpy()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a neural network on a Zarr dataset.")
    parser.add_argument('--zarr_path', type=str, required=True, help='Path to the Zarr file.')
    parser.add_argument('--aws', action='store_true', help='Flag to use extracted Zarr data if set.')
    parser.add_argument('--wandb_key', type=str, required=True, help='WandB API key for logging.')

    parser.add_argument('--learning_rate', type=float, default=0.001, help='Learning rate for the optimizer.')
    parser.add_argument('--epochs', type=int, default=50, help='Number of epochs to train.')
    parser.add_argument('--batch_size', type=int, default=32, help='Batch size for training.')
    parser.add_argument('--threshold', type=float, default=1.0, help='Fraction of non-NA values required to retain a row.')
    parser.add_argument('--calculate_importance', action='store_true', help='Whether to calculate feature importance.')
    parser.add_argument('--use_feature_selection', action='store_true', help='Whether to apply feature selection using the training data.')
    parser.add_argument('--importance_method', choices=['permutation', 'shap'], default='permutation', 
                            help="Method to use for calculating feature importance: 'permutation' or 'shap'.")
    
    parser.add_argument('--n_repeats', type=int, default=10, help='Number of permutations to perform for feature importance.')
    parser.add_argument('--num_p', type=int, default=multiprocessing.cpu_count(), help='Number of processes to use in multiprocessing.')
    parser.add_argument('--group_list', type=str, nargs='*', help='List of groups to filter features.')
    parser.add_argument('--feature_info_path', type=str, help='Path to the feature info file.')
    parser.add_argument('--gpu', action='store_true', help='Flag to use GPU if available, otherwise use CPU.')
    parser.add_argument('--seed', type=int, default=42, help='Seed for reproducibility.')

    parser.add_argument('--predict_prob_file', type=str, default=None, 
                        help="Output file name for predicted probabilities (e.g., 'predict_prob.tsv.gz' for a compressed file).")

    parser.add_argument('--log_file', type=str, default=None, 
                        help="Output file name for the training log (e.g., 'training_log.csv.gz' for a compressed file).")
    parser.add_argument('--feature_importance_file', type=str, default=None, 
                        help="Output file name for feature importance (e.g., 'feature_importance.tsv.gz' for a compressed file).")
    parser.add_argument('--feature_selection_file', type=str, default=None, 
                        help="Output file name for feature selection (e.g., 'feature_selection.tsv.gz' for a compressed file).")
    parser.add_argument('--stratify_col', type=str, help="Column to use for stratification; if provided, this column won't be included in training features.")

    # Adding new arguments for feature importance calculation on specific blocks
    parser.add_argument('--permutation_split', type=int, default=1, help="Number of splits to divide the features into for permutation.")
    parser.add_argument('--permutation_block', type=int, default=0, help="The block index (0-based) to perform permutation on for feature importance calculation.")

    parser.add_argument('--sample_meta_path', type=str, required=False,
                        help="Path to a TXT file containing sample meta information (e.g., Sex).")

    gender_group = parser.add_mutually_exclusive_group(required=False)
    gender_group.add_argument('--male', action='store_true', help='Filter for male (Sex == 0) samples.')
    gender_group.add_argument('--female', action='store_true', help='Filter for female (Sex == 1) samples.')

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--oversampling', action='store_true', help='Apply oversampling to balance classes in the training set.')
    group.add_argument('--undersampling', action='store_true', help='Apply undersampling to balance classes in the training set.')
    group.add_argument('--use_weight', action='store_true', help='Whether to apply class weights during training.')

    # Stratify by sex and label if --sex_fixed is provided
    parser.add_argument('--sex_fixed', action='store_true', help='Stratify based on both Sex and label.')

    parser.add_argument('--model_dir', type=str, default='/opt/ml/processing/output/data', help='Directory where the model will be saved.')
    

    args = parser.parse_args()

    # **Generate a unique timestamp before training starts**
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    random_number = random.randint(1000, 9999)  # Optional random number for added uniqueness
    model_filename = f"{timestamp}_{random_number}_best_model.pth"
    model_path = os.path.join(args.model_dir, model_filename)

    # Log in to WandB with the provided API key
    wandb.login(key=args.wandb_key)
    print("Successfully logged in to WandB.")

    # Set the seed for reproducibility
    set_seed(args.seed)

    sample_meta_df = None
    if args.sample_meta_path:
        sample_meta_df = load_sample_meta(args.sample_meta_path)
        print(f"Sample meta information loaded from {args.sample_meta_path}", flush=True)
        if 'SAMPLE' not in sample_meta_df.columns:
            raise ValueError("The sample meta file must contain a 'SAMPLE' column.")

        if args.sex_fixed and 'Sex' not in sample_meta_df.columns:
            raise ValueError("The sample meta file must contain a 'Sex' column for sex_fixed stratification.")
        
        if args.sex_fixed:
            print("Stratifying data using the 'Sex' column in the meta file.", flush=True)
        if args.male:
            print("Filtering for male samples.", flush=True)
        if args.female:
            print("Filtering for female samples.", flush=True)

    # Extract Zarr file if --aws flag is provided
    if args.aws:
        extract_zarr_file(args.zarr_path, "/opt/ml/processing/input/data/")
        zarr_file_name = os.path.basename(args.zarr_path).replace('.tar.gz', '')  # Remove '.tar.gz'
        extracted_zarr_path = os.path.join("/opt/ml/processing/input/data/", zarr_file_name)
        print(f"Using extracted Zarr data from {extracted_zarr_path}", flush=True)
    else:
        extracted_zarr_path = args.zarr_path  # Use the provided zarr_path if not using AWS
        print(f"Using original Zarr data from {extracted_zarr_path}", flush=True)
    
    device = torch.device('cuda' if args.gpu and torch.cuda.is_available() else 'cpu')

    if args.gpu and not torch.cuda.is_available():
        print("Warning: GPU flag is set but no GPU is available. Using CPU instead.", flush=True)
    print(f"Using device: {device}", flush=True)

    # Check if the zarr_path exists
    if not os.path.exists(args.zarr_path):
        raise FileNotFoundError(f"The specified Zarr file path does not exist: {args.zarr_path}")

    # Load feature_info_df if both group_list and feature_info_path are provided
    if args.group_list and args.feature_info_path:
        if not os.path.exists(args.feature_info_path):
            raise FileNotFoundError(f"The specified feature info file path does not exist: {args.feature_info_path}")
        
        feature_info_df = pd.read_csv(args.feature_info_path, sep='\t')
    else:
        feature_info_df = None

    # Determine the output paths for logs, feature importance, and feature selection files
    log_file = args.log_file if args.log_file else get_output_path(args.zarr_path, suffix="_training_log.csv")
    feature_importance_file = args.feature_importance_file if args.feature_importance_file else get_output_path(args.zarr_path, suffix="_feature_importance.tsv.gz")
    feature_selection_file = args.feature_selection_file if args.feature_selection_file else get_output_path(args.zarr_path, suffix="_feature_selection.tsv.gz")
    predict_prob_file = args.predict_prob_file if args.predict_prob_file else get_output_path(args.zarr_path, suffix="_predict_prob.tsv.gz")

    # Train the model and get the test tensors, and feature names
    X_train_tensor, y_train_tensor, X_test_tensor, y_test_tensor, feature_names, s_train, s_test, X_train, y_train = train_model(zarr_path=extracted_zarr_path,
                                                                                                                                learning_rate=args.learning_rate, 
                                                                                                                                epochs=args.epochs, 
                                                                                                                                batch_size=args.batch_size,
                                                                                                                                threshold=args.threshold,
                                                                                                                                group_list=args.group_list,
                                                                                                                                feature_info_df=feature_info_df,
                                                                                                                                log_file=log_file,
                                                                                                                                device=device,
                                                                                                                                male=args.male,
                                                                                                                                female=args.female,
                                                                                                                                sex_fixed=args.sex_fixed,
                                                                                                                                oversampling=args.oversampling,
                                                                                                                                undersampling=args.undersampling,
                                                                                                                                stratify_col=args.stratify_col,
                                                                                                                                model_path=model_path,
                                                                                                                                use_weight=args.use_weight,
                                                                                                                                sample_meta_df=sample_meta_df)
    
    print(f"Saving model as: {model_path}")

    model = SimpleNN(X_train_tensor.shape[1]).to(device)
    model.load_state_dict(torch.load(model_path, map_location=device, weights_only=True))  # Ensure the model is loaded to the right device
    model.eval()
    
    root = zarr.open(extracted_zarr_path, mode='r')
    sample_ids = root['metadata'].attrs['sample_ids']

    if args.predict_prob_file:
        with torch.no_grad():
            outputs = model(X_test_tensor).squeeze()
            y_pred_prob = torch.sigmoid(outputs).cpu().numpy()
            y_pred_test = (y_pred_prob >= 0.5).astype(int)

            results_df = pd.DataFrame({
                'SAMPLE': s_test,
                'TrueLabel': y_test_tensor.cpu().numpy(),
                'PredictedLabel': y_pred_test,
                'PredictedProbability': y_pred_prob
            })
            # Save the DataFrame to a TSV file
            results_df.to_csv(predict_prob_file, sep='\t', index=False)
            print(f"Predict probability saved to {predict_prob_file}", flush=True)

    if args.calculate_importance:
        start_time = time.time()
        X_train_unbalanced_tensor = torch.from_numpy(X_train).float().to(device)
        # Calculate feature importance using the selected method (SHAP or Permutation)
        feature_importance = calculate_importance(
            model, X_test_tensor, y_test_tensor, importance_method=args.importance_method,
            n_repeats=args.n_repeats, num_processes=args.num_p, 
            permutation_split=args.permutation_split, permutation_block=args.permutation_block,
            X_train_tensor=X_train_unbalanced_tensor
        )
        if args.importance_method=='shap':
            feature_importance.columns = feature_names
            feature_importance['SAMPLE'] = s_test
            feature_importance.to_csv(feature_importance_file, sep='\t', index=False, compression='gzip')
            
            mean_shap_values = feature_importance.drop(columns=['SAMPLE']).abs().mean().reset_index()
            mean_shap_df = mean_shap_values.rename(columns={'index': 'Feature', 0: 'MeanSHAP'})
            mean_feature_importance_file = feature_importance_file.replace(".tsv.gz", ".mean.tsv.gz")
            mean_shap_df.to_csv(mean_feature_importance_file, sep='\t', index=False, compression='gzip')
        else:
            # Convert to DataFrame and save as CSV
            df_feature_importance = pd.DataFrame({'Feature': feature_names, 'Importance': feature_importance})
            df_feature_importance.to_csv(feature_importance_file, sep='\t', index=False, compression='gzip')
        elapsed_time = time.time() - start_time
        print(f"Feature Importance saved to {feature_importance_file}", flush=True)
        print(f"Time taken for feature importance calculation: {elapsed_time:.2f} seconds", flush=True)

    
    if args.use_feature_selection:
        start_time = time.time()
        # Calculate feature selection importance using the selected method (SHAP or Permutation)
        X_train_unbalanced_tensor = torch.from_numpy(X_train).float().to(device)
        y_train_unbalanced_tensor = torch.from_numpy(y_train).float().to(device)
        feature_importance = calculate_importance(
            model, X_train_unbalanced_tensor, y_train_unbalanced_tensor, importance_method=args.importance_method,
            n_repeats=args.n_repeats, num_processes=args.num_p, 
            permutation_split=args.permutation_split, permutation_block=args.permutation_block,
            X_train_tensor=X_train_unbalanced_tensor
        )
        
        if args.importance_method=='shap':
            # Convert to DataFrame and save as CSV
            feature_importance.columns = feature_names
            feature_importance['SAMPLE'] = s_train
            feature_importance.to_csv(feature_selection_file, sep='\t', index=False, compression='gzip') 
            
            mean_shap_values = feature_importance.drop(columns=['SAMPLE']).abs().mean().reset_index()
            mean_shap_df = mean_shap_values.rename(columns={'index': 'Feature', 0: 'MeanSHAP'})
            mean_feature_selection_file = feature_selection_file.replace(".tsv.gz", ".mean.tsv.gz")
            mean_shap_df.to_csv(mean_feature_selection_file, sep='\t', index=False, compression='gzip')

        else:
            # Convert to DataFrame and save as CSV
            df_feature_selection = pd.DataFrame({'Feature': feature_names, 'Importance': feature_importance})
            df_feature_selection.to_csv(feature_selection_file, sep='\t', index=False, compression='gzip')
        elapsed_time = time.time() - start_time
        print(f"Feature Selection results saved to {feature_selection_file}", flush=True)
        print(f"Time taken for feature importance calculation: {elapsed_time:.2f} seconds", flush=True)
        

