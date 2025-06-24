import time
import subprocess
import os
import sys

# Function to run a single job
def run_job(i):
    current_group_list = ['all', 'noncoding', 'coding']

    # Define local paths for inputs and outputs
    input_path = '/data1/deepLN/train/data/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v27.modelA_train.zarr'
    output_dir = f'/data1/deepLN/train/train_sagemaker_scripts/output/new_agg_v27_modelA_train/'
    os.makedirs(output_dir, exist_ok=True)

    # Command to run the training script (equivalent to SageMaker's script_processor.run)
    command = [
        'python3', f'train_{i}L.py',
        '--zarr_path', input_path,
        '--model_dir', '/data1/deepLN/train/train_sagemaker_scripts/model',
        '--wandb_key', '26ed6c564d5f320f8879b1517bc3f19ae4b1e859',
        '--epochs', '500',
        '--learning_rate', '0.000001',
        '--batch_size', '32',
        '--threshold', '0.5',
        "--group_list", *current_group_list,
        "--feature_info_path", "/data1/deepLN/train/train_sagemaker_scripts/table.coding_noncoding.feature_info_v27.20240428.txt",
        '--log_file', os.path.join(output_dir, f'kor_sfari_mssng_training_log_{i}L.oversample.csv'),
        "--gpu",  # Use GPU
        # '--stratify_col', 'strat_col',
        # '--importance_method', 'shap',
        # '--use_weight',
        '--oversampling',
        # '--predict_prob_file', os.path.join(output_dir, f'kor_sfari_mssng_predict_prob_{i}L.oversample.tsv.gz'),
        # '--sample_meta_path', "/data1/deepLN/train/train_sagemaker_scripts/Korean_SFARI_MSSNG.14428samples.sample_list.20241213.txt",
        # "--female",
        # "--male",
        # '--undersampling',
        # '--num_p', '40',
        # '--feature_selection_file', os.path.join(output_dir, f'kor_sfari_mssng.feature_selection_{i}L.perm.tsv.gz'),
        # '--use_feature_selection',
        # '--permutation_split', '5',
        # '--permutation_block', '0'
    ]

    # Run the training job and print output to terminal
    process = subprocess.Popen(command, stdout=sys.stdout, stderr=sys.stderr)
    process.communicate()

    print(f"Job {i} completed.", flush=True)

# Range of jobs to run
job_indices = range(7, 8)  # Adjust job indices as needed

# Run jobs sequentially using a for loop
for i in job_indices:
    run_job(i)

print('Done', flush=True)
