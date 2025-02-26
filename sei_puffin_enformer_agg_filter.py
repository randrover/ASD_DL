import zarr
import numpy as np
import pandas as pd
from multiprocessing import Pool, Manager
from tqdm import tqdm
import sys
import argparse

def filter_data_by_variants(mat, var_ids, sample_ids, filter_variant):
    """
    Filters the data matrix (mat) and the variant ids (var_ids) based on a list of variants.
    
    Parameters:
    - mat: NumPy or Zarr array containing the data matrix.
    - var_ids: Array-like (list, NumPy array, or Zarr array) containing variant IDs.
    - filter_variant: A pandas DataFrame or similar with the 'variant' column containing the variants to filter by.
    
    Returns:
    - filtered_mat: The data matrix filtered based on the matching variant IDs.
    - filtered_variants: The variant IDs that match the filter.
    """
    
    # Convert var_ids to a NumPy array if it's a Zarr array or list
    var_ids = np.array(var_ids)  # Convert to NumPy array if it's not already
    sample_ids = np.array(sample_ids)
    
    # Ensure variant_set is a set for faster membership checking
    variant_set = set(filter_variant['variant'])  # Assuming filter_variant is a pandas DataFrame
    
    # Use np.isin to get a boolean mask
    variant_mask = np.isin(var_ids, list(variant_set))  # np.isin returns a boolean array
    
    # Filter the data matrix (mat) and variants using the boolean mask
    filtered_mat = mat[variant_mask]  # Assuming mat is a NumPy or Zarr array
    filtered_variants = var_ids[variant_mask]  # Now var_ids is a NumPy array
    filtered_sample_ids = sample_ids[variant_mask]
    
    return filtered_mat, filtered_variants, filtered_sample_ids

# Define functions for chunked processing
def compute_variant_aggregations(chunk):
    """Compute per-variant aggregations for a chunk of rows."""
    variant_sum = np.sum(chunk, axis=1)
    variant_mean = np.mean(chunk, axis=1)
    variant_median = np.median(chunk, axis=1)
    variant_extreme = chunk[np.arange(chunk.shape[0]), np.argmax(np.abs(chunk), axis=1)]
    return np.vstack((variant_sum, variant_mean, variant_median, variant_extreme)).T

def process_chunk_with_progress(args):
    """Wrapper function for multiprocessing with progress tracking."""
    chunk, func, progress_queue = args
    result = func(chunk)
    progress_queue.put(1)  # Notify the main process that one chunk is done
    return result

def process_in_chunks_multiprocessing(data, func, axis=0, chunk_size=10000, n_cpus=None):
    """
    Process data in parallel using multiprocessing with a visible progress bar.
    Args:
        data: The Zarr array.
        func: Function to apply to each chunk.
        axis: Axis along which to chunk (0 for rows, 1 for columns).
        chunk_size: Number of rows/columns per chunk.
        n_cpus: Number of CPUs to use (default: all available CPUs).
    Returns:
        Aggregated result as a numpy array.
    """
    # Create chunk indices
    if axis == 0:  # Row-based processing
        chunk_indices = [(i, min(i + chunk_size, data.shape[0])) for i in range(0, data.shape[0], chunk_size)]
        chunks = [(data[start:end, :], func) for start, end in chunk_indices]
    else:  # Column-based processing
        chunk_indices = [(i, min(i + chunk_size, data.shape[1])) for i in range(0, data.shape[1], chunk_size)]
        chunks = [(data[:, start:end].T, func) for start, end in chunk_indices]

    # Shared progress queue for tracking progress
    with Manager() as manager:
        progress_queue = manager.Queue()
        with Pool(processes=n_cpus) as pool:
            # Start the progress bar
            with tqdm(total=len(chunks), desc="Processing Chunks", file=sys.stdout, leave=True) as pbar:
                # Submit tasks with the progress queue
                results = pool.map(process_chunk_with_progress, [(chunk, func, progress_queue) for chunk, func in chunks])

                # Update progress bar as tasks complete
                for _ in range(len(chunks)):
                    progress_queue.get()
                    pbar.update(1)

        # Combine results
        return np.vstack(results)
    
def save_to_zarr(df, output_path, agg_type, filter_vcf, groups_to_keep):
    sample_ids = df['SAMPLE'].tolist()
    print(sample_ids[0:5], flush = True)
    
    if 'SAMPLE' in df.columns:
        df.drop(columns=['SAMPLE'], inplace=True)
        
    column_names = df.columns.tolist()
    
    print(column_names[0:5], flush = True)
    print(f"Final df.shape: {df.shape}", flush=True)
    
    root = zarr.open(output_path, mode='w')
    root.create_group('metadata')
    root['metadata'].attrs['agg_type'] = agg_type
    root['metadata'].attrs['filter_vcf'] = filter_vcf if filter_vcf else 'None'
    root['metadata'].attrs['keep_feature'] = groups_to_keep
    root['metadata'].attrs['columns'] = column_names
    root['metadata'].attrs['sample_ids'] = sample_ids
    root.create_dataset('data', data=df.values, chunks=(1000, 1000), dtype='float')
    print(f'Data saved to {output_path}', flush=True)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process data and save to Zarr")
    parser.add_argument('--filter_variant_path', type=str, required=True, help="Path to the filter variant file")
    parser.add_argument('--output_path', type=str, required=True, help="Output path for Zarr file")
    parser.add_argument('--sample_list_path', type=str, required=True, help="Path to the sample list file")
    parser.add_argument('--n_cpus', type=int, default=1, help="Number of CPUs to use for multiprocessing")
    return parser.parse_args()

def main():
    # Parse arguments from the terminal
    args = parse_arguments()
    
    # Use the parsed arguments
    output_path = args.output_path
    sample_list_path = args.sample_list_path
    # Read the filter variant file
    filter_variant = pd.read_table(args.filter_variant_path)
    n_cpus = args.n_cpus  # Number of CPUs for multiprocessing

    print('Load Enformer', flush=True)

    enformer_path = '/data2/deepLN/Enformer/enformer_result_kor_909_sfari_mssng.zarr'

    root = zarr.open(enformer_path, mode='r')
    mat = root['data'][:]
    column_names = root['metadata'].attrs['columns']
    var_ids = root['metadata'].attrs['variant']
    sample_ids = root['metadata'].attrs['samples']

    # Example usage:
    filtered_mat, filtered_variants, filtered_sample_ids = filter_data_by_variants(mat, var_ids, sample_ids, filter_variant)

    # Step 1: Per-Variant Aggregation (process rows)
    print("Processing Per-Variant Aggregations...")
    variant_aggregates = process_in_chunks_multiprocessing(filtered_mat, compute_variant_aggregations, axis=0, chunk_size=10000, n_cpus=n_cpus)

    # Convert variant_aggregates and sample_ids into a pandas DataFrame for easier grouping
    variant_df = pd.DataFrame(variant_aggregates, columns=["sum", "mean", "median", "directional_extreme"])
    variant_df["SAMPLE"] = filtered_sample_ids  # Add sample_id as a column

    # Define aggregation functions as a list
    aggregation_funcs = [
        ("normalized_sum", lambda x: x.sum() / len(x)),                # Normalized sum
        ("mean", "mean"),                                              # Mean
        ("median", "median"),                                          # Median
        ("directional_extreme", lambda x: x.dropna().loc[x.dropna().abs().idxmax()] if not x.dropna().empty else 0),    # Max absolute value preserving sign
        ("normalized_absolute_sum", lambda x: np.abs(x).sum() / len(x))  # Normalized absolute sum
    ]


    # Create a list to store the aggregated results for each feature
    aggregated_results = []

    # Iterate over each column in variant_df (except 'sample_id')
    for column in variant_df.columns[:-1]:  # Exclude the 'sample_id' column
        # Group by sample_id and apply the aggregations
        grouped = variant_df.groupby("SAMPLE")[column]
        aggregated = grouped.agg([func for _, func in aggregation_funcs])
        
        # Rename the columns to include the feature name (e.g., "sum_mean", "sum_median")
        aggregated.columns = [f"{column}_{name}" for name, _ in aggregation_funcs]
        
        # Append the aggregated DataFrame to the list
        aggregated_results.append(aggregated)
        
    # Combine all results into a single DataFrame
    final_aggregates = pd.concat(aggregated_results, axis=1)

    # Reset index to make 'sample_id' a column
    final_aggregates.reset_index(inplace=True)

    # Add 'Enformer_' prefix to all columns except 'sample_id'
    final_aggregates = final_aggregates.rename(
        columns={col: f"Enformer_{col}" for col in final_aggregates.columns if col != "SAMPLE"}
    )

    enformer_df = final_aggregates.copy()


    print('Load Sei', flush=True)
    
    root = zarr.open('/data2/deepLN/sei_result_kor_909_sfari_mssng.zarr', mode='r')
    column_names = root['metadata'].attrs['columns']
    sample_ids = root['metadata'].attrs['samples']
    var_ids = root['metadata'].attrs['variant']
    mat = root['data'][:]

    filtered_mat, filtered_variants, filtered_sample_ids = filter_data_by_variants(mat, var_ids, sample_ids, filter_variant)

    # Step 1: Per-Variant Aggregation (process rows)
    print("Processing Per-Variant Aggregations...")
    variant_aggregates = process_in_chunks_multiprocessing(filtered_mat, compute_variant_aggregations, axis=0, chunk_size=10000, n_cpus=n_cpus)

    # Convert variant_aggregates and sample_ids into a pandas DataFrame for easier grouping
    variant_df = pd.DataFrame(variant_aggregates, columns=["sum", "mean", "median", "directional_extreme"])
    variant_df["SAMPLE"] = filtered_sample_ids  # Add sample_id as a column

    # Create a list to store the aggregated results for each feature
    aggregated_results = []

    # Iterate over each column in variant_df (except 'sample_id')
    for column in variant_df.columns[:-1]:  # Exclude the 'sample_id' column
        # Group by sample_id and apply the aggregations
        grouped = variant_df.groupby("SAMPLE")[column]
        aggregated = grouped.agg([func for _, func in aggregation_funcs])
        
        # Rename the columns to include the feature name (e.g., "sum_mean", "sum_median")
        aggregated.columns = [f"{column}_{name}" for name, _ in aggregation_funcs]
        
        # Append the aggregated DataFrame to the list
        aggregated_results.append(aggregated)
        
    # Combine all results into a single DataFrame
    final_aggregates = pd.concat(aggregated_results, axis=1)

    # Reset index to make 'sample_id' a column
    final_aggregates.reset_index(inplace=True)

    # Add 'Enformer_' prefix to all columns except 'sample_id'
    final_aggregates = final_aggregates.rename(
        columns={col: f"Sei_{col}" for col in final_aggregates.columns if col != "SAMPLE"}
    )

    sei_df = final_aggregates.copy()


    print('Load Puffin', flush=True)

    puffin = pd.read_table('/data2/deepLN/Puffin/puffin_kor_sfari_mssng_output.14606samples.20250224.tsv.gz')

    sample_ids = puffin['SAMPLE'].tolist()
    var_ids = puffin['variant'].tolist()

    puffin = puffin.drop(columns=['SAMPLE', 'variant'])

    mat = puffin.values

    filtered_mat, filtered_variants, filtered_sample_ids = filter_data_by_variants(mat, var_ids, sample_ids, filter_variant)

    # Step 1: Per-Variant Aggregation (process rows)
    print("Processing Per-Variant Aggregations...")
    variant_aggregates = process_in_chunks_multiprocessing(filtered_mat, compute_variant_aggregations, axis=0, chunk_size=10000, n_cpus=n_cpus)

    # Convert variant_aggregates and sample_ids into a pandas DataFrame for easier grouping
    variant_df = pd.DataFrame(variant_aggregates, columns=["sum", "mean", "median", "directional_extreme"])
    variant_df["SAMPLE"] = filtered_sample_ids  # Add sample_id as a column

    # Create a list to store the aggregated results for each feature
    aggregated_results = []

    # Iterate over each column in variant_df (except 'sample_id')
    for column in variant_df.columns[:-1]:  # Exclude the 'sample_id' column
        # Group by sample_id and apply the aggregations
        grouped = variant_df.groupby("SAMPLE")[column]
        aggregated = grouped.agg([func for _, func in aggregation_funcs])
        
        # Rename the columns to include the feature name (e.g., "sum_mean", "sum_median")
        aggregated.columns = [f"{column}_{name}" for name, _ in aggregation_funcs]
        
        # Append the aggregated DataFrame to the list
        aggregated_results.append(aggregated)
        
    # Combine all results into a single DataFrame
    final_aggregates = pd.concat(aggregated_results, axis=1)

    # Reset index to make 'sample_id' a column
    final_aggregates.reset_index(inplace=True)

    # Add 'Enformer_' prefix to all columns except 'sample_id'
    final_aggregates = final_aggregates.rename(
        columns={col: f"Puffin_{col}" for col in final_aggregates.columns if col != "SAMPLE"}
    )

    s0 = pd.read_table(sample_list_path)

    s0 = s0[['SAMPLE', 'isASD']]

    puffin_df = pd.merge(
        final_aggregates,  # Left DataFrame
        s0,                # Right DataFrame
        left_on="SAMPLE",  # Column in final_aggregates
        right_on="SAMPLE",  # Column in s0
        how="inner"         # Inner join
    )

    # # Drop the 'vcf_iid' column
    # puffin_df = puffin_df.drop(columns=["vcf_iid"])

    # Rename 'isASD' to 'label'
    puffin_df = puffin_df.rename(columns={"isASD": "label"})

    # Perform sequential merging of the three DataFrames on 'SAMPLE'
    merged_df = (
        enformer_df
        .merge(sei_df, on='SAMPLE', how='inner')  # Merge enformer_df and sei_df
        .merge(puffin_df, on='SAMPLE', how='inner')  # Merge the result with puffin_df
    )

    # Save to Zarr with the provided arguments
    save_to_zarr(df=merged_df, 
                 output_path=args.output_path, 
                 agg_type='NA', 
                 filter_vcf='NA', 
                 groups_to_keep='NA')

if __name__ == "__main__":
    main()
