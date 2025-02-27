import argparse
import pandas as pd
import zarr

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

# Function to parse arguments from the terminal
def parse_args():
    parser = argparse.ArgumentParser(description="Process and merge Zarr data and filter variants.")
    
    # Define required arguments
    parser.add_argument('--zarr_input', type=str, required=True, help="Path to the input Zarr file")
    parser.add_argument('--zarr_output', type=str, required=True, help="Path to save the output Zarr file")
    parser.add_argument('--coding_combinations', type=str, required=True, help="Path to the coding_combinations file")
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    zarr_path = args.zarr_input

    root = zarr.open(zarr_path, mode='r')
    mat = root['data']
    column_names = root['metadata'].attrs['columns']
    sample_ids = root['metadata'].attrs['sample_ids']

    # Create a DataFrame
    df = pd.DataFrame(mat, columns=column_names, index=sample_ids)

    # Reset index to make SAMPLE a column
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'SAMPLE'}, inplace=True)

    s0 = pd.read_table(args.coding_combinations)

    merged_df = pd.merge(s0, df, on='SAMPLE', how='inner')

    save_to_zarr(df = merged_df, output_path=args.zarr_output, agg_type = 'NA', filter_vcf = 'NA', groups_to_keep = 'NA')

if __name__ == "__main__":
    main()
