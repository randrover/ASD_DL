import pandas as pd
import numpy as np
import zarr
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process and merge Zarr data and filter variants.")
    
    # Define arguments
    parser.add_argument('--zarr_input', type=str, help="Path to the input Zarr file")
    parser.add_argument('--zarr_output', type=str, help="Path to save the output Zarr file")
    parser.add_argument('--filter_variant', type=str, help="Path to the filter_variant file")
    
    return parser.parse_args()


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


def main():
    args = parse_args()
    
    def max_abs_value(series):
        return series.loc[series.abs().idxmax()] if not series.isna().all() else np.nan
    
    root = zarr.open(args.zarr_input, mode='r')
    mat = root['data']
    column_names = root['metadata'].attrs['columns']
    sample_ids = root['metadata'].attrs['sample_ids']

    # Create a DataFrame
    enformer_sei_puffin_df = pd.DataFrame(mat, columns=column_names, index=sample_ids)

    # Reset index to make SAMPLE a column
    enformer_sei_puffin_df.reset_index(inplace=True)
    enformer_sei_puffin_df.rename(columns={'index': 'SAMPLE'}, inplace=True)

    motif = pd.read_table('table.kor_sfari_mssng_14606.motif_diff_by_var.20241225.tsv.gz')

    filter_variant = pd.read_table(args.filter_variant)

    dnv = pd.read_table('table.kor_sfari_mssng.DNV_annotated.noncoding_raw_score.20250226.tsv.gz')

    dnv = dnv[['is_coding', 'variant', 'SAMPLE']]

    filter_variant_list = filter_variant['variant'].tolist()

    # Create a mask for is_coding == 0 and apply the filtering
    coding_0_mask = (dnv['is_coding'] == 0) & (dnv['variant'].isin(filter_variant_list))

    # Create the final DataFrame by combining the results
    dnv_filtered = pd.concat([
        dnv[coding_0_mask],  # Keep the rows where is_coding == 0 and variant is in the list
        dnv[dnv['is_coding'] == 1]  # Keep all rows where is_coding == 1
    ])

    motif_df = pd.merge(motif, dnv_filtered, on=['variant', 'SAMPLE'], how='inner')

    # List of columns to compute max values for (excluding 'variant', 'SAMPLE', and 'is_coding')
    columns_to_max = [col for col in motif_df.columns if col not in ['variant', 'SAMPLE', 'is_coding']]

    # (1) Max values for all rows by 'SAMPLE'
    max_all = motif_df.groupby('SAMPLE')[columns_to_max].apply(lambda x: x.apply(max_abs_value))
    max_all.columns = ['all_' + col for col in max_all.columns]

    # (2) Max values for rows where is_coding == 1
    max_is_coding_1 = motif_df[motif_df['is_coding'] == 1].groupby('SAMPLE')[columns_to_max].apply(lambda x: x.apply(max_abs_value))
    max_is_coding_1.columns = ['coding_' + col for col in max_is_coding_1.columns]

    # (3) Max values for rows where is_coding == 0
    max_is_coding_0 = motif_df[motif_df['is_coding'] == 0].groupby('SAMPLE')[columns_to_max].apply(lambda x: x.apply(max_abs_value))
    max_is_coding_0.columns = ['noncoding_' + col for col in max_is_coding_0.columns]

    # Merge the results into a single dataframe
    merged_df = pd.merge(max_all, max_is_coding_1, on='SAMPLE', how='outer')
    motif_merge_df = pd.merge(merged_df, max_is_coding_0, on='SAMPLE', how='outer')

    # Merge the three dataframes on 'SAMPLE'
    merged_df = pd.merge(motif_merge_df, enformer_sei_puffin_df, on='SAMPLE', how='inner')


    save_to_zarr(df = merged_df, output_path = args.zarr_output, agg_type = 'NA', filter_vcf = 'NA', groups_to_keep = 'NA')

if __name__ == "__main__":
    main()

