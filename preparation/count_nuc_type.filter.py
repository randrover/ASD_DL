import pandas as pd
import pysam
from collections import defaultdict
import argparse


pysam.set_verbosity(0)

# Function to process the VCF and count trinucleotide patterns for SNVs
def process_vcf_and_count_trinucs(vcf_file, fasta_file, filter_variants):
    fasta = pysam.FastaFile(fasta_file)
    vcf = pysam.VariantFile(vcf_file)
    sample_trinuc_counts = defaultdict(lambda: defaultdict(int))

    for record in vcf:
        ref = record.ref
        alt = record.alts[0]  # Assume only one alternate allele

        variant_id = f"{record.chrom}:{record.pos}:{ref}:{alt}"

        # If the variant is not in the filter list, skip it
        if filter_variants and variant_id not in filter_variants:
            continue

        if len(ref) != 1 or len(alt) != 1:
            continue

        chrom = record.chrom
        pos = record.pos
        start = max(0, pos - 2)
        end = pos + 1
        surrounding_ref = fasta.fetch(chrom, start, end).upper()

        mutated_trinuc = surrounding_ref[:1] + alt.upper() + surrounding_ref[2:]
        sample = record.info.get('SAMPLE', 'Unknown')
        sample_trinuc_counts[sample][f"{surrounding_ref}>{mutated_trinuc}"] += 1

    unique_samples = list(sample_trinuc_counts.keys())
    all_trinucleotides = list(set([trinuc for counts in sample_trinuc_counts.values() for trinuc in counts.keys()]))
    df_trinucs = pd.DataFrame(0, index=unique_samples, columns=all_trinucleotides)

    for sample, trinuc_dict in sample_trinuc_counts.items():
        for trinuc, count in trinuc_dict.items():
            df_trinucs.at[sample, trinuc] = count

    # Calculate proportions
    df_trinucs = df_trinucs.astype(float)  # Ensure entire DataFrame is float
    for sample in df_trinucs.index:
        total_count = df_trinucs.loc[sample].sum()
        if total_count > 0:
            df_trinucs.loc[sample] = df_trinucs.loc[sample] / total_count

    return df_trinucs

# Function to process the VCF and count single-nucleotide patterns for SNVs
def process_vcf_and_count_snvs(vcf_file, fasta_file, filter_variants):
    fasta = pysam.FastaFile(fasta_file)
    vcf = pysam.VariantFile(vcf_file)
    sample_snv_counts = defaultdict(lambda: defaultdict(int))

    for record in vcf:
        ref = record.ref
        alt = record.alts[0]
        
        variant_id = f"{record.chrom}:{record.pos}:{ref}:{alt}"
        
        # If the variant is not in the filter list, skip it
        if filter_variants and variant_id not in filter_variants:
            continue

        if len(ref) != 1 or len(alt) != 1:
            continue

        snv = f"{ref}>{alt}"
        sample = record.info.get('SAMPLE', 'Unknown')
        sample_snv_counts[sample][snv] += 1

    unique_samples = list(sample_snv_counts.keys())
    all_snvs = list(set([snv for counts in sample_snv_counts.values() for snv in counts.keys()]))
    df_snvs = pd.DataFrame(0, index=unique_samples, columns=all_snvs)

    for sample, snv_dict in sample_snv_counts.items():
        for snv, count in snv_dict.items():
            df_snvs.at[sample, snv] = count

    df_snvs = df_snvs.astype(float)  # Ensure entire DataFrame is float
    # Calculate proportions
    for sample in df_snvs.index:
        total_count = df_snvs.loc[sample].sum()
        if total_count > 0:
            df_snvs.loc[sample] = df_snvs.loc[sample] / total_count

    return df_snvs

# Main function to parse arguments and run the analysis
def main():
    parser = argparse.ArgumentParser(description="Count patterns for single-nucleotide variants (SNVs) in a VCF file.")
    parser.add_argument("--vcf_file", help="Path to the VCF file")
    parser.add_argument("--fasta_file", help="Path to the reference FASTA file")
    parser.add_argument("--output_file", help="Path to the output TSV file")
    parser.add_argument("--prefix", help="Prefix to add to column names", default="")
    parser.add_argument("--filter_variant_file", help="Path to the file containing variants to filter by")

    # Add mutually exclusive group for --single and --trin
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--single", action="store_true", help="Only count single-nucleotide patterns for SNVs (e.g., A>T, C>G)")
    group.add_argument("--trin", action="store_true", help="Only count trinucleotide patterns for SNVs (e.g., ACA>TCA)")

    args = parser.parse_args()

    # Read filter variants if provided
    if args.filter_variant_file:
        print(f"Reading filter variants from {args.filter_variant_file}", flush=True)
        filter_variants = pd.read_table(args.filter_variant_file)
        filter_variants = set(filter_variants)

    # Determine which DataFrames to create based on arguments
    df_trinucs = None
    df_snvs = None
    
    if args.single:
        print('Count single-nucleotide pattern', flush=True)
        df_snvs = process_vcf_and_count_snvs(args.vcf_file, args.fasta_file, filter_variants).reset_index(names='SAMPLE')
    elif args.trin:
        print('Count trinucleotide pattern', flush=True)
        df_trinucs = process_vcf_and_count_trinucs(args.vcf_file, args.fasta_file, filter_variants).reset_index(names='SAMPLE')
    else:  # Default: count both if neither is specified
        print('Count trinucleotide pattern', flush=True)
        df_trinucs = process_vcf_and_count_trinucs(args.vcf_file, args.fasta_file, filter_variants).reset_index(names='SAMPLE')
        print('Count single-nucleotide pattern', flush=True)
        df_snvs = process_vcf_and_count_snvs(args.vcf_file, args.fasta_file, filter_variants).reset_index(names='SAMPLE')

    # Merge and save DataFrames based on specified arguments
    if df_trinucs is not None and df_snvs is not None:
        merged_df = pd.merge(df_trinucs, df_snvs, on='SAMPLE', how='outer')
    elif df_trinucs is not None:
        merged_df = df_trinucs
    elif df_snvs is not None:
        merged_df = df_snvs
    else:
        print("Error: No data to save. Check your arguments.")
        return

    for col in merged_df.columns:
        if col != 'SAMPLE':
            merged_df.rename(columns={col: f"prop_{col}"}, inplace=True)

    # Apply the prefix to the column names (except for 'SAMPLE')
    for col in merged_df.columns:
        if col != 'SAMPLE':
            merged_df.rename(columns={col: f"{args.prefix}_{col}"}, inplace=True)

    # Save the merged DataFrame to a TSV file
    merged_df.to_csv(args.output_file, sep='\t', index=False)
    print(f"Merged DataFrame saved to {args.output_file}")

if __name__ == "__main__":
    main()
