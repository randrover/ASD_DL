#!/bin/bash

# Define the paths to your variables
FILTER_VARIANT_PATH="/data2/deepLN/variants/kor_sfari_mssng.noncoding_DNV.20250226.tsv.gz"
OUTPUT_PATH="/data2/deepLN/kor_sfari_mssng.enformer_sei_puffin.new_agg.v1.zarr"
OUTPUT_PATH2="/data2/deepLN/kor_sfari_mssng.noncoding_raw.new_agg.v1.zarr"
SAMPLE_LIST_PATH="/data2/deepLN/Korean_SFARI_MSSNG.14606samples.sample_list.20250224.txt"
NC_CPUS=10

# Step 1: Run sei_puffin_enformer_agg_filter.py
echo "Running sei_puffin_enformer_agg_filter.py..."
python sei_puffin_enformer_agg_filter.py \
  --filter_variant_path "$FILTER_VARIANT_PATH" \
  --output_path "$OUTPUT_PATH" \
  --sample_list_path "$SAMPLE_LIST_PATH" \
  --n_cpus "$NC_CPUS"

# Step 2: Run noncoding_annot_gather.py
echo "Running noncoding_annot_gather.py..."
python noncoding_annot_gather.py \
  --zarr_input "$OUTPUT_PATH" \
  --zarr_output "$OUTPUT_PATH2" \
  --filter_variant "$FILTER_VARIANT_PATH"

# Step 3: Run coding_noncoding_gather.py
echo "Running coding_noncoding_gather.py..."
python coding_noncoding_gather.py \
  --zarr_input "$OUTPUT_PATH2" \
  --zarr_output "/data2/deepLN/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v31_tmp.zarr" \
  --coding_combinations "/data2/deepLN/kor_sfari_mssng.coding_combinations.agg_by_sample.20250227.tsv.gz"

echo "All scripts completed successfully!"

FILTER_VARIANT_PATH="/data2/deepLN/variants/kor_sfari_mssng.noncoding_DNV.20250226.tsv.gz"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/noncoding.kor_sfari_mssng.ConstraintZ.20250226.tsv.gz \
-cutoffs 1 2 3 4 \
-output_path /data2/deepLN/annot_v7/noncoding_DNV.kor_sfari_mssng.ConstraintZ.20250226.cutoff_applied.tsv.gz \
-prefix "nc" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/noncoding.kor_sfari_mssng.JARVIS.20250226.tsv.gz \
-cutoffs 0.9 0.99 0.999 \
-output_path /data2/deepLN/annot_v7/noncoding_DNV.kor_sfari_mssng.JARVIS.20250226.cutoff_applied.tsv.gz \
-prefix "nc" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/noncoding.kor_sfari_mssng.SpliceAI.20250226.tsv.gz \
-cutoffs 0.2 0.5 0.7 \
-output_path /data2/deepLN/annot_v7/noncoding_DNV.kor_sfari_mssng.SpliceAI.20250226.cutoff_applied.tsv.gz \
-prefix "nc" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/noncoding.kor_sfari_mssng.phyloP447wayPrimates.20250226.tsv.gz \
-cutoffs 1 1.3 1.5 1.7 \
-output_path /data2/deepLN/annot_v7/noncoding_DNV.kor_sfari_mssng.phyloP447wayPrimates.20250226.cutoff_applied.tsv.gz \
-prefix "nc" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/noncoding.kor_sfari_mssng.phyloP46wayVt.20250226.tsv.gz \
-cutoffs 2 4 6 \
-output_path /data2/deepLN/annot_v7/noncoding_DNV.kor_sfari_mssng.phyloP46wayVt.20250226.cutoff_applied.tsv.gz \
-prefix "nc" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/noncoding.kor_sfari_mssng.phastCons46wayVt.20250226.tsv.gz \
-cutoffs 0.2 0.4 0.6 0.8 \
-output_path /data2/deepLN/annot_v7/noncoding_DNV.kor_sfari_mssng.phastCons46wayVt.20250226.cutoff_applied.tsv.gz \
-prefix "nc" \
-filter_variant_file "$FILTER_VARIANT_PATH"

FILTER_VARIANT_PATH="/data2/deepLN/variants/kor_sfari_mssng.coding_DNV.20250227.tsv.gz"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/coding.kor_sfari_mssng.ConstraintZ.20250226.tsv.gz \
-cutoffs 1 2 3 4 \
-output_path /data2/deepLN/annot_v7/coding_DNV.kor_sfari_mssng.ConstraintZ.20250226.cutoff_applied.tsv.gz \
-prefix "cd" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/coding.kor_sfari_mssng.JARVIS.20250226.tsv.gz \
-cutoffs 0.9 0.99 0.999 \
-output_path /data2/deepLN/annot_v7/coding_DNV.kor_sfari_mssng.JARVIS.20250226.cutoff_applied.tsv.gz \
-prefix "cd" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/coding.kor_sfari_mssng.SpliceAI.20250226.tsv.gz \
-cutoffs 0.2 0.5 0.7 \
-output_path /data2/deepLN/annot_v7/coding_DNV.kor_sfari_mssng.SpliceAI.20250226.cutoff_applied.tsv.gz \
-prefix "cd" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/coding.kor_sfari_mssng.phyloP447wayPrimates.20250226.tsv.gz \
-cutoffs 1 1.3 1.5 1.7 \
-output_path /data2/deepLN/annot_v7/coding_DNV.kor_sfari_mssng.phyloP447wayPrimates.20250226.cutoff_applied.tsv.gz \
-prefix "cd" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/coding.kor_sfari_mssng.phyloP46wayVt.20250226.tsv.gz \
-cutoffs 2 4 6 \
-output_path /data2/deepLN/annot_v7/coding_DNV.kor_sfari_mssng.phyloP46wayVt.20250226.cutoff_applied.tsv.gz \
-prefix "cd" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/coding.kor_sfari_mssng.phastCons46wayVt.20250226.tsv.gz \
-cutoffs 0.2 0.4 0.6 0.8 \
-output_path /data2/deepLN/annot_v7/coding_DNV.kor_sfari_mssng.phastCons46wayVt.20250226.cutoff_applied.tsv.gz \
-prefix "cd" \
-filter_variant_file "$FILTER_VARIANT_PATH"

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/kor_sfari_mssng.ConstraintZ.20250226.tsv.gz \
-cutoffs 1 2 3 4 \
-output_path /data2/deepLN/annot_v7/kor_sfari_mssng.ConstraintZ.20250226.cutoff_applied.tsv.gz

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/kor_sfari_mssng.JARVIS.20250226.tsv.gz \
-cutoffs 0.9 0.99 0.999 \
-output_path /data2/deepLN/annot_v7/kor_sfari_mssng.JARVIS.20250226.cutoff_applied.tsv.gz

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/kor_sfari_mssng.SpliceAI.20250226.tsv.gz \
-cutoffs 0.2 0.5 0.7 \
-output_path /data2/deepLN/annot_v7/kor_sfari_mssng.SpliceAI.20250226.cutoff_applied.tsv.gz

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/kor_sfari_mssng.phyloP447wayPrimates.20250226.tsv.gz \
-cutoffs 1 1.3 1.5 1.7 \
-output_path /data2/deepLN/annot_v7/kor_sfari_mssng.phyloP447wayPrimates.20250226.cutoff_applied.tsv.gz

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/kor_sfari_mssng.phyloP46wayVt.20250226.tsv.gz \
-cutoffs 2 4 6 \
-output_path /data2/deepLN/annot_v7/kor_sfari_mssng.phyloP46wayVt.20250226.cutoff_applied.tsv.gz

python /data2/deepLN/organize_cutoff_new_df.filter.py -file_path /data2/deepLN/annot_v7/kor_sfari_mssng.phastCons46wayVt.20250226.tsv.gz \
-cutoffs 0.2 0.4 0.6 0.8 \
-output_path /data2/deepLN/annot_v7/kor_sfari_mssng.phastCons46wayVt.20250226.cutoff_applied.tsv.gz

python merge_by_sample.ec2.py \
-file_with_label /data2/deepLN/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v31_tmp.zarr \
-other_files /data2/deepLN/annot_v7/by_sample/* \
-output_file /data2/deepLN/kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v31.zarra
