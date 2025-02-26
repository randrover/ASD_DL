# Define the paths to your variables
FILTER_VARIANT_PATH="/data2/deepLN/variants/kor_sfari_mssng.noncoding_DNV.20250226.tsv.gz"
OUTPUT_PATH="/data2/deepLN/kor_sfari_mssng.enformer_sei_puffin.new_agg.v1.zarr"
SAMPLE_LIST_PATH="/data2/deepLN/Korean_SFARI_MSSNG.14606samples.sample_list.20250224.txt"
NC_CPUS=10

# Step 1: Run sei_puffin_enformer_agg_filter.py
echo "Running sei_puffin_enformer_agg_filter.py..."
python sei_puffin_enformer_agg_filter.py \
  --filter_variant_path "$FILTER_VARIANT_PATH" \
  --output_path "$OUTPUT_PATH" \
  --sample_list_path "$SAMPLE_LIST_PATH" \
  --n_cpus "$NC_CPUS"
