# ASD_ML
ASD research with machine learning techniques

## 1. System Requirements
- Operating System: Ubuntu 24.04.2 LTS
- Python version: ≥ 3.9
- CUDA version: CUDA 12.6 (for GPU acceleration, optional)
- Hardware: CPU-only compatible; GPU (NVIDIA CUDA) recommended for large-scale training

## 2. Installation guide

## 3. Demo

## 4. Instructions for Use on Custom Data
- Prepare input data in Zarr format with the same structure as example_input/demo_dataset.zarr
- Optionally provide:
  - Feature group file (e.g., feature_info.txt)
  - Sample metadata (sample_meta.txt) for stratified splitting
- Run training using train_XL.py with desired parameters:

```
python train_7L.py \
  --zarr_path your_dataset.zarr \
  --log_file training_log.csv \
  --calculate_importance \
  --importance_method shap \
  --gpu
```

## 5. Instruction for reproduction of the manuscript

```
python train_7L.py \
  --zarr_path kor_sfari_mssng.coding_comb.noncoding_comb.new_agg.v27.modelA_train.zarr \
  --model_dir model \
  --wandb_key your_wandb_api_key \
  --log_file training_log.csv \
  --feature_info_path table.shap.feature_info_v28.20250121.txt \
  --group_list Sig03 Sig02 Sig01 \
  --gpu \
  --oversampling
```
