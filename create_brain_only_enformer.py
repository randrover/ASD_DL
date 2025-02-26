import zarr
import numpy as np
import pandas as pd

print('Load Enformer', flush=True)

enformer_path = '/data1/deepLN/Enformer/enformer_result_kor_909_sfari_mssng.zarr'

root = zarr.open(enformer_path, mode='r')
mat = root['data']
column_names1 = root['metadata'].attrs['columns']
var_ids = root['metadata'].attrs['variant']
samples = root['metadata'].attrs['samples']

print(mat.shape, flush=True)

# Read the Excel file
df = pd.read_excel("/data1/deepLN/Enformer/Enformer_training_features_Brain.xlsx", engine="openpyxl")

# Display the first few rows
print(df.head())

# 먼저 df["index"] 값들을 리스트로 저장
valid_indices = df["index"].astype(int).values  # 혹시라도 숫자가 문자열일 경우 변환

import numpy as np

# NumPy 배열로 변환 후 정렬
valid_indices_sorted = np.sort(valid_indices)

# valid_indices가 column_names1 리스트의 인덱스 범위 안에 있는지 확인
filtered_column_names1 = [column_names1[idx] for idx in valid_indices_sorted if 0 <= idx < len(column_names1)]

# 결과 출력
print(len(filtered_column_names1))

# valid_indices_sorted를 기준으로 data의 해당 컬럼들만 선택
mat_filtered = mat[:, valid_indices_sorted]

print("Filtered", flush=True)

print(mat_filtered.shape)

output_path = '/data1/deepLN/Enformer/enformer_result_kor_909_sfari_mssng_only_brain386.zarr'
root = zarr.open(output_path, mode='w')

# Create metadata group and add attributes
root.create_group('metadata')
root['metadata'].attrs['variant'] = var_ids
root['metadata'].attrs['samples'] = samples
root['metadata'].attrs['columns'] = filtered_column_names1

# Ensure data is in float format and create dataset
root.create_dataset('data', data=mat_filtered, chunks=(1000, 1000), dtype='float')
print(f'Data saved to {output_path}', flush=True)



