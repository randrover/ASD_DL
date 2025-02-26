#!/bin/sh
export R_HOME=/usr/lib/R
echo $(date)

input_dir=/data3/cwas_input/Kor_SFARI_MSSNG_20250226
output_dir=/data3/cwas_output/Kor_SFARI_MSSNG.v5.output.annotation_v8

cwas annotation -v $input_dir/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf.gz \
-o_dir $output_dir/annotation -p 40

cwas categorization -i $output_dir/annotation/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.annotated.vcf.gz \
-o_dir $output_dir/categorization/ \
-p 40

cwas binomial_test \
-i $output_dir/categorization/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.categorization_result.zarr \
-o_dir $output_dir/burden_test/ \
-s $input_dir/Korean_SFARI_MSSNG14606samples.sample_list.20250226.txt \
-a $input_dir/adjustFile_patAge_cohort.Korean_SFARI_MSSNG.WGS.autosomal_DNV.14606samples.20250226.txt

cwas extract_variant -i $output_dir/annotation/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.annotated.vcf.gz \
-o_dir $output_dir/extract_var -ai

echo $(date)
echo "Done"
