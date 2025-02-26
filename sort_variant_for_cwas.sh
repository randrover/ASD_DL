grep -a "^#" Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.20250226.vcf > Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf
grep -v "^#" Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.20250226.vcf | sort -k1,1V -k2n >> Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf

bgzip -c Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf > Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf.gz
tabix -p vcf Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf.gz
