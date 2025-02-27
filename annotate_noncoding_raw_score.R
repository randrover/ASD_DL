rm(list = ls())
library(tidyverse)
date = '20250226'

###### Coding
## Cell type gene
gm = read.delim('~/Dropbox/gene-mat-recipe/gene_matrix_20231221.txt')
ct_g = gm %>%
  dplyr::select(1:2, 10:22)

colnames(gm)
## ASD gene set
asd_g = gm %>%
  dplyr::select(1:2, 5:7, 9)

## PTV
dnv = data.table::fread('~/Dropbox/CWAS_outputs/Output_kor_sfari_mssng_20250226/Run_20250226/Kor_SFARI_MSSNG.v5.output.annotation_v8/extract_var/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.extracted_variants.txt.gz')
ptv = dnv %>%
  filter(is_PTVRegion==1)
# dnv$LOEUF37
ptv = ptv %>%
  dplyr::select(SAMPLE, variant = ID) %>%
  mutate(tid = paste(SAMPLE,
                     variant,
                     sep = ':'))

## Missense
dnv2 = data.table::fread('~/Dropbox/CWAS_outputs/Output_kor_sfari_mssng_20250226/Run_20250226/Kor_SFARI_MSSNG.v5.output.annotation_v8/extract_var/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.extracted_variants.txt.gz')
# dnv2$MisDb_MPC
dnv2$MPC = as.numeric(dnv2$MisDb_MPC)
dnv2 = dnv2 %>%
  dplyr::select(variant = ID, SAMPLE, MPC)

mg = merge(dnv %>%
             mutate(var_id = paste(ID, SAMPLE,
                                   sep = ':')),
           dnv2 %>%
             mutate(var_id = paste(variant, SAMPLE,
                                   sep = ':')) %>%
             dplyr::select(-variant, -SAMPLE),
           by.x = 'var_id',
           by.y = 'var_id')

mg = mg %>%
  dplyr::select(variant = ID, SAMPLE, MPC, is_MissenseRegion, is_PromoterRegion, is_UTRsRegion)
miss = mg %>%
  filter(is_MissenseRegion==1)

## 1~2
modMIS = miss %>%
  filter(!is.na(MPC)) %>%
  filter(MPC>=1 & MPC <2)
modMIS = modMIS %>%
  mutate(tid = paste(SAMPLE,
                     variant,
                     sep = ':'))
## 2~
dMIS = miss %>%
  filter(!is.na(MPC)) %>%
  filter(MPC>=2)
dMIS = dMIS %>%
  mutate(tid = paste(SAMPLE,
                     variant,
                     sep = ':'))

## promoter
prom = mg %>%
  filter(is_PromoterRegion == 1) %>%
  mutate(tid = paste(SAMPLE,
                     variant,
                     sep = ':'))

## UTR
UTR = mg %>%
  filter(is_UTRsRegion == 1) %>%
  mutate(tid = paste(SAMPLE,
                     variant,
                     sep = ':'))



## nearest gene
dnv = data.table::fread('~/Dropbox/AI_NDD/Data/variants/Korean_SFARI_MSSNG.DNV_list.all_nearest_transcript.20250226.tsv.gz') %>%
  mutate(tid = paste(SAMPLE, variant, sep = ':'))
mg1 = merge(dnv %>%
              dplyr::select(-gene_name),
            ct_g,
            by = 'gene_id')
mg2 = merge(mg1,
            asd_g %>%
              dplyr::select(-gene_name),
            by = 'gene_id')
mg2$is_PTV = ifelse(mg2$tid %in% ptv$tid,
                    1,
                    0)
mg2$is_modMIS = ifelse(mg2$tid %in% modMIS$tid,
                       1,
                       0)
mg2$is_dMIS = ifelse(mg2$tid %in% dMIS$tid,
                     1,
                     0)
mg2$is_promoter = ifelse(mg2$tid %in% prom$tid,
                         1,
                         0)
mg2$is_UTR = ifelse(mg2$tid %in% UTR$tid,
                    1,
                    0)

mg3 = merge(mg2,
            dnv2 %>%
              mutate(tid = paste(SAMPLE, variant, sep = ':')) %>%
              dplyr::select(-MPC, -SAMPLE, -variant),
            by = 'tid')
mg3 = mg3 %>%
  dplyr::select(SAMPLE, variant,
                gene_id, gene_name,
                12:33)

data.table::fwrite(mg3,
                   paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/table.kor_sfari_mssng.DNV_annotated.coding_combinations', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")


celltype_cols = c('MGE.dev', 'CGE.dev', 'AST', 'L2.3', 'L4', 'L5', 'L5.6.IT', 'L6', 'MG', 'OL', 'END', 'PER', 'SP')
asd_gene_cols = c('DDD285', 'CHD8Common', 'FMRPDarnell', 'ASD185')
varianttype_cols = c('is_PTV', 'is_modMIS', 'is_dMIS')

(length(celltype_cols) + 1) * (length(asd_gene_cols) + 1) * (length(varianttype_cols) + 1)








######## noncoding
rm(list = ls())
date = '20250226'

dnv = data.table::fread('~/Dropbox/AI_NDD/Data/variants/Korean_SFARI_MSSNG.DNV_list.annotated.20250226.tsv.gz')

## spliceAI
spliceai = dnv %>%
  dplyr::select(variant, SAMPLE, is_coding,
                5:8)
spliceai[is.na(spliceai)] <- 0
# Assuming 'spliceai' is a data.table
spliceai[, max_SpliceAI := pmax(SpliceAI_DS_AG, SpliceAI_DS_AL, 
                                SpliceAI_DS_DG, SpliceAI_DS_DL, na.rm = TRUE)]

## conservation
m1 = data.table::fread('~/Dropbox/Noncoding_kor_ASD_WD/Tables/phyloP447way.Korean_CSM.WGS.autosomal_DNV.178samples.sorted.20240426.vcf.gz')
# Use sub() to extract the value after 'ANNOT='
m1$phyloP447way <- sub(".*ANNOT=([^;]*).*", "\\1", m1$INFO)
# Extract value after 'SAMPLE=' and before the next ';'
m1$SAMPLE <- sub(".*SAMPLE=([^;]*).*", "\\1", m1$INFO)
m1 = m1 %>%
  dplyr::select(phyloP447way, variant = ID, SAMPLE)
m2 = data.table::fread('~/Dropbox/Noncoding_kor_ASD_WD/Tables/phastCons46way.Korean_CSM.WGS.autosomal_DNV.178samples.sorted.20240426.vcf.gz')
m2$phastCons46way <- sub(".*ANNOT=([^;]*).*", "\\1", m2$INFO)
m2$SAMPLE <- sub(".*SAMPLE=([^;]*).*", "\\1", m2$INFO)
m2 = m2 %>%
  dplyr::select(phastCons46way, variant = ID, SAMPLE)
m3 = data.table::fread('~/Dropbox/Noncoding_kor_ASD_WD/Tables/phyloP46way.Korean_CSM.WGS.autosomal_DNV.178samples.sorted.20240426.vcf.gz')
m3$phyloP46way <- sub(".*ANNOT=([^;]*).*", "\\1", m3$INFO)
m3$SAMPLE <- sub(".*SAMPLE=([^;]*).*", "\\1", m3$INFO)
m3 = m3 %>%
  dplyr::select(phyloP46way, variant = ID, SAMPLE)

conservation <- Reduce(function(x, y) merge(x, y, by = c("SAMPLE", "variant"), all = TRUE), list(m1, m2, m3))
# Convert all columns except 'SAMPLE' and 'variant' to numeric
numeric_cols <- setdiff(names(conservation), c("SAMPLE", "variant"))  # Select all except these two
conservation[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]  # Convert to numeric


## constraint
m1 = data.table::fread('~/Dropbox/Noncoding_kor_ASD_WD/Tables/ConstraintZ.Korean_CSM.WGS.autosomal_DNV.178samples.sorted.20240426.vcf.gz')
# Use sub() to extract the value after 'ANNOT='
m1$ConstraintZ <- sub(".*ANNOT=([^;]*).*", "\\1", m1$INFO)
# Extract value after 'SAMPLE=' and before the next ';'
m1$SAMPLE <- sub(".*SAMPLE=([^;]*).*", "\\1", m1$INFO)
m1 = m1 %>%
  dplyr::select(ConstraintZ, variant = ID, SAMPLE)
m2 = data.table::fread('~/Dropbox/Noncoding_kor_ASD_WD/Tables/JARVIS.Korean_CSM.WGS.autosomal_DNV.178samples.sorted.20240426.vcf.gz')
m2$JARVIS <- sub(".*ANNOT=([^;]*).*", "\\1", m2$INFO)
m2$SAMPLE <- sub(".*SAMPLE=([^;]*).*", "\\1", m2$INFO)
m2 = m2 %>%
  dplyr::select(JARVIS, variant = ID, SAMPLE)

constraint = merge(m1,
                   m2,
                   by = c('SAMPLE', 'variant'))
# Convert all columns except 'SAMPLE' and 'variant' to numeric
numeric_cols <- setdiff(names(constraint), c("SAMPLE", "variant"))  # Select all except these two
constraint[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]  # Convert to numeric


merged_df <- merge(conservation, constraint, by = c("SAMPLE", "variant"), all = FALSE)


dnv2 = data.table::fread('~/Dropbox/AI_NDD/Data/variants/Korean_SFARI_MSSNG.DNV_list.annotated.20241020.tsv.gz')
dnv2 = dnv2 %>%
  mutate(phyloP46way = phyloP46wayVt,
         phyloP447way = phyloP447wayPrimates,
         phastCons46way = phastCons46wayVt)
dnv3 = dnv2 %>%
  dplyr::select(colnames(merged_df))

merged_df2 = rbind(dnv3,
                   merged_df)

merged_df3 <- merge(merged_df2, spliceai, by = c("SAMPLE", "variant"), all = FALSE)


data.table::fwrite(merged_df3,
                   paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/table.kor_sfari_mssng.DNV_annotated.noncoding_raw_score', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")

