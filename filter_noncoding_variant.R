rm(list = ls())
library(tidyverse)
date = '20250226'

# Also VISTA v1

ss = read.delim('~/Dropbox/Noncoding_kor_ASD_WD/Data/Korean_SFARI_MSSNG.14606samples.sample_list.20250224.txt')

df = data.table::fread('~/Dropbox/CWAS_outputs/Output_kor_sfari_mssng_20250226/Run_20250226/Kor_SFARI_MSSNG.v5.output.annotation_v8/extract_var/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.extracted_variants.txt.gz')
df = df %>%
  filter(SAMPLE %in% ss$SAMPLE)

## noncoding variant
vcf1 = df %>%
  filter(is_NoncodingRegion==1)

data.table::fwrite(vcf1 %>%
                     dplyr::select(variant = ID,
                                   SAMPLE),
                   paste('~/Dropbox/Noncoding_kor_ASD_WD/Data/kor_sfari_mssng.noncoding_DNV', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")

df1 = df %>%
  filter(is_NoncodingRegion==1)

df2 = df1 %>% dplyr::select('ID', 'SAMPLE', 'PHENOTYPE', 'COHORT', 69:ncol(df1))
# df2 = df2 %>%
#   dplyr::select('ID', 'SAMPLE', 'PHENOTYPE', 'COHORT', 15:17, 47:62, 67:77, 78:97)
## Fetal
# Excitatory neuron
ex_col = c('DomckeCRE.ExN', 'HerringCRE.Fetal.L23', 'HerringCRE.Fetal.L4', 'HerringCRE.Fetal.L56', 'MannensCRE.GluN1', 'MannensCRE.GluN2', 'MannensCRE.GluN3', 'TrevinoCRE.EarlyGluN', 'TrevinoCRE.LateGluN', 'ZiffraCRE.dlEN', 'ZiffraCRE.earlyEN', 'ZiffraCRE.ulEN')
m1 = df2 %>%
  dplyr::select(all_of(ex_col))
df2$ExN = rowSums(m1) > 0

# Inhibitory neuron
inh_col = c('DomckeCRE.InN', 'HerringCRE.Fetal.CGE', 'HerringCRE.Fetal.MGE', 'MannensCRE.GABA', 'MannensCRE.LGEDGE', 'MannensCRE.MGE', 'TrevinoCRE.CGE', 'TrevinoCRE.MGE', 'ZiffraCRE.CGE', 'ZiffraCRE.MGE', 'pRE.CGE', 'pRE.MGE')
m2 = df2 %>%
  dplyr::select(all_of(inh_col))
df2$InN = rowSums(m2) > 0

# Oligodendrocyte
oligo_col = c('TrevinoCRE.Oligo', 'ZiffraCRE.AstroOligo')
m3 = df2 %>%
  dplyr::select(all_of(oligo_col))
df2$Oligo = rowSums(m3) > 0

# Microglia
micro_col = c('ZiffraCRE.Micro')
m4 = df2 %>%
  dplyr::select(all_of(micro_col))
df2$Micro = rowSums(m4) > 0

# Astrocyte
ast_col = c('DomckeCRE.Astro', 'HerringCRE.Fetal.Astro', 'ZiffraCRE.AstroOligo')
## ZiffraCRE.AstroOligo ?
m5 = df2 %>%
  dplyr::select(all_of(ast_col))
df2$Ast = rowSums(m5) > 0

## VISTA enhancer
vista_col = c("VISTA.Enh")
m6 = df2 %>%
  dplyr::select(all_of(vista_col))
df2$VISTA = rowSums(m6) > 0


df3 = df2 %>%
  dplyr::select('ID', 'SAMPLE', 'PHENOTYPE', 'COHORT', 75:ncol(df2))
df3 = as.data.frame(df3)

## Noncoding Regulatory DNV
r_idx = rowSums(df3[,5:10])>0
vcf2 = df3[r_idx,]
data.table::fwrite(vcf2 %>%
                     dplyr::select(variant = ID,
                                   SAMPLE),
                   paste('~/Dropbox/Noncoding_kor_ASD_WD/Data/kor_sfari_mssng.noncoding_reg_DNV', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")


## Noncoding Regulatory ExN DNV
r_idx = (df3$ExN==TRUE)
vcf3 = df3[r_idx,]
data.table::fwrite(vcf3 %>%
                     dplyr::select(variant = ID,
                                   SAMPLE),
                   paste('~/Dropbox/Noncoding_kor_ASD_WD/Data/kor_sfari_mssng.noncoding_reg_exN_DNV', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")

## Noncoding Regulatory InN DNV
r_idx = (df3$InN==TRUE)
vcf4 = df3[r_idx,]
data.table::fwrite(vcf4 %>%
                     dplyr::select(variant = ID,
                                   SAMPLE),
                   paste('~/Dropbox/Noncoding_kor_ASD_WD/Data/kor_sfari_mssng.noncoding_reg_inN_DNV', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")


## Noncoding Regulatory InN DNV
r_idx = rowSums(df3[,5:6]) > 0
vcf5 = df3[r_idx,]
data.table::fwrite(vcf5 %>%
                     dplyr::select(variant = ID,
                                   SAMPLE),
                   paste('~/Dropbox/Noncoding_kor_ASD_WD/Data/kor_sfari_mssng.noncoding_reg_neuron_DNV', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")

