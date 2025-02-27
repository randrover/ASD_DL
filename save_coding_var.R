rm(list = ls())
library(tidyverse)
date = '20250227'

dnv = data.table::fread('~/Dropbox/AI_NDD/Data/variants/Korean_SFARI_MSSNG.DNV_list.annotated.20250226.tsv.gz')
dnv = dnv %>%
  filter(is_coding==1)
data.table::fwrite(dnv %>%
                     dplyr::select(variant, SAMPLE),
                   paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/kor_sfari_mssng.coding_DNV', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")
