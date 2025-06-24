rm(list = ls())
library(tidyverse)
date = '20250226'

dt = data.table::fread('~/Dropbox/Noncoding_kor_ASD_WD/Tables/table.kor_sfari_mssng.DNV_annotated.noncoding_raw_score.20250226.tsv.gz')
is_coding=F

if(is_coding==T){
  df = dt %>%
    filter(is_coding==1) %>%
    dplyr::select(-is_coding)
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phastCons46way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/coding.kor_sfari_mssng.phastCons46wayVt', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phyloP46way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/coding.kor_sfari_mssng.phyloP46wayVt', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phyloP447way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/coding.kor_sfari_mssng.phyloP447wayPrimates', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, ConstraintZ)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/coding.kor_sfari_mssng.ConstraintZ', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, JARVIS)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/coding.kor_sfari_mssng.JARVIS', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, SpliceAI_DS_AG, SpliceAI_DS_AL,
                  SpliceAI_DS_DG, SpliceAI_DS_DL, max_SpliceAI)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/coding.kor_sfari_mssng.SpliceAI', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
}else{
  df = dt %>%
    filter(is_coding==0) %>%
    dplyr::select(-is_coding)
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phastCons46way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/noncoding.kor_sfari_mssng.phastCons46wayVt', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phyloP46way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/noncoding.kor_sfari_mssng.phyloP46wayVt', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phyloP447way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/noncoding.kor_sfari_mssng.phyloP447wayPrimates', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, ConstraintZ)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/noncoding.kor_sfari_mssng.ConstraintZ', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, JARVIS)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/noncoding.kor_sfari_mssng.JARVIS', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, SpliceAI_DS_AG, SpliceAI_DS_AL,
                  SpliceAI_DS_DG, SpliceAI_DS_DL, max_SpliceAI)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/noncoding.kor_sfari_mssng.SpliceAI', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
}


## total
df = dt

if(T){
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phastCons46way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/kor_sfari_mssng.phastCons46wayVt', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phyloP46way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/kor_sfari_mssng.phyloP46wayVt', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, phyloP447way)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/kor_sfari_mssng.phyloP447wayPrimates', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, ConstraintZ)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/kor_sfari_mssng.ConstraintZ', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, JARVIS)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/kor_sfari_mssng.JARVIS', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
  tmp = df %>%
    dplyr::select(variant, SAMPLE, SpliceAI_DS_AG, SpliceAI_DS_AL,
                  SpliceAI_DS_DG, SpliceAI_DS_DL, max_SpliceAI)
  
  data.table::fwrite(tmp,
                     paste('~/Dropbox/Noncoding_kor_ASD_WD/Tables/kor_sfari_mssng.SpliceAI', date, 'tsv.gz',
                           sep = '.'),
                     quote = F, row.names = F, col.names = T, sep = '\t',
                     compress="gzip")
  
}
