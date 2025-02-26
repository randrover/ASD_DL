library(tidyverse)
library(data.table)
rm(list = ls())
date = '20250226'

d2 = data.table::fread('~/Dropbox/Noncoding_kor_ASD_WD/Data/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf.gz')
d2$`#CHROM` = gsub(x = d2$`#CHROM`, pattern = 'chr', replacement = '', fixed = T)
d2$`#CHROM` = as.numeric(d2$`#CHROM`)

s0 = read.delim('~/Dropbox/Noncoding_kor_ASD_WD/Data/Korean_SFARI_MSSNG.14606samples.sample_list.20250224.txt')

d2 = d2 %>%
  dplyr::select(variant = ID,
                INFO,
                CHROM = `#CHROM`,
                POS)
d2$SAMPLE = do.call(rbind.data.frame, strsplit(x = d2$INFO, split = ';', fixed = T))[[1]]
d2$PHENOTYPE = do.call(rbind.data.frame, strsplit(x = d2$INFO, split = ';', fixed = T))[[2]]
d2$COHORT = do.call(rbind.data.frame, strsplit(x = d2$INFO, split = ';', fixed = T))[[3]]
d2$INFO = NULL
d2$SAMPLE = gsub(x = d2$SAMPLE, pattern = 'SAMPLE=', replacement = '', fixed = T)
d2$PHENOTYPE = gsub(x = d2$PHENOTYPE, pattern = 'PHENOTYPE=', replacement = '', fixed = T)
d2$COHORT = gsub(x = d2$COHORT, pattern = 'COHORT=', replacement = '', fixed = T)
d2 = d2 %>%
  mutate(Group = ifelse(PHENOTYPE=='case',
                        'Autism',
                        ifelse(PHENOTYPE=='ctrl',
                               'Control',
                               'Other'))) %>%
  dplyr::select(-PHENOTYPE)
table(d2$COHORT)
# KOR  MSSNG  SPARK    SSC 
# 61578 248026 378140 279462

d2 = d2 %>%
  arrange(CHROM, POS, SAMPLE)

d2 = d2 %>%
  filter(SAMPLE %in% s0$SAMPLE)
nrow(d2) # 967206 (14606 samples)
d2 = unique(d2)
nrow(d2) # 967206 (14606 samples)

t2 = fread('~/Dropbox/DeepLN/SpliceAI/spliceai_kor_sfari_mssng_output.14606samples.20250224.tsv.gz')

# Rename columns
setnames(t2, old = c("ID", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "MaxSpliceAI"),
         new = c("variant", "SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL", "SpliceAI_MaxSpliceAI"))

# Assuming your two data frames are named mg and df2

score <- merge(d2, t2, by = c("variant", "SAMPLE"), all = TRUE)

# score %>%
#   filter(variant=='chr10:108335424:GTCATAGCAA:G')

rm(list = ls()[!(ls()%in%c('score', 'date', 's0'))])




#### Annotate gene and TSS
dnv = data.table::fread('~/Dropbox/AI_NDD/Data/variants/Korean_SFARI_MSSNG.DNV_list.all_nearest_transcript.20250226.tsv.gz')
dnv$gene_id = do.call(rbind.data.frame, strsplit(x = dnv$gene_id, split = '.', fixed = T))[[1]]

lof = read.delim('~/Dropbox/Resources/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
lof = lof %>%
  dplyr::select(gene_name = gene,
                gene_id,
                exac_pLI,
                mis_z,
                oe_lof_upper)
tmp2 = dnv %>%
  dplyr::select(gene_id,
                gene_name) %>%
  unique() %>%
  filter(!is.na(gene_id) & !is.na(gene_name)) %>%
  unique()
table(tmp2$gene_id %in% lof$gene_id)
table(tmp2$gene_name %in% lof$gene_name)
if(T){
  notin = tmp2 %>%
    filter(!(gene_id %in% lof$gene_id))
}

mg = merge(dnv,
           lof %>%
             dplyr::select(-gene_name,
                           gene_id,
                           pLI = exac_pLI,
                           Mis_Z = mis_z,
                           LOEUF = oe_lof_upper),
           by = 'gene_id',
           all.x = T)

f_mg = merge(score,
             mg %>%
               dplyr::select(-transcript_id),
             by = c('SAMPLE', 'variant'))


################################################################################################
################################################################################################
table(f_mg$COHORT,
      is.na(f_mg$Mis_Z))

table(f_mg$COHORT,
      is.na(f_mg$LOEUF))



data.table::fwrite(f_mg %>%
                     dplyr::select(-CHROM, -POS),
                   paste('~/Dropbox/AI_NDD/Data/variants/Korean_SFARI_MSSNG.DNV_list.annotated', date, 'tsv.gz', sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")


