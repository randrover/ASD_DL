library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
rm(list = ls())
date = '20250226'

dnv = data.table::fread('Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.extracted_variants.txt.gz')

# table(is.na(dnv$Consequence))

## is coding or not
dnv$is_coding = dnv$is_CodingRegion

dnv$most_severe_consequence = dnv$Consequence
dnv$gene_id = dnv$Gene
dnv$gene_name = dnv$SYMBOL
dnv$transcript_id = dnv$Feature

dnv = dnv %>%
  dplyr::select(SAMPLE, variant = ID,
                most_severe_consequence, gene_id, gene_name, transcript_id,
                is_coding)

s0 = read.delim('Korean_SFARI_MSSNG.14606samples.sample_list.20250224.txt')

#### Annotate gene and TSS
# Load GTF file
d0 = d0 <- rtracklayer::import('gencode.v44.basic.annotation.gtf.gz') %>%
  as.data.frame()

######################## coding ########################
g1 = d0 %>% filter(type == 'gene')
g1 = g1 %>%
  filter(gene_type == 'protein_coding')
g1 = unique(g1)
coding_genes = g1$gene_id

d0 <- rtracklayer::import('gencode.v44.basic.annotation.gtf.gz') %>%
  as.data.frame() %>%
  filter(type == 'transcript')

# Define the TSS position based on the strand
d1 <- d0 %>%
  mutate(TSS = ifelse(strand == '+', start, end))

d1 = d1 %>%
  filter(gene_id %in% coding_genes)

hs.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
vcf = 'Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf.gz'
variants <- motifbreakR::variants.from.file(file = vcf,
                                            search.genome = hs.genome,
                                            format = "vcf")

meta = variants@ranges %>% as.data.frame()
df = data.table::fread(vcf)
meta$INFO = df$INFO
rm(df)
meta$SAMPLE = do.call(rbind.data.frame, strsplit(x = meta$INFO, split = ';', fixed = T))[[1]]
meta$SAMPLE = gsub(x = meta$SAMPLE, pattern = 'SAMPLE=', replacement = '', fixed = T)
meta$tid = paste(meta$SAMPLE, meta$names, sep = '|')

t2 = meta %>%
  filter(SAMPLE %in% s0$SAMPLE)
t2 = t2 %>%
  dplyr::select(-width, -names, -INFO)
t2 = unique(t2)

calculate_TSS <- function(TSS, start, end) {
  if(is.na(TSS)){
    return(NA)
  }
  
  if(TSS<start & TSS<end){
    return(min(abs(start-TSS), abs(end-TSS)))
  }else if(TSS>start & TSS<end){
    return(min(abs(TSS-start), abs(end-TSS)))
  }else if(TSS>start & TSS>end){
    return(min(abs(TSS-start), abs(TSS-end)))
  }else{
    return(NA)
  }
}

t2$start = as.numeric(t2$start)
t2$end = as.numeric(t2$end)

# Prepare GRanges objects for nearest TSS calculation
tss_gr <- GRanges(seqnames = d1$seqnames, ranges = IRanges(start = d1$TSS, end = d1$TSS),
                  strand = d1$strand, gene_id = d1$gene_id, gene_name = d1$gene_name,
                  transcript_id = d1$transcript_id,
                  transcript_type = d1$transcript_type,
                  hgnc_id = d1$hgnc_id,
                  transcript_start = d1$start,
                  transcript_end = d1$end)

# Find nearest TSS for each variant
nearest_idx <- nearest(variants, tss_gr)
nearest_tss <- tss_gr[nearest_idx, ]

# Calculate distances to nearest TSS
t2$TSS_distance <- abs(start(variants) - start(nearest_tss))

# Combine TSS information with the variants
t2$nearest_TSS_gene_id <- nearest_tss$gene_id
t2$nearest_TSS_gene_name <- nearest_tss$gene_name
t2$nearest_TSS_transcript_id <- nearest_tss$transcript_id
t2$nearest_TSS_transcript_type <- nearest_tss$transcript_type
t2$nearest_TSS_hgnc_id <- nearest_tss$hgnc_id
t2$nearest_TSS_start_site = nearest_tss$transcript_start
t2$nearest_TSS_end_site = nearest_tss$transcript_end
t2$nearest_TSS_strand = as.character(nearest_tss@strand)

head(t2)

t2 = t2 %>%
  dplyr::select(tid, nearest_TSS_gene_id, nearest_TSS_gene_name, nearest_TSS_transcript_id,
                nearest_TSS_start_site, nearest_TSS_end_site, nearest_TSS_strand)
dnv = dnv %>%
  mutate(tid = paste(SAMPLE, variant, sep = '|')) %>%
  dplyr::select(SAMPLE, variant, tid, gene_id, gene_name, most_severe_consequence, is_coding,
                transcript_id)

t2$nearest_TSS_transcript_id = do.call(rbind.data.frame, strsplit(x = t2$nearest_TSS_transcript_id, split = '.', fixed = T))[[1]]
t2$nearest_TSS_gene_id = do.call(rbind.data.frame, strsplit(x = t2$nearest_TSS_gene_id, split = '.', fixed = T))[[1]]

mg = merge(dnv, t2,
           by = 'tid')

mg$fin_gene_id = ifelse(mg$is_coding==1,
                        mg$gene_id,
                        mg$nearest_TSS_gene_id)
mg$fin_gene_name = ifelse(mg$is_coding==1,
                          mg$gene_name,
                          mg$nearest_TSS_gene_name)
d2 = d1 %>%
  dplyr::select(seqnames, start, end, strand, transcript_id, TSS, hgnc_id)
d2$transcript_id = do.call(rbind.data.frame, strsplit(x = d2$transcript_id, split = '.', fixed = T))[[1]]
mg2 = merge(mg,
            d2,
            by = 'transcript_id',
            all.x = T)
mg2$fin_transcript_id = ifelse(mg2$is_coding==1,
                               mg2$transcript_id,
                               mg2$nearest_TSS_transcript_id)
mg2$fin_transcrpt_start = ifelse(mg2$is_coding==1,
                                 mg2$start,
                                 mg2$nearest_TSS_start_site)
mg2$fin_transcrpt_end = ifelse(mg2$is_coding==1,
                               mg2$end,
                               mg2$nearest_TSS_end_site)
mg2$fin_strand = ifelse(mg2$is_coding==1,
                        mg2$strand,
                        mg2$nearest_TSS_strand)
mg2 = mg2 %>%
  dplyr::select(SAMPLE, variant,
                is_coding, most_severe_consequence,
                fin_gene_id, fin_gene_name,
                fin_transcript_id, fin_transcrpt_start, fin_transcrpt_end,
                fin_strand)
colnames(mg2) = gsub(x = colnames(mg2), pattern = 'fin_', replacement = '', fixed = T)

data.table::fwrite(mg2,
                   paste('Korean_SFARI_MSSNG.DNV_list.coding_gene_nearest_transcript', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")



######################## all ########################
g1 = d0 %>% filter(type == 'gene')
g1 = g1 %>%
  filter(gene_type == 'protein_coding')
g1 = unique(g1)
coding_genes = g1$gene_id

d0 <- rtracklayer::import('gencode.v44.basic.annotation.gtf.gz') %>%
  as.data.frame() %>%
  filter(type == 'transcript')

# Define the TSS position based on the strand
d1 <- d0 %>%
  mutate(TSS = ifelse(strand == '+', start, end))

hs.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
vcf = 'Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20250226.vcf.gz'
variants <- motifbreakR::variants.from.file(file = vcf,
                                            search.genome = hs.genome,
                                            format = "vcf")

meta = variants@ranges %>% as.data.frame()
df = data.table::fread(vcf)
meta$INFO = df$INFO
rm(df)
meta$SAMPLE = do.call(rbind.data.frame, strsplit(x = meta$INFO, split = ';', fixed = T))[[1]]
meta$SAMPLE = gsub(x = meta$SAMPLE, pattern = 'SAMPLE=', replacement = '', fixed = T)
meta$tid = paste(meta$SAMPLE, meta$names, sep = '|')

t2 = meta %>%
  filter(SAMPLE %in% s0$SAMPLE)
t2 = t2 %>%
  dplyr::select(-width, -names, -INFO)
t2 = unique(t2)

calculate_TSS <- function(TSS, start, end) {
  if(is.na(TSS)){
    return(NA)
  }
  
  if(TSS<start & TSS<end){
    return(min(abs(start-TSS), abs(end-TSS)))
  }else if(TSS>start & TSS<end){
    return(min(abs(TSS-start), abs(end-TSS)))
  }else if(TSS>start & TSS>end){
    return(min(abs(TSS-start), abs(TSS-end)))
  }else{
    return(NA)
  }
}

t2$start = as.numeric(t2$start)
t2$end = as.numeric(t2$end)

# Prepare GRanges objects for nearest TSS calculation
tss_gr <- GRanges(seqnames = d1$seqnames, ranges = IRanges(start = d1$TSS, end = d1$TSS),
                  strand = d1$strand, gene_id = d1$gene_id, gene_name = d1$gene_name,
                  transcript_id = d1$transcript_id,
                  transcript_type = d1$transcript_type,
                  hgnc_id = d1$hgnc_id,
                  transcript_start = d1$start,
                  transcript_end = d1$end)

# Find nearest TSS for each variant
nearest_idx <- nearest(variants, tss_gr)
nearest_tss <- tss_gr[nearest_idx, ]

# Calculate distances to nearest TSS
t2$TSS_distance <- abs(start(variants) - start(nearest_tss))

# Combine TSS information with the variants
t2$nearest_TSS_gene_id <- nearest_tss$gene_id
t2$nearest_TSS_gene_name <- nearest_tss$gene_name
t2$nearest_TSS_transcript_id <- nearest_tss$transcript_id
t2$nearest_TSS_transcript_type <- nearest_tss$transcript_type
t2$nearest_TSS_hgnc_id <- nearest_tss$hgnc_id
t2$nearest_TSS_start_site = nearest_tss$transcript_start
t2$nearest_TSS_end_site = nearest_tss$transcript_end
t2$nearest_TSS_strand = as.character(nearest_tss@strand)

head(t2)

t2 = t2 %>%
  dplyr::select(tid, nearest_TSS_gene_id, nearest_TSS_gene_name, nearest_TSS_transcript_id,
                nearest_TSS_start_site, nearest_TSS_end_site, nearest_TSS_strand)
dnv = dnv %>%
  mutate(tid = paste(SAMPLE, variant, sep = '|')) %>%
  dplyr::select(SAMPLE, variant, tid, gene_id, gene_name, most_severe_consequence, is_coding,
                transcript_id)

t2$nearest_TSS_transcript_id = do.call(rbind.data.frame, strsplit(x = t2$nearest_TSS_transcript_id, split = '.', fixed = T))[[1]]
t2$nearest_TSS_gene_id = do.call(rbind.data.frame, strsplit(x = t2$nearest_TSS_gene_id, split = '.', fixed = T))[[1]]

mg = merge(dnv, t2,
           by = 'tid')

mg$fin_gene_id = ifelse(mg$is_coding==1,
                        mg$gene_id,
                        mg$nearest_TSS_gene_id)
mg$fin_gene_name = ifelse(mg$is_coding==1,
                          mg$gene_name,
                          mg$nearest_TSS_gene_name)

d2 = d1 %>%
  dplyr::select(seqnames, start, end, strand, transcript_id, TSS, hgnc_id)
d2$transcript_id = do.call(rbind.data.frame, strsplit(x = d2$transcript_id, split = '.', fixed = T))[[1]]

mg2 = merge(mg,
            d2,
            by = 'transcript_id',
            all.x = T)
mg2$fin_transcript_id = ifelse(mg2$is_coding==1,
                               mg2$transcript_id,
                               mg2$nearest_TSS_transcript_id)
mg2$fin_transcrpt_start = ifelse(mg2$is_coding==1,
                                 mg2$start,
                                 mg2$nearest_TSS_start_site)
mg2$fin_transcrpt_end = ifelse(mg2$is_coding==1,
                               mg2$end,
                               mg2$nearest_TSS_end_site)
mg2$fin_strand = ifelse(mg2$is_coding==1,
                        mg2$strand,
                        mg2$nearest_TSS_strand)
mg2 = mg2 %>%
  dplyr::select(SAMPLE, variant,
                is_coding, most_severe_consequence,
                fin_gene_id, fin_gene_name,
                fin_transcript_id, fin_transcrpt_start, fin_transcrpt_end,
                fin_strand)
colnames(mg2) = gsub(x = colnames(mg2), pattern = 'fin_', replacement = '', fixed = T)

data.table::fwrite(mg2,
                   paste('Korean_SFARI_MSSNG.DNV_list.all_nearest_transcript', date, 'tsv.gz',
                         sep = '.'),
                   quote = F, row.names = F, col.names = T, sep = '\t',
                   compress="gzip")




