rm(list = ls())
library(tidyverse)
date = '20250226'

s0 = read.delim('~/Dropbox/Noncoding_kor_ASD_WD/Data/Korean_SFARI_MSSNG.14606samples.sample_list.20250224.txt')

dnv = data.table::fread('~/Dropbox/Noncoding_kor_ASD_WD/Data/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples.sorted.20241222.vcf.gz')
dnv$SAMPLE = do.call(rbind.data.frame, strsplit(x = dnv$INFO, split = ';', fixed = T))[[1]]
dnv$SAMPLE = gsub(x = dnv$SAMPLE, pattern = 'SAMPLE=', replacement = '', fixed = T)
dnv = dnv %>%
  filter(SAMPLE %in% s0$SAMPLE)
dnv2 = merge(dnv,
             s0 %>%
               dplyr::select(SAMPLE, isASD, cohort),
             by = 'SAMPLE')
dnv2$PHENOTYPE = ifelse(dnv2$isASD==1, 'case', 'ctrl')
dnv2$INFO = '.'
dnv2$INFO = paste('SAMPLE=', dnv2$SAMPLE, ';PHENOTYPE=', dnv2$PHENOTYPE, ';COHORT=', dnv2$cohort,
                  sep = '')


nrow(s0)
write.table(dnv2 %>%
              dplyr::select(-SAMPLE,
                            -isASD,
                            -cohort),
            paste('~/Dropbox/Noncoding_kor_ASD_WD/Data/Korean_SFARI_MSSNG_WGS_autosomal_DNV.14606samples', date, 'vcf',
                  sep = '.'),
            quote = F, row.names = F, col.names = T, sep = '\t')



sample_list = dnv2 %>%
  dplyr::select(SAMPLE, PHENOTYPE) %>%
  unique()
n_samples = nrow(sample_list)
## sample list
write.table(sample_list,
            paste('~/Dropbox/Noncoding_kor_ASD_WD/Data/Korean_SFARI_MSSNG', n_samples, 'samples.sample_list.', date, '.txt',
                  sep = ''),
            quote = F, row.names = F, col.names = T, sep = '\t')








## Adjust paternal age factor
if(T){
  t1 = read.delim('~/Dropbox/Noncoding_kor_ASD_WD/Data/adjustFile_patAge_cohort.Korean_SSC_SPARK.WGS.autosomal_DNV.10930samples.20240426.txt')
  t2 = read.delim('~/Dropbox/Noncoding_kor_ASD_WD/Data/adjustFile_patAge_cohort.MSSNG_WGS_autosomal_DNV.3696samples.20240819.txt')
  t3 = rbind(t1,t2)
  t3 = t3 %>%
    dplyr::select(SAMPLE, FatherAgeBirth, MotherAgeBirth)
  adj_info = t3 %>%
    filter(SAMPLE %in% sample_list$SAMPLE)
}

## NAs
table(is.na(adj_info$FatherAgeBirth)) # 2457


rates = dnv2 %>%
  group_by(SAMPLE, cohort) %>%
  dplyr::count(name = 'rate')
table(duplicated(rates$SAMPLE))
mg = merge(adj_info, rates, by = 'SAMPLE')

# class of ages
class(mg$FatherAgeBirth)

if(T){
  n_mg = mg
  
  # Separate out samples with missing mutation rate information
  table(is.na(n_mg$rate)) # no sample has missing mutation rate
  
  # Separate out samples with missing paternal age information
  naPat = n_mg[is.na(n_mg$FatherAgeBirth),] # 780 samples
  
  # Separate out samples with complete information
  all_mg <- n_mg[!is.na(n_mg$rate) & !is.na(n_mg$FatherAgeBirth),] # 9972 samples left
  
  # Run linear regression model on DNV rate vs paternal age
  # adjModel1 <- lm(rate ~ FatherAgeBirth, all_mg) # Adjusted R-squared:  0.4513
  # adjModel2 <- lm(rate ~ FatherAgeBirth + cohort, all_mg) # Adjusted R-squared:  0.4662
  # adjModel3 <- lm(rate ~ cohort, all_mg) # Adjusted R-squared:  0.004833
  
  adjModel <- lm(rate ~ FatherAgeBirth + cohort, all_mg) # Adjusted R-squared:  0.4702
  
  #----------------------------------------------------------------------------------------#
  # Incorporate shift so that mean of raw and unadjusted variant counts is the same
  meanRate <- mean(all_mg$rate) # 65.72282
  all_mg$N_dnv_adjRaw <- coef(adjModel)["(Intercept)"] + residuals(adjModel)
  meanAdjRaw <- mean(all_mg$N_dnv_adjRaw) # 14.36053
  shift <- meanRate - meanAdjRaw # 51.36229
  all_mg$N_dnv_adj <- all_mg$N_dnv_adjRaw + shift # Older father, lower adjusted rate
  #----------------------------------------------------------------------------------------#
  # Calculate the percent adjustment (adjustment factor)
  all_mg$AdjustFactor <- all_mg$N_dnv_adj / all_mg$rate # adjustment factor * rate = adjusted rate (Older father, lower adjusted rate, lower adjustment factor)
  
  # Format for output:
  allInfoOut <- all_mg[,c("SAMPLE","FatherAgeBirth", "MotherAgeBirth", "rate","N_dnv_adjRaw","N_dnv_adj","AdjustFactor")]
  colnames(allInfoOut) = c('SAMPLE', 'FatherAgeBirth', "MotherAgeBirth", 'N_dnv', 'N_dnv_adjRaw', 'N_dnv_adj', 'AdjustFactor')
  
  # Add the samples with mutation rates, but no paternal age information
  naPatAgeOut <- naPat[,c("SAMPLE","FatherAgeBirth", "MotherAgeBirth", "rate")]
  naPatAgeOut$N_dnv_adj <- naPatAgeOut$N_dnv_adjRaw <- NA
  naPatAgeOut$AdjustFactor <- 1 # Rates unchanged
  colnames(naPatAgeOut) = c('SAMPLE', 'FatherAgeBirth', 'MotherAgeBirth', 'N_dnv', 'N_dnv_adjRaw', 'N_dnv_adj', 'AdjustFactor')
  
  # Bind outputs
  allInfoOut <- rbind(allInfoOut, naPatAgeOut)
  
  # Write out to file
  outFile <- paste("~/Dropbox/Noncoding_kor_ASD_WD/Data/adjustFile_patAge_cohort.Korean_SFARI_MSSNG.WGS.autosomal_DNV.", nrow(allInfoOut), "samples.", date, ".txt", sep='')
  write.table(allInfoOut, file=outFile, sep="\t", row.names=F, col.names=T, quote=F)
}

## Done!

