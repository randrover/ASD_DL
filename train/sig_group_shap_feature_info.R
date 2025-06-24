rm(list = ls())
library(tidyverse)
date = '20250306'

mf1 = read.delim('kor_sfari_mssng.feature_selection_4L.shap.mean.tsv.gz')

# Define a function to extract top N% based on the MeanSHAP column
get_top_percent <- function(df, column_name, percent = 10) {
  # Calculate the threshold for the given percentile
  percent = percent*0.01
  threshold <- quantile(df[[column_name]], 1 - percent, na.rm=T)
  
  # Filter rows where the column value is greater than or equal to the threshold
  top_percent_df <- df[df[[column_name]] >= threshold, ]
  
  return(top_percent_df)
}

# Example usage: Get the top 10% of MeanSHAP
# top_10_percent <- get_top_percent(mf1, "MeanSHAP", 10)

if(T){
  # top 0.1%
  sig1a <- get_top_percent(mf1, "MeanSHAP", 0.1) %>% pull(Feature)
  # top 0.2%
  sig1b <- get_top_percent(mf1, "MeanSHAP", 0.2) %>% pull(Feature)
  # top 0.3%
  sig1c <- get_top_percent(mf1, "MeanSHAP", 0.3) %>% pull(Feature)
  # top 0.5%
  sig1d <- get_top_percent(mf1, "MeanSHAP", 0.5) %>% pull(Feature)
  # top 0.8%
  sig1e <- get_top_percent(mf1, "MeanSHAP", 0.8) %>% pull(Feature)
  # top 1%
  sig1 <- get_top_percent(mf1, "MeanSHAP", 1) %>% pull(Feature)
  # top 3%
  sig2 <- get_top_percent(mf1, "MeanSHAP", 3) %>% pull(Feature)
  # top 5%
  sig3 <- get_top_percent(mf1, "MeanSHAP", 5) %>% pull(Feature)
  # top 10%
  sig4 <- get_top_percent(mf1, "MeanSHAP", 10) %>% pull(Feature)
  # top 30%
  sig5 <- get_top_percent(mf1, "MeanSHAP", 30) %>% pull(Feature)
  # top 50%
  sig6 <- get_top_percent(mf1, "MeanSHAP", 50) %>% pull(Feature)
  # top 100% (All)
  sig7 <- get_top_percent(mf1, "MeanSHAP", 60) %>% pull(Feature)
  sig8 <- get_top_percent(mf1, "MeanSHAP", 70) %>% pull(Feature)
  sig9 <- get_top_percent(mf1, "MeanSHAP", 80) %>% pull(Feature)
  sig10 <- get_top_percent(mf1, "MeanSHAP", 90) %>% pull(Feature)
  sig11 <- get_top_percent(mf1, "MeanSHAP", 100) %>% pull(Feature)
  
  mf1 = mf1 %>%
    mutate(Group = ifelse(Feature %in% sig1a,
                          'Sig01',
                          ifelse(Feature %in% sig1b,
                                 'Sig02',
                                 ifelse(Feature %in% sig1c,
                                        'Sig03',
                                        ifelse(Feature %in% sig1d,
                                               'Sig05',
                                               ifelse(Feature %in% sig1e,
                                                      'Sig08',
                                                      ifelse(Feature %in% sig1,
                                                             'Sig1',
                                                             ifelse(Feature %in% sig2,
                                                                    'Sig3',
                                                                    ifelse(Feature %in% sig3,
                                                                           'Sig5',
                                                                           ifelse(Feature %in% sig4,
                                                                                  'Sig10',
                                                                                  ifelse(Feature %in% sig5,
                                                                                         'Sig30',
                                                                                         ifelse(Feature %in% sig6,
                                                                                                'Sig50',
                                                                                                ifelse(Feature %in% sig7,
                                                                                                       'Sig60',
                                                                                                       ifelse(Feature %in% sig8,
                                                                                                              'Sig70',
                                                                                                              ifelse(Feature %in% sig9,
                                                                                                                     'Sig80',
                                                                                                                     ifelse(Feature %in% sig10,
                                                                                                                            'Sig90',
                                                                                                                            'Other'))))))))))))))))
    

  print(table(mf1$Group))
  print(prop.table(table(mf1$Group)))
}


write.table(mf1,
            paste('table.shap.feature_info_v31', date, 'txt',
                  sep = '.'),
            quote = F, row.names = F, col.names = T, sep = '\t')
