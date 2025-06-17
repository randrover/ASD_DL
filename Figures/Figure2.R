## Figure 2. ASD cluster discovery and clustering

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(paletteer)
library(cowplot)
library(patchwork)
library(ggpubr)
library(data.table)
library(ggalluvial)

if(F){
  extrafont::font_import(pattern = "Arial", prompt = F)
  extrafont::loadfonts()
}

setwd('~/Dropbox/paper_ASD_ML/')

# Panel a: SHAP feature contribution heatmap + hierarchical clustering dendrogram
library(tidyverse)
library(pheatmap)
shap = data.table::fread('Data/model_result/new_agg_v27_modelA_train_after_fs/kor_sfari_mssng.feature_importance.test_sample.shap.tsv.gz')

s0 = read.delim('Data/Korean_SFARI_MSSNG.14428samples.sample_list.20241213.txt')

mg = merge(s0,
           shap,
           by = 'SAMPLE',
           all.y = T)
mg = mg %>%
  filter(isASD==1)

numeric_cols = 12:ncol(mg)
d0 = mg[, numeric_cols] %>% as.data.frame()
rownames(d0) = mg$SAMPLE

heatmap_matrix = d0 %>%
  as.matrix()

heatmap_matrix[is.na(heatmap_matrix)] <- 0  # NA를 0으로 대체
heatmap_matrix[is.infinite(heatmap_matrix)] <- 0  # Inf를 0으로 대체

sums = colSums(heatmap_matrix, na.rm = TRUE)
r_idx = which(sums==0)

heatmap_matrix <- heatmap_matrix[, sums != 0]

heatmap_matrix <- heatmap_matrix[rowSums(heatmap_matrix != 0) > 0, ]

annotation_data <- data.frame(
  row.names = c(mg$SAMPLE) # Feature 이름과 일치해야 함
)
annotation_data$SAMPLE = rownames(annotation_data)

annotation_data = merge(annotation_data, s0,
                        by = 'SAMPLE')
rownames(annotation_data) = annotation_data$SAMPLE
annotation_data = annotation_data %>%
  dplyr::select(-SAMPLE, -dmMIS_count, -dnPTV_count, -conPTV_count)

annotation_data = annotation_data %>%
  dplyr::select(-sequencing_type, -isASD, -Group)

# annotation_data에  cluster 정보 추가
s1 = read.delim('Data/table.data_cluster.20250125.txt') # Cluster데이터 가져오기

annotation_data$SAMPLE <- rownames(annotation_data)

annotation_data <- merge(annotation_data, s1[, c("SAMPLE", "Cluster")],
                         by = "SAMPLE", all.x = TRUE)
annotation_data$Cluster <- factor(annotation_data$Cluster, levels = c(1, 2, 3, 4))

rownames(annotation_data) <- annotation_data$SAMPLE
annotation_data$SAMPLE <- NULL



# isID에서 NA를 문자열 "NA"로 표시.
annotation_data$isID <- as.character(annotation_data$isID)
annotation_data$isID[is.na(annotation_data$isID)] <- "NA"

# annotation_data를 heatmap_matrix의 행 이름 순서로 정렬
annotation_data <- annotation_data[rownames(heatmap_matrix), , drop = FALSE]

# 확인
print(table(rownames(annotation_data) == rownames(heatmap_matrix)))


heatmap_matrix_t = t(heatmap_matrix)

table(colnames(heatmap_matrix_t) == rownames(annotation_data))

annotation_colors = list(
  isID = c(
    "0" = "lightgrey",
    "1" = "black",
    "NA" = "white"
  ),
  Sex = c(
    "Male" = "#1f77b4",
    "Female" = "lightpink"
  ),
  ancestry = c(
    "AFR" = "#8c564b",      # Deep brown
    "AMR" = "#e377c2",      # Light magenta
    "EAS" = "#17becf",      # Cyan-blue
    "EUR" = "#2ca02c",      # Green
    "SAS" = "#ff7f0e",      # Orange
    "UNKNOWN" = "#c7c7c7"   # Light grey for unknown
  ),
  cohort = c(
    "Korean" = "#bcbd22",   # Olive green
    "MSSNG" = "#9467bd",    # Purple
    "SPARK" = "#1f77b4",    # Blue (same as Male, adjust if needed)
    "SSC" = "#d62728"       # Red (same as Female, adjust if needed)
  ),
  Cluster = c(
    "1" = "#EE6677",  
    "2" = "#228833",  
    "3" = "#CCBB44",  
    "4" = "#66CCEE"
  )
)

set.seed(123)
date = '20250524'
pdf(paste('~/Dropbox/Noncoding_kor_ASD_WD/Figures/plot.heatmap.sample_cluster', date, 'pdf',
          sep = '.'),
    height = 5, width = 7)
pheatmap(
  mat = heatmap_matrix_t,
  scale = "column",
  clustering_method = "ward.D2", # Hierarchical clustering method
  color = colorRampPalette(c("#009E73", "white", "#D55E00"))(50), # 색상 그라데이션. # <- SHC changed it here
  # main = "Hierarchical Heatmap: ASD vs ASD+ID",
  breaks = seq(-3, 3, length.out = 51),
  fontsize = 6,                   
  annotation_col = annotation_data,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  angle_col=90
)
dev.off()


# Panel b: t-SNE or UMAP plot of samples colored by cluster
library(tidyverse)
library(class)  # For kNN
library(ggplot2)  # For visualization (optional)
library(umap)

if(T){
  df = data.table::fread('Data/model_result/new_agg_v27_modelA_train_after_fs/kor_sfari_mssng.feature_importance.test_sample.shap.tsv.gz')
  df2 = read.delim('Data/table.data_cluster.20250125.txt')
  df = df %>%
    filter(SAMPLE %in% df2$SAMPLE)
  df2 <- df2 %>%
    filter(SAMPLE %in% df$SAMPLE) %>%
    arrange(match(SAMPLE, df$SAMPLE))
  table(df$SAMPLE==df2$SAMPLE)
}
if(T){
  dt = data.table::fread('Data/model_result/new_agg_v27_modelA_train_after_fs/kor_sfari_mssng.feature_importance.modelB_sample.shap.tsv.gz')
  dt2 = read.delim('Data/table.modelB.data_cluster.20250125.txt')
  dt = dt %>%
    filter(SAMPLE %in% dt2$SAMPLE)
  dt2 <- dt2 %>%
    filter(SAMPLE %in% dt$SAMPLE) %>%
    arrange(match(SAMPLE, dt$SAMPLE))
  table(dt$SAMPLE==dt2$SAMPLE)
}

train_data = df %>% dplyr::select(-SAMPLE) %>% as.matrix()
test_data = dt %>% dplyr::select(-SAMPLE) %>% as.matrix()

scaled_train_data <- t(apply(train_data, 1, function(x) {
  (x - mean(x)) / sd(x)  # Z-score scaling: (값 - 평균) / 표준편차
}))
scaled_test_data <- t(apply(test_data, 1, function(x) {
  (x - mean(x)) / sd(x)  # Z-score scaling: (값 - 평균) / 표준편차
}))

set.seed(2025)
# Apply UMAP with adjusted parameters
umap_config <- umap.defaults
umap_config$n_neighbors <- 10  # Reduce neighbors to emphasize local structure
umap_config$min_dist <- 1   # Increase distance to spread points better
umap_config$metric <- "euclidean"  # Use cosine similarity for clustering
umap_config$spread = 1.5

umap_train <- umap(scaled_train_data, config = umap_config)
umap_test <- umap(scaled_test_data, config = umap_config)

umap_train_df <- data.frame(UMAP1 = umap_train$layout[,1], UMAP2 = umap_train$layout[,2], Cluster = as.factor(df2$Cluster),
                            SAMPLE = df2$SAMPLE)
umap_test_df <- data.frame(UMAP1 = umap_test$layout[,1], UMAP2 = umap_test$layout[,2], Cluster = as.factor(dt2$Cluster),
                           SAMPLE = dt2$SAMPLE)

cluster_colors <- c( # <- SHC changed it here 
  "1" = "#EE6677",  
  "2" = "#228833",  
  "3" = "#CCBB44",  
  "4" = "#66CCEE",  
  "Group A" = "#999933",
  "Group B" = "#AA4499",  
  "Control" = "#DDDDDD" 
)

# Save UMAP visualization
p3 = rbind(cbind(umap_train_df, Dataset = "Group A test set"), cbind(umap_test_df, Dataset = "Group B")) %>%
  mutate(Dataset = factor(x = Dataset, levels = c("Group A test set", 'Group B'))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = cluster_colors) +  # 🔥 Here!
  guides(colour = guide_legend(override.aes = list(size=3))) +
  facet_wrap(~Dataset) +
  labs(title = "Cluster assignment") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        plot.title = element_text(size = 10, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'),
        legend.title = element_text(size = 8, color = 'black'))

tm = rbind(cbind(umap_train_df, Dataset = "Group A test set"), cbind(umap_test_df, Dataset = "Group B"))
write.table(tm,
            '~/Dropbox/paper_ASD_ML/Data/table.figure2.UMAP.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')

cluster_colors <- c( # <- SHC changed it here
  "Cluster 1" = "#EE6677",  
  "Cluster 2" = "#228833",  
  "Cluster 3" = "#CCBB44",  
  "Cluster 4" = "#66CCEE",  
  "Group A" = "#999933",
  "Group B" = "#AA4499",  
  "Control" = "#DDDDDD" 
)


# panel d:
pred = read.delim('Data/model_result/new_agg_v27_modelA_train_after_fs/kor_sfari_mssng_predict_prob.modelB_sample.oversample.tsv.gz')

df = read.delim('Data/table.modelB.data_cluster.20250125.txt')
head(pred)

mg = merge(pred,
           df,
           by = 'SAMPLE')


########################
########################
cluster_colors <- c( # <- SHC changed it here
  "1" = "#EE6677",  
  "2" = "#228833",  
  "3" = "#CCBB44",  
  "4" = "#66CCEE",  
  "Group A" = "#999933",
  "Group B" = "#AA4499",  
  "Control" = "#DDDDDD" 
)


# Panel e: AUROC per cluster
d1 = read.delim('Data/model_result/seed_123/autoB_cluster1/evaluation_results.tsv') %>%
  mutate(Cluster = '1')
d2 = read.delim('Data/model_result/seed_123/autoB_cluster2/evaluation_results.tsv') %>%
  mutate(Cluster = '2')
d3 = read.delim('Data/model_result/seed_123/autoB_cluster3/evaluation_results.tsv') %>%
  mutate(Cluster = '3')
d4 = read.delim('Data/model_result/seed_123/autoB_cluster4/evaluation_results.tsv') %>%
  mutate(Cluster = '4')

# Bind all dataframes together
df_all = bind_rows(d1, d2, d3, d4)
p4 = df_all %>%
  mutate(Cluster = factor(Cluster, levels = seq(1,4))) %>%
  ggplot(aes(x = Cluster,
             y = roc_auc)) +
  geom_col(aes(fill = Cluster),
           color = 'black',
           linewidth = 0.3,
           show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = cluster_colors) +  # 🔥 Here!
  labs(y = 'AUROC',
       title = 'Cluster-specific model performance') +
  ylim(c(0,1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme(axis.title = element_text(size = 8, color ='black'),
        axis.text = element_text(size = 8, color = 'black'),
        plot.title = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'))



######################





df <- fread('Data/kor_sfari_mssng.all_features.tsv.gz')
# df <- fread('Data/model_result/kor_sfari_mssng.shap_top_fs.8L.tsv.gz')
d0 <- read.delim('Data/table.modelB.data_cluster.20250125.txt')
s0 <- read.delim('Data/Korean_SFARI_MSSNG.14428samples.sample_list.20241213.txt') %>%
  filter(isASD == 1)

# Merge and preprocess
mg <- merge(d0, s0, all.y = TRUE) %>%
  select(SAMPLE, Cluster) %>%
  merge(df, by = 'SAMPLE') %>%
  as.data.frame() %>%
  filter(!is.na(Cluster)) %>%
  mutate(Cluster = as.character(Cluster))

phenos = colnames(mg)[!(colnames(mg) %in% c('SAMPLE', 'Cluster'))]
target_phenos <- c('nc_prop_JARVIS_0.9', 'all_ZNF213', 'None_FMRPDarnell_None', 'all_SPIB_8.0')

kw_pvals <- sapply(phenos, function(pheno) {
  tryCatch({
    kruskal.test(reformulate("Cluster", response = pheno), data = mg)$p.value
  }, error = function(e) {
    return(NA)
  })
})

kw_pvals_bonf <- p.adjust(kw_pvals, method = "bonferroni")
bonferroni_table <- setNames(kw_pvals_bonf, phenos)

# STEP 3: Formatting function for adjusted p-values
format_adj_p <- function(p) {
  paste0("adj. p = ", formatC(p, format = "e", digits = 4))
}
bonferroni_table_fmt <- sapply(bonferroni_table, format_adj_p)

# STEP 4: Run and plot only for 4 phenos
fin_res <- data.frame()
p_list <- list()

for (m in target_phenos) {
  message("Analyzing: ", m)
  stats <- run_cluster_comparison(mg, m)
  fin_res <- rbind(fin_res, stats %>% mutate(phenotype = m))
  
  if (any(stats$Pvalue_FDR < 0.05, na.rm = TRUE)) {
    p_list[[m]] <- plot_phenotype(mg, m, stats, adj_pval = bonferroni_table_fmt[[m]])
  }
}

# Define function for pairwise cluster comparison
run_cluster_comparison <- function(data, phenotype) {
  data$value <- data[[phenotype]]
  data <- data %>% arrange(Cluster)
  
  cluster_pairs <- combn(unique(data$Cluster), 2, simplify = TRUE)
  fin_df <- data.frame(stringsAsFactors = FALSE)
  
  for (pair in seq_len(ncol(cluster_pairs))) {
    clust1 <- cluster_pairs[1, pair]
    clust2 <- cluster_pairs[2, pair]
    
    subset_data <- data %>%
      filter(Cluster %in% c(clust1, clust2)) %>%
      filter(!is.na(value))
    
    counts <- table(subset_data$Cluster)
    
    if (length(counts) != 2 || any(counts < 30)) {
      pval <- NA
    } else {
      pval <- wilcox.test(value ~ Cluster, data = subset_data)$p.value
    }
    
    fin_df <- rbind(fin_df, 
                    data.frame(Cluster1 = clust1,
                               Cluster2 = clust2,
                               Count_Cluster1 = counts[as.character(clust1)],
                               Count_Cluster2 = counts[as.character(clust2)],
                               Pvalue = pval))
  }
  
  fin_df$Pvalue_FDR <- p.adjust(fin_df$Pvalue, method = "fdr")
  fin_df$sigmark <- case_when(
    is.na(fin_df$Pvalue_FDR) ~ 'ns',
    fin_df$Pvalue_FDR < 0.0001 ~ '***',
    fin_df$Pvalue_FDR < 0.001 ~ '**',
    fin_df$Pvalue_FDR < 0.05 ~ '*',
    TRUE ~ 'ns'
  )
  
  return(fin_df)
}

cluster_colors <- c( # <- SHC changed it here 
  "1" = "#EE6677",  
  "2" = "#228833",  
  "3" = "#CCBB44",  
  "4" = "#66CCEE",  
  "Group A" = "#999933",
  "Group B" = "#AA4499",  
  "Control" = "#DDDDDD" 
)

# Define function for plotting
plot_phenotype <- function(data, phenotype, stats_df, adj_pval) {
  data$value <- data[[phenotype]]
  
  kw_p <- kruskal.test(value ~ Cluster, data = data)$p.value
  # kw_p_str <- paste0("p = ", formatC(kw_p, format = "f", digits = 4))
  kw_p_str <- paste0("p = ", format(signif(kw_p, digits = 3), scientific = TRUE))
  
  subtitle_text <- case_when(
    phenotype == 'nc_prop_JARVIS_0.9' ~ "JARVIS≥0.9 region",
    phenotype == 'all_ZNF213' ~ "ZNF213 motif overlap",
    phenotype == 'None_FMRPDarnell_None' ~ "FMRP target overlap",
    phenotype == 'all_SPIB_8.0' ~ "SPIB motif score ≥ 8.0",
    TRUE ~ ""
  )
  
  subtitle_full <- paste0(subtitle_text, "\n", adj_pval)
  
  plot <- data %>%
    ggplot(aes(x = Cluster, y = value)) +
    geom_boxplot(width = 0.4, color = "black", alpha = 0.7,
                 aes(fill = Cluster), linewidth = 0.2, outlier.size = 0.1,
                 show.legend = FALSE) +
    geom_violin(aes(fill = Cluster), alpha = 0.7,
                linewidth = 0.2, width = 0.55,
                show.legend = FALSE) +
    scale_fill_manual(values = cluster_colors) +
    labs(
      title = phenotype,
      subtitle = subtitle_full,
      x = "Cluster", y = "Proportion of noncoding variants"
    ) +
    theme_classic() +
    geom_hline(yintercept = mean(data$value, na.rm = TRUE),
               color = 'red', linetype = 'dashed') +
    theme(text = element_text(size = 8, color = 'black', family="Arial"),
          axis.text = element_text(size = 8, color = 'black'),
          axis.title = element_text(size = 8, color = 'black'),
          plot.title = element_text(size = 8, color = 'black', family="Arial"),
          plot.subtitle = element_text(size = 7, color = 'black', family="Arial"),
          plot.caption = element_text(size = 6, color = 'gray40', family="Arial", hjust = 1))
  
  return(plot)
}

# # Run analysis
# for (m in phenos) {
#   message("Analyzing: ", m)
#   stats <- run_cluster_comparison(mg, m)
#   fin_res <- rbind(fin_res, stats %>% mutate(phenotype = m))
#   
#   if (any(stats$Pvalue_FDR < 0.05, na.rm = TRUE)) {
#     p_list[[m]] <- plot_phenotype(mg, m, stats, adj_pval = bonferroni_table[[m]])
#   }
# }




######################




######################
p1a = NULL
p2a = NULL
p3a = p3 +
  theme(axis.text=element_text(size=8, color = 'black'),
        text = element_text(family="Arial"),
        axis.title=element_text(size=8)
        # plot.margin = ggplot2::margin(l = 0.35, b = 0.2, t = 0.3, unit = 'cm')
        )
p4a = p4 +
  theme(axis.text=element_text(size=8, color = 'black'),
        text = element_text(family="Arial"),
        axis.title=element_text(size=8)
        # plot.margin = ggplot2::margin(l = 0.35, b = 0.2, t = 0.3, unit = 'cm')
  )
p5 = plot_grid(plotlist = p_list, nrow = 1)
p5a = p5 +
  theme(axis.text=element_text(size=8, color = 'black'),
        text = element_text(family="Arial"),
        axis.title=element_text(size=8)
        # plot.margin = ggplot2::margin(l = 0.35, b = 0.2, t = 0.3, unit = 'cm')
        )

p_grid1 = plot_grid(p2a, p3a,
                    ncol = 1,
                    rel_heights = c(1.4, 1),
                    labels = c('b', 'c'),
                    label_size = 10)

p_grid2 = plot_grid(NULL, p_grid1,
                    nrow = 1,
                    rel_widths = c(1.2, 1),
                    labels = c('a', ''),
                    label_size = 10)
p_grid3 = plot_grid(p4a, p5a,
                    nrow = 1,
                    rel_widths = c(0.4, 1),
                    labels = c('d', 'e'),
                    label_size = 10)
final_p = plot_grid(p_grid2,
                    p_grid3,
                    ncol = 1,
                    rel_heights = c(2, 1.3),
                    label_size = 10) # <- SHC changed it here (label size 10)

# 21.61 × 27.9
ggsave('Figures/Figure2_v1.9.pdf', final_p, width = 21.61*1.1, height = 27.9*0.65, limitsize = F, units = "cm")

