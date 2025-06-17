rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dorothea)
library(igraph)
library(tidyverse)
library(cowplot)
library(tibble)
library(ggplot2)
library(ggraph)
library(decoupleR)
library(ggrepel)
library(ggrastr)
library(patchwork)

if(F){
  extrafont::font_import(pattern = "Arial", prompt = F)
  extrafont::loadfonts()
}

setwd('~/Dropbox/paper_ASD_ML/')

#####
# 3. disruption_mat: 샘플 x TF 행렬 (data.frame 또는 matrix)
dt = data.table::fread('Data/motif_dff/table.kor_sfari_mssng_14606.motif_diff_by_var.20241225.tsv.gz')
df = data.table::fread('Data/Korean_SFARI_MSSNG.DNV_list.annotated.20250226.tsv.gz')
df = df %>%
  dplyr::select(SAMPLE, variant, is_coding)
df = df %>%
  filter(is_coding==0)
dt = dt %>%
  filter(variant %in% df$variant)

# 2. TF motif columns만 선택
tf_cols <- setdiff(names(dt), c("variant", "SAMPLE"))

# 3. reshape 없이 non-zero 값만 long-format으로 추출
long_dt_sparse <- rbindlist(lapply(tf_cols, function(tf) {
  rows <- which(dt[[tf]] != 0)
  if (length(rows) == 0) return(NULL)
  data.table(
    variant = dt$variant[rows],
    SAMPLE = dt$SAMPLE[rows],
    TF = tf,
    score = dt[[tf]][rows]
  )
}))

if(F){
  # 3. 샘플-TF별 total score와 변이 개수 계산
  num_vars = dt %>%
    group_by(SAMPLE) %>%
    dplyr::count(name = 'num_variants')
  
  summary_dt <- long_dt_sparse %>%
    group_by(SAMPLE, TF) %>%
    summarise(
      total_score = sum(score),
      num_disrupt_var = n(),
      avg_score_per_variant = mean(score),
      .groups = "drop"
    )
  
  mg = merge(num_vars,
             summary_dt,
             by = 'SAMPLE')
  # TF별로: 샘플별 total_score vs. num_variants의 상관관계 계산
  cor_summary <- mg %>%
    group_by(TF) %>%
    summarise(
      pearson_r = cor(total_score, num_variants, method = "pearson"),
      spearman_r = cor(total_score, num_variants, method = "spearman"),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(pearson_r)))  # 절대값 기준 정렬
  
  # 결과 확인
  print(head(cor_summary))
  
  tf_low_corr <- cor_summary$TF[which.max(cor_summary$pearson_r)]
  
  plot_dt <- mg %>% filter(TF == tf_low_corr)
  
  ggplot(plot_dt, aes(x = num_variants, y = total_score)) +
    geom_point(alpha = 0.6, color = "tomato") +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    labs(
      title = paste0("TF = ", tf_low_corr, 
                     " | Pearson R = ", round(cor(plot_dt$num_variants, plot_dt$total_score), 3)),
      x = "Total Variant Count",
      y = paste0("Total Disruption Score (", tf_low_corr, ")")
    ) +
    theme_minimal(base_size = 13)
  
}

# 4. 샘플 × TF로 group_by → 총 disruption score 계산 (또는 평균도 가능)
agg_dt <- long_dt_sparse %>%
  group_by(SAMPLE, TF) %>%
  summarise(total_disruption = score[which.max(abs(score))], .groups = "drop")

# 5. wide-format으로 변환: sample × TF 행렬
disruption_mat <- pivot_wider(agg_dt, 
                              names_from = TF, 
                              values_from = total_disruption, 
                              values_fill = 0)

# 6. rownames로 샘플 지정
disruption_mat <- as.data.frame(disruption_mat)
rownames(disruption_mat) <- disruption_mat$SAMPLE
disruption_mat$SAMPLE <- NULL

# 결과 확인
print(dim(disruption_mat))
print(head(disruption_mat[,1:5]))
#####

true_cluster_df = read.delim('Data/table.modelB.data_cluster.20250125.txt')

# 1. long-format으로 변환
disruption_long <- disruption_mat %>%
  rownames_to_column("SAMPLE") %>%
  pivot_longer(-SAMPLE, names_to = "TF", values_to = "score") %>%
  dplyr::left_join(true_cluster_df, by = "SAMPLE")  # cluster column 필요

s0 = read.delim('Data/model_result/new_agg_v27_modelA_train_after_fs/kor_sfari_mssng_predict_prob.modelB_sample.oversample.tsv.gz')
s1 = read.delim('Data/Korean_SFARI_MSSNG.14606samples.sample_list.20250224.txt')
s1 = s1 %>%
  filter(SAMPLE %in% s0$SAMPLE)
controls = s1 %>%
  filter(isASD==0)

disruption_long = disruption_long %>%
  filter(SAMPLE %in% s0$SAMPLE)
disruption_long$Cluster = ifelse(disruption_long$SAMPLE %in% controls$SAMPLE,
                                 'Control',
                                 disruption_long$Cluster)

disruption_long = disruption_long %>%
  filter(!is.na(Cluster))

# control 샘플 필터링
control_scores <- disruption_long %>%
  filter(Cluster == "Control")

disruption_long <- disruption_long %>%
  mutate(is_high_disruption = ifelse(score>0 | score<0,
                                     TRUE,
                                     FALSE))

table(disruption_long$Cluster)
table(disruption_long$is_high_disruption)

# 2. Wilcoxon test: TF별로 cluster별 enrichment 비교
tf_enrichment_results <- list()

for (cluster_id in c(1, 2, 3, 4)) {
  
  df <- disruption_long %>%
    filter(Cluster %in% c(cluster_id, "Control")) %>%
    mutate(target = ifelse(Cluster == cluster_id, "IN", "CONTROL"))
  
  tf_stats <- df %>%
    group_by(TF) %>%
    summarise(
      wilcox_p = tryCatch(wilcox.test(score ~ target)$p.value, error = function(e) NA),
      effect = mean(score[target == "IN"]) - mean(score[target == "CONTROL"]),
      
      in_high  = sum(is_high_disruption & target == "IN"),
      in_low   = sum(!is_high_disruption & target == "IN"),
      ctrl_high = sum(is_high_disruption & target == "CONTROL"),
      ctrl_low  = sum(!is_high_disruption & target == "CONTROL"),
      
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      fisher_p = tryCatch(
        fisher.test(matrix(c(in_high, in_low, ctrl_high, ctrl_low), nrow = 2))$p.value,
        error = function(e) NA
      ),
      fisher_OR = tryCatch(
        fisher.test(matrix(c(in_high, in_low, ctrl_high, ctrl_low), nrow = 2))$estimate,
        error = function(e) NA
      )
    ) %>%
    ungroup() %>%
    mutate(
      wilcox_FDR = p.adjust(wilcox_p, method = "fdr"),
      fisher_FDR = p.adjust(fisher_p, method = "fdr"),
      Cluster = cluster_id
    )
  
  tf_enrichment_results[[as.character(cluster_id)]] <- tf_stats
}

# Combine all clusters into a single dataframe
enrichment_df <- bind_rows(tf_enrichment_results)

# plot_df: enrichment_df with Cluster, TF, effect, FDR
volcano_df <- enrichment_df %>%
  filter(!is.na(fisher_FDR), !is.na(fisher_OR)) %>%
  mutate(
    logFDR = -log10(fisher_FDR),
    logOR = log2(fisher_OR),
    Cluster = paste("Cluster", Cluster)
  )

cluster_colors <- c( # <- SHC changed it here 
  "Cluster 1" = "#EE6677",  
  "Cluster 2" = "#228833",  
  "Cluster 3" = "#CCBB44",  
  "Cluster 4" = "#66CCEE",  
  "Group A" = "#999933",
  "Group B" = "#AA4499",  
  "Control" = "#DDDDDD" 
)

volcano_df <- volcano_df %>%
  mutate(
    label_y_pos = ifelse(logFDR > 4, 4.2, logFDR)
  )

volcano_df = volcano_df %>%
  mutate(logOR = ifelse(Cluster=='Cluster 1' & logOR<(-1), -1,
                        logOR)) ## ZNF213
volcano_df = volcano_df %>%
  mutate(logOR = ifelse(Cluster=='Cluster 2' & logOR > (1.5), 1.5,
                        logOR)) ## ZNF213
volcano_df = volcano_df %>%
  mutate(logOR = ifelse(Cluster=='Cluster 4' & logOR < (-1), -1,
                        logOR)) ## ZNF213

xlims <- list(
  "Cluster 1" = c(-1, 1),
  "Cluster 2" = c(-1.5, 1.5),
  "Cluster 3" = c(-1, 1),
  "Cluster 4" = c(-1.0, 1.0)
)

xlim_dummy <- bind_rows(lapply(names(xlims), function(clu) {
  tibble(
    Cluster = clu,
    logOR = xlims[[clu]],  # min, max
    label_y_pos = 0  # y는 영향 없도록 더미 값
  )
}))


p1 <- ggplot(volcano_df, aes(x = logOR, y = label_y_pos, fill = Cluster)) +
  geom_blank(data = xlim_dummy) +
  geom_point_rast(alpha = 0.8, size = 2.4, shape = 21, color = "black", stroke = 0.3, show.legend = FALSE,
                  raster.dpi=600) +
  geom_hline(yintercept = -log10(0.1), linetype = "dotted", color = "red", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  facet_wrap(~ Cluster, nrow = 2, strip.position = "top", scales = "free") +
  ggrepel::geom_text_repel(
    data = subset(volcano_df, fisher_FDR < 0.1),
    aes(label = TF),
    size = 8*(5/14),
    max.overlaps = 28,
    box.padding = 0.25,
    min.segment.length = 0
  ) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "TF motif disruption enrichment per cluster",
    x = expression(log[2]~"(Odds Ratio)"),
    y = expression(-log[10]~"(FDR)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.spacing = unit(1.2, "lines"),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.2, -0.2, -0.2, -0.2), "cm")
  )
p1


######################
p1a = p1 +
  theme(axis.text=element_text(size=8, color = 'black'),
        text = element_text(family="Arial"),
        axis.title=element_text(size=8)
        # plot.margin = ggplot2::margin(l = 0.35, b = 0.2, t = 0.3, unit = 'cm')
  )

# p3a = p3 +
#   theme(axis.text=element_text(size=8, color = 'black'),
#         text = element_text(family="Arial"),
#         axis.title=element_text(size=8)
#         # plot.margin = ggplot2::margin(l = 0.35, b = 0.2, t = 0.3, unit = 'cm')
#   )

p_grid1 = plot_grid(p1a, NULL,
                    nrow = 2,
                    rel_heights = c(1, 0.35),
                    labels = c('a', 'c'),
                    label_size = 10)
p_grid2 = plot_grid(NULL,NULL,
                    ncol = 1)
final_p = plot_grid(p_grid1,
                    p_grid2,
                    ncol = 2,
                    rel_widths = c(1, 1),
                    labels = c('', 'b'),
                    label_size = 10)


# 21.61 × 27.9
ggsave('Figures/Figure3_v1.9.pdf', final_p, width = 21.61*1, height = 27.9*0.65, limitsize = F, units = "cm")



