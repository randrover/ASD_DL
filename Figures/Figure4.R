### Last updated: 2025.06.09. by SWK
### project: Figure 4 & TableS14

rm(list=ls())
sessionInfo()
setwd('~/Dropbox/paper_ASD_ML/')

options(stringsAsFactors = F)

library(tidyverse)
library(ggsignif)
library(paletteer)
library(ggpubr)
library(broom)
library(reshape2)
library(cowplot)
library(extrafont)
library(ggpubr)

#font_import()
loadfonts()

version <- 'v2.2'
`%notin%` <- Negate(`%in%`)
select <- dplyr::select

fig_path <- '~/Dropbox/paper_ASD_ML/Figures/Figure4_'

# Set colors
clr_1 = "#EE6677"
clr_2 = "#228833"
clr_3 = "#CCBB44"
clr_4 = "#66CCEE"



### Load data
# Cluster-wise
pheno_tb = read.delim("Data/PhenoOR_ClsWise.SWK.250610.txt")

pheno_tb$Cluster <- factor(pheno_tb$Cluster, levels = unique(pheno_tb$Cluster))
pheno_tb$Phenotype <- factor(pheno_tb$Phenotype, levels = unique(pheno_tb$Phenotype))

# TF-wise
mg2 = read.delim('Data/Phenotable_CaseOnly.TFannot.SWK.250610.txt')



### Plot
## Cluster-wise phenotype distribution (Fig 5A) -----------------------------------------------------------------
pheno_tb$Cluster = factor(pheno_tb$Cluster, levels = c(1,2,3,4))
pheno_tb = pheno_tb %>%
  arrange(Cluster)
p1 = ggplot(pheno_tb, aes(x = Phenotype, y = as.numeric(Cluster), fill = EffectSize)) +
  geom_tile(color = "grey30", linewidth=0.6) +
  geom_text(aes(label = Sigmark), color = "black", size = 3, family = 'Arial') +
  scale_fill_gradient2(
    low = "#5a5a83", mid = "white", high = "#8B3832",
    midpoint = 0, limits = c(-0.55, 0.55),  # Z-score 범위 조정 (원하시면 자동화도 가능)
    name = "Standardized\ncoefficient"
  ) +
  coord_polar(start = -0.15, direction = 1) +
  annotate(x="", y=1:4, label=c(paste0('Cluster ', 1:4)), size = 2.5, geom = "text", col = 'black', family = 'Arial') +
  ylim(c(-1, 4.5)) +
  xlim(c("", levels(pheno_tb$Phenotype))) +
  theme_minimal() +
  theme(
    text = element_text(size = 8, color = "black", family = "Arial"),  # 모든 기본 텍스트
    plot.title = element_text(size = 10, color = "black", family = "Arial", hjust = 0.5),  # 제목만 별도
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),  # text에서 기본 size=8 이지만 명시적으로 넣어도 OK
    legend.title = element_text(size = 8)
  ) +
  ggtitle("Cluster-wise phenotypic enrichment")

p1


## TF network-wise phenotype distribution (Fig 5B) -----------------------------------------------------------------
# VABS_dev (PPARD_carrier_ZNF148_carrier)
# Non_verbal_IQ (ZNF148_carrier_ZNF701_carrier)
# FSIQ (ZNF148_carrier_ZNF701_carrier)
# ADOS_total (ZNF281_carrier_ZNF701_carrier)
# ADOS_SA (ZNF281_carrier_ZNF701_carrier)
# first_word (ZNF701_carrier_ZNF384_carrier)

library(dplyr)
library(ggplot2)
library(purrr)

# 조합별 구성
combo_info <- list(
  list(tf_a = "ZNF148", tf_b = "ZNF701", phenotype = "FSIQ"),
  list(tf_a = "ZNF148",  tf_b = "ZNF701", phenotype = "Non_verbal_IQ"),
  list(tf_a = "PPARD",  tf_b = "ZNF148", phenotype = "VABS_dev"),
  list(tf_a = "ZNF281",  tf_b = "ZNF701", phenotype = "ADOS_total"),
  list(tf_a = "ZNF281",  tf_b = "ZNF701", phenotype = "ADOS_SA"),
  list(tf_a = "ZNF701",  tf_b = "ZNF384", phenotype = "first_word")
)

library(scales)
library(ggpubr)
library(cowplot)

plot_list <- list()

for (i in seq_along(combo_info)) {
  combo <- combo_info[[i]]
  tf_a <- combo$tf_a
  tf_b <- combo$tf_b
  pheno <- combo$phenotype
  var_a <- paste0(tf_a, "_carrier")
  var_b <- paste0(tf_b, "_carrier")
  type_col <- paste0(tf_a, "_", tf_b, "_Type")
  
  # Type 생성
  mg2[[type_col]] <- case_when(
    mg2[[var_a]] == 1 & mg2[[var_b]] == 1 ~ "Both",
    mg2[[var_a]] == 1 & mg2[[var_b]] == 0 ~ paste(tf_a),
    mg2[[var_a]] == 0 & mg2[[var_b]] == 1 ~ paste(tf_b),
    TRUE ~ "None"
  )
  
  df <- mg2 %>%
    filter(!is.na(.data[[pheno]]), !is.na(.data[[type_col]]))
  
  # 요약
  summary_df <- df %>%
    group_by(.data[[type_col]]) %>%
    summarise(
      mean = mean(.data[[pheno]], na.rm = TRUE),
      se = sd(.data[[pheno]], na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lower = mean - qt(0.975, df = n - 1) * se,
      ci_upper = mean + qt(0.975, df = n - 1) * se,
      label = paste0(.data[[type_col]], "\n(n = ", comma(n), ")")
    )
  
  desired_order <- c("Both", tf_a, tf_b, "None")
  type_levels <- desired_order
  
  summary_df[[type_col]] <- factor(summary_df[[type_col]], levels = type_levels)
  summary_df = summary_df %>%
    arrange(.data[[type_col]])
  summary_df <- summary_df %>%
    mutate(label = factor(label, levels = summary_df$label))
  
  
  label_map <- setNames(summary_df$label, summary_df[[type_col]])
  df$label <- factor(label_map[df[[type_col]]], levels = levels(summary_df$label))
  
  # 전체 그룹에서 최대 CI upper
  global_max_ci_upper <- max(summary_df$ci_upper, na.rm = TRUE)
  y_range <- max(summary_df$mean, na.rm = TRUE) - min(summary_df$mean, na.rm = TRUE)
  offset_unit <- y_range * 0.2
  
  # pairwise wilcox test
  combn_mat <- combn(type_levels, 2, simplify = FALSE)
  wilcox_df <- purrr::map_dfr(seq_along(combn_mat), function(k) {
    g1 <- combn_mat[[k]][1]
    g2 <- combn_mat[[k]][2]
    x1 <- df[df[[type_col]] == g1, pheno, drop = TRUE]
    x2 <- df[df[[type_col]] == g2, pheno, drop = TRUE]
    test_res <- suppressWarnings(wilcox.test(x1, x2))
    data.frame(
      group1 = g1,
      group2 = g2,
      group1_label = label_map[g1],
      group2_label = label_map[g2],
      p.value = test_res$p.value,
      stringsAsFactors = FALSE
    )
  }) %>%
    filter(p.value < 0.05) %>%
    arrange(p.value) %>%
    mutate(
      y.position = global_max_ci_upper + row_number() * offset_unit,
      p.adj = sprintf("p = %.3f", p.value)
    )
  
  # Plot
  p <- ggplot(summary_df, aes(x = label, y = mean)) +
    geom_point(size = 3, color = "black") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, color = "black") +
    stat_pvalue_manual(
      wilcox_df,
      label = "p.adj",
      xmin = "group1_label",
      xmax = "group2_label",
      y.position = "y.position",
      tip.length = 0.005,
      size = 2.5
    ) +
    theme_bw() +
    labs(
      x = "",
      y = pheno,
      title = paste(pheno, "by", tf_a, "+", tf_b, "carrier group")
    ) +
    theme(
      text = element_text(size = 8, family = "Arial"),
      # plot.title = element_blank(),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.text = element_text(size = 8, color = 'black')
    )
  
  plot_list[[i]] <- p
}


# 전체 subplot 구성
p2a <- plot_grid(plotlist = plot_list, ncol = 3)

# 예: 추가할 제목
title_text <- "Phenotypic differences by TF carrier combinations"

# 제목과 p2a를 함께 묶기
p2a_with_title <- plot_grid(
  ggdraw() + 
    draw_label(title_text, x = 0, hjust = 0, size = 10),
  p2a,
  ncol = 1,
  rel_heights = c(0.05, 1)  # 제목:그림 비율
)


# 최종 병합
plot_grid(p1, p2a_with_title, ncol = 2, rel_widths = c(1, 1.1),
          labels = c('a', 'b'), label_size = 10)

# 저장
ggsave(paste0(fig_path, paste(version, 'pdf', sep = '.')),
       last_plot(),
       width = 21.61 * 1.75, height = 27.9 * 0.6,
       limitsize = FALSE, units = "cm")
