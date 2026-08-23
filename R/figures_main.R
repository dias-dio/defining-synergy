#' Figure codes for the paper: A statistical framework for defining synergistic anticancer drug interactions 
#' Written by Diogo Dias <diogo.dias@helsinki.fi>, September 2025
#' 
#' @description Code for the main and supplementary figures using the synergy results from the Jaaks et al. dataset. The user can either choose: ZIP, Bliss, HSA, or Loewe scores.
#'
#' @import readxl
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ggpubr
#' @import stringr
#' @import viridisLite
#' @import grid
#' @import rlang
#' @import tidytext
#' @import forcats
#' 


#### Load all required libraries (PS: Install if any is missing).
### Required packages
pkgs <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggpubr",
  "stringr",
  "viridisLite",
  "grid",
  "rlang",
  "tidytext",
  "forcats"
)

### Install missing packages only
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  install.packages(to_install)
}

### Load packages
invisible(lapply(pkgs, library, character.only = TRUE))

#### Load data
Tissues <- c("Breast", "Colon", "Pancreas")

zip_bh_results <- bind_rows(
  lapply(Tissues, function(sheet) {
    read_excel("Data/ZIP_bh_pval_results.xlsx", sheet = sheet) %>%
      mutate(Tissue = sheet)
  })
)

zip_results   <- read_excel("Data/ZIP_results.xlsx")
bliss_results <- read_excel("Data/Bliss_results.xlsx")
hsa_results   <- read_excel("Data/HSA_results.xlsx")
loewe_results <- read_excel("Data/Loewe_results.xlsx")


#### Choose here which data you would like to use for the figures (i.e., zip_results, bliss_results, hsa_results, or loewe_results)
results_df <- zip_results ###### <------------ Change here the main synergy results file

### Define color palettes for the plots
col_trio_border <- c("Breast" = "#8815D4", "Colon" = "#D4B715", "Pancreas" = "#34857E")
col_volcano <- c("Breast" = "#8815D4", "Colon" = "#D4B715", "Pancreas" = "#34857E", "ns" = "#C1C1C1")
bar_colors <- c("Synergistic" = "#B2182B","Antagonistic" = "#2166AC","Non-significant" = "#C1C1C1")
role_colors <- c(ANCHOR  = "#D95F0E", LIBRARY = "#1C5C85")
synergistic_colors <- rev(rocket(30))[2:28]  
antagonistic_colors <- mako(30)[2:28]


######## Main Figures

####### Main Figure 2
#### Main Figure 2.c) - Reference (null hypothesis) distributions of the ZIP scores across the three cancer tissues
### Function to create a density plot for each tissue
plot_density_tissue <- function(df, tissue_name) {
  ggplot(df %>% filter(type == tissue_name), aes(x = Synergy.score)) +
    geom_density(fill = col_trio_border[[tissue_name]], color = "black", alpha = 0.7, linewidth = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey20", linewidth = 1) +
    coord_cartesian(xlim = c(-20, 20), ylim = c(0, 0.13)) +  # Anchor to y=0
    scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), expand = c(0, 0)) + # Remove vertical whitespace
    labs(title = NULL, x = "Synergy score", y = "Density") +
    theme_classic(base_family = "Arial") +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 24, color = "black"),
      axis.text = element_text(size = 22, color = "black"),
      axis.ticks = element_line(linewidth = 1.1),
      axis.ticks.length=unit(.40, "cm"),
      legend.position = "none"
    )
}

plot_breast <- plot_density_tissue(results_df, "Breast")
plot_colon  <- plot_density_tissue(results_df, "Colon")
plot_pancreas <- plot_density_tissue(results_df, "Pancreas")


plot_breast; plot_colon; plot_pancreas


#### Main Figure 2.d) - Violin distributions of log-transformed p-values across the tissue types
tissue_IQR <- results_df %>%
  group_by(type) %>%
  summarize(
    Q1 = quantile(log10_pval, 0.25, na.rm = TRUE),
    Median = quantile(log10_pval, 0.5, na.rm = TRUE),
    Q3 = quantile(log10_pval, 0.75, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(cols = c(Q1, Median, Q3), names_to = "stat", values_to = "log10_pval")

zip_violin <- ggplot(results_df, aes(x = type, y = log10_pval, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey20", linewidth = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black", size = 0.4) +
  geom_point(
    data = tissue_IQR,
    aes(x = type, y = log10_pval, fill = type),
    shape = 21, color = "black", size = 3
  ) +
  scale_fill_manual(values = col_trio_border) +  # same for violins and dots
  labs(x = NULL, y = "-log10 (P-value)") +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 1.1),
    axis.ticks.length=unit(.40, "cm"),
    legend.position = "none"
  )

zip_violin <- zip_violin +
  stat_compare_means(
    comparisons = list(
      c("Breast", "Colon"),
      c("Breast", "Pancreas"),
      c("Colon", "Pancreas")
    ),
    method = "wilcox.test", 
    label = "p.signif", 
    size = 7
  ) +
  geom_text(
    data = results_df %>%
      group_by(type) %>%
      summarise(n = n()),
    aes(x = type, y = max(results_df$log10_pval, na.rm = TRUE) + 0.5,
        label = paste0("n = ", n)),
    size = 7, family = "Arial", color = "black"
  )

print(zip_violin)


#### Main Figure 2.f - Proportion of synergistic, antagonistic, and non-significant drug combinations based on ZIP synergy scores and empirical p-value
### Categorize synergistic, antagonistic, and non-significant
category_tissue_df <- results_df %>%
  mutate(
    Category = case_when(
      Synergy.score >= 10 & log10_pval >= 2  ~ "Synergistic",
      Synergy.score <= -10 & log10_pval >= 2 ~ "Antagonistic",
      TRUE                                   ~ "Non-significant"
    )
  )

### Compute the proportions of the categories 
signif_bar_data <- category_tissue_df %>%
  count(type, Category) %>%
  group_by(type) %>%
  mutate(Percentage = 100 * n / sum(n)) %>%
  ungroup()

signif_bar_data$Category <- factor(signif_bar_data$Category, levels = c("Synergistic", "Antagonistic", "Non-significant"))

figure_2f <- ggplot(signif_bar_data, aes(x = type, y = Percentage, fill = Category)) +
  geom_col(position = position_dodge2(width = 0.7),width = 0.7) +
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.06)) 
  ) +
  labs(x = NULL, y = "Combinations (%)", fill = NULL) +
  theme_classic(base_family = "Arial") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.5, colour = "grey90"),
    axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 6)),
    axis.text.y = element_text(size = 22, colour = "black"),
    axis.title.y = element_text(size = 24, colour = "black"),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks.y = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.4, "cm"),
    legend.position = "none"
  ) 

print(figure_2f)


####### Main Figure 3
#### Main Figure 3.a) - Volcano plots of ZIP synergy scores and empirical p-values for breast, colon, and pancreatic cancer cell lines
volcano_df <- results_df %>%
  mutate(
    logs = dplyr::coalesce(log10_pval, -log10(pval)),
    zips = Synergy.score,
    plot_type = if_else(logs >= 2 & abs(zips) >= 10, type, "Not-significant")
  )

zip_volcano <- ggplot(volcano_df, aes(x = zips, y = logs, color = plot_type, fill = plot_type)) +
  geom_point(size = 3, stroke = 1) +
  geom_vline(xintercept = c(-10, 10), linetype = "dashed", color = "grey20", linewidth = 0.7) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "grey20", linewidth = 0.7) +
  scale_color_manual(values = col_volcano) +
  scale_fill_manual(values = col_volcano) +
  scale_x_continuous(n.breaks = 8) + 
  scale_y_continuous(limits = c(0, NA),             
                     n.breaks = 6,
                     expand = expansion(mult = c(0, .05))) +
  labs(x = "Synergy score", y = "-log10 (P-value)") +
  theme_classic(base_family = "Arial") +  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24, family = "Arial", color = "black"),
    axis.text = element_text(size = 22, family = "Arial", color = "black"),
    axis.ticks = element_line(linewidth = 1.1),
    axis.ticks.length = unit(.40, "cm"),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )

print(zip_volcano)


#### Main Figure 3.d) - Number of significant synergistic drug combinations (ZIP ≥ 10, p ≤ 0.01) grouped by biological pathway after combining anchor and library drugs
### Rename paths for better visualization
pathway_rename <- c(
  "Protein stability and degradation" = "Protein stability",
  "Chromatin histone acetylation"     = "Histone acetylation"
)

### Summarize per role within a tissue (ANCHOR or LIBRARY)
### Threshold: ZIP >= 10 and p <= 0.01  (i.e., -log10(p) >= 2)
summarize_pathways_role <- function(df, pathway_type) {
  
  pathway_sym <- rlang::sym(pathway_type)
  
  df %>%
    dplyr::filter(is.finite(Synergy.score), is.finite(pval)) %>%
    dplyr::filter(Synergy.score >= 10, pval <= 0.01) %>%
    dplyr::group_by(!!pathway_sym) %>%
    dplyr::summarise(
      mean_zip = mean(Synergy.score, na.rm = TRUE),
      n_hits   = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::rename(pathway = !!pathway_sym) %>%
    dplyr::mutate(
      role = sub("\\.PATHWAY$", "", pathway_type),
      role = toupper(role)
    ) %>%
    dplyr::filter(!is.na(pathway), pathway != "")
}

### Combine anchor + library per tissue (collapsed to pathway level)
prep_combined_tissue <- function(df) {
  
  dplyr::bind_rows(
    summarize_pathways_role(df, "ANCHOR.PATHWAY"),
    summarize_pathways_role(df, "LIBRARY.PATHWAY")
  ) %>%
    dplyr::mutate(pathway = dplyr::recode(pathway, !!!pathway_rename)) %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(
      mean_zip = mean(mean_zip, na.rm = TRUE),
      n_hits   = sum(n_hits, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    dplyr::filter(is.finite(mean_zip), is.finite(n_hits))
}

#### Build datasets
ZIP_breast <- results_df %>% filter(type == "Breast")
ZIP_colon <- results_df %>% filter(type == "Colon")
ZIP_pancreas <- results_df %>% filter(type == "Pancreas")

comb_B <- prep_combined_tissue(ZIP_breast)
comb_C <- prep_combined_tissue(ZIP_colon)
comb_P <- prep_combined_tissue(ZIP_pancreas)

### Global y max 
global_ymax_freq <- max(c(comb_B$n_hits, comb_C$n_hits, comb_P$n_hits), na.rm = TRUE)
if (!is.finite(global_ymax_freq)) global_ymax_freq <- 1

global_ymax_plot <- 120  
angle_x <- 45            

### Plot function: order by frequency (n_hits) desc; label mean ZIP above bars
plot_freq_with_meanzip_labels <- function(dat, y_max, angle = 45) {
  
  dat2 <- dat %>%
    dplyr::mutate(pathway = as.character(pathway)) %>%
    dplyr::arrange(dplyr::desc(n_hits), dplyr::desc(mean_zip)) %>%
    dplyr::mutate(pathway = factor(pathway, levels = unique(pathway)))
  
  ggplot(dat2, aes(x = pathway, y = n_hits)) +
    geom_col(
      width = 0.70,
      fill = "#1B6F6A",
      colour = "black",
      linewidth = 0.35
    ) +
    geom_label(
      aes(label = sprintf("%.1f", mean_zip)),
      vjust = -0.25,
      size = 11,
      family = "Arial",
      colour = "black",
      fill = "white",
      label.size = 0,
      label.padding = unit(0.16, "lines")
    ) +
    scale_x_discrete(expand = c(0.02, 0.02)) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, by = 20),
      expand = c(0, 2)   # prevents top clipping
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        size = 36,
        angle = angle,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.line  = element_line(linewidth = 1.0, colour = "black"),
      axis.ticks = element_line(linewidth = 1.0, colour = "black"),
      axis.ticks.length = unit(0.28, "cm"),
      panel.grid.major.y = element_line(colour = "grey75", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 18, r = 8, b = 22, l = 14)
    )
}

#### Plots
p_B_freq <- plot_freq_with_meanzip_labels(comb_B, y_max = global_ymax_plot, angle = angle_x)
p_C_freq <- plot_freq_with_meanzip_labels(comb_C, y_max = global_ymax_plot, angle = angle_x)
p_P_freq <- plot_freq_with_meanzip_labels(comb_P, y_max = global_ymax_plot, angle = angle_x)

p_B_freq; p_C_freq; p_P_freq



#### Main Figure 3.e) - Significant PI3K/mTOR-targeting drug combinations (ZIP ≥ 10, p ≤ 0.01) across breast, colon, and pancreatic cancer cell lines
pi3k_data <- results_df %>%
  mutate(logs = dplyr::coalesce(log10_pval, -log10(pval))) %>%
  filter(Synergy.score >= 10, logs >= 2) %>%
  filter(if_any(c(ANCHOR.PATHWAY, LIBRARY.PATHWAY),
                ~ .x == "PI3K/MTOR signaling")) %>%
  mutate(Tissue = type) 


pi3k_results <- ggplot(pi3k_data, aes(x = Synergy.score, y = log10_pval, color = type)) +
  geom_point(size = 8, alpha = 0.75) +
  scale_color_manual(values = col_trio_border) +
  labs(title = NULL, y = "-log10 (P-value)", x = "Synergy score") +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24, family = "Arial", color = "black"),
    axis.text = element_text(size = 22, family = "Arial", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 1.1),
    axis.ticks.length=unit(.40, "cm"),
    legend.position = "none"
  )

print(pi3k_results)



######## Supplementary Figures

#### Supplementary Figure 8 - Tissue-specific reference distributions of Bliss, HSA, and Loewe synergy scores
### Use the same function as in Main Figure. 2c)
plot_breast_bliss <- plot_density_tissue(bliss_results, "Breast")
plot_colon_bliss  <- plot_density_tissue(bliss_results, "Colon")
plot_pancreas_bliss <- plot_density_tissue(bliss_results, "Pancreas")

plot_breast_hsa <- plot_density_tissue(hsa_results, "Breast")
plot_colon_hsa  <- plot_density_tissue(hsa_results, "Colon")
plot_pancreas_hsa <- plot_density_tissue(hsa_results, "Pancreas")

### Please note: change x-limits of the function between -60 and 60 for Loewe
plot_breast_loewe <- plot_density_tissue(loewe_results, "Breast")
plot_colon_loewe  <- plot_density_tissue(loewe_results, "Colon")
plot_pancreas_loewe <- plot_density_tissue(loewe_results, "Pancreas")

plot_breast_bliss; plot_colon_bliss; plot_pancreas_bliss

plot_breast_hsa; plot_colon_hsa; plot_pancreas_hsa

### Again, please note: change x-limits of the function between -60 and 60 for Loewe
plot_breast_loewe; plot_colon_loewe; plot_pancreas_loewe 


#### TODO: align figure number with the final manuscript (this label collides
#### with a Supplementary Figure 10 in figures_external.R)
#### Supplementary Figure 10 - Target pathway prioritization under Benjamini-Hochberg (BH) FDR-controlled empirical significance
### Set threshold and data
q_main <- 0.20

df_bh_all <- zip_bh_results %>%
  filter(is.finite(Synergy.score), is.finite(p_bh_one)) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Breast","Colon","Pancreas")),
    sig = p_bh_one <= q_main,
    neglog10_q = -log10(pmax(p_bh_one, 1e-300)),
    plot_group = ifelse(sig, as.character(Tissue), "ns"),
    plot_group = factor(plot_group,
                        levels = c("Breast","Colon","Pancreas","ns"))
  )

# BH-adjusted p-values volcano plot for the key targeted pathways
p_volc_bh_adjusted_all <- ggplot(df_bh_all, aes(x = Synergy.score, y = neglog10_q)) +
  geom_point(aes(color = plot_group), size = 2.6, alpha = 0.80) +
  scale_color_manual(values = col_volcano, drop = FALSE) +
  geom_hline(
    yintercept = -log10(q_main),
    linetype = 2,
    linewidth = 1.1,
    color = "black"
  ) +
  scale_x_continuous(
    limits = c(-40, 40),
    breaks = seq(-40, 40, by = 10)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(title = NULL, y = "-log10 (BH-adjusted p-value)", x = "Synergy score") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, family = "Arial", color = "black"),
    axis.text = element_text(size = 22, family = "Arial", color = "black"),
    axis.ticks = element_line(linewidth = 1.1),
    axis.ticks.length = unit(.40, "cm"),
    plot.margin = margin(20, 30, 20, 30),
    legend.position = "none"
  )
  
p_volc_bh_adjusted_all


#### Supplementary Figure 12 - Cross-tissue enrichment of ERK/MAPK-targeting synergistic drug combinations
### Make data
erk_data <- results_df %>%
  mutate(logs = dplyr::coalesce(log10_pval, -log10(pval))) %>%
  filter(Synergy.score >= 10, logs >= 2) %>%
  filter(if_any(c(ANCHOR.PATHWAY, LIBRARY.PATHWAY),
                ~ .x == "ERK MAPK signaling")) %>%
  mutate(Tissue = type) 


erk_results <- ggplot(erk_data, aes(x = Synergy.score, y = log10_pval, color = type)) +
  geom_point(size = 8, alpha = 0.75) +
  scale_color_manual(values = col_trio_border) +
  labs(title = NULL, y = "-log10 (P-value)", x = "Synergy score") +
  theme_classic(base_family = "Arial") +
    theme(
      axis.title = element_text(size = 24, family = "Arial", color = "black"),
      axis.text = element_text(size = 22, family = "Arial", color = "black"),
      axis.ticks = element_line(linewidth = 1.1),
      axis.ticks.length = unit(.40, "cm"),
      legend.position = "none",
      panel.spacing = unit(0, "lines")
    )

print(erk_results)


#### Supplementary Figure 13 - Top-10 ranked synergistic and antagonistic drug combinations across cancer tissues
df_significant <- results_df %>%
  filter(log10_pval >= 2) %>%
  mutate(Drug.combination = sub(" - ", " + ", Drug.combination, fixed = TRUE))

### Top 10 synergists per tissue
top_synergisitc <- df_significant %>%
  filter(Synergy.score >= 10) %>%
  group_by(type, Drug.combination) %>%
  summarise(mean_synergy = mean(Synergy.score, na.rm = TRUE), .groups = "drop") %>%
  group_by(type) %>% slice_max(mean_synergy, n = 10, with_ties = FALSE) %>% ungroup()

### Top 10 antagonists per tissue
top_antagonistic <- df_significant %>%
  filter(Synergy.score <= -10) %>%
  group_by(type, Drug.combination) %>%
  summarise(mean_synergy = mean(Synergy.score, na.rm = TRUE), .groups = "drop") %>%
  group_by(type) %>% slice_min(mean_synergy, n = 10, with_ties = FALSE) %>% ungroup()

### Plot top 10
plot_top <- function(dat, cols, reverse_order = FALSE) {
  dat <- dat %>%
    mutate(x_lab = tidytext::reorder_within(Drug.combination, mean_synergy, type)) %>%
    { if (reverse_order) dplyr::mutate(., x_lab = forcats::fct_rev(x_lab)) else . }
  ggplot(dat, aes(x_lab, mean_synergy, fill = mean_synergy)) +
    geom_col(width = 0.8) +
    coord_flip() +
    facet_wrap(~ type, scales = "free_y") +
    tidytext::scale_x_reordered() +
    scale_fill_gradientn(colors = cols) +
    labs(title = NULL, y = "Average synergy", x = "Drug combination") +
    theme_classic(base_family = "Arial") +
    theme(
      axis.title = element_text(size = 24, family = "Arial", color = "black"),
      axis.text = element_text(size = 22, family = "Arial", color = "black"),
      axis.ticks = element_line(linewidth = 1.1),
      axis.ticks.length = unit(.40, "cm"),
      legend.position = "none",
      panel.spacing = unit(0, "lines")
    )
}

p_syn <- plot_top(top_synergisitc, synergistic_colors, reverse_order = FALSE)
p_ant <- plot_top(top_antagonistic, antagonistic_colors, reverse_order = TRUE)

p_syn; p_ant


##### End
