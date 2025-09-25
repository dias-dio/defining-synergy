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
#' 


#### Load all required libraries (PS: Install if any is missing).
pkgs <- c("readxl","dplyr","tidyr","ggplot2","ggpubr","stringr","viridisLite","grid","rlang")  
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

#### Load data
zip_results   <- read_xl("Data/ZIP_results.xlsx")
bliss_results <- read_xl("Data/Bliss_results.xlsx")
hsa_results   <- read_xl("Data/HSA_results.xlsx")
loewe_results <- read_xl("Data/Loewe_results.xlsx")

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
# Function to create a density plot for each tissue
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
      axis.ticks = element_line(size = 1.1),
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
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 1.1),
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


#### Main Figure 2.f - Proportion of synergistic, antagonistic, and non-significant drug combinations based on ZIP synergy scores and empirical p-values
## Categorize synergistic, antagonistic, and non-significant
category_tissue_df <- results_df %>%
  mutate(
    Category = case_when(
      Synergy.score >= 10 & log10_pval >= 2  ~ "Synergistic",
      Synergy.score <= -10 & log10_pval >= 2 ~ "Antagonistic",
      TRUE                                   ~ "Non-significant"
    )
  )
## Compute the proportions of the categories 
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
    axis.title = element_text(size = 24, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    axis.ticks = element_line(size = 1.1),
    axis.ticks.length = unit(.40, "cm"),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )

print(zip_volcano)


#### Main Figure 3.d) - Biological pathway enrichment analysis of the targets of synergistic drug combinations, stratified by anchor and library roles across the three tissues
# Biological pathways
all_pathways_levels <- c(
  "ABL signaling","Apoptosis regulation","Cell cycle","Chromatin histone acetylation",
  "Chromatin other","Cytoskeleton","DNA replication","EGFR signaling",
  "ERK/MAPK signaling","Genome integrity","IGF1R signaling","JNK and p38 signaling",
  "Mitosis","Other","Other, kinases","p53 pathway",
  "PI3K/MTOR signaling","Protein stability and degradation",
  "RTK signaling","WNT signaling"
)

# Summarize pathways for both anchor and library roles based on highly synergistic combinations (n) based on ZIP >= 10 and -log10(P-value) >= 2
path_summ <- results_df %>%
  mutate(
    logs   = dplyr::coalesce(log10_pval, -log10(pval)),
    Tissue = type
  ) %>%
  filter(Synergy.score >= 10, logs >= 2) %>%
  pivot_longer(
    c(ANCHOR.PATHWAY, LIBRARY.PATHWAY),
    names_to  = "RoleCol",
    values_to = "pathway"
  ) %>%
  mutate(
    Role = if_else(RoleCol == "ANCHOR.PATHWAY", "ANCHOR", "LIBRARY"),
    pathway = case_when(
      str_detect(pathway %||% "", regex("pi3k\\s*/?\\s*m/?tor", TRUE)) ~ "PI3K/MTOR signaling",
      str_detect(pathway %||% "", regex("erk.*mapk", TRUE))             ~ "ERK/MAPK signaling",
      TRUE ~ pathway
    )
  ) %>%
  filter(!is.na(pathway) & pathway != "") %>%
  mutate(pathway = factor(pathway, levels = rev(all_pathways_levels))) %>%
  group_by(Tissue, Role, pathway) %>%
  summarise(
    mean_synergy = mean(Synergy.score, na.rm = TRUE),
    count        = dplyr::n(),
    .groups      = "drop"
  )

# Global ranges of x and y-axes for optimal plotting visualization
global_range <- range(path_summ$mean_synergy, na.rm = TRUE)
global_min_synergy <- floor(global_range[1])
global_max_synergy <- ceiling(global_range[2]) + 0.2
global_max_count   <- max(path_summ$count, na.rm = TRUE)

# Function to provide a pathway drug role-based analysis per tissue (i.e., either the agent acting as an anchor or lirary), taking into account the global average synergy scores
plot_pathways <- function(tissue_name, role_colors) {
  df <- path_summ %>%
    filter(Tissue == tissue_name) %>%
    mutate(Role = factor(Role, levels = c("ANCHOR","LIBRARY")))
  
  ggplot(df, aes(x = pathway, y = mean_synergy)) +
    geom_point(
      aes(fill = Role, colour = Role, size = count, shape = Role),
      stroke = 1, alpha = 0.9, na.rm = TRUE
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(global_min_synergy, global_max_synergy),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_shape_manual(values = c(ANCHOR = 21, LIBRARY = 23), name = "Pathway role") +
    scale_fill_manual(values = role_colors, name = "Pathway role") +
    scale_colour_manual(values = c(ANCHOR = "#333333", LIBRARY = "#333333"), guide = "none") +
    scale_size_continuous(range = c(3, 8), limits = c(1, global_max_count)) +
    labs(x = "Biological pathway", y = "Average synergy (ZIP)") +
    theme_bw(base_family = "Arial") +
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.2),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 24, colour = "black"),
      axis.text = element_text(size = 22, colour = "black"),
      axis.ticks = element_line(linewidth = 1.1),
      axis.ticks.length = unit(.40, "cm"),
      panel.spacing = unit(0, "lines")
    ) +
    scale_x_discrete(drop = FALSE)   
}

breast_pathway_plot   <- plot_pathways("Breast", role_colors)
colon_pathway_plot    <- plot_pathways("Colon", role_colors)
pancreas_pathway_plot <- plot_pathways("Pancreas", role_colors)

breast_pathway_plot; colon_pathway_plot; pancreas_pathway_plot


#### Main Figure 3.e) - Synergistic drug combinations targeting the commonly enriched PI3K/mTOR pathway in breast, colon, and pancreatic cancers
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
    axis.ticks = element_line(size = 1.1),
    axis.ticks.length=unit(.40, "cm"),
    legend.position = "none"
  )

print(pi3k_results)


######## Supplementary Figures
#### Supplementary Figure 5 - Reference (null hypothesis) distributions of Bliss, HSA, and Loewe synergy scores across distinct cancer tissues
# Use the same function as in Main Figure. 2c)
plot_breast_bliss <- plot_density_tissue(bliss_results, "Breast")
plot_colon_bliss  <- plot_density_tissue(bliss_results, "Colon")
plot_pancreas_bliss <- plot_density_tissue(bliss_results, "Pancreas")

plot_breast_hsa <- plot_density_tissue(hsa_results, "Breast")
plot_colon_hsa  <- plot_density_tissue(hsa_results, "Colon")
plot_pancreas_hsa <- plot_density_tissue(hsa_results, "Pancreas")

# Please note: change x-limits of the function between -60 and 60 for Loewe
plot_breast_loewe <- plot_density_tissue(loewe_results, "Breast")
plot_colon_loewe  <- plot_density_tissue(loewe_results, "Colon")
plot_pancreas_loewe <- plot_density_tissue(loewe_results, "Pancreas")

plot_breast_bliss; plot_colon_bliss; plot_pancreas_bliss

plot_breast_hsa; plot_colon_hsa; plot_pancreas_hsa

# Again, please note: change x-limits of the function between -60 and 60 for Loewe
plot_breast_loewe; plot_colon_loewe; plot_pancreas_loewe 


#### Supplementary Figure 7 - Synergistic drug combinations targeting the ERK/MAPK signaling pathway, enriched across breast, colon, and pancreatic cancer cell lines
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
      axis.title = element_text(size = 24, color = "black"),
      axis.text = element_text(size = 22, color = "black"),
      axis.ticks = element_line(linewidth = 1.1),
      axis.ticks.length = unit(.40, "cm"),
      legend.position = "none",
      panel.spacing = unit(0, "lines")
    )

print(erk_results)


#### Supplementary Figure 8 - Top 10 synergistic and antagonistic drug combinations in breast, colon, and pancreatic cancer cell lines from the Jaaks et al. dataset, 
# ranked by average ZIP synergy scores across all tested cell lines.

# keep only significant rows; clean combo label
df_significant <- results_df %>%
  filter(log10_pval >= 2) %>%
  mutate(Drug.combination = sub(" - ", " + ", Drug.combination, fixed = TRUE))

# top 10 synergists per tissue
top_synergisitc <- df_significant %>%
  filter(Synergy.score >= 10) %>%
  group_by(type, Drug.combination) %>%
  summarise(mean_synergy = mean(Synergy.score, na.rm = TRUE), .groups = "drop") %>%
  group_by(type) %>% slice_max(mean_synergy, n = 10, with_ties = FALSE) %>% ungroup()

# top 10 antagonists per tissue
top_antagonistic <- df_significant %>%
  filter(Synergy.score <= -10) %>%
  group_by(type, Drug.combination) %>%
  summarise(mean_synergy = mean(Synergy.score, na.rm = TRUE), .groups = "drop") %>%
  group_by(type) %>% slice_min(mean_synergy, n = 10, with_ties = FALSE) %>% ungroup()

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
      axis.title = element_text(size = 24, color = "black"),
      axis.text = element_text(size = 22, color = "black"),
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
