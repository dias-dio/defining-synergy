#' Figure codes for the paper: A statistical framework for defining synergistic anticancer drug interactions 
#' Written by Diogo Dias <diogo.dias@helsinki.fi>, September 2025
#' 
#' @description Code for the main and supplementary figures using the synergy results from the Bashi et al. dataset using ZIP synergy scores
#' @import readxl
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' 


#### Load all required libraries (PS: Install if any is missing)
pkgs <- c("readxl", "dplyr", "ggplot2", "ggrepel")  
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

#### Load data (ZIP results)
jaaks_results = read_excel("Data/ZIP_results.xlsx") 
bashi_results = read_excel("Data/Bashi_results.xlsx") 

### Define color palettes for the plots
col_trio_border <- c("Breast" = "#8815D4", "Colon" = "#D4B715", "Pancreas" = "#34857E")
col_scatter <- c("Top 15 synergists" = "#B2182B","Bottom 15 antagonists" = "#2166AC", "Other" = "#C1C1C1")

######## Main Figures
####### Main Figure 4
#### Main Figure 4.b) - Comparison of the ZIP synergy score distributions between the Jaaks et al. and Bashi et al. datasets for each cancer type
# Results from the Jaaks et al. dataset
jaaks_results_2 <- jaaks_results %>%
  transmute(type, Synergy.score, screen = "Jaaks et al. (2x7)")
# Results from the Bashi et al. dataset
bashi_results_2 <- bashi_results %>%
  transmute(
    type = recode(type,
                  "Breast Carcinoma"     = "Breast",
                  "Pancreatic Carcinoma" = "Pancreas",
                  "Colorectal Carcinoma" = "Colon"
    ),
    Synergy.score,
    screen = "Bashi et al. (7x7)"
  )

df_all_zip <- bind_rows(jaaks_results_2, bashi_results_2) %>%
  filter(type %in% c("Breast","Colon","Pancreas")) %>%
  mutate(
    type   = factor(type, c("Breast","Colon","Pancreas")),
    screen = factor(screen, c("Jaaks et al. (2x7)", "Bashi et al. (7x7)"))
  )

# Color palette 
cols <- setNames(c("#b2171a", "#155a8a"), levels(df_all_zip$screen))

xlim <- ceiling(max(abs(df_all_zip$Synergy.score), na.rm = TRUE))

# Density plots comparison between the two datasets
density_bashi_plot <- ggplot(df_all_zip, aes(Synergy.score, fill = screen, color = screen)) +
  geom_density(colour = "black", alpha = 0.7, linewidth = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey20", linewidth = 1) +
  coord_cartesian(xlim = c(-xlim, xlim)) +
  scale_y_continuous(limits = c(0, 0.13), expand = expansion(mult = c(0, 0))) +  
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  labs(x = "Synergy score (ZIP)", y = "Density") +
  theme_classic(base_family = "Arial") +  
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    axis.ticks = element_line(size = 1.1),
    axis.ticks.length = unit(.40, "cm"),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )

print(density_bashi_plot)


# ECDF plots comparison between the two datasets
ecdf_bashi_plot <- ggplot(df_all_zip, aes(Synergy.score, color = screen)) +
  stat_ecdf(size =  1.15) +
  coord_cartesian(xlim = c(-xlim, xlim)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0, 0))) +                
  scale_color_manual(values = cols) +
  facet_wrap(~ type, ncol = 1) +
  theme_classic(base_family = "Helvetica Neue") +
  labs(x = "Synergy score (ZIP)", y = "Cumulative distribution (ECDF)") +
  theme_classic(base_family = "Arial") +  
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    axis.ticks = element_line(size = 1.1),
    axis.ticks.length = unit(.40, "cm"),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )

print(ecdf_bashi_plot)


# Kolmogorov-Smirnov normality test between the Jaaks et al. and Bashi et al. datasets
tissues <- c("Breast", "Pancreas", "Colon")

ks_results <- function(df, type) {
  d <- df[df$type == type, ]
  a <- d$Synergy.score[d$screen == "Jaaks et al. (2x7)"]
  b <- d$Synergy.score[d$screen == "Bashi et al. (7x7)"]
  if (length(a) && length(b)) suppressWarnings(stats::ks.test(a, b)$p.value) else NA_real_
}

ks_summary <- data.frame(
  Tissue     = tissues,
  KS_p_value = sapply(tissues, ks_results, df = df_all_zip)
)

print(ks_summary)


#### Main Figure 4.e) - Volcano plot of average ZIP scores against the average p-values across all cancer types
# summarize per tissue × combination
bashi_summary <- bashi_results %>%
  dplyr::filter(pval > 0) %>%                               
  dplyr::group_by(type, Drug.combination) %>%
  dplyr::summarise(
    mean_zip   = mean(Synergy.score, na.rm = TRUE),
    mean_log10 = mean(-log10(pval),  na.rm = TRUE),
    n_cells    = dplyr::n_distinct(cell),
    .groups    = "drop"
  )


# color categories (top/bottom 15 per tissue)
scatter_df <- bashi_summary %>%
  group_by(type) %>%
  mutate(cat = case_when(
    dense_rank(desc(mean_zip)) <= 15 ~ "Top 15 synergists",
    dense_rank(mean_zip) <= 15 ~ "Bottom 15 antagonists",
    TRUE ~ "Other"
  )) %>% ungroup()

# label only top 3 synergists + top 3 antagonists per tissue
label_df <- scatter_df %>%
  group_by(type) %>%
  filter(dense_rank(desc(mean_zip)) <= 3 | dense_rank(mean_zip) <= 3) %>%
  ungroup()

# Same axis ranges for all facets (i.e., tissues)
lims <- summarise(scatter_df,
                  xmin = floor(min(mean_zip,   na.rm = TRUE)),
                  xmax = ceiling(max(mean_zip, na.rm = TRUE)),
                  ymin = floor(min(mean_log10, na.rm = TRUE)),
                  ymax = ceiling(max(mean_log10, na.rm = TRUE))
)

scatter_plot_all_tissues <- ggplot(scatter_df, aes(mean_zip, mean_log10)) +
  geom_point(aes(fill = cat, size = n_cells), shape = 21, color = "black", alpha = 0.9) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "grey20", linewidth = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey20", linewidth = 1) +
  ggrepel::geom_text_repel(
    data = label_df, aes(label = Drug.combination),
    size = 6, family = "Arial", color = "black"
  ) +
  facet_wrap(~ type, nrow = 1, strip.position = "top") +                #
  scale_fill_manual(values = col_scatter) +
  scale_size(range = c(3, 13), guide = "none") +
  scale_x_continuous(limits = c(lims$xmin-3, lims$xmax)) +
  scale_y_continuous(limits = c(lims$ymin, lims$ymax)) +
  labs(x = "Average synergy score (ZIP)", y = "Average -log10 (P-value)") +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24, color = "black"),
    axis.text  = element_text(size = 22, color = "black"),
    strip.text.x    = element_text(hjust = 0.5)  ,
    axis.ticks = element_line(size = 1.1),
    axis.ticks.length = unit(.40, "cm"),
    legend.position = "none",
    panel.spacing.x = unit(1, "cm"),       
  )

print(scatter_plot_all_tissues)


######## Supplementary Figures
#### Supplementary Figure 10 - Boxplots showing the synergy score distributions across the shared breast, colon, and pancreatic cancer cell 
# lines shared between the Jaaks et al. and Bashi et al. datasets.
# Tissues in fixed order
tissues <- c("Breast","Colon","Pancreas")

# Colors
cols <- c("Jaaks et al." = "#b2171a",
          "Bashi et al." = "#155a8a")

# Build per-tissue shared-cell data (no purrr)
make_plot_df <- function(t) {
  j <- dplyr::filter(jaaks_results, type == t)
  b <- dplyr::filter(bashi_results, type == t)
  shared <- intersect(j$cell, b$cell)
  
  dplyr::bind_rows(
    dplyr::filter(j, cell %in% shared) %>%
      dplyr::transmute(type = factor(t, levels = tissues),
                       Source = "Jaaks et al.",
                       score = Synergy.score),
    dplyr::filter(b, cell %in% shared) %>%
      dplyr::transmute(type = factor(t, levels = tissues),
                       Source = "Bashi et al.",
                       score = Synergy.score)
  )
}

plot_df <- dplyr::bind_rows(lapply(tissues, make_plot_df)) %>%
  dplyr::mutate(Source = factor(Source, c("Jaaks et al.", "Bashi et al.")))

# KS p-values per tissue
ks_tbl <- plot_df %>%
  dplyr::group_by(type) %>%
  dplyr::reframe(KS_p_value = {
    s <- split(score, Source)
    if (length(s) == 2 && all(lengths(s) > 0))
      suppressWarnings(stats::ks.test(s[[1]], s[[2]])$p.value)
    else NA_real_
  })

# Select extreme points per Source × tissue as "outliers"
n_out <- 5
outliers_df <- plot_df %>%
  dplyr::group_by(type, Source) %>%
  dplyr::arrange(score, .by_group = TRUE) %>%
  dplyr::slice(c(seq_len(min(n_out, dplyr::n())),
                 seq(max(dplyr::n() - n_out + 1, 1), dplyr::n()))) %>%
  dplyr::ungroup()

# Boxplot (faceted by tissue)
boxplot_compare_all <- ggplot(plot_df, aes(Source, score, fill = Source)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "black", linewidth = 1.1, alpha = 0.75) +
  geom_point(data = outliers_df, aes(Source, score),
             colour = "black", size = 2, alpha = 0.85) +
  scale_fill_manual(values = cols) +
  facet_wrap(~ type, nrow = 1, scales = "free_y", strip.position = "top") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = "Helvetica Neue") +
  theme(
    axis.text        = element_text(size = 22, colour = "black"),
    axis.ticks       = element_line(linewidth = 1.1),
    axis.ticks.length= unit(0.4, "cm"),
    legend.position  = "none",
    panel.spacing.x  = unit(1, "cm"),
    strip.text.x     = element_text(hjust = 0.5, size = 24, face = "bold")
  )

# Show results
print(boxplot_compare_all)
print(ks_tbl)


##### End
