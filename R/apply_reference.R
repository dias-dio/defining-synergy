#' Figure codes for the paper: A statistical framework for defining synergistic anticancer drug interactions 
#' Written by Diogo Dias <diogo.dias@helsinki.fi>, September 2025
#'
#' @description Direct application of the reference null synergy score distributions per cancer type 
#'
#' Functions in this file:
#' compute_empirical_p: compute one-sided empirical p-value for a synergy score using a tissue-specific reference distribution
#'   
#' calculate_pval: derivation of empirical p-values retaining the original information for downstream analyses  
#'  
#' @param references_data Data frame or list with the reference null distributions
#' @param method Character, synergy metric ("ZIP","Bliss","HSA","Loewe").  
#' @param type Character, cancer type ("Breast","Colon","Pancreas").  
#' @param scores Numeric vector of synergy scores to evaluate.  
#' @return Data frame with columns: scores, type, method, p-value, and -log10(p-value). 
#'
#' @import readxl
#' @import dplyr
#' @import openxlsx
#' @import ggplot2
#'

#### Load all required libraries (PS: Install if any is missing).
pkgs <- c("dplyr", "readxl", "openxlsx","ggplot2")  
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))


#### Load dataset with the reference (null) synergy scores distributions based on the Jaaks et al. dataset
zip_results   <- readxl::read_excel("Data/ZIP_results.xlsx")
bliss_results <- readxl::read_excel("Data/Bliss_results.xlsx")
hsa_results   <- readxl::read_excel("Data/HSA_results.xlsx")
loewe_results <- readxl::read_excel("Data/Loewe_results.xlsx")

#### Load example synergy score data from the Bashi et al. breast cancer dataset
example_results <- readxl::read_excel("Data/example_results.xlsx") 

#### Empirical p-value function
compute_empirical_p <- function(scores, score) {
    # drop NA/Inf in reference scores
    s <- scores[is.finite(scores)]
    N <- length(s)
    if (!is.finite(score) || N == 0L) return(NA_real_)
    if (score >= 0) {
      p.val <- sum(sort(s) >= score) / N
    } else {
      p.val <- sum(sort(s) <= score) / N
    }
    # ensure non-zero p-value
    if (p.val == 0) p.val <- 1 / N
    p.val
}

# Reference datasets 
references_data <- list(ZIP = zip_results, BLISS = bliss_results, HSA = hsa_results, LOEWE = loewe_results)

# Helper function to calculate empirical p-values based on:
# i) Synergy model (ZIP, BLISS, HSA, LOEWE)
# ii) Tissue (Breast, Colon, Pancreas)
# iii) Synergy scores

# Main helper function to derive empirical p-values for any given synergy model, type (i.e., tissue), and synergy scores
calculate_pval <- function(refs, method, type, scores) {                 
  # Validate
  if (!is.list(refs) || length(refs) == 0L)
    stop("`refs` must be a non-empty named list of reference data frames!")
  if (!is.character(method) || length(method) != 1L)
    stop("`method` must be a single string (ZIP/BLISS/HSA/LOEWE).")
  if (!is.character(type) || length(type) != 1L)                         
    stop("`type` (cancer type) must be a single string (e.g., Breast/Colon/Pancreas).")
  if (!is.numeric(scores) || length(scores) == 0L)
    stop("`scores` must be a non-empty numeric vector.")
  
  m    <- toupper(method)
  have <- toupper(names(refs))
  if (!(m %in% have))
    stop("Unknown method: ", method, ". Available: ", paste(names(refs), collapse = ", "))
  
  allowed_types <- c("Breast","Colon","Pancreas")                        
  t <- trimws(type)                                                     
  if (!(t %in% allowed_types))
    stop("`type` must be one of: ", paste(allowed_types, collapse = ", "))
  
  ref_df <- refs[[match(m, have)]]
  if (!all(c("type","Synergy.score") %in% names(ref_df)))               
    stop("Reference for ", m, " must have columns `type` and `Synergy.score`.")
  
  ref_vec <- ref_df$Synergy.score[trimws(ref_df$type) == t]              
  ref_vec <- ref_vec[is.finite(ref_vec)]
  if (!length(ref_vec))
    stop("No reference scores for type '", t, "' in method ", m, ".")
  
  # Compute p-values
  emp <- vapply(scores, function(x) compute_empirical_p(ref_vec, x), numeric(1))
  
  data.frame(
    Synergy.score = scores,
    type   = t,                                                         
    Method = m,
    Pval = emp,
    Log10_pval = -log10(emp),
    row.names = NULL
  )
}

#### User input
# 1) Choose cancer type  
types <- c("Breast","Colon","Pancreas")
chosen_type   <- "Breast" ##### Change here <-----------------------------------

# 2) Choose synergy model
methods <- c("ZIP","BLISS","HSA","LOEWE")
chosen_method <- "ZIP" ###### Change here <-----------------------------------

# 3) Synergy scores to be evaluated under the reference (null) distribution per tissue
#scores_vec <- c(-21.1,-14.9, -9.6, -7,7, -4.3, -0.6, 0, 0.2, 2.4, 8.9, 10.1, 15.9, 25.3) # A simple vector for example purposes
scores_vec <- example_results %>% dplyr::filter(type == chosen_type) %>% dplyr::pull(Synergy.score) 

#### Function call
# Calculate empirical p-values for the user synergy scores
results <- calculate_pval(references_data, chosen_method, chosen_type, scores_vec)

#### User output
# Results
#print(results) # Optional

# View all the top synergistic and antagonistic drug combinations from the synergy scores and empirical p-values (no aggregated results)
col_volcano_synergy <- c("Synergistic" = "#B2182B", "Antagonistic" = "#2166AC", "Other" = "#C1C1C1")

results_synergy <- results %>%
  mutate(cat = case_when(
    Synergy.score >=  10 & Log10_pval >= 2 ~ "Synergistic",
    Synergy.score <= -10 & Log10_pval >= 2 ~ "Antagonistic",
    TRUE                                   ~ "Other"
  ))

volcano_plot_all <- ggplot(results_synergy, aes(Synergy.score, Log10_pval)) +
  geom_point(aes(color = cat), shape = 19, size = 4, alpha = 0.9) + 
  geom_hline(yintercept = 2,  linetype = "dashed", color = "grey20", linewidth = 1) +
  geom_vline(xintercept = -10,  linetype = "dashed", color = "grey20", linewidth = 1) +
  geom_vline(xintercept = 10,  linetype = "dashed", color = "grey20", linewidth = 1) +
  scale_color_manual(values = col_volcano_synergy, drop = FALSE, guide = "none") +
  labs(x = "Synergy score", y = "-log10 (P-value)") +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24, color = "black"),
    axis.text  = element_text(size = 22, color = "black"),
    axis.ticks = element_line(size = 1.1),
    axis.ticks.length = unit(.40, "cm"),
    legend.position = "none",
    panel.spacing.x = unit(1, "cm")
  )

print(volcano_plot_all)


#### Store the results
write.xlsx(results, "C:/Users/diogo/OneDrive/Ambiente de Trabalho/PhD/Tero Group/DATA/Supplementary_Files/results.xlsx", rowNames = FALSE)




