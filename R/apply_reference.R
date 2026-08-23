#' A statistical framework for defining synergistic anticancer drug interactions
#' Written by Diogo Dias <diogo.dias@helsinki.fi>
#'
#' @description Direct application of the tissue-specific reference (null)
#'   distributions of synergy scores to compute empirical p-values for
#'   user-provided drug combination synergy scores.
#'
#' Functions in this file:
#'   empirical_pvalues: vectorized one-sided empirical p-values for synergy
#'     scores against a reference distribution (external or leave-one-out
#'     formulation)
#'   calculate_pval: user-facing wrapper; selects the reference by synergy
#'     metric and cancer type and returns a results data frame
#'
#' @param refs Named list of reference data frames (ZIP/BLISS/HSA/LOEWE),
#'   each with columns `type` and `Synergy.score`
#' @param method Character, synergy metric ("ZIP","Bliss","HSA","Loewe")
#' @param type Character, cancer type ("Breast","Colon","Pancreas")
#' @param scores Numeric vector of synergy scores to evaluate
#' @param formulation "external" (default) for new/independent data,
#'   p = (b + 1) / (N + 1); "loo" for observations that are part of the
#'   reference itself, p = b / N with floor 1 / N
#' @return Data frame with columns: Synergy.score, type, Method, Pval,
#'   Log10_pval
#'
#' @import readxl
#' @import dplyr
#' @import ggplot2

#### Load all required libraries (PS: Install if any is missing)
pkgs <- c("dplyr", "readxl", "ggplot2")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

#### Load the reference (null) synergy score distributions (Jaaks et al. dataset,
#### self-combinations excluded)
zip_results   <- readxl::read_excel("Data/ZIP_results.xlsx")
bliss_results <- readxl::read_excel("Data/Bliss_results.xlsx")
hsa_results   <- readxl::read_excel("Data/HSA_results.xlsx")
loewe_results <- readxl::read_excel("Data/Loewe_results.xlsx")

#### Load example synergy score data (Bashi et al. breast cancer dataset)
example_results <- readxl::read_excel("Data/example_results.xlsx")


#### Reference datasets in list format for the main helper function
references_data <- list(ZIP = zip_results, BLISS = bliss_results,
                        HSA = hsa_results, LOEWE = loewe_results)

#### Vectorized one-sided empirical p-values
# For synergy scores (score >= 0), p is based on b = #{reference >= score};
# for antagonism (score < 0), b = #{reference <= score}.
# formulation = "external": p = (b + 1) / (N + 1)   (new, independent data)
# formulation = "loo":      p = b / N, floor 1 / N  (data within the reference)
empirical_pvalues <- function(scores, reference, formulation = c("external", "loo")) {
  formulation <- match.arg(formulation)
  ref <- sort(reference[is.finite(reference)])
  N <- length(ref)
  p <- rep(NA_real_, length(scores))
  ok <- is.finite(scores)
  if (N == 0L || !any(ok)) return(p)
  s <- scores[ok]
  b_ge <- N - findInterval(s, ref, left.open = TRUE)  # reference >= score
  b_le <- findInterval(s, ref)                        # reference <= score
  b <- ifelse(s >= 0, b_ge, b_le)
  p_ok <- if (formulation == "external") (b + 1) / (N + 1) else pmax(b, 1L) / N
  p[ok] <- pmin(p_ok, 1)
  p
}

#### Main helper function: empirical p-values for a given synergy metric,
#### cancer type and vector of synergy scores
calculate_pval <- function(refs, method, type, scores,
                           formulation = c("external", "loo")) {
  formulation <- match.arg(formulation)

  # Validate
  if (!is.list(refs) || length(refs) == 0L)
    stop("`refs` must be a non-empty named list of reference data frames.")
  if (!is.character(method) || length(method) != 1L)
    stop("`method` must be a single string (ZIP/Bliss/HSA/Loewe).")
  if (!is.character(type) || length(type) != 1L)
    stop("`type` (cancer type) must be a single string (Breast/Colon/Pancreas).")
  if (!is.numeric(scores) || length(scores) == 0L)
    stop("`scores` must be a non-empty numeric vector.")

  m    <- toupper(method)
  have <- toupper(names(refs))
  if (!(m %in% have))
    stop("Unknown method: ", method, ". Available: ",
         paste(names(refs), collapse = ", "))

  allowed_types <- c("Breast", "Colon", "Pancreas")
  t <- trimws(type)
  if (!(t %in% allowed_types))
    stop("`type` must be one of: ", paste(allowed_types, collapse = ", "))

  ref_df <- refs[[match(m, have)]]
  if (!all(c("type", "Synergy.score") %in% names(ref_df)))
    stop("Reference for ", m, " must have columns `type` and `Synergy.score`.")

  ref_vec <- ref_df$Synergy.score[trimws(ref_df$type) == t]
  ref_vec <- ref_vec[is.finite(ref_vec)]
  if (!length(ref_vec))
    stop("No reference scores for type '", t, "' in method ", m, ".")

  emp <- empirical_pvalues(scores, ref_vec, formulation = formulation)

  data.frame(
    Synergy.score = scores,
    type       = t,
    Method     = m,
    Pval       = emp,
    Log10_pval = -log10(emp),
    row.names  = NULL
  )
}

#### User input
# 1) Choose cancer type
types <- c("Breast", "Colon", "Pancreas")
chosen_type <- "Breast"    ###### Change here <-----------------------------------

# 2) Choose synergy model
methods <- c("ZIP", "Bliss", "HSA", "Loewe")
chosen_method <- "ZIP"     ###### Change here <-----------------------------------

# 3) Synergy scores to be evaluated against the tissue-specific reference
# scores_vec <- c(-21.1, -14.9, -9.6, -7.7, -4.3, -0.6, 0, 0.2, 2.4, 8.9, 10.1, 15.9, 25.3) # A simple vector for example purposes
scores_vec <- example_results %>%
  dplyr::filter(type == chosen_type) %>%
  dplyr::pull(Synergy.score)

#### Function call
# External formulation (default): use for any new dataset scored against the
# published reference, as in the example below
results <- calculate_pval(references_data, chosen_method, chosen_type, scores_vec)

#### User output
# print(results) # Optional

# Classify with the dual criterion (|synergy score| >= 10 and p <= 0.01)
col_volcano_synergy <- c("Synergistic" = "#B2182B", "Antagonistic" = "#2166AC",
                         "Other" = "#C1C1C1")

results_synergy <- results %>%
  mutate(cat = case_when(
    Synergy.score >=  10 & Log10_pval >= 2 ~ "Synergistic",
    Synergy.score <= -10 & Log10_pval >= 2 ~ "Antagonistic",
    TRUE                                   ~ "Other"
  ))

volcano_plot_all <- ggplot(results_synergy, aes(Synergy.score, Log10_pval)) +
  geom_point(aes(color = cat), shape = 19, size = 4, alpha = 0.9) +
  geom_hline(yintercept = 2,   linetype = "dashed", color = "grey20", linewidth = 1) +
  geom_vline(xintercept = -10, linetype = "dashed", color = "grey20", linewidth = 1) +
  geom_vline(xintercept = 10,  linetype = "dashed", color = "grey20", linewidth = 1) +
  scale_color_manual(values = col_volcano_synergy, drop = FALSE, guide = "none") +
  labs(x = "Synergy score", y = expression(-log[10]~(italic(p)*"-value"))) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 24, color = "black"),
    axis.text  = element_text(size = 22, color = "black"),
    axis.ticks = element_line(linewidth = 1.1),
    axis.ticks.length = unit(.40, "cm"),
    legend.position = "none"
  )

print(volcano_plot_all)

##### End
