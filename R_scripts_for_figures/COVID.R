suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(clinfun)
  library(mixOmics)
  library(vdjremix)
})

set.seed(1974)

files <- list(
  eigenvectors = file.path(base_dir, "Eigenvectors_BCR_COVID_BCR_with_uncorrelated_features.txt"),
  metadata     = file.path(base_dir, "COMBAT_BCR_IDs2.txt")
)

# ---- Plot settings ----
recovery_groups <- c("COVID_MILD", "COVID_HCW_MILD", "HV")
recovery_colors <- c(
  "COVID_MILD" = "#EAE7B2",
  "COVID_HCW_MILD" = "#C7C4E2",
  "HV" = "#9CCBC5"
)

pls_colors <- list(
  covid_vs_health = c("health" = "#F4A300", "Covid" = "#56AE57"),
  sepsis_vs_health = c("health" = "#F4A300", "Sepsis" = "#CC99FF")
)

additional_predictors <- c(
  "Percentage_max_cluster_size..IGHG1",
  "V_gene_replacement_clonal_expansion..d5_norm",
  "J_gene_freq_by_uniq_VDJ_IGHD.IGHM_unmutated..IGHJ1",
  "J_gene_freq_by_uniq_VDJ_IGHE..IGHJ4",
  "J_gene_freq_by_uniq_VDJ_IGHE..IGHJ6",
  "J_gene_freq_by_uniq_VDJ_IGHG4..IGHJ2",
  "J_gene_freq_by_uniq_VDJ_IGHG4..IGHJ4"
)

module_predictors <- paste0("Module_", 1:29)
all_predictors <- c(module_predictors, additional_predictors)

# ---- Data loading ----
load_combat_data <- function(eigenvector_file, metadata_file) {
  eig <- read.delim(eigenvector_file, row.names = 1, check.names = FALSE)
  eig <- eig[1:78, , drop = FALSE]
  eig[] <- lapply(eig, function(x) as.numeric(as.character(x)))
  eig <- eig[order(rownames(eig)), , drop = FALSE]
  eig$Sequencing.ID <- rownames(eig)

  meta <- read.delim(metadata_file, check.names = FALSE)

  full <- left_join(eig, meta, by = "Sequencing.ID")
  rownames(full) <- full$Sequencing.ID
  full <- full[order(rownames(full)), , drop = FALSE]

  # Standardise module column names if any are stored as bare numbers.
  bare_module_cols <- names(full)[grepl("^[0-9]+$", names(full))]
  if (length(bare_module_cols) > 0) {
    names(full)[match(bare_module_cols, names(full))] <- paste0("Module_", bare_module_cols)
  }

  full %>%
    mutate(
      age_numeric = as.numeric(as.character(Age)),
      age_group = case_when(
        dplyr::between(age_numeric, 20, 40.9) ~ "20-40",
        dplyr::between(age_numeric, 41, 60.9) ~ "41-60",
        dplyr::between(age_numeric, 61, 80.9) ~ "61-80",
        age_numeric >= 81                     ~ "81+",
        TRUE                                  ~ NA_character_
      ),
      covid_status_broad = case_when(
        Source %in% c("COVID_HCW_MILD", "COVID_MILD") ~ "mild",
        Source %in% c("COVID_SEV", "COVID_CRIT")      ~ "crit",
        TRUE                                             ~ Source
      )
    )
}

combat_df <- load_combat_data(files$eigenvectors, files$metadata)
print(table(combat_df$covid_status_broad, useNA = "ifany"))

# ---- Panel D: ordered recovery trend ----
run_recovery_trend <- function(data, group_levels = recovery_groups) {
  trend_df <- data %>%
    filter(Source %in% group_levels) %>%
    mutate(Source = factor(Source, levels = group_levels, ordered = TRUE))

  candidate_features <- names(trend_df)[grepl("^Module_", names(trend_df))]
  candidate_features <- unique(c(candidate_features, additional_predictors[additional_predictors %in% names(trend_df)]))

  trend_results <- lapply(candidate_features, function(feature_name) {
    feature_df <- trend_df %>%
      select(Source, score = all_of(feature_name)) %>%
      filter(!is.na(score))

    if (n_distinct(feature_df$Source) < 2) {
      return(NULL)
    }

    jt_result <- tryCatch(
      clinfun::jonckheere.test(x = feature_df$score, g = feature_df$Source, alternative = "decreasing"),
      error = function(e) NULL
    )

    if (is.null(jt_result)) {
      return(NULL)
    }

    data.frame(
      feature = feature_name,
      p_value = jt_result$p.value,
      statistic = unname(jt_result$statistic),
      stringsAsFactors = FALSE
    )
  })

  trend_results_df <- bind_rows(trend_results)
  if (nrow(trend_results_df) == 0) {
    stop("No recovery-trend results were generated.")
  }

  trend_results_df %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR)
}

plot_recovery_feature <- function(data, feature_name, stats_df, output_file = NULL) {
  feature_fdr <- stats_df %>%
    filter(feature == feature_name) %>%
    pull(FDR)

  if (length(feature_fdr) == 0) {
    stop(paste("Feature not found in trend results:", feature_name))
  }

  plot_df <- data %>%
    filter(Source %in% recovery_groups) %>%
    mutate(Source = factor(Source, levels = recovery_groups, ordered = TRUE))

  p <- ggplot(plot_df, aes(x = Source, y = .data[[feature_name]], fill = Source)) +
    geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.8, shape = 21, color = "black") +
    geom_smooth(
      method = "lm",
      aes(group = 1),
      se = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    scale_fill_manual(values = recovery_colors) +
    labs(
      title = gsub("_", " ", feature_name),
      subtitle = paste("Jonckheere-Terpstra trend test, FDR =", format.pval(feature_fdr, digits = 3)),
      x = NULL,
      y = "Module score"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black", angle = 45, hjust = 1)
    )

  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 4, height = 6)
  }

  p
}

recovery_results <- run_recovery_trend(combat_df)
print(recovery_results)
print(recovery_results %>% filter(FDR < 0.05))

# These are the two recovery features shown in panel D of the uploaded figure.
recovery_features_for_figure <- c(
  "V_gene_replacement_clonal_expansion..d5_norm",
  "Module_18"
)

for (feature_name in recovery_features_for_figure) {
  if (feature_name %in% recovery_results$feature) {
    print(plot_recovery_feature(
      data = combat_df,
      feature_name = feature_name,
      stats_df = recovery_results,
      output_file = paste0(feature_name, "_recovery_trend.pdf")
    ))
  }
}

# ---- Panel B: sPLS-DA publication plots ----
prepare_plsda_dataset <- function(data, comparison = c("covid_vs_health", "sepsis_vs_health")) {
  comparison <- match.arg(comparison)

  prepared <- data %>%
    filter(Source != "COVID_HCW_MILD", Source != "Batch control")

  if (comparison == "covid_vs_health") {
    prepared <- prepared %>%
      mutate(
        comparison_group = case_when(
          Source %in% c("COVID_MILD", "COVID_SEV", "COVID_CRIT") ~ "Covid",
          Source == "HV"                                             ~ "health",
          TRUE                                                        ~ NA_character_
        )
      ) %>%
      filter(comparison_group %in% c("Covid", "health"))
  }

  if (comparison == "sepsis_vs_health") {
    prepared <- prepared %>%
      mutate(
        comparison_group = case_when(
          Source == "Sepsis" ~ "Sepsis",
          Source == "HV"     ~ "health",
          TRUE                ~ NA_character_
        )
      ) %>%
      filter(comparison_group %in% c("Sepsis", "health"))
  }

  available_predictors <- intersect(all_predictors, names(prepared))

  model_df <- prepared %>%
    select(comparison_group, all_of(available_predictors)) %>%
    drop_na()

  list(
    data = model_df,
    X = model_df %>% select(-comparison_group) %>% as.matrix(),
    Y = factor(model_df$comparison_group),
    predictors = available_predictors,
    comparison = comparison
  )
}

fit_splsda_model <- function(X, Y, ncomp_tune = 4, folds = 5, nrepeat = 100, cpus = 2) {
  list_keepX <- c(1:10, seq(20, 300, 10))

  initial_model <- mixOmics::splsda(X, Y, ncomp = 10)
  perf_initial <- perf(
    initial_model,
    validation = "Mfold",
    folds = folds,
    nrepeat = nrepeat,
    progressBar = FALSE,
    auc = TRUE
  )

  tune_result <- tune.splsda(
    X, Y,
    ncomp = ncomp_tune,
    validation = "Mfold",
    folds = folds,
    nrepeat = nrepeat,
    dist = "max.dist",
    measure = "BER",
    test.keepX = list_keepX,
    cpus = cpus
  )

  optimal_ncomp <- tune_result$choice.ncomp$ncomp
  optimal_keepX <- tune_result$choice.keepX[1:optimal_ncomp]

  final_model <- mixOmics::splsda(
    X, Y,
    ncomp = optimal_ncomp,
    keepX = optimal_keepX
  )

  list(
    initial_model = initial_model,
    initial_perf = perf_initial,
    tuning = tune_result,
    final_model = final_model,
    optimal_ncomp = optimal_ncomp,
    optimal_keepX = optimal_keepX
  )
}

build_publication_pls_plot <- function(final_model, group_vector, palette, output_file = NULL) {
  plot_object <- plotIndiv(
    final_model,
    comp = c(1, 2),
    group = group_vector,
    ind.names = FALSE,
    ellipse = TRUE,
    legend = TRUE,
    title = ""
  )

  p <- ggplot(plot_object$df, aes(x = x, y = y, fill = group, color = group)) +
    stat_ellipse(geom = "polygon", alpha = 0.35) +
    geom_point(shape = 21, size = 3, stroke = 1) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    labs(x = "PLS1", y = "PLS2") +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(size = 14, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )

  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 5, height = 4)
  }

  p
}

# Covid vs health
covid_pls <- prepare_plsda_dataset(combat_df, comparison = "covid_vs_health")
cat("\nCovid vs health dimensions:\n")
print(dim(covid_pls$X))
print(summary(covid_pls$Y))

covid_model <- fit_splsda_model(covid_pls$X, covid_pls$Y)
print(covid_model$optimal_ncomp)
print(covid_model$optimal_keepX)

covid_plot <- build_publication_pls_plot(
  final_model = covid_model$final_model,
  group_vector = covid_pls$Y,
  palette = pls_colors$covid_vs_health,
  output_file = "health_vs_covid.pdf"
)
print(covid_plot)

# Sepsis vs health
sepsis_pls <- prepare_plsda_dataset(combat_df, comparison = "sepsis_vs_health")
cat("\nSepsis vs health dimensions:\n")
print(dim(sepsis_pls$X))
print(summary(sepsis_pls$Y))

sepsis_model <- fit_splsda_model(sepsis_pls$X, sepsis_pls$Y)
print(sepsis_model$optimal_ncomp)
print(sepsis_model$optimal_keepX)

sepsis_plot <- build_publication_pls_plot(
  final_model = sepsis_model$final_model,
  group_vector = sepsis_pls$Y,
  palette = pls_colors$sepsis_vs_health,
  output_file = "health_vs_sepsis.pdf"
)
print(sepsis_plot)

