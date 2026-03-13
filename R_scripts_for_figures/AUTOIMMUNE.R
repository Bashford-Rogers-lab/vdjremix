# Publication-ready analysis script: autoimmune module dynamics and heatmap

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(tibble)
})

files <- list(
  eigenvectors = file.path(base_dir, "Eigenvectors_AI_BCR.txt"),
  metadata     = file.path(base_dir, "Overall_autoimmune paper summary locations_share(Table S3).csv"),
  mmf          = "~/Downloads/Time_since_MMF.txt",
  rtx          = "~/Downloads/Time_since_RTX.txt"
)

# ---- Analysis constants ----
target_module <- "Module_9"
module_label  <- "Module 9"
module_cols   <- paste0("Module_", 1:20)
module_labels <- paste("Module", 1:20)

plot_colors <- c(AAV = "#228B22", SLE = "#FFB347", Healthy = "#F49595")

# ---- Load primary data ----
metadata <- read.csv(files$metadata, stringsAsFactors = FALSE)
eigenvectors <- read.delim(files$eigenvectors, check.names = FALSE)
eigenvectors <- as.data.frame(eigenvectors)

# Use sample IDs as row names so metadata and module scores can be aligned safely.
rownames(eigenvectors) <- rownames(eigenvectors)
metadata <- metadata[1:238, , drop = FALSE]
rownames(metadata) <- metadata$Sample

# Keep only samples present in both tables and enforce identical order.
common_ids <- intersect(rownames(metadata), rownames(eigenvectors))
metadata <- metadata[common_ids, , drop = FALSE]
eigenvectors <- eigenvectors[common_ids, , drop = FALSE]
stopifnot(identical(rownames(metadata), rownames(eigenvectors)))

# ---- Derive metadata fields and merge therapy tables ----
metadata <- metadata %>%
  mutate(
    age_numeric = as.numeric(Patient.age.at.time.of.sampling),
    age_group = case_when(
      dplyr::between(age_numeric, 20, 40) ~ "20-40",
      dplyr::between(age_numeric, 41, 60) ~ "41-60",
      dplyr::between(age_numeric, 61, 80) ~ "61-80",
      age_numeric >= 81                   ~ "81+",
      TRUE                                ~ NA_character_
    ),
    disease_broad = case_when(
      grepl("AAV|EGPA", Disease.group) ~ "AAV",
      grepl("SLE", Disease.group)      ~ "SLE",
      grepl("Healthy", Disease.group)  ~ "Healthy",
      TRUE                              ~ "Other"
    )
  ) %>%
  rownames_to_column("sample_id")

modules_df <- eigenvectors %>%
  rownames_to_column("sample_id")

mmf_df <- read.delim(files$mmf) %>%
  rename(sample_id = samples)

rtx_df <- read.delim(files$rtx) %>%
  rename(sample_id = samples)

analysis_df <- metadata %>%
  left_join(modules_df, by = "sample_id") %>%
  left_join(mmf_df, by = "sample_id") %>%
  left_join(rtx_df, by = "sample_id")

# Optional QC check for age-group counts.
print(table(analysis_df$age_group, useNA = "ifany"))

# ---- Subsets for module-dynamics figure ----
aav_diagnosis <- analysis_df %>%
  filter(
    disease_broad == "AAV",
    grepl("At diagnosis|At flare", Disease.and.flare.status)
  )

sle_diagnosis <- analysis_df %>%
  filter(
    disease_broad == "SLE",
    grepl("At diagnosis|At flare", Disease.and.flare.status)
  )

mmf_subset <- analysis_df %>%
  filter(
    therapies.x == "MMF",
    disease_broad %in% c("AAV", "SLE"),
    !is.na(time_since_MMFs),
    dplyr::between(time_since_MMFs, 0, 500)
  )

rtx_subset <- analysis_df %>%
  filter(
    therapies.y == "RTX",
    disease_broad %in% c("AAV", "SLE"),
    !is.na(time_since_RTXs),
    dplyr::between(time_since_RTXs, 0, 500)
  )

healthy_subset <- analysis_df %>%
  filter(disease_broad == "Healthy")

y_range <- range(
  c(
    aav_diagnosis[[target_module]],
    sle_diagnosis[[target_module]],
    mmf_subset[[target_module]],
    rtx_subset[[target_module]],
    healthy_subset[[target_module]]
  ),
  na.rm = TRUE
)

# ---- Small plotting helpers ----
make_single_boxplot <- function(data, x_label, fill_color, y_label = "", show_y = FALSE) {
  ggplot(data, aes(x = x_label, y = .data[[target_module]])) +
    geom_boxplot(fill = fill_color, width = 0.6, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(shape = 21, fill = "white", width = 0.15, size = 2) +
    coord_cartesian(ylim = y_range) +
    labs(x = NULL, y = y_label) +
    theme_classic2() +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = if (show_y) element_text() else element_blank(),
      axis.title.y = if (show_y) element_text() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      axis.line.y  = if (show_y) element_line() else element_blank()
    )
}

make_scatter_panel <- function(data, x_var, x_label) {
  ggscatter(
    data,
    x = x_var,
    y = target_module,
    color = "disease_broad",
    shape = "disease_broad",
    palette = plot_colors[c("AAV", "SLE")],
    add = "reg.line",
    add.params = list(color = "grey30"),
    conf.int = FALSE,
    cor.coef = TRUE,
    cor.method = "pearson",
    cor.coef.coord = c(50, y_range[2] * 0.9),
    xlab = x_label,
    ylab = "",
    ggtheme = theme_classic2()
  ) +
    coord_cartesian(ylim = y_range) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
}

# ---- Figure 1: module dynamics ----
p_aav_diag <- make_single_boxplot(
  aav_diagnosis,
  x_label = "AAV Pre-Treatment\n(day 0)",
  fill_color = plot_colors[["AAV"]],
  y_label = paste(module_label, "Score"),
  show_y = TRUE
)

p_sle_diag <- make_single_boxplot(
  sle_diagnosis,
  x_label = "SLE Pre-Treatment\n(day 0)",
  fill_color = plot_colors[["SLE"]]
)

p_mmf <- make_scatter_panel(
  mmf_subset,
  x_var = "time_since_MMFs",
  x_label = "Time since MMF therapy (d)"
)

p_rtx <- make_scatter_panel(
  rtx_subset,
  x_var = "time_since_RTXs",
  x_label = "Time since RTX therapy (d)"
)

p_healthy <- make_single_boxplot(
  healthy_subset,
  x_label = "Healthy",
  fill_color = plot_colors[["Healthy"]]
)

final_plot <- ggarrange(
  p_aav_diag, p_sle_diag, p_mmf, p_rtx, p_healthy,
  ncol = 5,
  nrow = 1,
  align = "h",
  widths = c(1.5, 1.5, 2.5, 2.5, 1.5),
  common.legend = TRUE,
  legend = "bottom"
)

print(final_plot)
ggsave("Module_9_dynamics.pdf", plot = final_plot, width = 10, height = 6)

# ---- Statistical comparisons reported to console ----
diagnosis_group <- analysis_df %>%
  filter(
    disease_broad %in% c("AAV", "SLE"),
    grepl("At diagnosis|At flare", Disease.and.flare.status)
  )

early_mmf_group <- analysis_df %>%
  filter(
    disease_broad %in% c("AAV", "SLE"),
    !is.na(time_since_MMFs),
    time_since_MMFs > 0,
    time_since_MMFs <= 90
  )

early_rtx_group <- analysis_df %>%
  filter(
    disease_broad %in% c("AAV", "SLE"),
    !is.na(time_since_RTXs),
    time_since_RTXs > 0,
    time_since_RTXs <= 90
  )

run_wilcoxon_summary <- function(group_1, group_2, label_1, label_2) {
  test <- wilcox.test(group_1[[target_module]], group_2[[target_module]])
  message("\n--- Statistical comparison ---")
  message("Group 1: ", label_1, " (n = ", nrow(group_1), ")")
  message("Group 2: ", label_2, " (n = ", nrow(group_2), ")")
  message("Wilcoxon p-value: ", format.pval(test$p.value, digits = 3))
}

run_wilcoxon_summary(diagnosis_group, early_mmf_group, "Diagnosis/flare", "1-3 months post-MMF")
run_wilcoxon_summary(diagnosis_group, early_rtx_group, "Diagnosis/flare", "1-3 months post-RTX")

# ---- Figure 2: heatmap of disease-vs-healthy significance ----
heatmap_df <- analysis_df %>%
  mutate(
    Disease.and.flare.status = as.character(Disease.and.flare.status),
    Disease.and.flare.status = na_if(Disease.and.flare.status, "NA")
  ) %>%
  filter(grepl("At diagnosis|Healthy", Disease.and.flare.status, ignore.case = TRUE))

healthy_data <- heatmap_df %>%
  filter(Disease.and.flare.status == "Healthy")

disease_groups <- setdiff(unique(heatmap_df$Disease.and.flare.status), "Healthy")
disease_groups <- disease_groups[!is.na(disease_groups)]

comparison_results <- list()
for (mod_idx in seq_along(module_cols)) {
  current_module <- module_cols[mod_idx]
  current_label <- module_labels[mod_idx]
  healthy_values <- healthy_data[[current_module]]
  healthy_values <- healthy_values[!is.na(healthy_values)]

  for (disease_group in disease_groups) {
    disease_values <- heatmap_df %>%
      filter(Disease.and.flare.status == disease_group) %>%
      pull(.data[[current_module]])
    disease_values <- disease_values[!is.na(disease_values)]

    if (length(healthy_values) >= 2 && length(disease_values) >= 2) {
      test_result <- wilcox.test(disease_values, healthy_values, exact = FALSE, conf.int = FALSE)
      p_value <- test_result$p.value
      direction <- if (mean(disease_values) > mean(healthy_values)) "higher" else "lower"
    } else {
      p_value <- 1
      direction <- "neutral"
    }

    comparison_results[[length(comparison_results) + 1]] <- data.frame(
      module = current_label,
      disease_group = disease_group,
      p_value = p_value,
      mean_disease = mean(disease_values, na.rm = TRUE),
      mean_healthy = mean(healthy_values, na.rm = TRUE),
      direction = direction,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- bind_rows(comparison_results) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    category = case_when(
      FDR < 0.005 & direction == "higher" ~ "Significantly higher\nFDR < 0.005",
      FDR < 0.05  & direction == "higher" ~ "Significantly higher\nFDR < 0.05",
      FDR < 0.005 & direction == "lower"  ~ "Significantly lower\nFDR < 0.005",
      FDR < 0.05  & direction == "lower"  ~ "Significantly lower\nFDR < 0.05",
      TRUE                                  ~ "Not significant"
    )
  )

# Asterisks mark modules with a significant difference across disease groups.
disease_only_long <- heatmap_df %>%
  filter(Disease.and.flare.status != "Healthy", !is.na(Disease.and.flare.status)) %>%
  pivot_longer(
    cols = all_of(module_cols),
    names_to = "module_original",
    values_to = "score"
  )

kw_results <- lapply(module_cols, function(current_module) {
  current_data <- disease_only_long %>%
    filter(module_original == current_module) %>%
    mutate(Disease.and.flare.status = factor(Disease.and.flare.status))

  if (length(unique(current_data$Disease.and.flare.status)) >= 2) {
    kw_test <- tryCatch(
      kruskal.test(score ~ Disease.and.flare.status, data = current_data),
      error = function(e) NULL
    )
    p_value <- if (!is.null(kw_test)) kw_test$p.value else 1
  } else {
    p_value <- 1
  }

  data.frame(module_original = current_module, p_value = p_value)
})

asterisk_column <- " "
asterisk_data <- bind_rows(kw_results) %>%
  mutate(
    FDR_across_diseases = p.adjust(p_value, method = "BH"),
    module = factor(gsub("_", " ", module_original), levels = rev(module_labels)),
    label = if_else(FDR_across_diseases < 0.05, "*", ""),
    disease_group = asterisk_column
  ) %>%
  select(module, disease_group, label)

category_order <- c(
  "Significantly higher\nFDR < 0.005",
  "Significantly higher\nFDR < 0.05",
  "Significantly lower\nFDR < 0.05",
  "Significantly lower\nFDR < 0.005",
  "Not significant"
)

category_colors <- c(
  "Significantly higher\nFDR < 0.005" = "#8B0000",
  "Significantly higher\nFDR < 0.05"  = "#F8BB96",
  "Significantly lower\nFDR < 0.005"  = "#00008B",
  "Significantly lower\nFDR < 0.05"   = "#8EC4E5",
  "Not significant"                    = "#FFFFFF"
)

results_df <- results_df %>%
  mutate(
    category = factor(category, levels = category_order),
    module = factor(module, levels = rev(module_labels))
  )

desired_disease_order <- c(
  "Behcets At diagnosis",
  "SLE At diagnosis",
  "AAV MPO+ At diagnosis",
  "AAV PR3+ At diagnosis",
  "EGPA ANCA- At diagnosis",
  "EGPA MPO+ At diagnosis",
  "CD At diagnosis",
  "IgAV At diagnosis"
)

final_disease_groups <- intersect(desired_disease_order, unique(results_df$disease_group))
if (length(final_disease_groups) == 0) {
  final_disease_groups <- sort(unique(results_df$disease_group))
}

all_x_levels <- c(final_disease_groups, asterisk_column)
results_df$disease_group <- factor(results_df$disease_group, levels = all_x_levels)

heatmap_plot <- ggplot(results_df, aes(x = disease_group, y = module)) +
  geom_tile(aes(fill = category), color = "grey50", linewidth = 0.38) +
  geom_text(
    data = asterisk_data,
    aes(label = label),
    size = 8,
    color = "black",
    fontface = "bold",
    vjust = 0.75
  ) +
  scale_fill_manual(
    values = category_colors,
    name = "Significance compared to health",
    labels = function(x) gsub("\\n", " ", x),
    breaks = category_order,
    drop = FALSE
  ) +
  scale_x_discrete(drop = FALSE, labels = c(final_disease_groups, "")) +
  labs(x = "Diseases and flare status", y = "Module scores") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(face = "bold", size = 10, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 10, margin = margin(r = 10)),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "grey80"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

print(heatmap_plot)
ggsave(
  "heatmap_significance_with_inter-disease_asterisks.pdf",
  plot = heatmap_plot,
  width = 6.5,
  height = 5
)
