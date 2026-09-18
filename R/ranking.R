rank_patient_dependencies <- function(cfg, model_csv = NULL, gene_effect_csv = NULL,
                                      depmap_release = NULL,
                                      track = NULL, group_1 = NULL, group_2 = NULL,
                                      test = "anova", model_regex = NULL,
                                      r2_top_n = 2000L, result_top_n = 50L,
                                      p_threshold = 0.05, min_patient_effect = 0,
                                      dependency_cutoff = -0.5,
                                      expression_weight = 0.45,
                                      dependency_weight = 0.40,
                                      selectivity_weight = 0.15,
                                      refresh_r2 = FALSE) {
  model_csv <- normalize_local_path(model_csv, "DEPMAP_MODEL_CSV")
  gene_effect_csv <- normalize_local_path(gene_effect_csv, "DEPMAP_GENE_EFFECT_CSV")
  depmap_release <- depmap_release %||% Sys.getenv("DEPMAP_RELEASE", "unspecified")
  model_regex <- model_regex %||% cfg$depmap_default_model_regex
  result_top_n <- as_integer_scalar(result_top_n, 50L, 1L, 500L)

  weight_sum <- expression_weight + dependency_weight + selectivity_weight
  if (!is.finite(weight_sum) || weight_sum <= 0) stop("Ranking weights must sum to a positive value", call. = FALSE)
  expression_weight <- expression_weight / weight_sum
  dependency_weight <- dependency_weight / weight_sum
  selectivity_weight <- selectivity_weight / weight_sum

  r2 <- r2_compare_groups(
    cfg = cfg,
    track = track,
    group_1 = group_1,
    group_2 = group_2,
    test = test,
    top_n = r2_top_n,
    p_threshold = p_threshold,
    refresh = refresh_r2
  )
  patient <- r2$data
  patient <- patient[
    is.finite(patient$adjusted_p) & patient$adjusted_p <= p_threshold &
      patient$patient_effect_group1_minus_group2 >= min_patient_effect,
    , drop = FALSE
  ]
  if (!nrow(patient)) stop("No R2 genes passed the patient-expression filters", call. = FALSE)

  dep <- summarize_depmap_for_genes(
    patient$gene,
    model_csv,
    gene_effect_csv,
    model_regex,
    dependency_cutoff
  )
  combined <- merge(patient, dep, by = "gene", all = FALSE, sort = FALSE)
  combined <- combined[is.finite(combined$relevant_mean_gene_effect), , drop = FALSE]
  if (!nrow(combined)) stop("No genes overlapped between the R2 result and DepMap", call. = FALSE)

  combined$patient_expression_z <- safe_zscore(combined$patient_effect_group1_minus_group2)
  combined$dependency_strength_z <- safe_zscore(-combined$relevant_mean_gene_effect)
  selectivity_for_score <- combined$dependency_selectivity
  selectivity_for_score[!is.finite(selectivity_for_score)] <- 0
  combined$dependency_selectivity_z <- safe_zscore(selectivity_for_score)
  combined$priority_score <-
    expression_weight * combined$patient_expression_z +
    dependency_weight * combined$dependency_strength_z +
    selectivity_weight * combined$dependency_selectivity_z

  combined <- combined[order(-combined$priority_score, combined$adjusted_p), , drop = FALSE]
  combined$rank <- seq_len(nrow(combined))
  combined <- combined[, c(
    "rank", "gene", "priority_score", "patient_effect_group1_minus_group2",
    "adjusted_p", "present_n", "relevant_mean_gene_effect",
    "relevant_median_gene_effect", "relevant_min_gene_effect",
    "relevant_fraction_below_cutoff", "background_mean_gene_effect",
    "dependency_selectivity", "relevant_models_n", "patient_expression_z",
    "dependency_strength_z", "dependency_selectivity_z"
  )]

  list(
    ranking = utils::head(combined, result_top_n),
    counts = list(
      r2_rows_returned = nrow(r2$data),
      patient_filter_passed = nrow(patient),
      r2_depmap_overlap = nrow(combined),
      rows_returned = min(result_top_n, nrow(combined))
    ),
    scoring = list(
      expression_weight = expression_weight,
      dependency_weight = dependency_weight,
      selectivity_weight = selectivity_weight,
      dependency_cutoff = dependency_cutoff,
      patient_effect_definition = "positive means group_1 is higher than group_2",
      dependency_effect_definition = "more-negative DepMap gene effect means stronger dependency",
      selectivity_definition = "background mean minus relevant-model mean; positive is more selective"
    ),
    provenance = list(
      r2 = r2$provenance,
      depmap_release = depmap_release,
      depmap_model_filter_regex = model_regex,
      depmap_model_csv = model_csv,
      depmap_gene_effect_csv = gene_effect_csv,
      generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
  )
}
