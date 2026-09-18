evidence_section <- function(code) {
  tryCatch(list(status = "ok", result = force(code)), error = function(e) list(status = "unavailable", reason = conditionMessage(e)))
}

query_gene_context <- function(cfg, gene, disease, analysis_id = NULL, cohort_id = NULL, subgroup = NULL,
                                model_csv = NULL, gene_effect_csv = NULL, model_regex = NULL, depmap_release = NULL,
                                covariates = character(), libraries = character(), gmt_path = NULL,
                                include_r2 = FALSE, r2_group_1 = NULL, r2_group_2 = NULL) {
  gene <- study_symbols(gene)
  if (length(gene) != 1) stop("Provide a single gene", call. = FALSE)
  result <- list(gene = gene, disease = disease)
  if (!is.null(analysis_id)) {
    a <- load_study_object(cfg, analysis_id, "analysis")
    if (tolower(a$disease) != tolower(disease)) stop("Analysis disease does not match requested context", call. = FALSE)
    result$differential_expression <- list(summary = a$summary[a$summary$gene == gene, , drop = FALSE],
      contrasts = a$data[a$data$gene == gene, , drop = FALSE], provenance = a$provenance)
    cohort_id <- cohort_id %||% a$cohort_id
    subgroup <- subgroup %||% a$target
  }
  if (!is.null(cohort_id)) {
    c <- load_study_object(cfg, cohort_id, "cohort")
    if (tolower(c$disease) != tolower(disease)) stop("Cohort disease does not match requested context", call. = FALSE)
    result$expression <- evidence_section({
      expr <- clinical_expression(c)
      if (!gene %in% rownames(expr)) stop("Gene absent from cohort")
      stats <- lapply(split(as.numeric(expr[gene, ]), c$metadata$subgroup), function(v) list(n = length(v), mean = mean(v), median = stats::median(v), sd = stats::sd(v)))
      list(by_subgroup = stats, provenance = c$provenance)
    })
    result$clinical <- evidence_section(clinical_associations(cfg, cohort_id, gene, subgroup, covariates))
  } else result$clinical <- list(status = "unavailable", reason = "Supply a registered cohort with sample-level expression and matched clinical outcomes")
  model_regex <- model_regex %||% paste0("(?i)", disease)
  result$depmap <- evidence_section(depmap_gene_dependency(gene, model_csv, gene_effect_csv, model_regex, depmap_release))
  result$depmap_context_note <- "Disease-level model matches are not subgroup-specific models unless the selected metadata explicitly support that subgroup. Dependencies do not establish a therapeutic window."
  if (include_r2) result$r2 <- evidence_section({
    if (tolower(disease) != tolower(cfg$r2_disease %||% "medulloblastoma")) stop("Configured R2 dataset does not match disease; set r2_disease and dataset configuration")
    r <- r2_compare_groups(cfg, group_1 = r2_group_1, group_2 = r2_group_2, top_n = 5000L, p_threshold = 1)
    rows <- r$data[toupper(r$data$gene) == gene, , drop = FALSE]
    list(rows = rows, provenance = r$provenance, note = "Absent genes are unobserved in a capped result, not evidence of no differential expression")
  })
  libraries <- unlist(libraries, use.names = FALSE)
  if (length(libraries) || !is.null(gmt_path)) result$pathway_membership <- evidence_section({
    inputs <- if (!is.null(gmt_path)) list(load_gene_sets(cfg, gmt_path = gmt_path)) else lapply(libraries, function(l) load_gene_sets(cfg, l))
    lapply(inputs, function(x) list(terms = as.list(names(x$sets)[vapply(x$sets, function(s) gene %in% s, logical(1))]), provenance = x$provenance))
  })
  result
}

prioritize_candidates <- function(cfg, analysis_id, model_csv = NULL, gene_effect_csv = NULL, model_regex = NULL,
                                   depmap_release = NULL, dependency_cutoff = -0.5, direction = "up", limit = 50L,
                                   reference_analysis_id = NULL, clinical_cohort_id = NULL, clinical_covariates = character()) {
  a <- load_study_object(cfg, analysis_id, "analysis")
  rows <- a$summary[a$summary$passes_all_comparators, , drop = FALSE]
  if (direction != "both") rows <- rows[rows$direction == direction, , drop = FALSE]
  if (!nrow(rows)) return(list(analysis_id = analysis_id, candidates = rows, total = 0L, note = "No genes pass every comparator at the stored thresholds"))
  model_regex <- model_regex %||% paste0("(?i)", a$disease)
  dep <- evidence_section({
    mp <- normalize_local_path(model_csv, "DEPMAP_MODEL_CSV")
    gp <- normalize_local_path(gene_effect_csv, "DEPMAP_GENE_EFFECT_CSV")
    list(data = summarize_depmap_for_genes(rows$gene, mp, gp, model_regex, dependency_cutoff),
      provenance = list(model_file = file_provenance(mp), gene_effect_file = file_provenance(gp),
        release = depmap_release %||% Sys.getenv("DEPMAP_RELEASE", "unspecified"), model_regex = model_regex))
  })
  if (dep$status == "ok") {
    dep$result$data$gene <- toupper(dep$result$data$gene)
    rows <- merge(rows, dep$result$data, by = "gene", all.x = TRUE, sort = FALSE)
  } else {
    rows$relevant_mean_gene_effect <- NA_real_
    rows$dependency_selectivity <- NA_real_
    rows$relevant_models_n <- NA_integer_
  }
  rows$dependency_supported <- is.finite(rows$relevant_mean_gene_effect) & rows$relevant_mean_gene_effect < dependency_cutoff
  # Missing dependency/selectivity produces no integrated score; it is not imputed as evidence.
  rows$priority_score <- rows$min_abs_log2fc + pmax(0, -rows$relevant_mean_gene_effect) + pmax(0, rows$dependency_selectivity)
  rows <- rows[order(!rows$dependency_supported, -rows$priority_score, rows$max_adjusted_p, na.last = TRUE), , drop = FALSE]
  reference <- NULL
  if (!is.null(reference_analysis_id)) {
    ref <- load_study_object(cfg, reference_analysis_id, "analysis")
    if (tolower(ref$disease) != tolower(a$disease) || ref$target != a$target) stop("Reference must match disease and target subgroup", call. = FALSE)
    ri <- match(rows$gene, ref$summary$gene)
    rows$reference_passes_all <- ref$summary$passes_all_comparators[ri]
    rows$reference_direction <- ref$summary$direction[ri]
    rows$reference_concordant <- rows$reference_passes_all & rows$direction == rows$reference_direction
    reference <- study_overview(ref)
  }
  rows$rank <- seq_len(nrow(rows))
  total <- nrow(rows)
  rows <- utils::head(rows, limit)
  clinical <- NULL
  if (!is.null(clinical_cohort_id)) {
    c <- load_study_object(cfg, clinical_cohort_id, "cohort")
    if (tolower(c$disease) != tolower(a$disease)) stop("Clinical cohort disease does not match analysis", call. = FALSE)
    clinical <- clinical_associations(cfg, clinical_cohort_id, rows$gene, a$target, clinical_covariates)
  }
  list(analysis_id = analysis_id, candidates = rows, total = total, depmap_status = dep$status,
    depmap_reason = dep$reason, depmap_provenance = dep$result$provenance, clinical = clinical, reference = reference,
    scoring = "Order by dependency support, then min(abs(log2FC)) + max(0,-mean dependency) + max(0,background mean - disease mean), then max adjusted P. Missing dependency/selectivity => null score. Heuristic for follow-up, not a probability.",
    limitations = "Specificity is relative to supplied comparators. Disease model selection does not imply molecular subgroup specificity. Clinical associations are exploratory within the returned candidate panel and are not used in ranking.",
    provenance = a$provenance)
}

r2_subgroup_analysis <- function(cfg, target, comparators, disease = "medulloblastoma", track = "subgroup",
                                  fdr = 0.05, min_abs_log2fc = 1, background_path = NULL, refresh = FALSE) {
  if (tolower(disease) != tolower(cfg$r2_disease %||% "medulloblastoma")) stop("Disease does not match configured R2 cohort", call. = FALSE)
  comparators <- unique(tolower(unlist(comparators, use.names = FALSE)))
  target <- tolower(target)
  if (!length(comparators) || target %in% comparators) stop("Supply comparators different from target", call. = FALSE)
  contrasts <- lapply(comparators, function(g) r2_compare_groups(cfg, track, target, g, top_n = 5000L, p_threshold = 1, refresh = refresh))
  data <- do.call(rbind, lapply(seq_along(contrasts), function(i) {
    d <- contrasts[[i]]$data
    data.frame(gene = toupper(d$gene), log2fc = d$patient_effect_group1_minus_group2, adjusted_p = d$adjusted_p, comparator = comparators[i])
  }))
  if (anyDuplicated(data[, c("gene", "comparator")])) stop("R2 returned ambiguous duplicate gene mappings", call. = FALSE)
  universe <- character()
  if (!is.null(background_path)) {
    bg <- read_study_table(background_path)
    require_columns(bg, "gene")
    universe <- unique(study_symbols(bg$gene))
    if (!all(data$gene %in% universe)) stop("Background must include returned R2 genes", call. = FALSE)
  }
  a <- list(label = cfg$r2_dataset_label, disease = disease, target = target, comparators = comparators, data = data,
    universe = universe, summary = summarize_specificity(data, comparators, fdr, min_abs_log2fc),
    thresholds = list(fdr = fdr, min_abs_log2fc = min_abs_log2fc),
    provenance = list(source = "R2", contrasts = lapply(contrasts, `[[`, "provenance"),
      limitations = "R2 results are capped at 5000 genes per comparator. Missing genes are unobserved. Pathway analysis requires an explicit full tested-gene background.",
      background = if (!is.null(background_path)) file_provenance(background_path) else NULL))
  a$id <- save_study_object(cfg, "analysis", a)
  study_overview(a)
}
