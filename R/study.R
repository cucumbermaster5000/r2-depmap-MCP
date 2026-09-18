# Local, versioned study objects. No patient data are sent to external services.
utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

need_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install ", package, " using scripts/install-dependencies.R --analysis", call. = FALSE)
}

file_provenance <- function(path) {
  path <- normalize_local_path(path)
  list(path = path, md5 = unname(tools::md5sum(path)), bytes = unname(file.info(path)$size))
}

object_digest <- function(value) {
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(value, path, version = 3)
  unname(tools::md5sum(path))
}

save_study_object <- function(cfg, kind, value) {
  id <- paste0(kind, "_", object_digest(value))
  value$id <- id
  cache_put(cfg, id, value)
  id
}

load_study_object <- function(cfg, id, kind) {
  if (length(id) != 1L || !grepl(paste0("^", kind, "_[a-f0-9]{32}$"), id))
    stop("Invalid ", kind, " identifier", call. = FALSE)
  value <- cache_get(cfg, id)
  if (is.null(value)) stop("Unknown ", kind, " identifier: ", id, call. = FALSE)
  value
}

read_study_table <- function(path) {
  data.table::fread(normalize_local_path(path), data.table = FALSE, check.names = FALSE,
                   showProgress = FALSE, na.strings = c("", "NA", "NaN"))
}

require_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(names(data))) stop("Duplicate column names are not supported", call. = FALSE)
}

study_symbols <- function(x) {
  out <- toupper(trimws(as.character(x)))
  if (anyNA(out) || any(!nzchar(out)) || any(grepl("[;|,\\s]", out, perl = TRUE)))
    stop("Use one non-empty human gene symbol per row; resolve probes, aliases and multiple mappings before import", call. = FALSE)
  if (any(grepl("^ENS[A-Z]*G[0-9]|^[0-9]+$", out)))
    stop("Map Ensembl/Entrez identifiers to human gene symbols before import", call. = FALSE)
  out
}

numeric_column <- function(x, label, allow_na = FALSE) {
  y <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & !is.finite(y)) || (!allow_na && anyNA(y)))
    stop(label, " must contain finite numeric values", call. = FALSE)
  y
}

study_overview <- function(x) {
  list(id = x$id, disease = x$disease, label = x$label, target = x$target,
       comparators = as.list(x$comparators), genes_n = if (!is.null(x$expression)) nrow(x$expression) else length(x$universe),
       samples_n = if (!is.null(x$metadata)) nrow(x$metadata) else NULL,
       subgroups = if (!is.null(x$metadata)) as.list(table(x$metadata$subgroup)) else NULL,
       provenance = x$provenance)
}

list_studies <- function(cfg) {
  paths <- list.files(cfg$cache_dir, "^(cohort|analysis)_[a-f0-9]{32}\\.rds$", full.names = TRUE)
  list(studies = lapply(paths, function(p) study_overview(readRDS(p))))
}

summarize_specificity <- function(data, comparators, fdr = 0.05, min_abs_log2fc = 1) {
  if (!is.finite(fdr) || fdr <= 0 || fdr > 1 || !is.finite(min_abs_log2fc) || min_abs_log2fc < 0)
    stop("Require 0 < fdr <= 1 and min_abs_log2fc >= 0", call. = FALSE)
  rows <- lapply(split(data, data$gene), function(d) {
    tested <- is.finite(d$adjusted_p) & is.finite(d$log2fc)
    significant <- tested & d$adjusted_p <= fdr & abs(d$log2fc) >= min_abs_log2fc & d$log2fc != 0
    complete <- setequal(d$comparator[tested], comparators)
    same_direction <- all(d$log2fc[tested] > 0) || all(d$log2fc[tested] < 0)
    specific <- complete && all(significant) && same_direction
    data.frame(gene = d$gene[1], comparisons_tested = sum(tested), comparisons_required = length(comparators),
      comparisons_passed = sum(significant), consistent_direction = same_direction && any(tested),
      subgroup_specific = specific && length(comparators) > 1L,
      passes_all_comparators = specific,
      direction = if (!any(tested) || !same_direction) "mixed_or_unknown" else if (d$log2fc[which(tested)[1]] > 0) "up" else "down",
      min_abs_log2fc = if (any(tested)) min(abs(d$log2fc[tested])) else NA_real_,
      max_adjusted_p = if (any(tested)) max(d$adjusted_p[tested]) else NA_real_, stringsAsFactors = FALSE)
  })
  if (!length(rows)) stop("No differential-expression rows available", call. = FALSE)
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(!result$passes_all_comparators, result$max_adjusted_p, -result$min_abs_log2fc), ]
}

import_de_results <- function(cfg, path, disease, target, label = "Imported differential expression",
                              gene_column = "gene", effect_column = "log2fc", adjusted_p_column = "adjusted_p",
                              comparator_column = "comparator", comparator = NULL, background_path = NULL,
                              fdr = 0.05, min_abs_log2fc = 1) {
  d <- read_study_table(path)
  require_columns(d, c(gene_column, effect_column, adjusted_p_column))
  if (is.null(comparator)) require_columns(d, comparator_column)
  comparisons <- if (is.null(comparator)) as.character(d[[comparator_column]]) else rep(comparator, nrow(d))
  if (!nrow(d) || anyNA(comparisons) || any(!nzchar(comparisons)) || any(comparisons == target))
    stop("Provide non-empty comparator labels different from target", call. = FALSE)
  data <- data.frame(gene = study_symbols(d[[gene_column]]),
    log2fc = numeric_column(d[[effect_column]], effect_column, TRUE),
    adjusted_p = numeric_column(d[[adjusted_p_column]], adjusted_p_column, TRUE),
    comparator = comparisons, stringsAsFactors = FALSE)
  if (any(data$adjusted_p < 0 | data$adjusted_p > 1, na.rm = TRUE)) stop("Adjusted P values must be between 0 and 1", call. = FALSE)
  if (anyDuplicated(data[, c("gene", "comparator")])) stop("Duplicate gene/comparator rows; resolve mappings before import", call. = FALSE)
  comparators <- sort(unique(comparisons))
  # With no explicit universe, all finite tests from every comparator form the tested intersection.
  universe <- Reduce(intersect, lapply(comparators, function(g) data$gene[data$comparator == g & is.finite(data$adjusted_p) & is.finite(data$log2fc)]))
  if (!is.null(background_path)) {
    bg <- read_study_table(background_path)
    require_columns(bg, "gene")
    universe <- unique(study_symbols(bg$gene))
    if (!all(data$gene %in% universe)) stop("Background must contain every imported gene", call. = FALSE)
  }
  result <- list(label = label, disease = disease, target = target, comparators = comparators,
    data = data, universe = universe, summary = summarize_specificity(data, comparators, fdr, min_abs_log2fc),
    thresholds = list(fdr = fdr, min_abs_log2fc = min_abs_log2fc),
    provenance = list(source = "user differential-expression results", file = file_provenance(path),
      background_file = if (!is.null(background_path)) file_provenance(background_path) else NULL,
      universe_definition = if (is.null(background_path)) "Intersection of finite tests; input MUST include all tested genes, not only significant hits" else "Explicit user-supplied tested-gene universe",
      effect_definition = "log2 fold change: target minus comparator; retained without recalculation",
      imported_at = utc_now(), limitations = "Specificity is relative to supplied comparators; a pooled rest comparison cannot establish subgroup specificity. Input adjusted P values are retained; no joint FDR claim is made."))
  result$id <- save_study_object(cfg, "analysis", result)
  study_overview(result)
}

register_cohort <- function(cfg, expression_path, metadata_path, disease, label = "Local cohort",
                            scale = "log2", gene_column = "gene", sample_column = "sample_id",
                            subgroup_column = "subgroup", source = "user") {
  if (!scale %in% c("log2", "tpm", "counts")) stop("scale must be log2, tpm or counts", call. = FALSE)
  e <- read_study_table(expression_path)
  m <- read_study_table(metadata_path)
  require_columns(e, gene_column)
  require_columns(m, c(sample_column, subgroup_column))
  genes <- study_symbols(e[[gene_column]])
  if (anyDuplicated(genes)) stop("Duplicate gene symbols; collapse or resolve mappings before import", call. = FALSE)
  samples <- as.character(m[[sample_column]])
  if (anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) stop("Sample IDs must be non-empty and unique", call. = FALSE)
  if (anyNA(m[[subgroup_column]]) || any(!nzchar(as.character(m[[subgroup_column]])))) stop("Subgroup labels must not be missing", call. = FALSE)
  expression_samples <- setdiff(names(e), gene_column)
  if (!setequal(samples, expression_samples)) stop("Expression columns and metadata sample IDs must match exactly", call. = FALSE)
  mat <- vapply(samples, function(s) numeric_column(e[[s]], s), numeric(nrow(e)))
  dimnames(mat) <- list(genes, samples)
  if (nrow(mat) < 2 || ncol(mat) < 3) stop("Provide at least two genes and three biological samples", call. = FALSE)
  if (scale != "log2" && any(mat < 0)) stop("Counts/TPM cannot be negative", call. = FALSE)
  if (scale == "counts" && (any(abs(mat - round(mat)) > 1e-7) || any(colSums(mat) == 0))) stop("Counts must be non-negative integers with non-empty libraries", call. = FALSE)
  if (scale == "tpm") mat <- log2(mat + 1)
  m$sample_id <- samples
  m$subgroup <- as.character(m[[subgroup_column]])
  value <- list(label = label, disease = disease, expression = mat, metadata = m, scale = scale,
    provenance = list(source = source, expression = file_provenance(expression_path), metadata = file_provenance(metadata_path),
      scale = scale, transformation = if (scale == "tpm") "log2(TPM + 1)" else "none", imported_at = utc_now()))
  value$id <- save_study_object(cfg, "cohort", value)
  study_overview(value)
}

study_design <- function(metadata, covariates = character()) {
  covariates <- unlist(covariates, use.names = FALSE)
  require_columns(metadata, c("subgroup", covariates))
  design_data <- data.frame(group = factor(metadata$subgroup))
  if (length(covariates)) for (i in seq_along(covariates)) design_data[[paste0("cov", i)]] <- metadata[[covariates[i]]]
  if (anyNA(design_data)) stop("Subgroup/covariate values must be complete for differential expression", call. = FALSE)
  design <- stats::model.matrix(~ 0 + ., design_data)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop("Design is confounded or lacks residual degrees of freedom", call. = FALSE)
  list(matrix = design, groups = levels(design_data$group))
}

analyze_subgroup <- function(cfg, cohort_id, target, covariates = character(), fdr = 0.05, min_abs_log2fc = 1) {
  need_package("limma")
  cohort <- load_study_object(cfg, cohort_id, "cohort")
  counts <- table(cohort$metadata$subgroup)
  if (!target %in% names(counts) || length(counts) < 2 || any(counts < 3)) stop("Need target and comparators with at least three biological replicates per subgroup", call. = FALSE)
  design <- study_design(cohort$metadata, covariates)
  expr <- cohort$expression
  if (cohort$scale == "counts") {
    need_package("edgeR")
    y <- edgeR::DGEList(expr)
    keep <- edgeR::filterByExpr(y, design = design$matrix)
    if (sum(keep) < 2) stop("Too few genes pass count filtering", call. = FALSE)
    y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
    fit <- limma::lmFit(limma::voom(y, design$matrix, plot = FALSE), design$matrix)
  } else {
    keep <- apply(expr, 1, stats::sd) > 0
    if (sum(keep) < 2) stop("Too few variable genes", call. = FALSE)
    fit <- limma::lmFit(expr[keep, , drop = FALSE], design$matrix)
  }
  comparators <- setdiff(design$groups, target)
  contrast <- matrix(0, ncol(design$matrix), length(comparators), dimnames = list(colnames(design$matrix), comparators))
  contrast[match(target, design$groups), ] <- 1
  for (i in seq_along(comparators)) contrast[match(comparators[i], design$groups), i] <- -1
  fitted <- limma::eBayes(limma::contrasts.fit(fit, contrast), trend = cohort$scale != "counts")
  data <- do.call(rbind, lapply(seq_along(comparators), function(i) {
    tt <- limma::topTable(fitted, coef = i, number = Inf, sort.by = "none", adjust.method = "BH")
    data.frame(gene = rownames(tt), log2fc = tt$logFC, adjusted_p = tt$adj.P.Val, p_value = tt$P.Value, comparator = comparators[i])
  }))
  result <- list(label = paste(cohort$label, target), disease = cohort$disease, target = target, comparators = comparators,
    cohort_id = cohort_id, data = data, universe = unique(data$gene), summary = summarize_specificity(data, comparators, fdr, min_abs_log2fc),
    thresholds = list(fdr = fdr, min_abs_log2fc = min_abs_log2fc),
    provenance = list(source = cohort$provenance, method = if (cohort$scale == "counts") "edgeR TMM + limma voom" else "limma trend",
      covariates = as.list(covariates), samples_per_group = as.list(counts), multiple_testing = "BH separately over all genes per contrast",
      effect_definition = "target minus comparator", generated_at = utc_now(), limma_version = as.character(utils::packageVersion("limma"))))
  result$id <- save_study_object(cfg, "analysis", result)
  study_overview(result)
}

analysis_results <- function(cfg, analysis_id, only_candidates = FALSE, direction = "both", offset = 0L, limit = 100L) {
  a <- load_study_object(cfg, analysis_id, "analysis")
  rows <- a$summary
  if (only_candidates) rows <- rows[rows$passes_all_comparators, , drop = FALSE]
  if (direction != "both") rows <- rows[rows$direction == direction, , drop = FALSE]
  n <- nrow(rows)
  rows <- utils::head(rows[seq_len(n) > offset, , drop = FALSE], limit)
  list(analysis = study_overview(a), total = n, offset = offset, rows = rows,
       contrasts = a$data[a$data$gene %in% rows$gene, , drop = FALSE], thresholds = a$thresholds)
}
