clinical_expression <- function(cohort) {
  if (cohort$scale != "counts") return(cohort$expression)
  need_package("edgeR")
  y <- edgeR::calcNormFactors(edgeR::DGEList(cohort$expression))
  edgeR::cpm(y, log = TRUE, prior.count = 2)
}

fit_clinical_gene <- function(gene, expression, metadata, endpoint, covariates,
                               time_column, event_column, metastasis_column) {
  row <- data.frame(gene = gene, endpoint = endpoint, status = "unavailable", n = 0L, events = NA_integer_,
    effect_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
    ph_test_p = NA_real_, detail = "", stringsAsFactors = FALSE)
  tryCatch({
    if (!gene %in% rownames(expression)) stop("Gene absent from expression matrix")
    fields <- if (endpoint == "survival") c(time_column, event_column) else metastasis_column
    require_columns(metadata, c(fields, covariates))
    d <- data.frame(expression = as.numeric(expression[gene, ]))
    if (endpoint == "survival") {
      d$time <- numeric_column(metadata[[time_column]], time_column, TRUE)
      d$event <- numeric_column(metadata[[event_column]], event_column, TRUE)
      if (any(d$time < 0, na.rm = TRUE)) stop("Survival times must be non-negative")
    } else d$event <- numeric_column(metadata[[metastasis_column]], metastasis_column, TRUE)
    if (any(!is.na(d$event) & !d$event %in% c(0, 1))) stop("Event/metastasis must be explicitly coded 0 or 1 (1=event/metastasis)")
    if (length(covariates)) for (i in seq_along(covariates)) d[[paste0("cov", i)]] <- metadata[[covariates[i]]]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    row$n <- nrow(d)
    row$events <- sum(d$event)
    if (nrow(d) < 10 || row$events < 5 || (endpoint == "metastasis" && sum(d$event == 0) < 5))
      stop("Insufficient complete cases or events (need >=10 samples, >=5 events; metastasis also >=5 controls)")
    s <- stats::sd(d$expression)
    if (!is.finite(s) || s == 0) stop("Expression has zero variance")
    d$expression <- as.numeric(scale(d$expression))
    rhs <- paste(c("expression", paste0("cov", seq_along(covariates))[seq_along(covariates)]), collapse = " + ")
    design <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
    if (qr(design)$rank != ncol(design)) stop("Clinical design is confounded")
    notes <- character()
    fit <- withCallingHandlers({
      if (endpoint == "survival") {
        need_package("survival")
        survival::coxph(stats::as.formula(paste("survival::Surv(time, event) ~", rhs)), data = d, x = TRUE)
      } else stats::glm(stats::as.formula(paste("event ~", rhs)), data = d, family = stats::binomial())
    }, warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (length(notes) || (endpoint == "metastasis" && (!isTRUE(fit$converged) || any(fit$fitted.values < 1e-8 | fit$fitted.values > 1 - 1e-8))))
      stop(paste("Unstable clinical fit:", paste(notes, collapse = "; "), "possible separation/non-convergence"))
    cf <- summary(fit)$coefficients
    estimate <- cf["expression", 1]
    se <- cf["expression", if (endpoint == "survival") "se(coef)" else "Std. Error"]
    if (!is.finite(estimate) || !is.finite(se) || se <= 0 || se > 10) stop("Unstable clinical coefficient")
    row$effect_ratio <- exp(estimate)
    row$ci_lower <- exp(estimate - 1.96 * se)
    row$ci_upper <- exp(estimate + 1.96 * se)
    row$p_value <- cf["expression", if (endpoint == "survival") "Pr(>|z|)" else "Pr(>|z|)"]
    if (endpoint == "survival") row$ph_test_p <- tryCatch(survival::cox.zph(fit)$table["expression", "p"], error = function(e) NA_real_)
    row$status <- "ok"
    row$detail <- if (endpoint == "survival") "Hazard ratio per SD higher expression" else "Odds ratio per SD higher expression"
    row
  }, error = function(e) { row$detail <- conditionMessage(e); row })
}

clinical_associations <- function(cfg, cohort_id, genes, subgroup = NULL, covariates = character(),
                                   time_column = "survival_time", event_column = "survival_event",
                                   metastasis_column = "metastasis", endpoints = c("survival", "metastasis")) {
  c <- load_study_object(cfg, cohort_id, "cohort")
  genes <- unique(study_symbols(unlist(genes)))
  covariates <- unlist(covariates, use.names = FALSE)
  endpoints <- unlist(endpoints, use.names = FALSE)
  if (!length(genes) || length(genes) > 500) stop("Provide between 1 and 500 genes", call. = FALSE)
  if (!length(endpoints) || any(!endpoints %in% c("survival", "metastasis"))) stop("Unknown endpoint", call. = FALSE)
  expression <- clinical_expression(c)
  metadata <- c$metadata
  if (!is.null(subgroup)) {
    keep <- metadata$subgroup == subgroup
    if (!any(keep)) stop("Subgroup absent from cohort", call. = FALSE)
    metadata <- metadata[keep, , drop = FALSE]
    expression <- expression[, keep, drop = FALSE]
  }
  rows <- do.call(rbind, lapply(endpoints, function(endpoint) do.call(rbind, lapply(genes, function(gene)
    fit_clinical_gene(gene, expression, metadata, endpoint, covariates, time_column, event_column, metastasis_column)))))
  rows$adjusted_p <- NA_real_
  for (endpoint in endpoints) {
    i <- which(rows$endpoint == endpoint)
    rows$adjusted_p[i] <- stats::p.adjust(rows$p_value[i], "BH", n = length(genes))
  }
  list(rows = rows, provenance = list(cohort = study_overview(c), subgroup = subgroup,
    covariates = as.list(covariates), time_column = time_column, event_column = event_column, metastasis_column = metastasis_column,
    multiple_testing = "BH over the requested gene panel separately per endpoint, including failed fits in family size",
    limitations = "Exploratory association, not causation. Continuous expression; no optimal-cutpoint search. Check PH test and confidence intervals. Pooled subgroups require appropriate covariate adjustment.",
    generated_at = utc_now()))
}
