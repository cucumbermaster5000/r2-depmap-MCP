# Interface-independent services. Collection may use network/cache; analysis is
# in-memory and returns ordinary R lists, data frames and ggplot objects.
profile_attempt <- function(expr) tryCatch(expr,
  error = function(e) list(status = "unavailable", reason = conditionMessage(e)))

collect_gene_sources <- function(cfg, gene, reporter = NULL, publication_limit = 10L,
    publication_disease = NULL, aliases = character(), include_aliases = TRUE,
    refresh_publications = FALSE, pfister_reporter = NULL) {
  retrieval <- get_r2_expression(cfg, gene, reporter)
  list(retrieval = retrieval,
    survival = profile_attempt(get_cavalli_survival(cfg, retrieval)),
    pfister = profile_attempt(get_pfister_expression(cfg, gene, pfister_reporter)),
    enrichr = profile_attempt(profile_enrichr(cfg, gene)),
    depmap = profile_attempt(get_public_depmap_gene(cfg, gene)),
    publications = profile_attempt(get_gene_publications(cfg, gene,
      publication_limit, publication_disease, aliases, include_aliases, refresh_publications)))
}

analyze_gene_sources <- function(sources, cutoff = "median") {
  cutoff <- match.arg(cutoff, c("median", "mean", "optimal"))
  retrieval <- sources$retrieval
  gene <- retrieval$provenance$queried_gene
  unavailable <- function(x) is.null(x) || identical(x$status, "unavailable")
  optional <- function(x) if (is.null(x)) list(status = "unavailable", reason = "Source not supplied") else x
  source <- optional(sources$survival)
  if (!unavailable(source)) {
    retrieval$audit$survival$event_coding_verified <- TRUE
    retrieval$audit$survival$note <- paste(source$provenance$event, source$provenance$time)
  }
  validation <- analyze_r2_subgroups(retrieval)
  subtypes <- profile_attempt(profile_group_comparison(retrieval$data, "subtype",
    paste(toupper(gene), "Cavalli molecular subtypes")))
  pf <- optional(sources$pfister)
  pfister <- if (unavailable(pf)) pf else profile_attempt({
    x <- profile_group_comparison(pf$data[pf$data$include_cns, , drop = FALSE],
      "cancer_type", paste(toupper(gene), "across primary CNS tumours"), "log2(1 + FPKM)")
    x$excluded_samples <- pf$data[!pf$data$include_cns,
      c("sample_id", "cancer_type", "sample_type", "exclusion_reason"), drop = FALSE]
    x$all_counts <- as.data.frame(table(pf$data$cancer_type, pf$data$sample_type))
    names(x$all_counts) <- c("cancer_type", "sample_type", "n")
    x$provenance <- pf$provenance
    x$all_samples <- pf$data
    x
  })
  en <- optional(sources$enrichr)
  dep <- optional(sources$depmap)
  if (!unavailable(dep)) dep <- profile_attempt(plot_public_depmap(dep))
  clinical <- profile_clinical(retrieval, if (unavailable(source)) NULL else source)
  pubs <- optional(sources$publications)
  pub_error <- if (unavailable(pubs)) pubs$reason else NULL
  if (unavailable(pubs)) pubs <- NULL
  modules <- list(cavalli_subtypes = subtypes, pfister_cns = pfister,
    survival_source = source, enrichr = en, depmap = dep)
  status <- data.frame(module = names(modules),
    status = vapply(modules, function(x) if (unavailable(x)) "unavailable" else "ok", character(1)),
    detail = vapply(modules, function(x) x$reason %||% "", character(1)))
  status <- rbind(status, data.frame(module = "pubmed",
    status = if (is.null(pubs)) "unavailable" else "ok", detail = pub_error %||% ""))
  profile <- list(subtypes = subtypes, pfister = pfister, clinical = clinical,
    enrichr = en, depmap = dep, status = status,
    survival_provenance = source$provenance %||% source)
  list(retrieval = retrieval, validation = validation, publications = pubs,
    publication_error = pub_error, profile = profile, cutoff = cutoff)
}

analyze_gene <- function(cfg, gene, reporter = NULL, publication_limit = 10L,
    publication_disease = NULL, aliases = character(), include_aliases = TRUE,
    refresh_publications = FALSE, cutoff = "median", pfister_reporter = NULL) {
  cutoff <- match.arg(cutoff, c("median", "mean", "optimal"))
  sources <- collect_gene_sources(cfg, gene, reporter, publication_limit,
    publication_disease, aliases, include_aliases, refresh_publications, pfister_reporter)
  analyze_gene_sources(sources, cutoff)
}
