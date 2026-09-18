study_tool_definitions <- function() {
  s <- function(description) list(type = "string", minLength = 1L, description = description)
  arr <- function(description) list(type = "array", items = list(type = "string", minLength = 1L), uniqueItems = TRUE, description = description)
  n <- function(description, default, minimum = NULL, maximum = NULL) {
    p <- list(type = "number", description = description, default = default)
    if (!is.null(minimum)) p$minimum <- minimum
    if (!is.null(maximum)) p$maximum <- maximum
    p
  }
  int <- function(description, default, minimum, maximum) { x <- n(description, default, minimum, maximum); x$type <- "integer"; x }
  b <- function(description) list(type = "boolean", default = FALSE, description = description)
  choice <- function(values, default) list(type = "string", enum = as.list(values), default = default)
  define <- function(name, description, properties = list(), required = character())
    list(name = name, description = description, inputSchema = object_schema(properties, required),
      annotations = list(readOnlyHint = FALSE, destructiveHint = FALSE, openWorldHint = name %in% c("list_enrichr_libraries", "enrich_pathways", "query_gene_context", "r2_subgroup_analysis")))
  paths <- list(model_csv = s("Pinned DepMap Model.csv path, or DEPMAP_MODEL_CSV"), gene_effect_csv = s("Pinned CRISPRGeneEffect.csv path, or DEPMAP_GENE_EFFECT_CSV"),
    model_regex = s("Regex selecting relevant models; default is case-insensitive disease name. To assert subgroup dependency, select verified subgroup models."), depmap_release = s("Pinned release label, or DEPMAP_RELEASE"))
  thresholds <- list(fdr = n("Maximum adjusted P per supplied comparator", 0.05, 0.0000001, 1), min_abs_log2fc = n("Minimum absolute log2 fold change in every comparator", 1, 0))
  list(
    define("import_de_results", "Import existing human bulk or biological-replicate pseudobulk differential-expression results (CSV/TSV). Positive log2FC MUST mean target minus comparator. Supply ALL tested genes or an explicit tested-gene background; significant-only lists are not valid backgrounds.",
      c(list(path = s("CSV/TSV path"), disease = s("Disease label, e.g. medulloblastoma"), target = s("Target subgroup, e.g. group3"), label = s("Study label"),
        gene_column = s("Gene symbol column; default gene"), effect_column = s("Log2FC column; default log2fc; e.g. log2FoldChange"),
        adjusted_p_column = s("Adjusted P column; default adjusted_p; e.g. padj"), comparator_column = s("Comparator column; default comparator"),
        comparator = s("Use this label for every row in a single-contrast file"), background_path = s("Optional CSV/TSV with gene column listing all tested genes")), thresholds), c("path", "disease", "target")),
    define("register_cohort", "Register local expression and matched clinical metadata, including exports from R2. Stores a local snapshot. One independent biological sample per column; no cells-as-replicates analysis.",
      list(expression_path = s("Gene by sample CSV/TSV; human gene symbols"), metadata_path = s("CSV/TSV with sample_id and subgroup; optional survival_time, survival_event and metastasis (0/1)"),
        disease = s("Disease label"), label = s("Cohort label"), scale = choice(c("log2", "tpm", "counts"), "log2"), gene_column = s("Default gene"),
        sample_column = s("Default sample_id"), subgroup_column = s("Default subgroup"), source = s("Source and accession, e.g. R2 Cavalli GSE85217 export")), c("expression_path", "metadata_path", "disease")),
    define("list_studies", "List imported cohorts and differential-expression analyses with identifiers and provenance."),
    define("analyze_subgroup", "Generate target-versus-every-other-subgroup contrasts from a local matrix using limma (TMM/voom for counts). At least three biological samples per group. Batch/clinical covariates are supported.",
      c(list(cohort_id = s("Registered cohort identifier"), target = s("Exact subgroup label"), covariates = arr("Metadata covariates, e.g. batch and age")), thresholds), c("cohort_id", "target")),
    define("analysis_results", "Retrieve paginated gene candidates and their individual contrasts. A single supplied comparator is not called subgroup-specific.",
      list(analysis_id = s("Analysis identifier"), only_candidates = b("Keep genes passing every comparator"), direction = choice(c("both", "up", "down"), "both"), offset = int("Rows to skip", 0L, 0L, 1000000L), limit = int("Page size", 100L, 1L, 500L)), "analysis_id"),
    define("clinical_associations", "Fit continuous-expression Cox survival and/or logistic metastasis associations locally. Report effects per SD, confidence intervals, model failures and BH correction within requested panel. Requires patient-level data.",
      list(cohort_id = s("Registered cohort identifier"), genes = arr("1-500 gene symbols"), subgroup = s("Optional exact subgroup filter"), covariates = arr("Metadata adjustment covariates; e.g. age, sex, subgroup"),
        time_column = s("Survival time; default survival_time; keep a consistent unit"), event_column = s("0=censored, 1=event; default survival_event"), metastasis_column = s("0=absent, 1=present; default metastasis"), endpoints = arr("survival and/or metastasis")), c("cohort_id", "genes")),
    define("list_enrichr_libraries", "Discover exact Enrichr library names from the cached/live catalogue.", list(pattern = list(type = "string", description = "Optional name regex, e.g. Reactome|KEGG|Hallmark"), refresh = b("Refresh catalogue"))),
    define("enrich_pathways", "Download Enrichr annotations and perform local, measured-background over-representation separately on up/down genes passing every comparator. Does not upload patient data or gene lists; does not infer pathway activation.",
      list(analysis_id = s("Analysis identifier"), libraries = arr("Exact Enrichr library names; choose with list_enrichr_libraries"), gmt_path = s("Local GMT alternative for offline analysis"),
        min_set_size = int("Minimum term size inside tested background", 5L, 1L, 100000L), max_set_size = int("Maximum term size", 2000L, 1L, 100000L),
        fdr = n("Maximum term BH adjusted P", 0.05, 0, 1), limit = int("Terms per direction/library", 50L, 1L, 500L), refresh = b("Redownload library")), "analysis_id"),
    define("query_gene_context", "Combine gene expression, imported differential expression, clinical associations, DepMap dependency, optional R2 validation and Enrichr pathway memberships. Missing sources are explicit unavailable sections.",
      c(list(gene = s("Human gene symbol"), disease = s("Disease label"), analysis_id = s("Optional imported analysis"), cohort_id = s("Optional patient-level cohort"), subgroup = s("Optional clinical subgroup"),
        covariates = arr("Clinical adjustment covariates"), libraries = arr("Enrichr libraries for membership (not single-gene enrichment)"), gmt_path = s("Local GMT alternative"),
        include_r2 = b("Fetch gene's contrast from configured R2 cohort"), r2_group_1 = s("R2 target group"), r2_group_2 = s("R2 comparator")), paths), c("gene", "disease")),
    define("prioritize_candidates", "Rank consistently deregulated genes with DepMap dependency and selectivity evidence, optionally cross-validate against an R2/imported reference analysis and add clinical associations for the returned panel.",
      c(list(analysis_id = s("Analysis identifier"), dependency_cutoff = n("Mean gene effect threshold", -0.5), direction = choice(c("up", "down", "both"), "up"),
        limit = int("Candidates returned", 50L, 1L, 500L), reference_analysis_id = s("Independent reference with same disease and target label"), clinical_cohort_id = s("Optional registered clinical cohort"), clinical_covariates = arr("Clinical adjustment covariates")), paths), "analysis_id"),
    define("r2_subgroup_analysis", "Build a reusable R2 analysis from target-versus-comparator contrasts for ranking and reference validation. Capped R2 results require explicit full tested-gene background for pathway enrichment.",
      c(list(target = s("R2 target, e.g. group3"), comparators = arr("R2 groups, e.g. group4, shh, wnt"), disease = s("Must match configured R2 disease; default medulloblastoma"), track = s("Default subgroup"), background_path = s("Full tested-gene CSV/TSV, column gene"), refresh = b("Refresh R2 results")), thresholds), c("target", "comparators"))
  )
}

tool_definitions <- function() c(base_tool_definitions(), study_tool_definitions())

validate_tool_arguments <- function(name, args) {
  definitions <- tool_definitions()
  names <- vapply(definitions, `[[`, character(1), "name")
  if (!name %in% names) stop("Unknown tool: ", name, call. = FALSE)
  schema <- definitions[[match(name, names)]]$inputSchema
  if (!is.list(args) || (length(args) && is.null(names(args)))) stop("Tool arguments must be an object", call. = FALSE)
  extra <- setdiff(names(args), names(schema$properties))
  missing <- setdiff(unlist(schema$required), names(args))
  if (length(extra)) stop("Unknown arguments: ", paste(extra, collapse = ", "), call. = FALSE)
  if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  for (key in names(args)) {
    value <- args[[key]]
    prop <- schema$properties[[key]]
    valid <- switch(prop$type,
      string = is.character(value) && length(value) == 1L && !is.na(value),
      number = is.numeric(value) && length(value) == 1L && is.finite(value),
      integer = is.numeric(value) && length(value) == 1L && is.finite(value) && value == floor(value),
      boolean = is.logical(value) && length(value) == 1L && !is.na(value),
      array = (is.list(value) || is.character(value)) && all(vapply(as.list(value), function(v) is.character(v) && length(v) == 1L && !is.na(v) && nzchar(v), logical(1))), FALSE)
    if (!valid) stop("Invalid type for ", key, ": expected ", prop$type, call. = FALSE)
    if (!is.null(prop$enum) && !value %in% unlist(prop$enum)) stop("Invalid choice for ", key, call. = FALSE)
    if (!is.null(prop$minimum) && value < prop$minimum || !is.null(prop$maximum) && value > prop$maximum) stop("Out of range: ", key, call. = FALSE)
    if (!is.null(prop$minLength) && nchar(value) < prop$minLength) stop("Empty value: ", key, call. = FALSE)
  }
  invisible(TRUE)
}

dispatch_tool <- function(name, args, cfg) {
  validate_tool_arguments(name, args)
  if (name %in% vapply(study_tool_definitions(), `[[`, character(1), "name")) {
    fn <- get(name, mode = "function")
    return(do.call(fn, c(list(cfg = cfg), args)))
  }
  base_dispatch_tool(name, args, cfg)
}
