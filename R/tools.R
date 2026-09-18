object_schema <- function(properties = list(), required = character()) {
  out <- list(type = "object", additionalProperties = FALSE)
  if (length(properties)) out$properties <- properties
  if (length(required)) out$required <- as.list(required)
  out
}

base_tool_definitions <- function() {
  string_prop <- function(description, default = NULL) {
    out <- list(type = "string", description = description)
    if (!is.null(default)) out$default <- default
    out
  }
  number_prop <- function(description, default = NULL, minimum = NULL, maximum = NULL) {
    out <- list(type = "number", description = description)
    if (!is.null(default)) out$default <- default
    if (!is.null(minimum)) out$minimum <- minimum
    if (!is.null(maximum)) out$maximum <- maximum
    out
  }
  integer_prop <- function(description, default = NULL, minimum = NULL, maximum = NULL) {
    out <- number_prop(description, default, minimum, maximum)
    out$type <- "integer"
    out
  }
  boolean_prop <- function(description, default = FALSE) {
    list(type = "boolean", description = description, default = default)
  }

  list(
    list(
      name = "describe_r2_dataset",
      description = "Describe the configured public R2 Cavalli medulloblastoma cohort and analysis conventions.",
      inputSchema = object_schema()
    ),
    list(
      name = "r2_compare_groups",
      description = "Run or retrieve a cached two-group R2 differential-expression analysis. Positive patient_effect_group1_minus_group2 means group 1 is higher.",
      inputSchema = object_schema(list(
        track = string_prop("R2 categorical annotation track.", "subgroup"),
        group_1 = string_prop("First R2 group.", "group3"),
        group_2 = string_prop("Second R2 group.", "group4"),
        test = list(type = "string", enum = c("anova", "kruskal", "log2fc", "limma"), default = "anova"),
        top_n = integer_prop("Maximum R2 rows to calculate/cache.", 1000L, 1L, 5000L),
        return_n = integer_prop("Maximum rows returned to the MCP client.", 100L, 1L, 500L),
        p_threshold = number_prop("R2 FDR threshold.", 0.05, 0, 1),
        min_present = integer_prop("Minimum present-call count used by R2.", 1L, 0L),
        refresh = boolean_prop("Ignore a cached result and query R2 again.", FALSE)
      ))
    ),
    list(
      name = "depmap_gene_dependency",
      description = "Return gene-effect values for one gene across DepMap models matching a disease or lineage regular expression.",
      inputSchema = object_schema(list(
        gene = string_prop("HGNC gene symbol."),
        model_csv = string_prop("Path to the pinned DepMap model metadata CSV; may instead be set with DEPMAP_MODEL_CSV."),
        gene_effect_csv = string_prop("Path to the pinned DepMap CRISPR gene-effect CSV; may instead be set with DEPMAP_GENE_EFFECT_CSV."),
        model_regex = string_prop("Case-sensitive or inline-modified regular expression applied to model metadata.", "(?i)medulloblastoma"),
        depmap_release = string_prop("Human-readable DepMap release label."),
        dependency_cutoff = number_prop("Gene-effect threshold counted as dependent.", -0.5),
        return_n = integer_prop("Maximum model rows returned.", 100L, 1L, 500L)
      ), required = "gene")
    ),
    list(
      name = "rank_patient_dependencies",
      description = "Rank genes that are enriched in an R2 patient subgroup and dependencies in selected DepMap models.",
      inputSchema = object_schema(list(
        model_csv = string_prop("Path to DepMap model metadata CSV; may instead be set with DEPMAP_MODEL_CSV."),
        gene_effect_csv = string_prop("Path to DepMap CRISPR gene-effect CSV; may instead be set with DEPMAP_GENE_EFFECT_CSV."),
        depmap_release = string_prop("Pinned DepMap release label."),
        track = string_prop("R2 categorical annotation track.", "subgroup"),
        group_1 = string_prop("Patient group expected to be enriched.", "group3"),
        group_2 = string_prop("Patient comparator group.", "group4"),
        test = list(type = "string", enum = c("anova", "kruskal", "log2fc", "limma"), default = "anova"),
        model_regex = string_prop("Regex selecting relevant DepMap models.", "(?i)medulloblastoma"),
        r2_top_n = integer_prop("Number of R2 genes retained before joining.", 2000L, 1L, 5000L),
        result_top_n = integer_prop("Number of ranked targets returned.", 50L, 1L, 500L),
        p_threshold = number_prop("Maximum R2 adjusted P value.", 0.05, 0, 1),
        min_patient_effect = number_prop("Minimum group-1-minus-group-2 expression effect.", 0),
        dependency_cutoff = number_prop("DepMap gene-effect threshold counted as dependent.", -0.5),
        expression_weight = number_prop("Weight for patient expression enrichment.", 0.45, 0),
        dependency_weight = number_prop("Weight for dependency strength.", 0.40, 0),
        selectivity_weight = number_prop("Weight for dependency selectivity.", 0.15, 0),
        refresh_r2 = boolean_prop("Ignore cached R2 results.", FALSE)
      ))
    ),
    list(
      name = "cache_status",
      description = "List cached R2 differential-expression comparisons and retrieval dates.",
      inputSchema = object_schema()
    )
  )
}

compact_r2_result <- function(result, return_n) {
  list(
    rows = utils::head(result$data, return_n),
    returned_n = min(return_n, nrow(result$data)),
    available_n = nrow(result$data),
    details = result$details,
    provenance = result$provenance
  )
}

base_dispatch_tool <- function(name, args, cfg) {
  args <- args %||% list()
  if (name == "describe_r2_dataset") {
    return(list(
      dataset_table = cfg$r2_dataset_table,
      dataset_label = cfg$r2_dataset_label,
      public_url = paste0(cfg$r2_base_url, "?table=", cfg$r2_dataset_table),
      default_track = cfg$r2_default_track,
      supported_workflow = "Cached differential expression and multi-comparator analysis; import sample-level exports for survival/metastasis",
      important_note = "R2 Log2FC is converted so positive patient_effect_group1_minus_group2 means group 1 is higher.",
      connector = "unofficial read-only HTML workflow; no R2 credentials are used"
    ))
  }

  if (name == "r2_compare_groups") {
    result <- r2_compare_groups(
      cfg = cfg,
      track = args$track,
      group_1 = args$group_1,
      group_2 = args$group_2,
      test = args$test %||% "anova",
      top_n = args$top_n %||% 1000L,
      p_threshold = args$p_threshold %||% 0.05,
      min_present = args$min_present %||% 1L,
      refresh = as_logical_scalar(args$refresh, FALSE)
    )
    return(compact_r2_result(result, as_integer_scalar(args$return_n, 100L, 1L, 500L)))
  }

  if (name == "depmap_gene_dependency") {
    result <- depmap_gene_dependency(
      gene = args$gene,
      model_csv = args$model_csv,
      gene_effect_csv = args$gene_effect_csv,
      model_regex = args$model_regex %||% cfg$depmap_default_model_regex,
      release = args$depmap_release,
      dependency_cutoff = args$dependency_cutoff %||% -0.5
    )
    return_n <- as_integer_scalar(args$return_n, 100L, 1L, 500L)
    result$models <- utils::head(result$models, return_n)
    return(result)
  }

  if (name == "rank_patient_dependencies") {
    return(rank_patient_dependencies(
      cfg = cfg,
      model_csv = args$model_csv,
      gene_effect_csv = args$gene_effect_csv,
      depmap_release = args$depmap_release,
      track = args$track,
      group_1 = args$group_1,
      group_2 = args$group_2,
      test = args$test %||% "anova",
      model_regex = args$model_regex,
      r2_top_n = args$r2_top_n %||% 2000L,
      result_top_n = args$result_top_n %||% 50L,
      p_threshold = args$p_threshold %||% 0.05,
      min_patient_effect = args$min_patient_effect %||% 0,
      dependency_cutoff = args$dependency_cutoff %||% -0.5,
      expression_weight = args$expression_weight %||% 0.45,
      dependency_weight = args$dependency_weight %||% 0.40,
      selectivity_weight = args$selectivity_weight %||% 0.15,
      refresh_r2 = as_logical_scalar(args$refresh_r2, FALSE)
    ))
  }

  if (name == "cache_status") {
    return(list(cache_dir = cfg$cache_dir, entries = cache_inventory(cfg)))
  }

  stop("Unknown tool: ", name, call. = FALSE)
}
