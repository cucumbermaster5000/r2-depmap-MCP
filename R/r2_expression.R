# Version 1: public R2 patient-level retrieval. No downstream clinical inference.
r2_missing <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "na", "NA")

r2_numeric <- function(x, label) {
  missing <- r2_missing(x)
  value <- suppressWarnings(as.numeric(x))
  if (any(!missing & !is.finite(value))) stop("Invalid numeric values in ", label, call. = FALSE)
  value[missing] <- NA_real_
  value
}

r2_json <- function(text, context) {
  tryCatch(jsonlite::fromJSON(text, simplifyVector = TRUE), error = function(e)
    stop("R2 returned an unexpected response for ", context, "; inspect saved source files. ", conditionMessage(e), call. = FALSE))
}

r2_script_payload <- function(html, module) {
  doc <- xml2::read_html(charToRaw(html), encoding = "UTF-8")
  scripts <- xml2::xml_text(xml2::xml_find_all(doc, ".//script"))
  candidates <- scripts[grepl(module, scripts, fixed = TRUE)]
  if (length(candidates) != 1L) stop("Expected exactly one R2 data payload for ", module, call. = FALSE)
  marker <- ".then(m => m['default']("
  start <- regexpr(marker, candidates, fixed = TRUE)[1]
  if (start < 0) stop("R2 script wrapper changed: ", module, call. = FALSE)
  payload <- trimws(substring(candidates, start + nchar(marker)))
  if (!endsWith(payload, "));")) stop("Unexpected R2 script suffix", call. = FALSE)
  r2_json(substr(payload, 1, nchar(payload) - 3L), module)
}

parse_r2_patient_expression <- function(html, gene, reporter, dataset, transformation = "transform_log2") {
  table <- r2_script_payload(html, "/loadDataTableModal.js")$tableData
  plot <- r2_script_payload(html, "/d3/plots/plot.js")
  if (!identical(plot$table, dataset) || !identical(as.character(plot$plotData$reporter), as.character(reporter)) ||
      !identical(toupper(plot$plotData$reporterSymbol), toupper(gene))) stop("R2 returned a different dataset/gene/reporter", call. = FALSE)
  if (!identical(plot$corType, transformation)) stop("R2 did not confirm the requested expression transformation", call. = FALSE)
  if (!is.data.frame(table$data) || !all(c("samplenames", "x") %in% names(table$data))) stop("R2 expression table schema changed", call. = FALSE)
  d <- table$data
  ids <- as.character(d$samplenames)
  if (any(r2_missing(ids)) || anyDuplicated(ids)) stop("R2 expression sample IDs are missing or duplicated", call. = FALSE)
  data.frame(queried_gene = toupper(gene), reporter = as.character(reporter), sample_id = ids,
    expression = r2_numeric(d$x, "expression"), expression_source_value = as.character(d$x), stringsAsFactors = FALSE)
}

r2_metadata_audit <- function(metadata, dataset_info, track_info) {
  rows <- lapply(names(metadata), function(field) {
    value <- as.character(metadata[[field]])
    list(field = field, source_type = "string", n = length(value), missing_n = sum(r2_missing(value)),
      unique_values = if (length(unique(value)) <= 30L) as.list(sort(unique(value))) else NULL,
      value_counts = if (length(unique(value)) <= 30L) as.list(table(value, useNA = "ifany")) else NULL)
  })
  subtypes <- if (all(c("subgroup", "subtype") %in% names(metadata))) {
    as.data.frame(table(subgroup = metadata$subgroup, subtype = metadata$subtype), stringsAsFactors = FALSE)
  } else NULL
  if (!is.null(subtypes)) subtypes <- subtypes[subtypes$Freq > 0, , drop = FALSE]
  met_field <- "met_status_(1_met__0_m0)"
  verified <- met_field %in% names(metadata) && any(grepl(paste0(met_field, ":"), track_info$trackDescriptions, fixed = TRUE))
  if (verified && any(!r2_missing(metadata[[met_field]]) & !metadata[[met_field]] %in% c("0", "1")))
    stop("Unexpected metastasis codes despite declared 0/1 mapping; no labels derived", call. = FALSE)
  list(fields = rows, dataset = dataset_info, track_descriptions = as.list(track_info$trackDescriptions),
    subtype_by_subgroup = subtypes,
    classification = list(field = if ("subtype" %in% names(metadata)) "subtype" else NULL,
      source_publication = dataset_info$pubmedId, note = "Original R2 labels and observed subgroup/subtype cross-tabulation. No classifier, relabelling, or alternative classification inferred."),
    metastasis = list(field = if (met_field %in% names(metadata)) met_field else NULL, coding_verified = verified,
      code_0 = if (verified) "M0 (no metastasis)" else NULL, code_1 = if (verified) "Metastatic" else NULL,
      evidence = if (verified) "R2's exact source field name met_status_(1_met__0_m0), present in both annotation data and dataset track descriptions" else "Coding unresolved; labels not generated",
      missing_definition = "Source token na is retained in raw data and treated as missing only in derived fields"),
    survival = list(time_field = if ("os_(years)" %in% names(metadata)) "os_(years)" else NULL,
      time_unit = if ("os_(years)" %in% names(metadata)) "years (explicit in source field name)" else NULL,
      status_field = if ("dead" %in% names(metadata)) "dead" else NULL,
      event_coding_verified = FALSE, note = "Observed dead codes are retained verbatim. An independent codebook for event semantics has not been retrieved; no survival model or event recoding is performed."),
    identifier_note = dataset_info$design,
    missing_tokens = as.list(c("", "na", "NA")))
}

assemble_r2_patient_data <- function(expression, annotation_payload, dataset_info, track_info) {
  metadata <- annotation_payload$data
  if (!is.data.frame(metadata) || !"samplenames" %in% names(metadata)) stop("R2 metadata schema changed", call. = FALSE)
  if (any(r2_missing(metadata$samplenames)) || anyDuplicated(metadata$samplenames)) stop("R2 metadata IDs are missing or duplicated", call. = FALSE)
  # Exact identifiers only. GEO expression accession is NOT the R2 sample key.
  if (!setequal(expression$sample_id, metadata$samplenames)) stop("Expression/metadata sample ID sets differ; refusing a partial or positional join", call. = FALSE)
  if (nrow(metadata) != as.integer(dataset_info$sampleSize)) stop("Retrieved samples differ from R2 dataset's declared sample count", call. = FALSE)
  aligned <- metadata[match(expression$sample_id, metadata$samplenames), , drop = FALSE]
  audit <- r2_metadata_audit(metadata, dataset_info, track_info)
  derived <- expression
  get_raw <- function(field) if (field %in% names(aligned)) as.character(aligned[[field]]) else rep(NA_character_, nrow(aligned))
  derived$subgroup <- get_raw("subgroup")
  derived$subtype <- get_raw("subtype")
  derived$metastasis_raw <- get_raw(audit$metastasis$field %||% "")
  derived$metastasis_label <- NA_character_
  if (isTRUE(audit$metastasis$coding_verified)) {
    derived$metastasis_label[which(derived$metastasis_raw == "0")] <- audit$metastasis$code_0
    derived$metastasis_label[which(derived$metastasis_raw == "1")] <- audit$metastasis$code_1
  }
  derived$survival_time_source_value <- get_raw("os_(years)")
  derived$survival_time_years <- r2_numeric(derived$survival_time_source_value, "os_(years)")
  derived$survival_status_raw <- get_raw("dead")
  derived$survival_event <- NA_integer_
  list(data = derived, metadata = metadata, audit = audit)
}

get_r2_expression <- function(cfg, gene, reporter = NULL, refresh = FALSE) {
  if (length(gene) != 1L || !grepl("^[A-Za-z][A-Za-z0-9._-]*$", gene)) stop("Provide one gene symbol", call. = FALSE)
  gene <- toupper(gene)
  transformation <- cfg$r2_transformation %||% 'transform_log2'
  key <- cache_key(if(transformation=='transform_log2') 'r2_patient_v1' else 'r2_patient_v2', list(dataset = cfg$r2_dataset_table, gene = gene, reporter = reporter %||% "unique", transformation = if(transformation=='transform_log2') 'log2' else transformation))
  if (!refresh) {
    cached <- cache_get(cfg, key)
    if (!is.null(cached)) { cached$provenance$cache_hit <- TRUE; return(cached) }
  }
  root <- cfg$output_dir %||% file.path(dirname(cfg$cache_dir), "..", "outputs")
  run_id <- basename(tempfile(paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-")))
  out <- file.path(root, gene, run_id)
  raw_dir <- file.path(out, "source")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(raw_dir)) stop("Cannot create output directory", call. = FALSE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)
  h <- r2_session(120)
  requests <- list()
  capture <- function(name, params = NULL, post = FALSE) {
    url <- cfg$r2_base_url
    if (!is.null(params) && !post) url <- paste0(url, "?", encode_form(params))
    text <- r2_request(h, url, if (post) params else NULL)
    path <- file.path(raw_dir, name)
    writeBin(charToRaw(text), path)
    requests[[length(requests) + 1L]] <<- list(file = name, url = url, method = if (post) "POST" else "GET", fields = if (post) params else NULL,
      retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), md5 = unname(tools::md5sum(path)))
    text
  }
  landing <- capture("landing.html", list(table = cfg$r2_dataset_table))
  if (grepl("Use R2 without an account", landing, fixed = TRUE)) capture("guest.html", list(method = "guest", open_page = "auth"), TRUE)
  capture("dataset.html", list(table = cfg$r2_dataset_table))
  info <- r2_json(capture("dataset.json", list(json_option = "json_dataset_info", table = cfg$r2_dataset_table)), "dataset")
  if (!identical(info$internalIdentifier, cfg$r2_dataset_table)) stop("Dataset identity mismatch", call. = FALSE)
  tracks <- r2_json(capture("tracks.json", list(json_option = "json_dataset_info_tracks", table = cfg$r2_dataset_table)), "tracks")
  annotations <- r2_json(capture("annotations.json", list(json_option = "json_cg_sampleannotation_v1", ctable = cfg$r2_dataset_table, subset = "")), "annotations")
  reporters <- r2_json(capture("reporters.json", list(json_option = "json_find_reporter_in_dataset_v1", dataset = cfg$r2_dataset_table, search_by = "gene_symbol", query = gene)), "reporters")
  if (!is.data.frame(reporters) || !all(c("gene_symbol", "reporter") %in% names(reporters))) stop("Gene absent from R2 reporter search", call. = FALSE)
  reporters <- reporters[toupper(reporters$gene_symbol) == gene, , drop = FALSE]
  if (!is.null(reporter)) reporters <- reporters[as.character(reporters$reporter) == reporter, , drop = FALSE]
  if (nrow(reporters) != 1L) stop("Expected one exact reporter for ", gene, "; found ", nrow(reporters), ". Specify reporter explicitly if ambiguous.", call. = FALSE)
  reporter <- as.character(reporters$reporter[1])
  html <- capture("expression.html", list(table = cfg$r2_dataset_table, option = "display2", analysis_type = "single_reporter",
    factor = reporter, cortype = transformation, subset = "", grouping_track = cfg$r2_grouping_track %||% 'subgroup', graphtype = "yy_i"), TRUE)
  expression <- parse_r2_patient_expression(html, gene, reporter, cfg$r2_dataset_table, transformation)
  result <- assemble_r2_patient_data(expression, annotations, info, tracks)
  result$provenance <- list(dataset_table = cfg$r2_dataset_table, accession = info$accession, dataset_label = info$datasetTitle,
    queried_gene = gene, reporter = reporter, transformation = if(transformation=='transform_log2') 'R2 transform_log2 (log2 of source expression)' else 'R2 transform_none (untransformed source expression)',
    sample_join = "Exact equality: expression table samplenames to annotation samplenames; no positional or GEO-expression-ID join",
    cache_hit = FALSE, retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    output_dir = out, requests = requests, limitations = "Unofficial read-only R2 interface; validated against observed responses. No clinical/subtype inference. All samples retained.")
  result$files <- list(patient_data = file.path(out, paste0(gene, "_R2_patient_data.csv")), metadata = file.path(out, paste0(gene, "_R2_metadata.csv")),
    metadata_audit = file.path(out, "metadata_audit.json"), provenance = file.path(out, "provenance.json"), r_object = file.path(out, "retrieval.rds"))
  data.table::fwrite(result$data, result$files$patient_data, na = "NA")
  data.table::fwrite(result$metadata, result$files$metadata, na = "NA")
  writeLines(json_text(result$audit), result$files$metadata_audit)
  writeLines(json_text(result$provenance), result$files$provenance)
  saveRDS(result, result$files$r_object)
  writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
  cache_put(cfg, key, result)
  result
}

analyze_r2_subgroups <- function(result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 for the Version 1 plot", call. = FALSE)
  d <- result$data
  eligible <- is.finite(d$expression) & !r2_missing(d$subgroup)
  exclusions <- d[!eligible, c("sample_id", "expression", "subgroup"), drop = FALSE]
  used <- d[eligible, , drop = FALSE]
  sizes <- table(used$subgroup)
  if (length(sizes) < 2 || any(sizes < 3)) stop("Need at least two annotated subgroups with three samples each", call. = FALSE)
  if (anyDuplicated(used$sample_id)) stop("Repeated sample IDs; unpaired comparison is invalid", call. = FALSE)
  # Predefined descriptive comparison: heterogeneous independent tumour cohorts,
  # non-parametric ranks; no normality-test-driven selection or post-hoc gene claims.
  test <- stats::kruskal.test(expression ~ subgroup, used)
  summaries <- do.call(rbind, lapply(split(used$expression, used$subgroup), function(x)
    data.frame(n = length(x), mean = mean(x), median = stats::median(x), sd = stats::sd(x), iqr = stats::IQR(x), min = min(x), max = max(x))))
  summaries$subgroup <- rownames(summaries)
  rownames(summaries) <- NULL
  statistics <- data.frame(test = test$method, statistic = unname(test$statistic), df = unname(test$parameter), p_value = test$p.value,
    n = nrow(used), groups = length(sizes), excluded_n = nrow(exclusions),
    epsilon_squared = max(0, (unname(test$statistic) - length(sizes) + 1) / (nrow(used) - length(sizes))),
    multiple_testing = "One predefined omnibus test for one gene; no pairwise tests or multiplicity adjustment in Version 1")
  labels <- setNames(paste0(names(sizes), "\n(n = ", as.integer(sizes), ")"), names(sizes))
  p <- ggplot2::ggplot(used, ggplot2::aes(x = subgroup, y = expression, fill = subgroup)) +
    ggplot2::geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.35, linewidth = 0.45) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.16, height = 0, seed = 85217), size = 1.1, alpha = 0.5, shape = 16) +
    ggplot2::scale_x_discrete(labels = labels) + ggplot2::guides(fill = "none") +
    ggplot2::labs(title = paste(unique(used$queried_gene), "expression across medulloblastoma subgroups"),
      subtitle = sprintf("Kruskal-Wallis H(%d) = %.2f; p = %.3g", test$parameter, test$statistic, test$p.value),
      x = "Original R2 subgroup labels", y = "Expression (R2 log2 transformation)",
      caption = paste0("GSE85217 / Cavalli; reporter ", result$provenance$reporter, ". Each point is one R2 sample.\n",
        "Overall distribution comparison; this test does not establish subgroup specificity. Excluded: ", nrow(exclusions), ".")) +
    ggplot2::theme_classic(base_size = 12) + ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 9))
  list(statistics = statistics, summaries = summaries, exclusions = exclusions, plot = p)
}

# Compatibility/export service; statistical calculations are in the core above.
plot_r2_subgroup_test <- function(result, output_dir = result$provenance$output_dir, export_images = FALSE) {
  analysis <- analyze_r2_subgroups(result)
  statistics <- analysis$statistics
  summaries <- analysis$summaries
  exclusions <- analysis$exclusions
  p <- analysis$plot
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  gene <- result$provenance$queried_gene
  files <- list(statistics = file.path(output_dir, paste0(gene, "_subgroup_test.csv")), summaries = file.path(output_dir, paste0(gene, "_subgroup_summaries.csv")), exclusions = file.path(output_dir, "subgroup_test_exclusions.csv"))
  if (isTRUE(export_images)) {
    files$pdf <- file.path(output_dir, paste0(gene, "_subgroup_test.pdf"))
    files$png <- file.path(output_dir, paste0(gene, "_subgroup_test.png"))
    ggplot2::ggsave(files$pdf, p, width = 8, height = 5.8, device = grDevices::pdf)
    ggplot2::ggsave(files$png, p, width = 8, height = 5.8, dpi = 300)
  }
  data.table::fwrite(statistics, files$statistics)
  data.table::fwrite(summaries, files$summaries)
  data.table::fwrite(exclusions, files$exclusions)
  writeLines("Predefined unpaired Kruskal-Wallis comparison of distributions. Samples are treated as independent primary tumours as described by the source. Batch, age and other confounding are not adjusted in this Version 1 smoke test. No pairwise or subtype/clinical analyses were performed. Missing annotations/expression are reported in the exclusions file. Epsilon-squared is the non-negative rank-based omnibus effect estimate (H-k+1)/(n-k).", file.path(output_dir, "subgroup_test_method.txt"))
  dump(c("r2_missing", "analyze_r2_subgroups", "plot_r2_subgroup_test"), file = file.path(output_dir, "analysis_functions.R"), envir = environment(plot_r2_subgroup_test))
  writeLines(c("# Open this output directory as the working directory in RStudio.",
    "# Uses the saved patient data; no network request and no new clinical analysis.",
    "source('analysis_functions.R')", "retrieval <- readRDS('retrieval.rds')",
    "validation <- plot_r2_subgroup_test(retrieval, output_dir = 'rerun')", "print(validation$plot)", "validation$statistics"), file.path(output_dir, "reproduce_plot.R"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "analysis_session-info.txt"))
  list(statistics = statistics, summaries = summaries, files = files, plot = p)
}
