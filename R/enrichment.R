enrichr_get <- function(url) {
  h <- curl::new_handle(timeout = 90, connecttimeout = 20, followlocation = TRUE, useragent = "r2-depmap-mcp/0.2.0")
  response <- curl::curl_fetch_memory(url, handle = h)
  if (response$status_code != 200) stop("Enrichr returned HTTP ", response$status_code, call. = FALSE)
  rawToChar(response$content)
}

list_enrichr_libraries <- function(cfg, pattern = "", refresh = FALSE) {
  key <- "enrichr_library_catalog"
  cached <- if (!refresh) cache_get(cfg, key) else NULL
  if (is.null(cached)) {
    url <- "https://maayanlab.cloud/Enrichr/datasetStatistics"
    payload <- jsonlite::fromJSON(enrichr_get(url))
    if (!is.data.frame(payload$statistics) || !"libraryName" %in% names(payload$statistics)) stop("Unexpected Enrichr library catalogue", call. = FALSE)
    cached <- list(data = payload$statistics, provenance = list(source = url, retrieved_at = utc_now()))
    cache_put(cfg, key, cached)
  }
  list(libraries = cached$data[grepl(pattern, cached$data$libraryName, ignore.case = TRUE), , drop = FALSE], provenance = cached$provenance)
}

parse_gmt <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  if (!length(fields) || any(lengths(fields) < 3)) stop("Invalid GMT: expected term, description and gene symbols", call. = FALSE)
  result <- list()
  for (f in fields) {
    term <- f[1]
    genes <- toupper(trimws(sub(",.*$", "", f[-c(1, 2)])))
    genes <- genes[nzchar(genes)]
    result[[term]] <- unique(c(result[[term]], genes))
  }
  result
}

load_gene_sets <- function(cfg, library = NULL, gmt_path = NULL, refresh = FALSE) {
  if (!is.null(gmt_path)) {
    raw <- paste(readLines(normalize_local_path(gmt_path), warn = FALSE), collapse = "\n")
    return(list(sets = parse_gmt(raw), provenance = list(source = "local GMT", file = file_provenance(gmt_path))))
  }
  if (is.null(library) || !nzchar(library)) stop("Specify an Enrichr library or local GMT path", call. = FALSE)
  key <- paste0("enrichr_gmt_", object_digest(library))
  cached <- if (!refresh) cache_get(cfg, key) else NULL
  if (!is.null(cached)) return(cached)
  catalog <- list_enrichr_libraries(cfg, refresh = refresh)
  if (!library %in% catalog$libraries$libraryName) stop("Unknown Enrichr library; use list_enrichr_libraries with refresh=true", call. = FALSE)
  url <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", utils::URLencode(library, reserved = TRUE))
  raw <- enrichr_get(url)
  result <- list(sets = parse_gmt(raw), provenance = list(source = "Enrichr", library = library, url = url,
    content_md5 = object_digest(raw), retrieved_at = utc_now(), computation = "Local over-representation using downloaded annotations; not Enrichr combined scores"))
  cache_put(cfg, key, result)
  result
}

overrepresentation <- function(genes, universe, sets, min_set_size = 5L, max_set_size = 2000L) {
  universe <- unique(universe)
  genes <- intersect(unique(genes), universe)
  empty <- data.frame(term = character(), overlap_n = integer(), set_n = integer(), query_n = integer(),
    universe_n = integer(), fold_enrichment = numeric(), p_value = numeric(), adjusted_p = numeric(), overlap_genes = character())
  if (!length(genes) || !length(universe)) return(empty)
  sets <- lapply(sets, intersect, universe)
  sets <- sets[lengths(sets) >= min_set_size & lengths(sets) <= max_set_size]
  if (!length(sets)) return(empty)
  rows <- lapply(names(sets), function(term) {
    members <- sets[[term]]
    overlap <- intersect(genes, members)
    data.frame(term = term, overlap_n = length(overlap), set_n = length(members), query_n = length(genes), universe_n = length(universe),
      fold_enrichment = (length(overlap) / length(genes)) / (length(members) / length(universe)),
      p_value = stats::phyper(length(overlap) - 1, length(members), length(universe) - length(members), length(genes), lower.tail = FALSE),
      overlap_genes = paste(overlap, collapse = ";"), stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  result$adjusted_p <- stats::p.adjust(result$p_value, "BH")
  result[order(result$adjusted_p, result$p_value, -result$fold_enrichment), ]
}

enrich_pathways <- function(cfg, analysis_id, libraries = character(), gmt_path = NULL,
                            min_set_size = 5L, max_set_size = 2000L, fdr = 0.05, limit = 50L, refresh = FALSE) {
  a <- load_study_object(cfg, analysis_id, "analysis")
  libraries <- unlist(libraries, use.names = FALSE)
  if (!length(libraries) && is.null(gmt_path)) stop("Specify libraries or gmt_path", call. = FALSE)
  if (!length(a$universe)) stop("Analysis has no tested-gene background", call. = FALSE)
  inputs <- if (!is.null(gmt_path)) list(load_gene_sets(cfg, gmt_path = gmt_path)) else lapply(libraries, function(l) load_gene_sets(cfg, l, refresh = refresh))
  labels <- if (!is.null(gmt_path)) basename(gmt_path) else libraries
  results <- list()
  for (i in seq_along(inputs)) for (direction in c("up", "down")) {
    genes <- a$summary$gene[a$summary$passes_all_comparators & a$summary$direction == direction]
    table <- overrepresentation(genes, a$universe, inputs[[i]]$sets, min_set_size, max_set_size)
    results[[paste(labels[i], direction, sep = ":")]] <- list(direction = direction, library = labels[i],
      input_genes_n = length(genes), tested_terms_n = nrow(table), significant_terms_n = sum(table$adjusted_p <= fdr),
      rows = utils::head(table[table$adjusted_p <= fdr, , drop = FALSE], limit), provenance = inputs[[i]]$provenance)
  }
  list(analysis_id = analysis_id, results = results, background_n = length(a$universe),
    method = "One-sided hypergeometric over-representation; BH across ALL size-eligible terms separately per library and direction, including zero overlaps",
    interpretation = "Pathways enriched among consistently up/down deregulated genes; enrichment does not establish pathway activation or inhibition.",
    privacy = "Only library names are requested from Enrichr. Gene lists, expression and clinical data remain local.",
    analysis_provenance = a$provenance)
}
