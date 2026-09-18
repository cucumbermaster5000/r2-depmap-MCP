v1_tool_definitions <- function() {
  gene <- list(type = "string", minLength = 1, description = "One human gene symbol, e.g. HLX")
  reporter <- list(type = "string", description = "Optional exact R2 reporter; required if a gene maps to multiple reporters")
  aliases <- list(type = "array", items = list(type = "string"), description = "Additional explicitly supplied gene aliases; NCBI human Gene aliases are also included by default")
  schema <- function(properties, required = character()) {
    x <- list(type = "object", additionalProperties = FALSE)
    if (length(properties)) x$properties <- properties
    if (length(required)) x$required <- as.list(required)
    x
  }
  list(
    list(name = "describe_r2_dataset", description = "Describe Version 1 scope and configured R2 source. Clinical coding and subtype labels are verified during retrieval, not assumed.", inputSchema = schema(list())),
    list(name = "get_r2_expression", description = "Retrieve patient-level log2 expression and original metadata from the configured R2 medulloblastoma cohort. Preserves source responses, identifiers, clinical raw codes and provenance locally. Returns paginated rows and paths to the full data. No survival, metastasis, subtype or enrichment analysis.",
      inputSchema = schema(list(gene = gene, reporter = reporter,
        refresh = list(type = "boolean", default = FALSE, description = "Retrieve a new source snapshot instead of using the cached result"),
        offset = list(type = "integer", minimum = 0, default = 0), limit = list(type = "integer", minimum = 1, maximum = 1000, default = 20)), "gene"),
      annotations = list(readOnlyHint = FALSE, destructiveHint = FALSE, openWorldHint = TRUE)),
    list(name = "plot_r2_subgroup_test", description = "Generate a subgroup plot object, summaries, a Kruskal-Wallis test and exclusions. Use build_gene_page to view the plot in your browser; PDF/PNG exports are optional in R.",
      inputSchema = schema(list(gene = gene, reporter = reporter), "gene"), annotations = list(readOnlyHint = FALSE, destructiveHint = FALSE, openWorldHint = TRUE)),
    list(name = "get_gene_publications", description = "Get the most recent PubMed gene papers sorted by publication date: title, journal, citation/electronic date, DOI, PMID and one relevant verbatim abstract sentence (title-only fallback when no abstract). Verified human NCBI Gene aliases included by default. Searches all contexts unless disease is provided. Cache expires in 24 hours; refresh bypasses it.",
      inputSchema = schema(list(gene = gene, limit = list(type = "integer", minimum = 1, maximum = 30, default = 10),
        disease = list(type = "string", description = "Optional explicit title/abstract disease restriction, e.g. medulloblastoma"), aliases = aliases,
        include_aliases = list(type = "boolean", default = TRUE), refresh = list(type = "boolean", default = FALSE)), "gene"),
      annotations = list(readOnlyHint = FALSE, destructiveHint = FALSE, openWorldHint = TRUE)),
    list(name = "get_pfister_expression", description = "Retrieve all 272 Pfister FPKM samples with original cancer labels and explicit primary CNS selection; values for plotting are log2(1+FPKM).", inputSchema = schema(list(gene=gene,reporter=reporter), 'gene')),
    list(name = "get_depmap_dependencies", description = "Query the official pinned DepMap 24Q4 v1 archive for a gene's Chronos dependencies, highlighting explicitly annotated medulloblastoma models. First use downloads 429 MB once; later queries use the local cache. This is not the latest release.", inputSchema = schema(list(gene=gene), 'gene')),
    list(name = "build_gene_page", description = "Build the integrated browser gene report: Cavalli subgroups and subtypes, Pfister CNS expression, verified survival and metastasis (median/mean/optimised cutoffs), Enrichr functional memberships and TF target enrichment, public DepMap dependencies and PubMed. Missing sources are explicit; no pathway activity or therapeutic score is inferred.",
      inputSchema = schema(list(gene = gene, reporter = reporter, publication_limit = list(type = "integer", minimum = 1, maximum = 30, default = 10),
        publication_disease = list(type = "string", description = "Optional publication disease restriction; by default gene literature across all contexts"), aliases = aliases,
        include_aliases = list(type = "boolean", default = TRUE), refresh_publications = list(type = "boolean", default = FALSE),
        cutoff=list(type='string',enum=as.list(c('median','mean','optimal')),default='median',description='Initially selected cutoff in the HTML; all three modes are included'),pfister_reporter=reporter), "gene"),
      annotations = list(readOnlyHint = FALSE, destructiveHint = FALSE, openWorldHint = TRUE))
  )
}

dispatch_v1_tool <- function(name, args, cfg) {
  definitions <- v1_tool_definitions()
  tool_names <- vapply(definitions, `[[`, character(1), "name")
  if (!name %in% tool_names) stop("Tool is not enabled in Version 1: ", name, call. = FALSE)
  schema <- definitions[[match(name, tool_names)]]$inputSchema
  if (!is.list(args) || (length(args) && is.null(names(args)))) stop("Arguments must be a named object", call. = FALSE)
  if (length(setdiff(names(args), names(schema$properties)))) stop("Unknown argument", call. = FALSE)
  if (length(setdiff(unlist(schema$required), names(args)))) stop("Missing required gene argument", call. = FALSE)
  for (key in names(args)) {
    v <- args[[key]]; p <- schema$properties[[key]]
    valid <- switch(p$type, string = is.character(v) && length(v) == 1 && !is.na(v) && nzchar(v),
      boolean = is.logical(v) && length(v) == 1 && !is.na(v),
      integer = is.numeric(v) && length(v) == 1 && is.finite(v) && v == floor(v),
      array = (is.list(v) || is.character(v)) && all(vapply(as.list(v), function(x) is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x), logical(1))), FALSE)
    if (!valid) stop("Invalid argument type: ", key, call. = FALSE)
    if(!is.null(p$enum)&&!v %in% unlist(p$enum))stop('Invalid choice: ',key,call.=FALSE)
    if ((!is.null(p$minimum) && v < p$minimum) || (!is.null(p$maximum) && v > p$maximum)) stop("Argument out of range: ", key, call. = FALSE)
  }
  if (name == "describe_r2_dataset") return(list(version = "2", dataset_table = cfg$r2_dataset_table, dataset_label = cfg$r2_dataset_label,
    endpoint = cfg$r2_base_url, source_accession = "GSE85217", scope = "Cavalli subgroups/subtypes, Pfister CNS expression, verified survival/metastasis, Enrichr/TF targets, archived public DepMap and recent PubMed in one HTML report",
    enabled_tools = as.list(tool_names), scope_note = "Integrated single-gene report; no composite scores or inferred pathway activity"))
  if(name=='get_pfister_expression') {x<-do.call(get_pfister_expression,c(list(cfg=cfg),args));return(list(rows=x$data,provenance=x$provenance))}
  if(name=='get_depmap_dependencies')return(get_public_depmap_gene(cfg,args$gene))
  if (name == "get_gene_publications") {
    result <- do.call(get_gene_publications, c(list(cfg = cfg), args))
    result$articles <- lapply(result$articles, function(x) { x$abstract <- NULL; x })
    return(result)
  }
  if (name == "build_gene_page") return(do.call(build_gene_page, c(list(cfg = cfg), args)))
  r <- get_r2_expression(cfg, args$gene, args$reporter, args$refresh %||% FALSE)
  if (name == "plot_r2_subgroup_test") {
    test <- plot_r2_subgroup_test(r)
    test$plot <- NULL
    return(test)
  }
  offset <- args$offset %||% 0L; limit <- args$limit %||% 20L
  list(total_samples = nrow(r$data), offset = offset,
    rows = utils::head(r$data[seq_len(nrow(r$data)) > offset, , drop = FALSE], limit),
    full_patient_data_file = r$files$patient_data, files = r$files, metadata_audit = r$audit, provenance = r$provenance)
}
