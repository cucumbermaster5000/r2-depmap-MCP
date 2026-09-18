`%||%` <-
function (x, y) 
{
    if (is.null(x) || length(x) == 0L || (length(x) == 1L && 
        is.na(x))) 
        y
    else x
}
analysis_results <-
function (cfg, analysis_id, only_candidates = FALSE, direction = "both", 
    offset = 0L, limit = 100L) 
{
    a <- load_study_object(cfg, analysis_id, "analysis")
    rows <- a$summary
    if (only_candidates) 
        rows <- rows[rows$passes_all_comparators, , drop = FALSE]
    if (direction != "both") 
        rows <- rows[rows$direction == direction, , drop = FALSE]
    n <- nrow(rows)
    rows <- utils::head(rows[seq_len(n) > offset, , drop = FALSE], 
        limit)
    list(analysis = study_overview(a), total = n, offset = offset, 
        rows = rows, contrasts = a$data[a$data$gene %in% rows$gene, 
            , drop = FALSE], thresholds = a$thresholds)
}
analyze_subgroup <-
function (cfg, cohort_id, target, covariates = character(), fdr = 0.050000000000000003, 
    min_abs_log2fc = 1) 
{
    need_package("limma")
    cohort <- load_study_object(cfg, cohort_id, "cohort")
    counts <- table(cohort$metadata$subgroup)
    if (!target %in% names(counts) || length(counts) < 2 || any(counts < 
        3)) 
        stop("Need target and comparators with at least three biological replicates per subgroup", 
            call. = FALSE)
    design <- study_design(cohort$metadata, covariates)
    expr <- cohort$expression
    if (cohort$scale == "counts") {
        need_package("edgeR")
        y <- edgeR::DGEList(expr)
        keep <- edgeR::filterByExpr(y, design = design$matrix)
        if (sum(keep) < 2) 
            stop("Too few genes pass count filtering", call. = FALSE)
        y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
        fit <- limma::lmFit(limma::voom(y, design$matrix, plot = FALSE), 
            design$matrix)
    }
    else {
        keep <- apply(expr, 1, stats::sd) > 0
        if (sum(keep) < 2) 
            stop("Too few variable genes", call. = FALSE)
        fit <- limma::lmFit(expr[keep, , drop = FALSE], design$matrix)
    }
    comparators <- setdiff(design$groups, target)
    contrast <- matrix(0, ncol(design$matrix), length(comparators), 
        dimnames = list(colnames(design$matrix), comparators))
    contrast[match(target, design$groups), ] <- 1
    for (i in seq_along(comparators)) contrast[match(comparators[i], 
        design$groups), i] <- -1
    fitted <- limma::eBayes(limma::contrasts.fit(fit, contrast), 
        trend = cohort$scale != "counts")
    data <- do.call(rbind, lapply(seq_along(comparators), function(i) {
        tt <- limma::topTable(fitted, coef = i, number = Inf, 
            sort.by = "none", adjust.method = "BH")
        data.frame(gene = rownames(tt), log2fc = tt$logFC, adjusted_p = tt$adj.P.Val, 
            p_value = tt$P.Value, comparator = comparators[i])
    }))
    result <- list(label = paste(cohort$label, target), disease = cohort$disease, 
        target = target, comparators = comparators, cohort_id = cohort_id, 
        data = data, universe = unique(data$gene), summary = summarize_specificity(data, 
            comparators, fdr, min_abs_log2fc), thresholds = list(fdr = fdr, 
            min_abs_log2fc = min_abs_log2fc), provenance = list(source = cohort$provenance, 
            method = if (cohort$scale == "counts") "edgeR TMM + limma voom" else "limma trend", 
            covariates = as.list(covariates), samples_per_group = as.list(counts), 
            multiple_testing = "BH separately over all genes per contrast", 
            effect_definition = "target minus comparator", generated_at = utc_now(), 
            limma_version = as.character(utils::packageVersion("limma"))))
    result$id <- save_study_object(cfg, "analysis", result)
    study_overview(result)
}
as_integer_scalar <-
function (x, default, minimum = NULL, maximum = NULL) 
{
    out <- suppressWarnings(as.integer(x %||% default))
    if (is.na(out)) 
        out <- as.integer(default)
    if (!is.null(minimum)) 
        out <- max(out, minimum)
    if (!is.null(maximum)) 
        out <- min(out, maximum)
    out
}
as_logical_scalar <-
function (x, default = FALSE) 
{
    if (is.null(x) || !length(x)) 
        return(default)
    if (is.logical(x)) 
        return(isTRUE(x[[1]]))
    tolower(as.character(x[[1]])) %in% c("1", "true", "yes", 
        "y")
}
as_number <-
function (x, default = NA_real_) 
{
    out <- suppressWarnings(as.numeric(x))
    if (!length(out) || is.na(out[[1]])) 
        default
    else out[[1]]
}
assemble_r2_patient_data <-
function (expression, annotation_payload, dataset_info, track_info) 
{
    metadata <- annotation_payload$data
    if (!is.data.frame(metadata) || !"samplenames" %in% names(metadata)) 
        stop("R2 metadata schema changed", call. = FALSE)
    if (any(r2_missing(metadata$samplenames)) || anyDuplicated(metadata$samplenames)) 
        stop("R2 metadata IDs are missing or duplicated", call. = FALSE)
    if (!setequal(expression$sample_id, metadata$samplenames)) 
        stop("Expression/metadata sample ID sets differ; refusing a partial or positional join", 
            call. = FALSE)
    if (nrow(metadata) != as.integer(dataset_info$sampleSize)) 
        stop("Retrieved samples differ from R2 dataset's declared sample count", 
            call. = FALSE)
    aligned <- metadata[match(expression$sample_id, metadata$samplenames), 
        , drop = FALSE]
    audit <- r2_metadata_audit(metadata, dataset_info, track_info)
    derived <- expression
    get_raw <- function(field) if (field %in% names(aligned)) 
        as.character(aligned[[field]])
    else rep(NA_character_, nrow(aligned))
    derived$subgroup <- get_raw("subgroup")
    derived$subtype <- get_raw("subtype")
    derived$metastasis_raw <- get_raw(audit$metastasis$field %||% 
        "")
    derived$metastasis_label <- NA_character_
    if (isTRUE(audit$metastasis$coding_verified)) {
        derived$metastasis_label[which(derived$metastasis_raw == 
            "0")] <- audit$metastasis$code_0
        derived$metastasis_label[which(derived$metastasis_raw == 
            "1")] <- audit$metastasis$code_1
    }
    derived$survival_time_source_value <- get_raw("os_(years)")
    derived$survival_time_years <- r2_numeric(derived$survival_time_source_value, 
        "os_(years)")
    derived$survival_status_raw <- get_raw("dead")
    derived$survival_event <- NA_integer_
    list(data = derived, metadata = metadata, audit = audit)
}
assert_r2_page <-
function (html, expected, stage) 
{
    plain <- gsub("<[^>]+>", " ", html)
    if (!grepl(expected, plain, ignore.case = TRUE, perl = TRUE)) {
        stop("R2 workflow changed or failed at ", stage, ": expected ", 
            expected, call. = FALSE)
    }
    invisible(TRUE)
}
build_gene_page <-
function (cfg, gene, reporter = NULL, publication_limit = 10L, 
    publication_disease = NULL, aliases = character(), include_aliases = TRUE, 
    refresh_publications = FALSE, cutoff = "median", pfister_reporter = NULL) 
{
    cutoff <- match.arg(cutoff, c("median", "mean", "optimal"))
    retrieval <- get_r2_expression(cfg, gene, reporter)
    validation <- plot_r2_subgroup_test(retrieval, export_images = FALSE)
    source <- profile_attempt(get_cavalli_survival(cfg, retrieval))
    if (!identical(source$status, "unavailable")) {
        retrieval$audit$survival$event_coding_verified <- TRUE
        retrieval$audit$survival$note <- paste(source$provenance$event, 
            source$provenance$time)
    }
    subtypes <- profile_attempt(profile_group_comparison(retrieval$data, 
        "subtype", paste(toupper(gene), "Cavalli molecular subtypes")))
    pfister <- profile_attempt({
        pf <- get_pfister_expression(cfg, gene, pfister_reporter)
        x <- profile_group_comparison(pf$data[pf$data$include_cns, 
            , drop = FALSE], "cancer_type", paste(toupper(gene), 
            "across primary CNS tumours"), "log2(1 + FPKM)")
        x$excluded_samples <- pf$data[!pf$data$include_cns, c("sample_id", 
            "cancer_type", "sample_type", "exclusion_reason"), 
            drop = FALSE]
        x$all_counts <- as.data.frame(table(pf$data$cancer_type, 
            pf$data$sample_type))
        names(x$all_counts) <- c("cancer_type", "sample_type", 
            "n")
        x$provenance <- pf$provenance
        x$all_samples <- pf$data
        x
    })
    en <- profile_attempt(profile_enrichr(cfg, gene))
    dep <- profile_attempt(plot_public_depmap(get_public_depmap_gene(cfg, 
        gene)))
    clinical <- profile_clinical(retrieval, if (identical(source$status, 
        "unavailable")) 
        NULL
    else source)
    pub_error <- NULL
    pubs <- tryCatch(get_gene_publications(cfg, gene, publication_limit, 
        publication_disease, aliases, include_aliases, refresh_publications), 
        error = function(e) {
            pub_error <<- conditionMessage(e)
            NULL
        })
    modules <- list(cavalli_subtypes = subtypes, pfister_cns = pfister, 
        survival_source = source, enrichr = en, depmap = dep)
    status <- data.frame(module = names(modules), status = vapply(modules, 
        function(x) if (identical(x$status, "unavailable")) 
            "unavailable"
        else "ok", character(1)), detail = vapply(modules, function(x) x$reason %||% 
        "", character(1)))
    status <- rbind(status, data.frame(module = "pubmed", status = if (is.null(pubs)) 
        "unavailable"
    else "ok", detail = pub_error %||% ""))
    profile <- list(subtypes = subtypes, pfister = pfister, clinical = clinical, 
        enrichr = en, depmap = dep, status = status, survival_provenance = source$provenance %||% 
            source)
    root <- cfg$output_dir %||% file.path(dirname(cfg$cache_dir), 
        "..", "outputs")
    out <- file.path(root, toupper(gene), basename(tempfile(paste0("report-", 
        format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-"))))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    out <- normalizePath(out, winslash = "/")
    saveRDS(list(retrieval = retrieval, validation = validation, 
        publications = pubs, publication_error = pub_error, profile = profile, 
        cutoff = cutoff), file.path(out, "report_data.rds"))
    save_profile_tables(profile, file.path(out, "tables"))
    data.table::fwrite(retrieval$data, file.path(out, "Cavalli_patient_data.csv"), 
        na = "NA")
    writeLines(json_text(list(generated_at = utc_now(), cavalli = retrieval$provenance, 
        survival = source$provenance, pfister = pfister$provenance, 
        enrichr = en$provenance, depmap = dep$provenance, pubmed = pubs$provenance)), 
        file.path(out, "provenance.json"))
    if (!is.null(pubs)) 
        writeLines(json_text(pubs), file.path(out, "pubmed_results.json"))
    path <- write_gene_page(retrieval, validation, pubs, pub_error, 
        file.path(out, paste0(toupper(gene), "_gene_page.html")), 
        profile_extra_html(profile, cutoff))
    functions <- ls(envir = globalenv())
    functions <- functions[vapply(functions, function(n) is.function(get(n, 
        envir = globalenv())) && !startsWith(n, "."), logical(1))]
    dump(functions, file.path(out, "analysis_functions.R"), envir = globalenv())
    writeLines(c("source('analysis_functions.R')", "x <- readRDS('report_data.rds')", 
        "write_gene_page(x$retrieval, x$validation, x$publications, x$publication_error, 'regenerated_report.html', profile_extra_html(x$profile, x$cutoff))"), 
        file.path(out, "reproduce_report.R"))
    writeLines(capture.output(sessionInfo()), file.path(out, 
        "session-info.txt"))
    list(gene = toupper(gene), html_path = path, samples = nrow(retrieval$data), 
        publication_count = length(pubs$articles), publication_status = if (is.null(pubs)) "unavailable" else "ok", 
        modules = status, clinical_cutoff_default = cutoff, available_cutoffs = as.list(c("median", 
            "mean", "optimal")), report_data = file.path(out, 
            "report_data.rds"), instructions = "Open html_path in a browser; use cohort and cutoff selectors. Source data, tables, R functions and offline reproduction script are saved alongside the report.")
}
cache_get <-
function (cfg, key) 
{
    path <- cache_path(cfg, key)
    if (!file.exists(path)) 
        return(NULL)
    tryCatch(readRDS(path), error = function(e) NULL)
}
cache_inventory <-
function (cfg) 
{
    files <- list.files(cfg$cache_dir, pattern = "\\.rds$", full.names = TRUE)
    if (!length(files)) 
        return(data.frame())
    rows <- lapply(files, function(path) {
        obj <- tryCatch(readRDS(path), error = function(e) NULL)
        data.frame(file = basename(path), created_at = obj$provenance$retrieved_at %||% 
            format(file.info(path)$mtime, tz = "UTC"), dataset = obj$provenance$dataset_table %||% 
            NA_character_, contrast = obj$provenance$contrast %||% 
            NA_character_, rows = if (is.data.frame(obj$data)) 
            nrow(obj$data)
        else NA_integer_, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}
cache_key <-
function (prefix, parts) 
{
    raw <- paste(c(prefix, unlist(parts, use.names = TRUE)), 
        collapse = "__")
    safe <- gsub("[^A-Za-z0-9._-]+", "-", raw)
    if (nchar(safe) > 180L) {
        tmp <- tempfile()
        writeBin(charToRaw(raw), tmp)
        on.exit(unlink(tmp), add = TRUE)
        safe <- paste0(substr(safe, 1L, 120L), "-", unname(tools::md5sum(tmp)))
    }
    safe
}
cache_path <-
function (cfg, key) 
{
    file.path(cfg$cache_dir, paste0(key, ".rds"))
}
cache_put <-
function (cfg, key, value) 
{
    path <- cache_path(cfg, key)
    tmp <- paste0(path, ".tmp")
    saveRDS(value, tmp, version = 3)
    if (!file.rename(tmp, path)) {
        unlink(tmp)
        stop("Could not write cache file: ", path, call. = FALSE)
    }
    invisible(path)
}
clean_symbol <-
function (x) 
{
    trimws(sub("\\s*\\([^)]*\\)\\s*$", "", as.character(x)))
}
clean_table_names <-
function (x) 
{
    x <- trimws(gsub("[\r\n]+", " ", x))
    x <- gsub("[^A-Za-z0-9]+", "", x)
    tolower(x)
}
column_means <-
function (x) 
{
    out <- colMeans(x, na.rm = TRUE)
    out[!is.finite(out)] <- NA_real_
    out
}
column_medians <-
function (x) 
{
    out <- apply(x, 2L, stats::median, na.rm = TRUE)
    out[!is.finite(out)] <- NA_real_
    out
}
depmap_gene_columns <-
function (path) 
{
    columns <- depmap_header(path)
    id_col <- depmap_id_column(columns)
    gene_cols <- setdiff(columns, id_col)
    data.frame(column = gene_cols, gene = clean_symbol(gene_cols), 
        stringsAsFactors = FALSE)
}
depmap_gene_dependency <-
function (gene, model_csv = NULL, gene_effect_csv = NULL, model_regex = "(?i)medulloblastoma", 
    release = NULL, dependency_cutoff = -0.5) 
{
    model_csv <- normalize_local_path(model_csv, "DEPMAP_MODEL_CSV")
    gene_effect_csv <- normalize_local_path(gene_effect_csv, 
        "DEPMAP_GENE_EFFECT_CSV")
    release <- release %||% Sys.getenv("DEPMAP_RELEASE", "unspecified")
    gene <- toupper(clean_symbol(gene))
    if (!nzchar(gene)) 
        stop("gene must not be empty", call. = FALSE)
    models <- read_model_metadata(model_csv)
    filtered <- filter_depmap_models(models, model_regex)
    effect <- read_gene_effect_columns(gene_effect_csv, gene)
    hit <- effect$mapping[toupper(effect$mapping$gene) == gene, 
        , drop = FALSE]
    if (!nrow(hit)) 
        stop("Gene not found in DepMap gene-effect file: ", gene, 
            call. = FALSE)
    effect_col <- hit$column[[1]]
    joined <- merge(filtered$data, effect$data[, c(effect$id_col, 
        effect_col), drop = FALSE], by.x = filtered$id_col, by.y = effect$id_col, 
        all.x = TRUE, sort = FALSE)
    values <- suppressWarnings(as.numeric(joined[[effect_col]]))
    rows <- data.frame(model_id = joined[[filtered$id_col]], 
        model_name = preferred_model_label(joined), gene = gene, 
        gene_effect = values, stringsAsFactors = FALSE)
    rows <- rows[order(rows$gene_effect, na.last = TRUE), , drop = FALSE]
    list(summary = list(gene = gene, matched_models = nrow(rows), 
        models_with_data = sum(is.finite(rows$gene_effect)), 
        mean_gene_effect = mean(rows$gene_effect, na.rm = TRUE), 
        median_gene_effect = stats::median(rows$gene_effect, 
            na.rm = TRUE), minimum_gene_effect = if (any(is.finite(rows$gene_effect))) min(rows$gene_effect, 
            na.rm = TRUE) else NA_real_, fraction_below_cutoff = mean(rows$gene_effect < 
            dependency_cutoff, na.rm = TRUE), dependency_cutoff = dependency_cutoff), 
        models = rows, provenance = list(depmap_release = release, 
            model_filter_regex = model_regex, model_csv = model_csv, 
            gene_effect_csv = gene_effect_csv, effect_definition = "More-negative DepMap gene effect indicates stronger dependency"))
}
depmap_header <-
function (path) 
{
    names(data.table::fread(path, nrows = 0L, check.names = FALSE, 
        showProgress = FALSE))
}
depmap_id_column <-
function (columns) 
{
    candidates <- c("ModelID", "DepMap_ID", "DepMapID", "model_id", 
        "depmap_id")
    hit <- candidates[candidates %in% columns]
    if (!length(hit) && length(columns) && columns[1] %in% c("", 
        "V1")) 
        return(columns[1])
    if (!length(hit)) 
        stop("Could not find a DepMap model identifier column", 
            call. = FALSE)
    hit[[1]]
}
depmap_public_files <-
function (cfg) 
{
    folder <- file.path(cfg$cache_dir, "depmap-24Q4-v1")
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    manifest_path <- file.path(folder, "manifest.json")
    if (!file.exists(manifest_path)) {
        txt <- enrichr_get("https://api.figshare.com/v2/articles/27993248/versions/1")
        manifest <- jsonlite::fromJSON(txt)
        if (!identical(manifest$title, "DepMap 24Q4 Public") || 
            manifest$version != 1L) 
            stop("DepMap archive identity mismatch")
        writeBin(charToRaw(txt), manifest_path)
    }
    manifest <- jsonlite::fromJSON(manifest_path)
    wanted <- c("Model.csv", "CRISPRGeneEffect.csv")
    paths <- setNames(character(length(wanted)), wanted)
    for (name in wanted) {
        f <- manifest$files[manifest$files$name == name, , drop = FALSE]
        if (nrow(f) != 1) 
            stop("Ambiguous DepMap archive file")
        path <- file.path(folder, name)
        if (!file.exists(path)) {
            message("Downloading pinned DepMap ", name, " (", 
                round(f$size/1000000), " MB), once per installation")
            part <- paste0(path, ".part")
            response <- curl::curl_fetch_disk(f$download_url, 
                part, curl::new_handle(timeout = 900, connecttimeout = 30, 
                  followlocation = TRUE))
            if (response$status_code != 200 || unname(file.info(part)$size) != 
                f$size || unname(tools::md5sum(part)) != f$computed_md5) 
                stop("DepMap download failed integrity validation: ", 
                  name)
            if (!file.rename(part, path)) 
                stop("Could not finalize DepMap download")
        }
        if (unname(file.info(path)$size) != f$size || unname(tools::md5sum(path)) != 
            f$computed_md5) 
            stop("Cached DepMap file integrity mismatch")
        paths[name] <- normalizePath(path, winslash = "/")
    }
    list(paths = paths, provenance = list(release = manifest$title, 
        doi = manifest$doi, manifest = normalizePath(manifest_path, 
            winslash = "/"), note = "Pinned official archived 24Q4 v1 release, not the latest release. Portal automated access currently requires human verification.", 
        files = manifest$files[manifest$files$name %in% wanted, 
            c("name", "download_url", "computed_md5", "size")]))
}
encode_form <-
function (fields) 
{
    encode_one <- function(x) utils::URLencode(as.character(x %||% 
        ""), reserved = TRUE)
    paste0(vapply(names(fields), encode_one, character(1)), "=", 
        vapply(fields, encode_one, character(1)), collapse = "&")
}
enrich_pathways <-
function (cfg, analysis_id, libraries = character(), gmt_path = NULL, 
    min_set_size = 5L, max_set_size = 2000L, fdr = 0.050000000000000003, 
    limit = 50L, refresh = FALSE) 
{
    a <- load_study_object(cfg, analysis_id, "analysis")
    libraries <- unlist(libraries, use.names = FALSE)
    if (!length(libraries) && is.null(gmt_path)) 
        stop("Specify libraries or gmt_path", call. = FALSE)
    if (!length(a$universe)) 
        stop("Analysis has no tested-gene background", call. = FALSE)
    inputs <- if (!is.null(gmt_path)) 
        list(load_gene_sets(cfg, gmt_path = gmt_path))
    else lapply(libraries, function(l) load_gene_sets(cfg, l, 
        refresh = refresh))
    labels <- if (!is.null(gmt_path)) 
        basename(gmt_path)
    else libraries
    results <- list()
    for (i in seq_along(inputs)) for (direction in c("up", "down")) {
        genes <- a$summary$gene[a$summary$passes_all_comparators & 
            a$summary$direction == direction]
        table <- overrepresentation(genes, a$universe, inputs[[i]]$sets, 
            min_set_size, max_set_size)
        results[[paste(labels[i], direction, sep = ":")]] <- list(direction = direction, 
            library = labels[i], input_genes_n = length(genes), 
            tested_terms_n = nrow(table), significant_terms_n = sum(table$adjusted_p <= 
                fdr), rows = utils::head(table[table$adjusted_p <= 
                fdr, , drop = FALSE], limit), provenance = inputs[[i]]$provenance)
    }
    list(analysis_id = analysis_id, results = results, background_n = length(a$universe), 
        method = "One-sided hypergeometric over-representation; BH across ALL size-eligible terms separately per library and direction, including zero overlaps", 
        interpretation = "Pathways enriched among consistently up/down deregulated genes; enrichment does not establish pathway activation or inhibition.", 
        privacy = "Only library names are requested from Enrichr. Gene lists, expression and clinical data remain local.", 
        analysis_provenance = a$provenance)
}
enrichr_get <-
function (url) 
{
    h <- curl::new_handle(timeout = 90, connecttimeout = 20, 
        followlocation = TRUE, useragent = "r2-depmap-mcp/0.2.0")
    response <- curl::curl_fetch_memory(url, handle = h)
    if (response$status_code != 200) 
        stop("Enrichr returned HTTP ", response$status_code, 
            call. = FALSE)
    rawToChar(response$content)
}
file_provenance <-
function (path) 
{
    path <- normalize_local_path(path)
    list(path = path, md5 = unname(tools::md5sum(path)), bytes = unname(file.info(path)$size))
}
filter_depmap_models <-
function (models, regex) 
{
    id_col <- depmap_id_column(names(models))
    keep <- grepl(regex, model_filter_text(models), perl = TRUE)
    selected <- models[keep & !is.na(models[[id_col]]), , drop = FALSE]
    if (!nrow(selected)) 
        stop("No DepMap models matched regex: ", regex, call. = FALSE)
    list(data = selected, id_col = id_col)
}
get_cavalli_survival <-
function (cfg, retrieval, refresh = FALSE) 
{
    key <- cache_key("cavalli_survival_v2", list(dataset = cfg$r2_dataset_table, 
        reporter = retrieval$provenance$reporter))
    cached <- if (!refresh) 
        cache_get(cfg, key)
    else NULL
    if (!is.null(cached)) 
        return(cached)
    h <- r2_session(120)
    landing <- r2_request(h, cfg$r2_base_url)
    if (grepl("Use R2 without an account", landing, fixed = TRUE)) 
        r2_request(h, cfg$r2_base_url, list(method = "guest", 
            open_page = "auth"))
    html <- r2_request(h, cfg$r2_base_url, list(table = cfg$r2_dataset_table, 
        option = "kaplanscan", inputprobeset = retrieval$provenance$reporter, 
        scanmodus = "median", survival = "overall", mingrpsize = "8", 
        subset = ""))
    doc <- xml2::read_html(charToRaw(html), encoding = "UTF-8")
    scripts <- xml2::xml_text(xml2::xml_find_all(doc, ".//script"))
    line <- scripts[grepl("/d3/plots/kaplan.js", scripts, fixed = TRUE)]
    if (length(line) != 1) 
        stop("R2 Kaplan source module unavailable")
    url <- sub(".*import\\('([^']+)'.*", "\\1", line)
    if (!grepl("^https://hgserver1[.]amc[.]nl/r2/assets/[^/]+/js/src/d3/plots/kaplan[.]js$", 
        url)) 
        stop("Unexpected Kaplan renderer URL")
    js <- r2_request(h, url)
    result <- verify_r2_survival_payload(html, js, retrieval)
    folder <- file.path(retrieval$provenance$output_dir, "survival_source")
    dir.create(folder, showWarnings = FALSE)
    writeBin(charToRaw(html), file.path(folder, "r2-kaplan.html"))
    writeBin(charToRaw(js), file.path(folder, "kaplan.js"))
    data.table::fwrite(result$data, file.path(folder, "verified_survival.csv"))
    result$provenance$source_dir <- folder
    result$provenance$retrieved_at <- utc_now()
    result$provenance$renderer_url <- url
    cache_put(cfg, key, result)
    result
}
get_gene_publications <-
function (cfg, gene, limit = 10L, disease = NULL, aliases = character(), 
    include_aliases = TRUE, refresh = FALSE) 
{
    terms <- pubmed_gene_terms(gene, aliases)
    gene <- terms[1]
    if (length(limit) != 1L || !is.numeric(limit) || !is.finite(limit) || 
        limit != floor(limit) || limit < 1L || limit > 30L) 
        stop("limit must be 1-30")
    pubmed_query(terms, disease)
    key <- cache_key("pubmed_v2", list(gene = gene, disease = disease %||% 
        "all", aliases = paste(terms, collapse = "|"), expand = include_aliases, 
        limit = limit))
    cached <- if (!refresh) 
        cache_get(cfg, key)
    else NULL
    if (!is.null(cached)) {
        age <- as.numeric(difftime(Sys.time(), as.POSIXct(cached$provenance$retrieved_at, 
            format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), units = "hours"))
        if (is.finite(age) && age >= 0 && age < 24) {
            cached$provenance$cache_hit <- TRUE
            return(cached)
        }
    }
    resolved <- if (include_aliases) 
        tryCatch(resolve_pubmed_aliases(gene), error = function(e) list(aliases = character(), 
            status = paste("Alias lookup unavailable:", conditionMessage(e))))
    else list(aliases = character(), status = "Automatic aliases disabled")
    terms <- pubmed_gene_terms(gene, unique(c(aliases, resolved$aliases)))
    query <- pubmed_query(terms, disease)
    screening_limit <- min(200L, max(50L, limit * 5L))
    raw_search <- ncbi_get("esearch.fcgi", list(db = "pubmed", 
        term = query, retmode = "json", retmax = screening_limit, 
        sort = "pub_date"))
    search <- r2_json(raw_search, "PubMed search")
    if (!is.null(search$error) || !is.null(search$esearchresult$errorlist)) 
        stop("PubMed rejected the query: ", json_text(search$error %||% 
            search$esearchresult$errorlist))
    if (is.null(search$esearchresult$count) || is.null(search$esearchresult$idlist)) 
        stop("Unexpected PubMed search result")
    ids <- as.character(search$esearchresult$idlist)
    raw_xml <- NULL
    articles <- list()
    if (length(ids)) {
        raw_xml <- ncbi_get("efetch.fcgi", list(db = "pubmed", 
            id = paste(ids, collapse = ","), retmode = "xml"))
        articles <- parse_pubmed_articles(raw_xml, terms)
        returned <- vapply(articles, `[[`, character(1), "pmid")
        if (anyDuplicated(returned) || !setequal(ids, returned)) 
            stop("PubMed search/fetch record mismatch; refusing partial results")
        articles <- articles[match(ids, returned)]
    }
    relevant <- vapply(articles, function(a) pubmed_mentions(paste(a$title, 
        a$abstract), terms), logical(1))
    excluded_ids <- vapply(articles[!relevant], `[[`, character(1), 
        "pmid")
    articles <- head(articles[relevant], limit)
    raw_dir <- file.path(cfg$cache_dir, "pubmed_sources", basename(tempfile("retrieval-")))
    dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
    writeBin(charToRaw(raw_search), file.path(raw_dir, "search.json"))
    if (!is.null(raw_xml)) 
        writeBin(charToRaw(raw_xml), file.path(raw_dir, "articles.xml"))
    if (!is.null(resolved$raw_search)) 
        writeBin(charToRaw(resolved$raw_search), file.path(raw_dir, 
            "gene-search.json"))
    if (!is.null(resolved$raw_summary)) 
        writeBin(charToRaw(resolved$raw_summary), file.path(raw_dir, 
            "gene-summary.json"))
    resolved$raw_search <- resolved$raw_summary <- NULL
    result <- list(gene = gene, disease_filter = disease, total_matches = as.integer(search$esearchresult$count), 
        articles = articles, provenance = list(source = "NCBI PubMed E-utilities", 
            query = query, query_translation = search$esearchresult$querytranslation, 
            searched_terms = as.list(terms), alias_resolution = resolved, 
            sort = "PubMed publication date descending (pub_date), not indexing date", 
            date_definition = "Journal citation date displayed verbatim, with electronic publication date separately when provided; partial dates are not invented", 
            summary_method = "One verbatim abstract sentence, preferring gene mentions and conclusions/results; title-only fallback when no abstract is supplied", 
            screening = list(records_checked = length(ids), excluded_pmids = as.list(excluded_ids), 
                rule = "Require a standalone gene symbol or alias in title/abstract; hyphenated drug codes such as HLX-02 do not qualify", 
                maximum_records = screening_limit), limitations = "Keyword/alias matches can include unrelated uses or non-human studies; no implication that a paper studies medulloblastoma unless separately checked", 
            search_warnings = search$esearchresult$warninglist, 
            retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", 
                tz = "UTC"), cache_hit = FALSE, raw_source_dir = normalizePath(raw_dir, 
                winslash = "/")))
    cache_put(cfg, key, result)
    result
}
get_pfister_expression <-
function (cfg, gene, reporter = NULL, refresh = FALSE) 
{
    pcfg <- cfg
    pcfg$r2_dataset_table <- "ps_avgpres_pfisterb272_informp3"
    pcfg$r2_transformation <- "transform_none"
    pcfg$r2_grouping_track <- "cancer_type"
    pcfg$output_dir <- file.path(cfg$output_dir %||% file.path(dirname(cfg$cache_dir), 
        "..", "outputs"), "pfister")
    r <- get_r2_expression(pcfg, gene, reporter, refresh)
    if (nrow(r$data) != 272L || !grepl("Pfister", r$provenance$dataset_label, 
        fixed = TRUE)) 
        stop("Pfister dataset identity/count mismatch")
    d <- r$data
    a <- r$metadata[match(d$sample_id, r$metadata$samplenames), 
        ]
    require_columns(a, c("cancer_type", "type"))
    if (any(d$expression < 0, na.rm = TRUE)) 
        stop("Negative source FPKM values")
    d$fpkm <- d$expression
    d$expression <- log2(1 + d$fpkm)
    d$cancer_type <- a$cancer_type
    d$sample_type <- a$type
    cns <- c("atrt", "epd_it", "etmr", "hgg_k27m", "hggother", 
        "mb_group3", "mb_group4", "mb_shh", "mb_wnt", "pa")
    known <- c(cns, "ews", "nb", "os", "rms", "t-all")
    if (any(!a$cancer_type %in% known)) 
        stop("Unrecognised Pfister cancer label; update documented CNS mapping before filtering")
    d$include_cns <- a$cancer_type %in% cns & a$type == "primary"
    d$exclusion_reason <- ifelse(!a$cancer_type %in% cns, "Non-CNS cancer label", 
        ifelse(a$type != "primary", "Non-primary sample; avoid paired relapse inference", 
            ""))
    r$data <- d
    r$provenance$analysis_transform <- "log2(1 + FPKM), applied locally to R2 transform_none values"
    r$provenance$cns_labels <- as.list(cns)
    r$provenance$selection <- "Explicit cancer_type CNS label whitelist and type=primary; full 272 samples retained in source files"
    r
}
get_public_depmap_gene <-
function (cfg, gene) 
{
    key <- cache_key("public_depmap_24q4_v1", list(gene = toupper(gene)))
    cached <- cache_get(cfg, key)
    if (!is.null(cached)) 
        return(cached)
    files <- depmap_public_files(cfg)
    models <- read_model_metadata(files$paths[["Model.csv"]])
    require_columns(models, c("ModelID", "CellLineName", "OncotreePrimaryDisease"))
    effects <- read_gene_effect_columns(files$paths[["CRISPRGeneEffect.csv"]], 
        gene)
    effect_col <- effects$mapping$column[1]
    ids <- as.character(effects$data[[effects$id_col]])
    if (any(!ids %in% models$ModelID)) 
        stop("DepMap effect IDs absent from same-release metadata")
    values <- as.numeric(effects$data[[effect_col]])
    idx <- match(models$ModelID, ids)
    data <- models[, intersect(c("ModelID", "CellLineName", "OncotreeLineage", 
        "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode"), 
        names(models)), drop = FALSE]
    data$gene_effect <- values[idx]
    disease_fields <- intersect(c("OncotreePrimaryDisease", "OncotreeSubtype"), 
        names(data))
    data$is_medulloblastoma <- Reduce(`|`, lapply(disease_fields, 
        function(f) !is.na(data[[f]]) & grepl("^medulloblastoma(,|$)", 
            data[[f]], ignore.case = TRUE)))
    if (!any(data$is_medulloblastoma)) 
        stop("No explicit medulloblastoma disease annotation found")
    data$gene <- toupper(gene)
    result <- list(data = data, provenance = c(files$provenance, 
        list(gene_column = effect_col, model_selection = "Explicit OncotreePrimaryDisease or OncotreeSubtype label Medulloblastoma (including comma-qualified subtypes); no subgroup inferred from model names", 
            retrieved_at = utc_now(), interpretation = "Chronos gene effect: more negative indicates greater in-vitro dependency; 0 and -1 are reference values, not clinical thresholds")))
    cache_put(cfg, key, result)
    result
}
get_r2_expression <-
function (cfg, gene, reporter = NULL, refresh = FALSE) 
{
    if (length(gene) != 1L || !grepl("^[A-Za-z][A-Za-z0-9._-]*$", 
        gene)) 
        stop("Provide one gene symbol", call. = FALSE)
    gene <- toupper(gene)
    transformation <- cfg$r2_transformation %||% "transform_log2"
    key <- cache_key(if (transformation == "transform_log2") 
        "r2_patient_v1"
    else "r2_patient_v2", list(dataset = cfg$r2_dataset_table, 
        gene = gene, reporter = reporter %||% "unique", transformation = if (transformation == 
            "transform_log2") "log2" else transformation))
    if (!refresh) {
        cached <- cache_get(cfg, key)
        if (!is.null(cached)) {
            cached$provenance$cache_hit <- TRUE
            return(cached)
        }
    }
    root <- cfg$output_dir %||% file.path(dirname(cfg$cache_dir), 
        "..", "outputs")
    run_id <- basename(tempfile(paste0(format(Sys.time(), "%Y%m%dT%H%M%S", 
        tz = "UTC"), "-")))
    out <- file.path(root, gene, run_id)
    raw_dir <- file.path(out, "source")
    dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(raw_dir)) 
        stop("Cannot create output directory", call. = FALSE)
    out <- normalizePath(out, winslash = "/", mustWork = TRUE)
    h <- r2_session(120)
    requests <- list()
    capture <- function(name, params = NULL, post = FALSE) {
        url <- cfg$r2_base_url
        if (!is.null(params) && !post) 
            url <- paste0(url, "?", encode_form(params))
        text <- r2_request(h, url, if (post) 
            params
        else NULL)
        path <- file.path(raw_dir, name)
        writeBin(charToRaw(text), path)
        requests[[length(requests) + 1L]] <<- list(file = name, 
            url = url, method = if (post) "POST" else "GET", 
            fields = if (post) params else NULL, retrieved_at = format(Sys.time(), 
                "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), md5 = unname(tools::md5sum(path)))
        text
    }
    landing <- capture("landing.html", list(table = cfg$r2_dataset_table))
    if (grepl("Use R2 without an account", landing, fixed = TRUE)) 
        capture("guest.html", list(method = "guest", open_page = "auth"), 
            TRUE)
    capture("dataset.html", list(table = cfg$r2_dataset_table))
    info <- r2_json(capture("dataset.json", list(json_option = "json_dataset_info", 
        table = cfg$r2_dataset_table)), "dataset")
    if (!identical(info$internalIdentifier, cfg$r2_dataset_table)) 
        stop("Dataset identity mismatch", call. = FALSE)
    tracks <- r2_json(capture("tracks.json", list(json_option = "json_dataset_info_tracks", 
        table = cfg$r2_dataset_table)), "tracks")
    annotations <- r2_json(capture("annotations.json", list(json_option = "json_cg_sampleannotation_v1", 
        ctable = cfg$r2_dataset_table, subset = "")), "annotations")
    reporters <- r2_json(capture("reporters.json", list(json_option = "json_find_reporter_in_dataset_v1", 
        dataset = cfg$r2_dataset_table, search_by = "gene_symbol", 
        query = gene)), "reporters")
    if (!is.data.frame(reporters) || !all(c("gene_symbol", "reporter") %in% 
        names(reporters))) 
        stop("Gene absent from R2 reporter search", call. = FALSE)
    reporters <- reporters[toupper(reporters$gene_symbol) == 
        gene, , drop = FALSE]
    if (!is.null(reporter)) 
        reporters <- reporters[as.character(reporters$reporter) == 
            reporter, , drop = FALSE]
    if (nrow(reporters) != 1L) 
        stop("Expected one exact reporter for ", gene, "; found ", 
            nrow(reporters), ". Specify reporter explicitly if ambiguous.", 
            call. = FALSE)
    reporter <- as.character(reporters$reporter[1])
    html <- capture("expression.html", list(table = cfg$r2_dataset_table, 
        option = "display2", analysis_type = "single_reporter", 
        factor = reporter, cortype = transformation, subset = "", 
        grouping_track = cfg$r2_grouping_track %||% "subgroup", 
        graphtype = "yy_i"), TRUE)
    expression <- parse_r2_patient_expression(html, gene, reporter, 
        cfg$r2_dataset_table, transformation)
    result <- assemble_r2_patient_data(expression, annotations, 
        info, tracks)
    result$provenance <- list(dataset_table = cfg$r2_dataset_table, 
        accession = info$accession, dataset_label = info$datasetTitle, 
        queried_gene = gene, reporter = reporter, transformation = if (transformation == 
            "transform_log2") "R2 transform_log2 (log2 of source expression)" else "R2 transform_none (untransformed source expression)", 
        sample_join = "Exact equality: expression table samplenames to annotation samplenames; no positional or GEO-expression-ID join", 
        cache_hit = FALSE, retrieved_at = format(Sys.time(), 
            "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), output_dir = out, 
        requests = requests, limitations = "Unofficial read-only R2 interface; validated against observed responses. No clinical/subtype inference. All samples retained.")
    result$files <- list(patient_data = file.path(out, paste0(gene, 
        "_R2_patient_data.csv")), metadata = file.path(out, paste0(gene, 
        "_R2_metadata.csv")), metadata_audit = file.path(out, 
        "metadata_audit.json"), provenance = file.path(out, "provenance.json"), 
        r_object = file.path(out, "retrieval.rds"))
    data.table::fwrite(result$data, result$files$patient_data, 
        na = "NA")
    data.table::fwrite(result$metadata, result$files$metadata, 
        na = "NA")
    writeLines(json_text(result$audit), result$files$metadata_audit)
    writeLines(json_text(result$provenance), result$files$provenance)
    saveRDS(result, result$files$r_object)
    writeLines(capture.output(sessionInfo()), file.path(out, 
        "session-info.txt"))
    cache_put(cfg, key, result)
    result
}
get_tf_annotation <-
function (cfg, gene) 
{
    key <- "human_tfs_lambert_v101"
    x <- cache_get(cfg, key)
    if (is.null(x)) {
        url <- "https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
        raw <- enrichr_get(url)
        genes <- trimws(strsplit(raw, "\n", fixed = TRUE)[[1]])
        if (length(genes) < 1000 || any(grepl("<html", genes, 
            ignore.case = TRUE))) 
            stop("Unexpected Human TFs response")
        path <- file.path(cfg$cache_dir, "TF_names_v_1.01.txt")
        writeBin(charToRaw(raw), path)
        x <- list(genes = genes, source = url, version = "Lambert human TF list v1.01", 
            retrieved_at = utc_now(), md5 = unname(tools::md5sum(path)), 
            raw_file = path)
        cache_put(cfg, key, x)
    }
    list(gene = gene, classification = if (toupper(gene) %in% 
        toupper(x$genes)) "Listed human transcription factor" else "Not listed in this version; absence is not proof of non-TF function", 
        source = x$source, version = x$version, raw_file = x$raw_file, 
        md5 = x$md5)
}
html_escape <-
function (x) 
{
    x <- as.character(x)
    x[is.na(x)] <- "Not reported"
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    gsub("'", "&#39;", x, fixed = TRUE)
}
html_table <-
function (d, id = NULL) 
{
    if (is.null(d) || !nrow(d)) 
        return("<p>No rows available.</p>")
    display <- d
    for (key in names(display)) if (is.numeric(display[[key]])) 
        display[[key]] <- format(signif(display[[key]], 6), trim = TRUE)
    header <- paste0("<th scope='col'>", html_escape(names(display)), 
        "</th>", collapse = "")
    rows <- apply(display, 1, function(row) paste0("<tr>", paste0("<td>", 
        html_escape(row), "</td>", collapse = ""), "</tr>"))
    paste0("<div class='table-wrap'><table", if (!is.null(id)) 
        paste0(" id='", html_escape(id), "'"), "><thead><tr>", 
        header, "</tr></thead><tbody>", paste(rows, collapse = ""), 
        "</tbody></table></div>")
}
import_de_results <-
function (cfg, path, disease, target, label = "Imported differential expression", 
    gene_column = "gene", effect_column = "log2fc", adjusted_p_column = "adjusted_p", 
    comparator_column = "comparator", comparator = NULL, background_path = NULL, 
    fdr = 0.050000000000000003, min_abs_log2fc = 1) 
{
    d <- read_study_table(path)
    require_columns(d, c(gene_column, effect_column, adjusted_p_column))
    if (is.null(comparator)) 
        require_columns(d, comparator_column)
    comparisons <- if (is.null(comparator)) 
        as.character(d[[comparator_column]])
    else rep(comparator, nrow(d))
    if (!nrow(d) || anyNA(comparisons) || any(!nzchar(comparisons)) || 
        any(comparisons == target)) 
        stop("Provide non-empty comparator labels different from target", 
            call. = FALSE)
    data <- data.frame(gene = study_symbols(d[[gene_column]]), 
        log2fc = numeric_column(d[[effect_column]], effect_column, 
            TRUE), adjusted_p = numeric_column(d[[adjusted_p_column]], 
            adjusted_p_column, TRUE), comparator = comparisons, 
        stringsAsFactors = FALSE)
    if (any(data$adjusted_p < 0 | data$adjusted_p > 1, na.rm = TRUE)) 
        stop("Adjusted P values must be between 0 and 1", call. = FALSE)
    if (anyDuplicated(data[, c("gene", "comparator")])) 
        stop("Duplicate gene/comparator rows; resolve mappings before import", 
            call. = FALSE)
    comparators <- sort(unique(comparisons))
    universe <- Reduce(intersect, lapply(comparators, function(g) data$gene[data$comparator == 
        g & is.finite(data$adjusted_p) & is.finite(data$log2fc)]))
    if (!is.null(background_path)) {
        bg <- read_study_table(background_path)
        require_columns(bg, "gene")
        universe <- unique(study_symbols(bg$gene))
        if (!all(data$gene %in% universe)) 
            stop("Background must contain every imported gene", 
                call. = FALSE)
    }
    result <- list(label = label, disease = disease, target = target, 
        comparators = comparators, data = data, universe = universe, 
        summary = summarize_specificity(data, comparators, fdr, 
            min_abs_log2fc), thresholds = list(fdr = fdr, min_abs_log2fc = min_abs_log2fc), 
        provenance = list(source = "user differential-expression results", 
            file = file_provenance(path), background_file = if (!is.null(background_path)) file_provenance(background_path) else NULL, 
            universe_definition = if (is.null(background_path)) "Intersection of finite tests; input MUST include all tested genes, not only significant hits" else "Explicit user-supplied tested-gene universe", 
            effect_definition = "log2 fold change: target minus comparator; retained without recalculation", 
            imported_at = utc_now(), limitations = "Specificity is relative to supplied comparators; a pooled rest comparison cannot establish subgroup specificity. Input adjusted P values are retained; no joint FDR claim is made."))
    result$id <- save_study_object(cfg, "analysis", result)
    study_overview(result)
}
json_text <-
function (x, pretty = TRUE) 
{
    jsonlite::toJSON(x, auto_unbox = TRUE, dataframe = "rows", 
        null = "null", na = "null", digits = 8, pretty = pretty)
}
list_enrichr_libraries <-
function (cfg, pattern = "", refresh = FALSE) 
{
    key <- "enrichr_library_catalog"
    cached <- if (!refresh) 
        cache_get(cfg, key)
    else NULL
    if (is.null(cached)) {
        url <- "https://maayanlab.cloud/Enrichr/datasetStatistics"
        payload <- jsonlite::fromJSON(enrichr_get(url))
        if (!is.data.frame(payload$statistics) || !"libraryName" %in% 
            names(payload$statistics)) 
            stop("Unexpected Enrichr library catalogue", call. = FALSE)
        cached <- list(data = payload$statistics, provenance = list(source = url, 
            retrieved_at = utc_now()))
        cache_put(cfg, key, cached)
    }
    list(libraries = cached$data[grepl(pattern, cached$data$libraryName, 
        ignore.case = TRUE), , drop = FALSE], provenance = cached$provenance)
}
list_studies <-
function (cfg) 
{
    paths <- list.files(cfg$cache_dir, "^(cohort|analysis)_[a-f0-9]{32}\\.rds$", 
        full.names = TRUE)
    list(studies = lapply(paths, function(p) study_overview(readRDS(p))))
}
load_gene_sets <-
function (cfg, library = NULL, gmt_path = NULL, refresh = FALSE) 
{
    if (!is.null(gmt_path)) {
        raw <- paste(readLines(normalize_local_path(gmt_path), 
            warn = FALSE), collapse = "\n")
        return(list(sets = parse_gmt(raw), provenance = list(source = "local GMT", 
            file = file_provenance(gmt_path))))
    }
    if (is.null(library) || !nzchar(library)) 
        stop("Specify an Enrichr library or local GMT path", 
            call. = FALSE)
    key <- paste0("enrichr_gmt_", object_digest(library))
    cached <- if (!refresh) 
        cache_get(cfg, key)
    else NULL
    if (!is.null(cached)) 
        return(cached)
    catalog <- list_enrichr_libraries(cfg, refresh = refresh)
    if (!library %in% catalog$libraries$libraryName) 
        stop("Unknown Enrichr library; use list_enrichr_libraries with refresh=true", 
            call. = FALSE)
    url <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", 
        utils::URLencode(library, reserved = TRUE))
    raw <- enrichr_get(url)
    result <- list(sets = parse_gmt(raw), provenance = list(source = "Enrichr", 
        library = library, url = url, content_md5 = object_digest(raw), 
        retrieved_at = utc_now(), computation = "Local over-representation using downloaded annotations; not Enrichr combined scores"))
    cache_put(cfg, key, result)
    result
}
load_study_object <-
function (cfg, id, kind) 
{
    if (length(id) != 1L || !grepl(paste0("^", kind, "_[a-f0-9]{32}$"), 
        id)) 
        stop("Invalid ", kind, " identifier", call. = FALSE)
    value <- cache_get(cfg, id)
    if (is.null(value)) 
        stop("Unknown ", kind, " identifier: ", id, call. = FALSE)
    value
}
log_stderr <-
function (...) 
{
    cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), 
        ..., "\n", file = stderr())
}
model_filter_text <-
function (models) 
{
    candidates <- c("ModelID", "CellLineName", "CCLEName", "StrippedCellLineName", 
        "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", 
        "ModelType", "Tissue", "lineage", "disease", "primary_disease")
    cols <- intersect(candidates, names(models))
    if (!length(cols)) 
        cols <- names(models)[vapply(models, is.character, logical(1))]
    apply(models[, cols, drop = FALSE], 1L, function(row) paste(row, 
        collapse = " | "))
}
ncbi_get <-
function (endpoint, params) 
{
    if (!endpoint %in% c("esearch.fcgi", "esummary.fcgi", "efetch.fcgi")) 
        stop("Unsupported NCBI endpoint")
    params$tool <- "r2_mb_gene_profile"
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
        endpoint, "?", encode_form(params))
    for (attempt in 1:3) {
        delay <- 0.35999999999999999 - as.numeric(difftime(Sys.time(), 
            last_request, units = "secs"))
        if (delay > 0) 
            Sys.sleep(delay)
        last_request <<- Sys.time()
        response <- curl::curl_fetch_memory(url, curl::new_handle(timeout = 60, 
            connecttimeout = 15, followlocation = TRUE, useragent = "r2-mb-gene-profile/0.3 (+NCBI E-utilities)"))
        if (response$status_code == 200L) 
            return(rawToChar(response$content))
        if (!response$status_code %in% c(429L, 500L, 502L, 503L, 
            504L) || attempt == 3L) 
            stop("NCBI returned HTTP ", response$status_code, 
                call. = FALSE)
        Sys.sleep(attempt)
    }
}
need_package <-
function (package) 
{
    if (!requireNamespace(package, quietly = TRUE)) 
        stop("Install ", package, " using scripts/install-dependencies.R --analysis", 
            call. = FALSE)
}
normalize_local_path <-
function (path, env_name = NULL) 
{
    value <- path %||% if (!is.null(env_name)) 
        Sys.getenv(env_name, "")
    else ""
    if (!nzchar(value)) {
        stop("A file path is required", if (!is.null(env_name)) 
            paste0(" (or set ", env_name, ")")
        else "", call. = FALSE)
    }
    value <- path.expand(value)
    if (!file.exists(value)) 
        stop("File not found: ", value, call. = FALSE)
    normalizePath(value, winslash = "/", mustWork = TRUE)
}
numeric_column <-
function (x, label, allow_na = FALSE) 
{
    y <- suppressWarnings(as.numeric(x))
    if (any(!is.na(x) & !is.finite(y)) || (!allow_na && anyNA(y))) 
        stop(label, " must contain finite numeric values", call. = FALSE)
    y
}
object_digest <-
function (value) 
{
    path <- tempfile()
    on.exit(unlink(path))
    saveRDS(value, path, version = 3)
    unname(tools::md5sum(path))
}
overrepresentation <-
function (genes, universe, sets, min_set_size = 5L, max_set_size = 2000L) 
{
    universe <- unique(universe)
    genes <- intersect(unique(genes), universe)
    empty <- data.frame(term = character(), overlap_n = integer(), 
        set_n = integer(), query_n = integer(), universe_n = integer(), 
        fold_enrichment = numeric(), p_value = numeric(), adjusted_p = numeric(), 
        overlap_genes = character())
    if (!length(genes) || !length(universe)) 
        return(empty)
    sets <- lapply(sets, intersect, universe)
    sets <- sets[lengths(sets) >= min_set_size & lengths(sets) <= 
        max_set_size]
    if (!length(sets)) 
        return(empty)
    rows <- lapply(names(sets), function(term) {
        members <- sets[[term]]
        overlap <- intersect(genes, members)
        data.frame(term = term, overlap_n = length(overlap), 
            set_n = length(members), query_n = length(genes), 
            universe_n = length(universe), fold_enrichment = (length(overlap)/length(genes))/(length(members)/length(universe)), 
            p_value = stats::phyper(length(overlap) - 1, length(members), 
                length(universe) - length(members), length(genes), 
                lower.tail = FALSE), overlap_genes = paste(overlap, 
                collapse = ";"), stringsAsFactors = FALSE)
    })
    result <- do.call(rbind, rows)
    result$adjusted_p <- stats::p.adjust(result$p_value, "BH")
    result[order(result$adjusted_p, result$p_value, -result$fold_enrichment), 
        ]
}
parse_gmt <-
function (text) 
{
    lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
    lines <- lines[nzchar(trimws(lines))]
    fields <- strsplit(lines, "\t", fixed = TRUE)
    if (!length(fields) || any(lengths(fields) < 3)) 
        stop("Invalid GMT: expected term, description and gene symbols", 
            call. = FALSE)
    result <- list()
    for (f in fields) {
        term <- f[1]
        genes <- toupper(trimws(sub(",.*$", "", f[-c(1, 2)])))
        genes <- genes[nzchar(genes)]
        result[[term]] <- unique(c(result[[term]], genes))
    }
    result
}
parse_pubmed_articles <-
function (xml, terms) 
{
    doc <- xml2::read_xml(charToRaw(xml), options = "NONET")
    errors <- xml2::xml_find_all(doc, "//ERROR")
    if (length(errors)) 
        stop("NCBI EFetch error: ", paste(xml2::xml_text(errors), 
            collapse = "; "))
    if (!identical(xml2::xml_name(doc), "PubmedArticleSet")) 
        stop("Unexpected PubMed XML root", call. = FALSE)
    articles <- xml2::xml_find_all(doc, "/PubmedArticleSet/PubmedArticle")
    lapply(articles, function(a) {
        pmid <- xml_value(a, "./MedlineCitation/PMID")
        title <- xml_value(a, "./MedlineCitation/Article/ArticleTitle")
        if (is.na(pmid) || !grepl("^[0-9]+$", pmid) || is.na(title)) 
            stop("Incomplete PubMed record")
        doi <- xml_value(a, "./PubmedData/ArticleIdList/ArticleId[@IdType=\"doi\"]")
        if (is.na(doi)) 
            doi <- xml_value(a, "./MedlineCitation/Article/ELocationID[@EIdType=\"doi\"]")
        if (!is.na(doi) && !grepl("^10\\.[0-9]{4,9}/[^[:space:]]+$", 
            doi)) 
            doi <- NA_character_
        description <- pubmed_description(a, terms, title)
        types <- xml2::xml_text(xml2::xml_find_all(a, "./MedlineCitation/Article/PublicationTypeList/PublicationType"))
        notices <- xml2::xml_find_all(a, "./MedlineCitation/CommentsCorrectionsList/CommentsCorrections[@RefType=\"RetractionIn\" or @RefType=\"ErratumIn\" or @RefType=\"ExpressionOfConcernIn\"]")
        list(pmid = pmid, title = title, journal = xml_value(a, 
            "./MedlineCitation/Article/Journal/Title"), publication_date = pubmed_date(xml2::xml_find_first(a, 
            "./MedlineCitation/Article/Journal/JournalIssue/PubDate")), 
            electronic_date = pubmed_date(xml2::xml_find_first(a, 
                "./MedlineCitation/Article/ArticleDate[@DateType=\"Electronic\"]")), 
            doi = doi, doi_url = if (!is.na(doi)) paste0("https://doi.org/", 
                utils::URLencode(doi, reserved = TRUE)) else NULL, 
            pubmed_url = paste0("https://pubmed.ncbi.nlm.nih.gov/", 
                pmid, "/"), description = description$text, description_source = description$source, 
            abstract = paste(xml2::xml_text(xml2::xml_find_all(a, 
                "./MedlineCitation/Article/Abstract/AbstractText")), 
                collapse = "\n"), publication_types = as.list(types), 
            notices = lapply(notices, function(n) list(type = xml2::xml_attr(n, 
                "RefType"), pmid = xml_value(n, "./PMID"))), 
            retracted = "Retracted Publication" %in% types || 
                any(xml2::xml_attr(notices, "RefType") == "RetractionIn"))
    })
}
parse_r2_diff_html <-
function (html, group_1, group_2, track) 
{
    doc <- xml2::read_html(html)
    tabs <- rvest::html_table(doc, fill = TRUE, trim = TRUE)
    required <- c("gene", "p", "log2fc", "group", "present")
    selected <- NULL
    for (tab in tabs) {
        names(tab) <- clean_table_names(names(tab))
        if (all(required %in% names(tab))) {
            selected <- tab
            break
        }
    }
    if (is.null(selected)) {
        stop("R2 result did not contain the expected Gene/P/Log2FC/Group/Present table", 
            call. = FALSE)
    }
    out <- data.frame(gene = clean_symbol(selected$gene), adjusted_p = suppressWarnings(as.numeric(selected$p)), 
        r2_log2fc = suppressWarnings(as.numeric(selected$log2fc)), 
        r2_direction = as.character(selected$group), present_n = suppressWarnings(as.integer(selected$present)), 
        stringsAsFactors = FALSE)
    out <- out[nzchar(out$gene) & is.finite(out$r2_log2fc), , 
        drop = FALSE]
    out$patient_effect_group1_minus_group2 <- -out$r2_log2fc
    out$group1_higher <- out$patient_effect_group1_minus_group2 > 
        0
    option_nodes <- xml2::xml_find_all(doc, ".//select[@name='group_1']/option")
    option_text <- trimws(xml2::xml_text(option_nodes))
    parse_count <- function(group) {
        hit <- option_text[grepl(paste0("^", group, " \\("), 
            option_text)]
        if (!length(hit)) 
            return(NA_integer_)
        suppressWarnings(as.integer(sub(".*\\(([0-9]+)\\).*", 
            "\\1", hit[[1]])))
    }
    list(data = out, details = list(track = track, group_1 = group_1, 
        group_2 = group_2, group_1_n = parse_count(group_1), 
        group_2_n = parse_count(group_2), fold_change_definition = "patient_effect_group1_minus_group2 = -R2 Log2FC"))
}
parse_r2_patient_expression <-
function (html, gene, reporter, dataset, transformation = "transform_log2") 
{
    table <- r2_script_payload(html, "/loadDataTableModal.js")$tableData
    plot <- r2_script_payload(html, "/d3/plots/plot.js")
    if (!identical(plot$table, dataset) || !identical(as.character(plot$plotData$reporter), 
        as.character(reporter)) || !identical(toupper(plot$plotData$reporterSymbol), 
        toupper(gene))) 
        stop("R2 returned a different dataset/gene/reporter", 
            call. = FALSE)
    if (!identical(plot$corType, transformation)) 
        stop("R2 did not confirm the requested expression transformation", 
            call. = FALSE)
    if (!is.data.frame(table$data) || !all(c("samplenames", "x") %in% 
        names(table$data))) 
        stop("R2 expression table schema changed", call. = FALSE)
    d <- table$data
    ids <- as.character(d$samplenames)
    if (any(r2_missing(ids)) || anyDuplicated(ids)) 
        stop("R2 expression sample IDs are missing or duplicated", 
            call. = FALSE)
    data.frame(queried_gene = toupper(gene), reporter = as.character(reporter), 
        sample_id = ids, expression = r2_numeric(d$x, "expression"), 
        expression_source_value = as.character(d$x), stringsAsFactors = FALSE)
}
plot_inline_svg <-
function (plot) 
{
    path <- tempfile(fileext = ".svg")
    on.exit(unlink(path))
    grDevices::svg(path, width = 9, height = 6.2999999999999998, 
        onefile = TRUE)
    tryCatch(print(plot), finally = grDevices::dev.off())
    svg <- paste(readLines(path, warn = FALSE), collapse = "\n")
    start <- regexpr("<svg", svg, fixed = TRUE)[1]
    if (start < 1) 
        stop("Plot rendering did not produce SVG")
    svg <- substring(svg, start)
    prefix <- paste0("svg", basename(tempfile()), "-")
    svg <- gsub("id=\"", paste0("id=\"", prefix), svg, fixed = TRUE)
    svg <- gsub("href=\"#", paste0("href=\"#", prefix), svg, 
        fixed = TRUE)
    gsub("url(#", paste0("url(#", prefix), svg, fixed = TRUE)
}
plot_public_depmap <-
function (result) 
{
    d <- result$data
    d <- d[is.finite(d$gene_effect), , drop = FALSE]
    d <- d[order(d$gene_effect, d$ModelID), ]
    d$rank <- seq_len(nrow(d))
    d$context <- ifelse(d$is_medulloblastoma, "Medulloblastoma", 
        "Other models")
    p <- ggplot2::ggplot(d, ggplot2::aes(rank, gene_effect, color = context)) + 
        ggplot2::geom_point(alpha = 0.65000000000000002, size = 1.5) + 
        ggplot2::geom_hline(yintercept = c(0, -1), linetype = 3, 
            color = "grey50") + ggplot2::geom_text(data = d[d$is_medulloblastoma, 
        ], ggplot2::aes(label = CellLineName), size = 2.6000000000000001, 
        angle = 45, hjust = 0, vjust = -0.40000000000000002, 
        check_overlap = FALSE) + ggplot2::scale_color_manual(values = c(Medulloblastoma = "#b43b47", 
        `Other models` = "#b6c0c7")) + ggplot2::theme_classic() + 
        ggplot2::labs(title = paste(unique(d$gene), "dependency across DepMap models"), 
            subtitle = result$provenance$release, x = "Models ordered by gene effect", 
            y = "Chronos gene effect", color = "Disease context")
    mb <- result$data[result$data$is_medulloblastoma, , drop = FALSE]
    mp <- ggplot2::ggplot(mb[is.finite(mb$gene_effect), ], ggplot2::aes(reorder(CellLineName, 
        gene_effect), gene_effect)) + ggplot2::geom_point(color = "#b43b47", 
        size = 3) + ggplot2::geom_hline(yintercept = c(0, -1), 
        linetype = 3, color = "grey50") + ggplot2::coord_flip() + 
        ggplot2::theme_classic() + ggplot2::labs(title = "Medulloblastoma models", 
        x = "DepMap cell-line name", y = "Chronos gene effect")
    list(plot = p, mb_plot = mp, models = mb, all_models = result$data, 
        provenance = result$provenance)
}
plot_r2_subgroup_test <-
function (result, output_dir = result$provenance$output_dir, 
    export_images = FALSE) 
{
    if (!requireNamespace("ggplot2", quietly = TRUE)) 
        stop("Install ggplot2 for the Version 1 plot", call. = FALSE)
    d <- result$data
    eligible <- is.finite(d$expression) & !r2_missing(d$subgroup)
    exclusions <- d[!eligible, c("sample_id", "expression", "subgroup"), 
        drop = FALSE]
    used <- d[eligible, , drop = FALSE]
    sizes <- table(used$subgroup)
    if (length(sizes) < 2 || any(sizes < 3)) 
        stop("Need at least two annotated subgroups with three samples each", 
            call. = FALSE)
    if (anyDuplicated(used$sample_id)) 
        stop("Repeated sample IDs; unpaired comparison is invalid", 
            call. = FALSE)
    test <- stats::kruskal.test(expression ~ subgroup, used)
    summaries <- do.call(rbind, lapply(split(used$expression, 
        used$subgroup), function(x) data.frame(n = length(x), 
        mean = mean(x), median = stats::median(x), sd = stats::sd(x), 
        iqr = stats::IQR(x), min = min(x), max = max(x))))
    summaries$subgroup <- rownames(summaries)
    rownames(summaries) <- NULL
    statistics <- data.frame(test = test$method, statistic = unname(test$statistic), 
        df = unname(test$parameter), p_value = test$p.value, 
        n = nrow(used), groups = length(sizes), excluded_n = nrow(exclusions), 
        epsilon_squared = max(0, (unname(test$statistic) - length(sizes) + 
            1)/(nrow(used) - length(sizes))), multiple_testing = "One predefined omnibus test for one gene; no pairwise tests or multiplicity adjustment in Version 1")
    labels <- setNames(paste0(names(sizes), "\n(n = ", as.integer(sizes), 
        ")"), names(sizes))
    p <- ggplot2::ggplot(used, ggplot2::aes(x = subgroup, y = expression, 
        fill = subgroup)) + ggplot2::geom_boxplot(width = 0.5, 
        outlier.shape = NA, alpha = 0.34999999999999998, linewidth = 0.45000000000000001) + 
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.16, 
            height = 0, seed = 85217), size = 1.1000000000000001, 
            alpha = 0.5, shape = 16) + ggplot2::scale_x_discrete(labels = labels) + 
        ggplot2::guides(fill = "none") + ggplot2::labs(title = paste(unique(used$queried_gene), 
        "expression across medulloblastoma subgroups"), subtitle = sprintf("Kruskal-Wallis H(%d) = %.2f; p = %.3g", 
        test$parameter, test$statistic, test$p.value), x = "Original R2 subgroup labels", 
        y = "Expression (R2 log2 transformation)", caption = paste0("GSE85217 / Cavalli; reporter ", 
            result$provenance$reporter, ". Each point is one R2 sample.\n", 
            "Overall distribution comparison; this test does not establish subgroup specificity. Excluded: ", 
            nrow(exclusions), ".")) + ggplot2::theme_classic(base_size = 12) + 
        ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, 
            size = 9))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    gene <- result$provenance$queried_gene
    files <- list(statistics = file.path(output_dir, paste0(gene, 
        "_subgroup_test.csv")), summaries = file.path(output_dir, 
        paste0(gene, "_subgroup_summaries.csv")), exclusions = file.path(output_dir, 
        "subgroup_test_exclusions.csv"))
    if (isTRUE(export_images)) {
        files$pdf <- file.path(output_dir, paste0(gene, "_subgroup_test.pdf"))
        files$png <- file.path(output_dir, paste0(gene, "_subgroup_test.png"))
        ggplot2::ggsave(files$pdf, p, width = 8, height = 5.7999999999999998, 
            device = grDevices::pdf)
        ggplot2::ggsave(files$png, p, width = 8, height = 5.7999999999999998, 
            dpi = 300)
    }
    data.table::fwrite(statistics, files$statistics)
    data.table::fwrite(summaries, files$summaries)
    data.table::fwrite(exclusions, files$exclusions)
    writeLines("Predefined unpaired Kruskal-Wallis comparison of distributions. Samples are treated as independent primary tumours as described by the source. Batch, age and other confounding are not adjusted in this Version 1 smoke test. No pairwise or subtype/clinical analyses were performed. Missing annotations/expression are reported in the exclusions file. Epsilon-squared is the non-negative rank-based omnibus effect estimate (H-k+1)/(n-k).", 
        file.path(output_dir, "subgroup_test_method.txt"))
    dump(c("r2_missing", "plot_r2_subgroup_test"), file = file.path(output_dir, 
        "analysis_functions.R"))
    writeLines(c("# Open this output directory as the working directory in RStudio.", 
        "# Uses the saved patient data; no network request and no new clinical analysis.", 
        "source('analysis_functions.R')", "retrieval <- readRDS('retrieval.rds')", 
        "validation <- plot_r2_subgroup_test(retrieval, output_dir = 'rerun')", 
        "print(validation$plot)", "validation$statistics"), file.path(output_dir, 
        "reproduce_plot.R"))
    writeLines(capture.output(sessionInfo()), file.path(output_dir, 
        "analysis_session-info.txt"))
    list(statistics = statistics, summaries = summaries, files = files, 
        plot = p)
}
preferred_model_label <-
function (models) 
{
    candidates <- c("CellLineName", "StrippedCellLineName", "CCLEName", 
        "ModelID")
    hit <- candidates[candidates %in% names(models)]
    if (!length(hit)) 
        return(rep(NA_character_, nrow(models)))
    as.character(models[[hit[[1]]]])
}
profile_attempt <-
function (expr) 
tryCatch(expr, error = function(e) list(status = "unavailable", 
    reason = conditionMessage(e)))
profile_clinical <-
function (retrieval, survival_source, methods = c("median", "mean", 
    "optimal")) 
{
    results <- list()
    d <- retrieval$data
    d$time <- d$event <- NA_real_
    if (!is.null(survival_source)) {
        i <- match(d$sample_id, survival_source$data$sample_id)
        d$time <- survival_source$data$time[i]
        d$event <- survival_source$data$event[i]
    }
    for (cohort in c("all", sort(unique(d$subgroup[!r2_missing(d$subgroup)])))) for (method in methods) {
        subset <- if (cohort == "all") 
            d
        else d[which(d$subgroup == cohort), , drop = FALSE]
        complete <- subset[is.finite(subset$expression) & is.finite(subset$time) & 
            !is.na(subset$event), , drop = FALSE]
        attempt <- function(expr) tryCatch(expr, error = function(e) list(status = "unavailable", 
            reason = conditionMessage(e)))
        cut <- attempt(profile_cutoff(subset$expression, method, 
            complete))
        survival <- if (!is.null(cut$value)) 
            attempt(profile_survival(subset, cut, cohort))
        else cut
        met <- if (!isTRUE(retrieval$audit$metastasis$coding_verified)) 
            list(status = "unavailable", reason = "Metastasis coding is not verified")
        else if (!is.null(cut$value)) 
            attempt(profile_metastasis(subset, cut, cohort))
        else cut
        results[[paste(cohort, method, sep = ":")]] <- list(cohort = cohort, 
            method = method, cutoff = cut, survival = survival, 
            metastasis = met)
    }
    for (endpoint in c("survival", "metastasis")) {
        ps <- vapply(results, function(r) {
            s <- r[[endpoint]]$statistics
            if (is.null(s)) 
                return(NA_real_)
            if (endpoint == "survival") 
                s$cutoff_adjusted_p
            else s$p_value
        }, numeric(1))
        adjusted <- p.adjust(ps, "BH", n = length(results))
        for (i in seq_along(results)) if (!is.null(results[[i]][[endpoint]]$statistics)) 
            results[[i]][[endpoint]]$statistics$family_adjusted_p <- adjusted[i]
    }
    results
}
profile_cutoff <-
function (expression, method, survival_data = NULL) 
{
    method <- match.arg(method, c("median", "mean", "optimal"))
    x <- expression[is.finite(expression)]
    if (length(x) < 10 || length(unique(x)) < 2) 
        stop("Insufficient variable expression for a high/low split")
    if (method != "optimal") 
        return(list(value = if (method == "median") median(x) else mean(x), 
            method = method, candidates = 1L, scan = data.frame(), 
            reference_n = length(x)))
    if (is.null(survival_data)) 
        stop("Optimised cutoff requires verified survival data")
    d <- survival_data
    values <- sort(unique(d$expression))
    cuts <- head(values, -1)
    minimum <- max(10L, ceiling(0.10000000000000001 * nrow(d)))
    rows <- lapply(cuts, function(cut) {
        high <- d$expression > cut
        if (min(sum(high), sum(!high)) < minimum) 
            return(NULL)
        test <- tryCatch(survival::survdiff(survival::Surv(time, 
            event) ~ high, data = d), error = function(e) NULL)
        if (is.null(test) || !is.finite(test$chisq)) 
            return(NULL)
        data.frame(cutoff = cut, low_n = sum(!high), high_n = sum(high), 
            p_value = pchisq(test$chisq, 1, lower.tail = FALSE))
    })
    scan <- do.call(rbind, rows)
    if (is.null(scan) || !nrow(scan)) 
        stop("No eligible optimal cutoffs (minimum 10 patients and 10% per group)")
    scan$adjusted_p <- p.adjust(scan$p_value, "bonferroni")
    best <- which.min(scan$p_value)
    list(value = scan$cutoff[best], method = method, candidates = nrow(scan), 
        scan = scan, reference_n = nrow(d))
}
profile_enrichr <-
function (cfg, gene, libraries = c("GO_Biological_Process_2026", 
    "Reactome_Pathways_2024", "KEGG_2021_Human", "WikiPathways_2024_Human")) 
{
    gene <- toupper(gene)
    inputs <- setNames(lapply(libraries, function(l) load_gene_sets(cfg, 
        l)), libraries)
    rows <- lapply(names(inputs), function(l) {
        sets <- inputs[[l]]$sets
        terms <- names(sets)[vapply(sets, function(s) gene %in% 
            s, logical(1))]
        data.frame(library = rep(l, length(terms)), term = terms, 
            set_size = as.integer(lengths(sets[terms])), evidence = rep("Annotated gene membership; not enrichment or pathway activity", 
                length(terms)))
    })
    membership <- do.call(rbind, rows)
    chea <- load_gene_sets(cfg, "ChEA_2022")
    labels <- names(chea$sets)
    selected <- labels[startsWith(labels, paste0(gene, "_")) & 
        grepl("HUMAN", labels, ignore.case = TRUE)]
    target_sets <- chea$sets[selected]
    targets <- if (length(target_sets)) 
        do.call(rbind, lapply(names(target_sets), function(term) data.frame(source_term = term, 
            gene = target_sets[[term]], evidence = "ChEA ChIP binding-associated target; not proof of regulation in MB")))
    else data.frame()
    enrichments <- list()
    human <- chea$sets[grepl("HUMAN", names(chea$sets), ignore.case = TRUE)]
    universe <- unique(unlist(human, use.names = FALSE))
    for (term in names(target_sets)) for (l in names(inputs)) {
        tab <- overrepresentation(target_sets[[term]], universe, 
            inputs[[l]]$sets)
        if (nrow(tab)) {
            tab$source_term <- term
            tab$library <- l
            enrichments[[paste(term, l, sep = ":")]] <- tab
        }
    }
    tf <- tryCatch(get_tf_annotation(cfg, gene), error = function(e) list(classification = paste("Unavailable:", 
        conditionMessage(e))))
    list(membership = membership, targets = targets, tf = tf, 
        enrichment = if (length(enrichments)) do.call(rbind, 
            enrichments) else data.frame(), provenance = list(libraries = lapply(inputs, 
            `[[`, "provenance"), tf_source = chea$provenance, 
            tf_target_sets = as.list(selected), background_n = length(universe), 
            method = "TF target over-representation: hypergeometric; BH over all eligible terms separately per target experiment and library; human ChEA annotated genes as background. This is an annotation background, not a measured ChIP assay universe.", 
            interpretation = "Single-gene membership does not establish pathway deregulation. TF ChIP-associated targets are kept separate from perturbation, predictions and correlations; no MB-specific regulation is inferred. No user transcriptomics data were supplied."))
}
profile_extra_html <-
function (profile, default_cutoff = "median") 
{
    clinical <- profile$clinical
    panels <- paste(vapply(clinical, function(x) {
        paste0("<div class='clinical-panel' data-cohort='", html_escape(x$cohort), 
            "' data-cutoff='", html_escape(x$method), "'", if (x$cohort != 
                "all" || x$method != default_cutoff) 
                " hidden"
            else "", ">", "<h3>Overall survival</h3>", profile_plot_html(x$survival), 
            if (x$method == "optimal") 
                "<p class='notice'>Optimised in these same patients: the cutoff-adjusted p-value uses Bonferroni correction over all eligible cutoffs. Hazard ratios and confidence intervals are descriptive after selection and may be optimistic. This is not independent validation.</p>", 
            "<h3>Metastasis status</h3>", profile_plot_html(x$metastasis), 
            "<details><summary>Group counts and cutoff search</summary>", 
            html_table(x$metastasis$counts), html_table(x$cutoff$scan), 
            "</details></div>")
    }, character(1)), collapse = "")
    cohorts <- unique(vapply(clinical, `[[`, character(1), "cohort"))
    controls <- paste0("<label for='clinical-cohort'>Cohort</label> <select id='clinical-cohort'>", 
        paste0("<option value=\"", html_escape(cohorts), "\">", 
            html_escape(cohorts), "</option>", collapse = ""), 
        "</select> ", "<label for='clinical-cutoff'>Expression cutoff</label> <select id='clinical-cutoff'>", 
        paste(vapply(c("median", "mean", "optimal"), function(m) paste0("<option value=\"", 
            m, "\"", if (m == default_cutoff) 
                " selected"
            else "", ">", m, "</option>"), character(1)), collapse = ""), 
        "</select>")
    pf <- profile$pfister
    en <- profile$enrichr
    dep <- profile$depmap
    en_html <- if (identical(en$status, "unavailable")) 
        paste0("<p class=\"notice\">", html_escape(en$reason), 
            "</p>")
    else {
        sig <- en$enrichment
        if (nrow(sig)) 
            sig <- sig[sig$adjusted_p <= 0.050000000000000003, 
                , drop = FALSE]
        if (nrow(sig)) 
            sig <- sig[order(sig$adjusted_p), , drop = FALSE]
        paste0("<p>Known gene-set memberships, not evidence that a pathway is activated, inhibited or deregulated.</p>", 
            html_table(en$membership), "<h3>Transcription factor annotation</h3><p>", 
            html_escape(en$tf$classification), " (", html_escape(en$tf$version %||% 
                "source unavailable"), ")</p>", "<h3>Binding-associated TF targets</h3>", 
            if (!nrow(en$targets)) 
                "<p>No exact human TF-target experiment for this gene was found in ChEA 2022. This does not mean the gene has no targets.</p>"
            else paste0("<details><summary>Inspect target genes and source experiments</summary>", 
                html_table(en$targets), "</details>"), "<h3>TF-target enrichment</h3>", 
            if (!nrow(sig)) 
                "<p>No eligible/significant TF-target enrichment results at BH FDR 0.05.</p>"
            else html_table(head(sig, 100)), "<p>", html_escape(en$provenance$method), 
            "</p><p>Top 100 significant rows are shown; all tested terms are retained in the R object and CSV files.</p>")
    }
    dep_html <- if (identical(dep$status, "unavailable")) 
        paste0("<p class=\"notice\">", html_escape(dep$reason), 
            "</p>")
    else paste0("<p>", html_escape(dep$provenance$note), "</p>", 
        profile_plot_html(dep), "<div class='figure'>", plot_inline_svg(dep$mb_plot), 
        "</div>", html_table(dep$models), "<p>", html_escape(dep$provenance$model_selection), 
        ". ", html_escape(dep$provenance$interpretation), "</p>")
    paste0("<section class='card' id='profile-summary'><h2>Report coverage</h2>", 
        html_table(profile$status), "<p>Separate evidence layers: expression distributions, unadjusted clinical associations, functional annotations, in-vitro dependency and literature. No composite therapeutic score.</p>", 
        "<nav><a href=\"#subtypes\">Subtypes</a><a href=\"#pfister\">CNS tumours</a><a href=\"#clinical\">Survival and metastasis</a><a href=\"#enrichr\">Enrichr / TF targets</a><a href=\"#depmap\">DepMap</a></nav></section>", 
        "<section class='card' id='subtypes'><h2>Cavalli molecular subtypes</h2>", 
        profile_plot_html(profile$subtypes), "<p>All original subtype labels are preserved. Pairwise rank tests compare distributions; they do not prove subtype specificity.</p></section>", 
        "<section class='card' id='pfister'><h2>Expression across CNS tumours: Pfister 272</h2>", 
        profile_plot_html(pf), "<p>FPKM is transformed locally as log2(1 + FPKM). This panel is separate from Cavalli microarray expression. CNS labels are explicitly listed in provenance; only primary samples enter the comparison.</p>", 
        "<details><summary>Source group counts and excluded samples</summary>", 
        html_table(pf$all_counts), html_table(pf$excluded_samples), 
        "</details></section>", "<section class='card' id='clinical'><h2>Survival and metastasis</h2>", 
        controls, "<p>Mean and median are calculated on all finite log2 expression values within the selected cohort, before clinical missing-data exclusions; the same threshold is used for survival and metastasis. High &gt; cutoff; Low &lt;= cutoff. The optimal threshold uses survival-complete cases. Subgroup cutoffs are recalculated within each subgroup.</p>", 
        "<p>Family-adjusted p-values use BH across all five cohorts and all three cutoff methods, separately for survival and metastasis. Survival first corrects the optimised cutoff search. Pooled results are not adjusted for subgroup, age or treatment. Check cohort-specific estimates and the Cox proportional-hazards diagnostic.</p>", 
        panels, "<details><summary>Survival source verification</summary><pre>", 
        html_escape(json_text(profile$survival_provenance)), 
        "</pre></details></section>", "<section class='card' id='enrichr'><h2>Enrichr functional evidence</h2>", 
        en_html, "</section>", "<section class='card' id='depmap'><h2>DepMap dependencies</h2>", 
        dep_html, "</section>", "<script>function updateClinical(){const c=document.getElementById(\"clinical-cohort\").value;const m=document.getElementById(\"clinical-cutoff\").value;document.querySelectorAll(\".clinical-panel\").forEach(x=>{x.hidden=x.dataset.cohort!==c||x.dataset.cutoff!==m;});}document.getElementById(\"clinical-cohort\").addEventListener(\"change\",updateClinical);document.getElementById(\"clinical-cutoff\").addEventListener(\"change\",updateClinical);updateClinical();</script>")
}
profile_group_comparison <-
function (data, field, title, y_label = "Expression (R2 log2)", 
    pairwise = TRUE) 
{
    d <- data.frame(sample_id = data$sample_id, expression = data$expression, 
        group = as.character(data[[field]]))
    ok <- is.finite(d$expression) & !r2_missing(d$group)
    excluded <- d[!ok, , drop = FALSE]
    d <- d[ok, , drop = FALSE]
    sizes <- table(d$group)
    if (length(sizes) < 2 || anyDuplicated(d$sample_id)) 
        stop("Need independent samples in at least two groups")
    summaries <- do.call(rbind, lapply(names(sizes), function(g) {
        x <- d$expression[d$group == g]
        data.frame(group = g, n = length(x), mean = mean(x), 
            median = median(x), iqr = IQR(x))
    }))
    kw <- kruskal.test(expression ~ group, d)
    comparisons <- data.frame()
    if (pairwise) {
        pairs <- combn(names(sizes), 2, simplify = FALSE)
        comparisons <- do.call(rbind, lapply(pairs, function(g) {
            x <- d$expression[d$group == g[1]]
            y <- d$expression[d$group == g[2]]
            test <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
            data.frame(group1 = g[1], group2 = g[2], n1 = length(x), 
                n2 = length(y), median_difference = median(x) - 
                  median(y), rank_biserial = 2 * as.numeric(test$statistic)/(length(x) * 
                  length(y)) - 1, p_value = test$p.value)
        }))
        comparisons$adjusted_p <- p.adjust(comparisons$p_value, 
            "BH")
    }
    statistics <- data.frame(test = "Kruskal-Wallis", n = nrow(d), 
        excluded_n = nrow(excluded), groups = length(sizes), 
        statistic = unname(kw$statistic), p_value = kw$p.value, 
        epsilon_squared = max(0, (unname(kw$statistic) - length(sizes) + 
            1)/(nrow(d) - length(sizes))))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = expression, 
        fill = group)) + ggplot2::geom_boxplot(outlier.shape = NA, 
        alpha = 0.40000000000000002) + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.17000000000000001, 
        height = 0, seed = 85217), size = 0.80000000000000004, 
        alpha = 0.5) + ggplot2::scale_x_discrete(labels = setNames(paste0(names(sizes), 
        " (n=", sizes, ")"), names(sizes))) + ggplot2::coord_flip() + 
        ggplot2::guides(fill = "none") + ggplot2::theme_classic(base_size = 11) + 
        ggplot2::labs(title = title, subtitle = sprintf("Kruskal-Wallis p = %.3g; %d samples", 
            kw$p.value, nrow(d)), x = "Original source labels", 
            y = y_label, caption = "Exploratory unadjusted distribution comparison; small groups remain visible. Pairwise tests use BH correction.")
    list(plot = p, statistics = statistics, summaries = summaries, 
        pairwise = comparisons, exclusions = excluded, data = d)
}
profile_metastasis <-
function (data, cutoff, cohort) 
{
    if (!all(na.omit(data$metastasis_raw[!r2_missing(data$metastasis_raw)]) %in% 
        c("0", "1"))) 
        stop("Unexpected metastasis codes")
    eligible <- is.finite(data$expression) & !r2_missing(data$metastasis_raw)
    d <- data[eligible, , drop = FALSE]
    d$group <- factor(ifelse(d$expression > cutoff$value, "High", 
        "Low"), levels = c("Low", "High"))
    d$metastasis <- factor(ifelse(d$metastasis_raw == "1", "Metastatic", 
        "M0"), levels = c("M0", "Metastatic"))
    tab <- table(d$group, d$metastasis)
    if (any(rowSums(tab) < 5) || any(colSums(tab) == 0)) 
        stop("Metastasis comparison needs >=5 per expression group and both metastasis states")
    fisher <- fisher.test(tab)
    counts <- as.data.frame(tab)
    names(counts) <- c("expression_group", "metastasis", "n")
    counts$fraction <- counts$n/as.numeric(rowSums(tab)[as.character(counts$expression_group)])
    statistics <- data.frame(cohort = cohort, cutoff_method = cutoff$method, 
        cutoff = cutoff$value, n = nrow(d), excluded_n = sum(!eligible), 
        odds_ratio_high_vs_low = unname(fisher$estimate), ci_lower = fisher$conf.int[1], 
        ci_upper = fisher$conf.int[2], p_value = fisher$p.value, 
        interpretation = if (cutoff$method == "optimal") 
            "Exploratory: threshold selected using survival in overlapping patients; not an independent validation"
        else "Fisher exact test; unadjusted for clinical confounding")
    p <- ggplot2::ggplot(counts, ggplot2::aes(expression_group, 
        fraction, fill = metastasis)) + ggplot2::geom_col(width = 0.59999999999999998) + 
        ggplot2::geom_text(ggplot2::aes(label = paste0("n=", 
            n)), position = ggplot2::position_stack(vjust = 0.5), 
            color = "white") + ggplot2::scale_fill_manual(values = c(M0 = "#477b9e", 
        Metastatic = "#b43b47")) + ggplot2::theme_classic(base_size = 12) + 
        ggplot2::labs(title = paste("Metastasis by expression:", 
            cohort), subtitle = sprintf("%s cutoff %.4g; Fisher p %.3g", 
            cutoff$method, cutoff$value, fisher$p.value), x = "Expression group", 
            y = "Fraction within expression group", fill = "Verified status")
    list(plot = p, statistics = statistics, counts = counts, 
        data = d, exclusions = data[!eligible, , drop = FALSE])
}
profile_plot_html <-
function (x) 
{
    if (is.null(x$plot)) 
        return(paste0("<p class='notice'>Unavailable: ", html_escape(x$reason %||% 
            "No eligible observations"), "</p>"))
    paste0("<div class='figure'>", plot_inline_svg(x$plot), "</div>", 
        html_table(x$statistics), if (!is.null(x$risk_table)) 
            paste0("<details><summary>Number at risk</summary>", 
                html_table(x$risk_table), "</details>"), if (!is.null(x$pairwise)) 
            paste0("<details><summary>Pairwise comparisons (BH-adjusted)</summary>", 
                html_table(x$pairwise), "</details>"), if (!is.null(x$summaries)) 
            paste0("<details><summary>Group summaries</summary>", 
                html_table(x$summaries), "</details>"))
}
profile_survival <-
function (data, cutoff, cohort) 
{
    d <- data
    eligible <- is.finite(d$expression) & is.finite(d$time) & 
        d$time >= 0 & !is.na(d$event) & d$event %in% c(0, 1)
    excluded <- d[!eligible, , drop = FALSE]
    d <- d[eligible, , drop = FALSE]
    if (nrow(d) < 10 || sum(d$event) < 5) 
        stop("Survival needs >=10 complete samples and >=5 deaths")
    d$group <- factor(ifelse(d$expression > cutoff$value, "High", 
        "Low"), levels = c("Low", "High"))
    if (any(table(d$group) < 5)) 
        stop("Survival split needs >=5 samples in each expression group")
    fit <- survival::survfit(survival::Surv(time, event) ~ group, 
        data = d, conf.type = "log-log")
    test <- survival::survdiff(survival::Surv(time, event) ~ 
        group, data = d)
    rawp <- pchisq(test$chisq, 1, lower.tail = FALSE)
    notes <- character()
    cox <- withCallingHandlers(tryCatch(survival::coxph(survival::Surv(time, 
        event) ~ group, data = d, x = TRUE), error = function(e) {
        notes <<- c(notes, conditionMessage(e))
        NULL
    }), warning = function(w) {
        notes <<- c(notes, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    hr <- lo <- hi <- ph <- NA_real_
    if (!is.null(cox) && !length(notes)) {
        cf <- summary(cox)$conf.int
        hr <- cf[1, "exp(coef)"]
        lo <- cf[1, "lower .95"]
        hi <- cf[1, "upper .95"]
        ph <- tryCatch(survival::cox.zph(cox)$table[1, "p"], 
            error = function(e) NA_real_)
    }
    curves <- data.frame(time = fit$time, survival = fit$surv, 
        lower = fit$lower, upper = fit$upper, censored = fit$n.censor, 
        group = rep(sub("group=", "", names(fit$strata), fixed = TRUE), 
            as.integer(fit$strata)))
    curves <- rbind(data.frame(time = 0, survival = 1, lower = 1, 
        upper = 1, censored = 0, group = c("Low", "High")), curves)
    curves <- curves[order(curves$group, curves$time), ]
    times <- pretty(c(0, max(d$time)), n = 6)
    times <- times[times >= 0 & times <= max(d$time)]
    risk <- do.call(rbind, lapply(c("Low", "High"), function(g) data.frame(group = g, 
        time_years = times, n_at_risk = vapply(times, function(t) sum(d$group == 
            g & d$time >= t), integer(1)))))
    statistics <- data.frame(cohort = cohort, cutoff_method = cutoff$method, 
        cutoff = cutoff$value, n = nrow(d), deaths = sum(d$event), 
        low_n = sum(d$group == "Low"), high_n = sum(d$group == 
            "High"), excluded_n = nrow(excluded), logrank_p = rawp, 
        cutoff_adjusted_p = min(1, rawp * cutoff$candidates), 
        candidate_cutoffs = cutoff$candidates, hazard_ratio_high_vs_low = hr, 
        ci_lower = lo, ci_upper = hi, ph_test_p = ph, cox_note = paste(notes, 
            collapse = "; "))
    p <- ggplot2::ggplot(curves, ggplot2::aes(time, survival, 
        color = group, group = group)) + ggplot2::geom_step(linewidth = 0.80000000000000004) + 
        ggplot2::geom_step(ggplot2::aes(y = lower), linetype = 3, 
            alpha = 0.40000000000000002, na.rm = TRUE) + ggplot2::geom_step(ggplot2::aes(y = upper), 
        linetype = 3, alpha = 0.40000000000000002, na.rm = TRUE) + 
        ggplot2::geom_point(data = curves[curves$censored > 0, 
            ], shape = 3, size = 1.8) + ggplot2::scale_color_manual(values = c(High = "#b43b47", 
        Low = "#236fa1")) + ggplot2::coord_cartesian(ylim = c(0, 
        1)) + ggplot2::theme_classic(base_size = 12) + ggplot2::labs(title = paste("Overall survival:", 
        cohort), subtitle = sprintf("%s cutoff %.4g; low %d / high %d; log-rank p %.3g", 
        cutoff$method, cutoff$value, statistics$low_n, statistics$high_n, 
        rawp), x = "Follow-up (years)", y = "Survival probability", 
        color = "Expression", caption = "High > cutoff; Low <= cutoff. Dashed lines: pointwise 95% CI; + censoring. Unadjusted exploratory association.")
    list(plot = p, statistics = statistics, risk_table = risk, 
        curves = curves, data = d, exclusions = excluded, cutoff = cutoff)
}
publication_cards <-
function (publications) 
{
    if (!length(publications$articles)) 
        return("<p>No publications matched this search. The disease filter, alias choices and search date are shown above.</p>")
    paste(vapply(publications$articles, function(p) {
        notices <- c(if (isTRUE(p$retracted)) "Retracted publication", 
            vapply(p$notices, function(n) n$type, character(1)))
        paste0("<article class='paper'>", if (length(notices)) 
            paste0("<p class='notice'>", html_escape(paste(unique(notices), 
                collapse = "; ")), " &mdash; inspect the PubMed record.</p>"), 
            "<h3><a href='", html_escape(p$pubmed_url), "' target='_blank' rel='noopener noreferrer'>", 
            html_escape(p$title), "</a></h3>", "<p class='paper-meta'><strong>", 
            html_escape(p$journal), "</strong><br>Publication date: ", 
            html_escape(p$publication_date), if (!is.na(p$electronic_date)) 
                paste0(" &middot; Electronic publication: ", 
                  html_escape(p$electronic_date)), "</p>", "<p class='description'>", 
            html_escape(p$description), "</p>", "<p class='small'>", 
            if (identical(p$description_source, "abstract_sentence_excerpt")) 
                "One-sentence abstract excerpt; the authors' wording."
            else "No abstract supplied; title-only description.", 
            "</p>", "<p class='paper-links'>PMID ", html_escape(p$pmid), 
            " &middot; DOI: ", if (is.null(p$doi_url)) 
                "Not reported in PubMed"
            else paste0("<a href='", html_escape(p$doi_url), 
                "' target='_blank' rel='noopener noreferrer'>", 
                html_escape(p$doi), "</a>"), "</p></article>")
    }, character(1)), collapse = "\n")
}
pubmed_date <-
function (node) 
{
    if (inherits(node, "xml_missing")) 
        return(NA_character_)
    medline <- xml_value(node, "./MedlineDate")
    if (!is.na(medline)) 
        return(medline)
    parts <- vapply(c("Year", "Month", "Day", "Season"), function(x) xml_value(node, 
        paste0("./", x)), character(1))
    if (!any(!is.na(parts))) 
        return(NA_character_)
    paste(parts[!is.na(parts)], collapse = " ")
}
pubmed_description <-
function (article, terms, title) 
{
    blocks <- xml2::xml_find_all(article, "./MedlineCitation/Article/Abstract/AbstractText")
    candidates <- character()
    scores <- numeric()
    for (b in blocks) {
        text <- trimws(gsub("[\\s\\x{00A0}]+", " ", xml2::xml_text(b), 
            perl = TRUE))
        sentences <- unlist(strsplit(text, "(?<=[.!?])\\s+(?=[A-Z0-9])", 
            perl = TRUE))
        sentences <- sentences[nzchar(sentences)]
        label <- paste(xml2::xml_attr(b, "Label"), xml2::xml_attr(b, 
            "NlmCategory"))
        matches <- vapply(sentences, pubmed_mentions, logical(1), 
            terms = terms)
        score <- 10 * matches + 3 * grepl("conclu", label, ignore.case = TRUE) + 
            2 * grepl("result", label, ignore.case = TRUE)
        candidates <- c(candidates, sentences)
        scores <- c(scores, score)
    }
    if (!length(candidates)) 
        return(list(text = paste0("No abstract is available in PubMed; the publication is titled <U+201C>", 
            sub("[.]$", "", title), "<U+201D>."), source = "title_only_no_abstract"))
    list(text = candidates[which.max(scores)], source = "abstract_sentence_excerpt")
}
pubmed_gene_terms <-
function (gene, aliases = character()) 
{
    if (length(gene) != 1L || !grepl("^[A-Za-z][A-Za-z0-9._-]*$", 
        gene)) 
        stop("Provide one gene symbol", call. = FALSE)
    aliases <- unlist(aliases, use.names = FALSE)
    terms <- unique(c(toupper(gene), aliases))
    if (anyNA(terms) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9 ._()-]{0,100}$", 
        terms))) 
        stop("Invalid gene alias", call. = FALSE)
    terms
}
pubmed_mentions <-
function (text, terms) 
{
    text <- gsub("(*UTF)[\\x{2010}-\\x{2015}\\x{2212}]", "-", 
        enc2utf8(text), perl = TRUE)
    any(vapply(terms, function(term) grepl(paste0("(?i)(?<![A-Za-z0-9_-])\\Q", 
        term, "\\E(?![A-Za-z0-9_-])"), text, perl = TRUE), logical(1)))
}
pubmed_query <-
function (terms, disease = NULL) 
{
    q <- paste0("(", paste0("\"", terms, "\"[Title/Abstract]", 
        collapse = " OR "), ")")
    if (!is.null(disease)) {
        if (length(disease) != 1 || !grepl("^[A-Za-z0-9][A-Za-z0-9 -]{0,100}$", 
            disease)) 
            stop("Invalid disease text", call. = FALSE)
        q <- paste0(q, " AND \"", disease, "\"[Title/Abstract]")
    }
    q
}
r2_compare_groups <-
function (cfg, track = NULL, group_1 = NULL, group_2 = NULL, 
    test = "anova", top_n = 1000L, p_threshold = 0.050000000000000003, 
    transformation = "transform_log2", min_present = 1L, refresh = FALSE) 
{
    track <- track %||% cfg$r2_default_track
    group_1 <- tolower(group_1 %||% cfg$r2_default_group_1)
    group_2 <- tolower(group_2 %||% cfg$r2_default_group_2)
    top_n <- as_integer_scalar(top_n, 1000L, 1L, 5000L)
    min_present <- as_integer_scalar(min_present, 1L, 0L, 100000L)
    p_threshold <- as_number(p_threshold, 0.050000000000000003)
    parts <- list(dataset = cfg$r2_dataset_table, track = track, 
        group1 = group_1, group2 = group_2, test = test, top = top_n, 
        p = p_threshold, transform = transformation, present = min_present)
    key <- cache_key("r2diff", parts)
    if (!isTRUE(refresh)) {
        cached <- cache_get(cfg, key)
        if (!is.null(cached)) {
            cached$provenance$cache_hit <- TRUE
            return(cached)
        }
    }
    parsed <- r2_fetch_diff(cfg, track, group_1, group_2, test, 
        top_n, p_threshold, transformation, min_present)
    result <- list(data = parsed$data, details = parsed$details, 
        provenance = list(source = "R2 Genomics Analysis and Visualization Platform", 
            endpoint = cfg$r2_base_url, dataset_table = cfg$r2_dataset_table, 
            dataset_label = cfg$r2_dataset_label, contrast = paste0(track, 
                ":", group_1, "_vs_", group_2), test = test, 
            multiple_testing = "False Discovery Rate (R2)", transformation = transformation, 
            p_threshold = p_threshold, top_n_requested = top_n, 
            retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", 
                tz = "UTC"), cache_hit = FALSE, connector = "unofficial read-only R2 HTML workflow"))
    cache_put(cfg, key, result)
    result
}
r2_fetch_diff <-
function (cfg, track, group_1, group_2, test = "anova", top_n = 1000L, 
    p_threshold = 0.050000000000000003, transformation = "transform_log2", 
    min_present = 1L) 
{
    if (identical(group_1, group_2)) 
        stop("group_1 and group_2 must differ", call. = FALSE)
    allowed_tests <- c("anova", "kruskal", "log2fc", "limma")
    if (!test %in% allowed_tests) 
        stop("Unsupported R2 test: ", test, call. = FALSE)
    handle <- r2_session()
    dataset_url <- paste0(cfg$r2_base_url, "?table=", utils::URLencode(cfg$r2_dataset_table, 
        reserved = TRUE))
    html_1 <- r2_request(handle, dataset_url)
    assert_r2_page(html_1, cfg$r2_dataset_table, "dataset selection")
    html_2 <- r2_request(handle, cfg$r2_base_url, list(perspective = "singleds", 
        table = cfg$r2_dataset_table, option = "displaygene_two_group_diff", 
        button1 = "Next"))
    assert_r2_page(html_2, "Two-group differential expression|Select a test", 
        "analysis selection")
    html_3 <- r2_request(handle, cfg$r2_base_url, list(table = "", 
        option = "", test = test, grouping_track = track, subsettracksubset = "", 
        subset = ""))
    assert_r2_page(html_3, "group_1|Group 1", "group selection")
    html_4 <- r2_request(handle, cfg$r2_base_url, list(test = "", 
        subset = "", grouping_track = "", table = "", option = "", 
        test_mode = "", display = "list", group_1 = group_1, 
        group_2 = group_2, floor = "", cortype = transformation, 
        mtc = "fdr", minpval = format(p_threshold, scientific = FALSE, 
            trim = TRUE), top_x = as.character(top_n), hugoonce = "yes", 
        minpres = as.character(min_present), minmax = "", mindif = "0", 
        gopath = "", goid = "", geneset = ""))
    assert_r2_page(html_4, "Scan result for track|combinations meet your criteria", 
        "differential-expression result")
    parse_r2_diff_html(html_4, group_1, group_2, track)
}
r2_json <-
function (text, context) 
{
    tryCatch(jsonlite::fromJSON(text, simplifyVector = TRUE), 
        error = function(e) stop("R2 returned an unexpected response for ", 
            context, "; inspect saved source files. ", conditionMessage(e), 
            call. = FALSE))
}
r2_metadata_audit <-
function (metadata, dataset_info, track_info) 
{
    rows <- lapply(names(metadata), function(field) {
        value <- as.character(metadata[[field]])
        list(field = field, source_type = "string", n = length(value), 
            missing_n = sum(r2_missing(value)), unique_values = if (length(unique(value)) <= 
                30L) as.list(sort(unique(value))) else NULL, 
            value_counts = if (length(unique(value)) <= 30L) as.list(table(value, 
                useNA = "ifany")) else NULL)
    })
    subtypes <- if (all(c("subgroup", "subtype") %in% names(metadata))) {
        as.data.frame(table(subgroup = metadata$subgroup, subtype = metadata$subtype), 
            stringsAsFactors = FALSE)
    }
    else NULL
    if (!is.null(subtypes)) 
        subtypes <- subtypes[subtypes$Freq > 0, , drop = FALSE]
    met_field <- "met_status_(1_met__0_m0)"
    verified <- met_field %in% names(metadata) && any(grepl(paste0(met_field, 
        ":"), track_info$trackDescriptions, fixed = TRUE))
    if (verified && any(!r2_missing(metadata[[met_field]]) & 
        !metadata[[met_field]] %in% c("0", "1"))) 
        stop("Unexpected metastasis codes despite declared 0/1 mapping; no labels derived", 
            call. = FALSE)
    list(fields = rows, dataset = dataset_info, track_descriptions = as.list(track_info$trackDescriptions), 
        subtype_by_subgroup = subtypes, classification = list(field = if ("subtype" %in% 
            names(metadata)) "subtype" else NULL, source_publication = dataset_info$pubmedId, 
            note = "Original R2 labels and observed subgroup/subtype cross-tabulation. No classifier, relabelling, or alternative classification inferred."), 
        metastasis = list(field = if (met_field %in% names(metadata)) met_field else NULL, 
            coding_verified = verified, code_0 = if (verified) "M0 (no metastasis)" else NULL, 
            code_1 = if (verified) "Metastatic" else NULL, evidence = if (verified) "R2's exact source field name met_status_(1_met__0_m0), present in both annotation data and dataset track descriptions" else "Coding unresolved; labels not generated", 
            missing_definition = "Source token na is retained in raw data and treated as missing only in derived fields"), 
        survival = list(time_field = if ("os_(years)" %in% names(metadata)) "os_(years)" else NULL, 
            time_unit = if ("os_(years)" %in% names(metadata)) "years (explicit in source field name)" else NULL, 
            status_field = if ("dead" %in% names(metadata)) "dead" else NULL, 
            event_coding_verified = FALSE, note = "Observed dead codes are retained verbatim. An independent codebook for event semantics has not been retrieved; no survival model or event recoding is performed."), 
        identifier_note = dataset_info$design, missing_tokens = as.list(c("", 
            "na", "NA")))
}
r2_missing <-
function (x) 
is.na(x) | trimws(as.character(x)) %in% c("", "na", "NA")
r2_numeric <-
function (x, label) 
{
    missing <- r2_missing(x)
    value <- suppressWarnings(as.numeric(x))
    if (any(!missing & !is.finite(value))) 
        stop("Invalid numeric values in ", label, call. = FALSE)
    value[missing] <- NA_real_
    value
}
r2_request <-
function (handle, url, fields = NULL) 
{
    if (is.null(fields)) {
        curl::handle_setopt(handle, httpget = TRUE)
    }
    else {
        curl::handle_setheaders(handle, `Content-Type` = "application/x-www-form-urlencoded")
        curl::handle_setopt(handle, post = TRUE, postfields = encode_form(fields))
    }
    response <- curl::curl_fetch_memory(url, handle = handle)
    if (response$status_code < 200L || response$status_code >= 
        300L) {
        stop("R2 returned HTTP ", response$status_code, call. = FALSE)
    }
    rawToChar(response$content)
}
r2_script_payload <-
function (html, module) 
{
    doc <- xml2::read_html(charToRaw(html), encoding = "UTF-8")
    scripts <- xml2::xml_text(xml2::xml_find_all(doc, ".//script"))
    candidates <- scripts[grepl(module, scripts, fixed = TRUE)]
    if (length(candidates) != 1L) 
        stop("Expected exactly one R2 data payload for ", module, 
            call. = FALSE)
    marker <- ".then(m => m['default']("
    start <- regexpr(marker, candidates, fixed = TRUE)[1]
    if (start < 0) 
        stop("R2 script wrapper changed: ", module, call. = FALSE)
    payload <- trimws(substring(candidates, start + nchar(marker)))
    if (!endsWith(payload, "));")) 
        stop("Unexpected R2 script suffix", call. = FALSE)
    r2_json(substr(payload, 1, nchar(payload) - 3L), module)
}
r2_session <-
function (timeout_seconds = 240) 
{
    h <- curl::new_handle()
    curl::handle_setopt(h, followlocation = TRUE, timeout = timeout_seconds, 
        connecttimeout = 20, cookiefile = "", useragent = "r2-depmap-mcp/0.1.0 (+read-only research connector)")
    h
}
read_gene_effect_columns <-
function (path, genes) 
{
    mapping <- depmap_gene_columns(path)
    requested <- unique(toupper(clean_symbol(genes)))
    hits <- mapping[toupper(mapping$gene) %in% requested, , drop = FALSE]
    if (anyDuplicated(toupper(hits$gene))) 
        stop("Ambiguous DepMap symbol maps to multiple gene columns; resolve Entrez IDs first", 
            call. = FALSE)
    if (!nrow(hits)) 
        stop("None of the requested genes were found in the DepMap gene-effect file", 
            call. = FALSE)
    header <- depmap_header(path)
    id_col <- depmap_id_column(header)
    data <- data.table::fread(path, select = c(id_col, hits$column), 
        data.table = FALSE, check.names = FALSE, showProgress = FALSE)
    if (anyNA(data[[id_col]]) || anyDuplicated(data[[id_col]])) 
        stop("DepMap gene-effect IDs must be present and unique", 
            call. = FALSE)
    if (any(!vapply(data[hits$column], is.numeric, logical(1)))) 
        stop("DepMap gene effects must be numeric", call. = FALSE)
    list(data = data, mapping = hits, id_col = id_col)
}
read_model_metadata <-
function (path) 
{
    data <- data.table::fread(path, data.table = FALSE, check.names = FALSE, 
        showProgress = FALSE)
    ids <- data[[depmap_id_column(names(data))]]
    if (anyNA(ids) || anyDuplicated(ids)) 
        stop("DepMap metadata IDs must be present and unique", 
            call. = FALSE)
    data
}
read_server_config <-
function (path) 
{
    if (!file.exists(path)) 
        stop("Configuration file not found: ", path, call. = FALSE)
    cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
    cfg$cache_dir <- path.expand(Sys.getenv("R2_DEPMAP_CACHE_DIR", 
        cfg$cache_dir))
    if (!grepl("^(/|[A-Za-z]:|\\\\)", cfg$cache_dir)) 
        cfg$cache_dir <- file.path(dirname(path), "..", cfg$cache_dir)
    cfg$r2_base_url <- Sys.getenv("R2_BASE_URL", cfg$r2_base_url)
    cfg$r2_dataset_table <- Sys.getenv("R2_DATASET_TABLE", cfg$r2_dataset_table)
    dir.create(cfg$cache_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(cfg$cache_dir)) 
        stop("Cannot create cache directory: ", cfg$cache_dir, 
            call. = FALSE)
    cfg$cache_dir <- normalizePath(cfg$cache_dir, winslash = "/", 
        mustWork = TRUE)
    cfg
}
read_study_table <-
function (path) 
{
    data.table::fread(normalize_local_path(path), data.table = FALSE, 
        check.names = FALSE, showProgress = FALSE, na.strings = c("", 
            "NA", "NaN"))
}
register_cohort <-
function (cfg, expression_path, metadata_path, disease, label = "Local cohort", 
    scale = "log2", gene_column = "gene", sample_column = "sample_id", 
    subgroup_column = "subgroup", source = "user") 
{
    if (!scale %in% c("log2", "tpm", "counts")) 
        stop("scale must be log2, tpm or counts", call. = FALSE)
    e <- read_study_table(expression_path)
    m <- read_study_table(metadata_path)
    require_columns(e, gene_column)
    require_columns(m, c(sample_column, subgroup_column))
    genes <- study_symbols(e[[gene_column]])
    if (anyDuplicated(genes)) 
        stop("Duplicate gene symbols; collapse or resolve mappings before import", 
            call. = FALSE)
    samples <- as.character(m[[sample_column]])
    if (anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) 
        stop("Sample IDs must be non-empty and unique", call. = FALSE)
    if (anyNA(m[[subgroup_column]]) || any(!nzchar(as.character(m[[subgroup_column]])))) 
        stop("Subgroup labels must not be missing", call. = FALSE)
    expression_samples <- setdiff(names(e), gene_column)
    if (!setequal(samples, expression_samples)) 
        stop("Expression columns and metadata sample IDs must match exactly", 
            call. = FALSE)
    mat <- vapply(samples, function(s) numeric_column(e[[s]], 
        s), numeric(nrow(e)))
    dimnames(mat) <- list(genes, samples)
    if (nrow(mat) < 2 || ncol(mat) < 3) 
        stop("Provide at least two genes and three biological samples", 
            call. = FALSE)
    if (scale != "log2" && any(mat < 0)) 
        stop("Counts/TPM cannot be negative", call. = FALSE)
    if (scale == "counts" && (any(abs(mat - round(mat)) > 9.9999999999999995e-08) || 
        any(colSums(mat) == 0))) 
        stop("Counts must be non-negative integers with non-empty libraries", 
            call. = FALSE)
    if (scale == "tpm") 
        mat <- log2(mat + 1)
    m$sample_id <- samples
    m$subgroup <- as.character(m[[subgroup_column]])
    value <- list(label = label, disease = disease, expression = mat, 
        metadata = m, scale = scale, provenance = list(source = source, 
            expression = file_provenance(expression_path), metadata = file_provenance(metadata_path), 
            scale = scale, transformation = if (scale == "tpm") "log2(TPM + 1)" else "none", 
            imported_at = utc_now()))
    value$id <- save_study_object(cfg, "cohort", value)
    study_overview(value)
}
require_columns <-
function (data, columns) 
{
    missing <- setdiff(columns, names(data))
    if (length(missing)) 
        stop("Missing columns: ", paste(missing, collapse = ", "), 
            call. = FALSE)
    if (anyDuplicated(names(data))) 
        stop("Duplicate column names are not supported", call. = FALSE)
}
require_runtime_packages <-
function () 
{
    required <- c("jsonlite", "curl", "xml2", "rvest", "data.table")
    missing <- required[!vapply(required, requireNamespace, logical(1), 
        quietly = TRUE)]
    if (length(missing)) {
        stop("Missing R packages: ", paste(missing, collapse = ", "), 
            ". Run scripts/install-dependencies.R.", call. = FALSE)
    }
}
resolve_pubmed_aliases <-
function (gene, fetch = ncbi_get) 
{
    raw_search <- fetch("esearch.fcgi", list(db = "gene", term = paste0(gene, 
        "[Gene Name] AND Homo sapiens[Organism]"), retmode = "json", 
        retmax = 20))
    search <- r2_json(raw_search, "NCBI Gene search")
    ids <- search$esearchresult$idlist
    if (!length(ids)) 
        return(list(aliases = character(), status = "No matching NCBI human Gene record; using supplied symbol", 
            raw_search = raw_search))
    raw_summary <- fetch("esummary.fcgi", list(db = "gene", id = paste(ids, 
        collapse = ","), retmode = "json"))
    summary <- r2_json(raw_summary, "NCBI Gene summary")$result
    hits <- lapply(as.character(ids), function(id) summary[[id]])
    hits <- Filter(function(x) !is.null(x$name) && toupper(x$name) == 
        gene && as.character(x$organism$taxid) == "9606", hits)
    if (length(hits) != 1L) 
        return(list(aliases = character(), status = "Official symbol could not be uniquely verified; using supplied symbol", 
            raw_search = raw_search, raw_summary = raw_summary))
    hit <- hits[[1]]
    aliases <- trimws(unlist(strsplit(hit$otheraliases %||% "", 
        ",", fixed = TRUE)))
    aliases <- aliases[nzchar(aliases) & grepl("^[A-Za-z0-9][A-Za-z0-9 ._()-]{0,100}$", 
        aliases)]
    list(aliases = unique(aliases), gene_id = hit$uid, official_symbol = hit$name, 
        description = hit$description, status = "Human gene symbol verified against NCBI Gene; documented aliases included", 
        raw_search = raw_search, raw_summary = raw_summary)
}
safe_zscore <-
function (x) 
{
    x <- as.numeric(x)
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) 
        return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE))/s
}
save_profile_tables <-
function (x, path) 
{
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (is.data.frame(x)) {
        if (ncol(x)) 
            data.table::fwrite(x, file.path(path, "data.csv"), 
                na = "NA")
        return(invisible(NULL))
    }
    if (!is.list(x) || inherits(x, "ggplot")) 
        return(invisible(NULL))
    for (n in names(x)) if (n != "plot" && n != "mb_plot") {
        child <- file.path(path, gsub("[^a-zA-Z0-9_-]", "_", 
            n))
        if (is.data.frame(x[[n]]) || is.list(x[[n]])) 
            save_profile_tables(x[[n]], child)
    }
}
save_study_object <-
function (cfg, kind, value) 
{
    id <- paste0(kind, "_", object_digest(value))
    value$id <- id
    cache_put(cfg, id, value)
    id
}
study_design <-
function (metadata, covariates = character()) 
{
    covariates <- unlist(covariates, use.names = FALSE)
    require_columns(metadata, c("subgroup", covariates))
    design_data <- data.frame(group = factor(metadata$subgroup))
    if (length(covariates)) 
        for (i in seq_along(covariates)) design_data[[paste0("cov", 
            i)]] <- metadata[[covariates[i]]]
    if (anyNA(design_data)) 
        stop("Subgroup/covariate values must be complete for differential expression", 
            call. = FALSE)
    design <- stats::model.matrix(~0 + ., design_data)
    if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) 
        stop("Design is confounded or lacks residual degrees of freedom", 
            call. = FALSE)
    list(matrix = design, groups = levels(design_data$group))
}
study_overview <-
function (x) 
{
    list(id = x$id, disease = x$disease, label = x$label, target = x$target, 
        comparators = as.list(x$comparators), genes_n = if (!is.null(x$expression)) nrow(x$expression) else length(x$universe), 
        samples_n = if (!is.null(x$metadata)) nrow(x$metadata) else NULL, 
        subgroups = if (!is.null(x$metadata)) as.list(table(x$metadata$subgroup)) else NULL, 
        provenance = x$provenance)
}
study_symbols <-
function (x) 
{
    out <- toupper(trimws(as.character(x)))
    if (anyNA(out) || any(!nzchar(out)) || any(grepl("[;|,\\s]", 
        out, perl = TRUE))) 
        stop("Use one non-empty human gene symbol per row; resolve probes, aliases and multiple mappings before import", 
            call. = FALSE)
    if (any(grepl("^ENS[A-Z]*G[0-9]|^[0-9]+$", out))) 
        stop("Map Ensembl/Entrez identifiers to human gene symbols before import", 
            call. = FALSE)
    out
}
summarize_depmap_for_genes <-
function (genes, model_csv, gene_effect_csv, model_regex, dependency_cutoff = -0.5) 
{
    models <- read_model_metadata(model_csv)
    model_id_col <- depmap_id_column(names(models))
    filtered <- filter_depmap_models(models, model_regex)
    effect <- read_gene_effect_columns(gene_effect_csv, genes)
    selected_ids <- as.character(filtered$data[[filtered$id_col]])
    effect_ids <- as.character(effect$data[[effect$id_col]])
    relevant_index <- which(effect_ids %in% selected_ids)
    background_index <- which(!effect_ids %in% selected_ids)
    if (!length(relevant_index)) 
        stop("Selected DepMap models were absent from the gene-effect file", 
            call. = FALSE)
    gene_columns <- effect$mapping$column
    relevant <- as.matrix(effect$data[relevant_index, gene_columns, 
        drop = FALSE])
    storage.mode(relevant) <- "numeric"
    background <- as.matrix(effect$data[background_index, gene_columns, 
        drop = FALSE])
    storage.mode(background) <- "numeric"
    relevant_mean <- column_means(relevant)
    relevant_median <- column_medians(relevant)
    relevant_min <- apply(relevant, 2L, min, na.rm = TRUE)
    relevant_min[!is.finite(relevant_min)] <- NA_real_
    fraction_dependent <- colMeans(relevant < dependency_cutoff, 
        na.rm = TRUE)
    fraction_dependent[!is.finite(fraction_dependent)] <- NA_real_
    background_mean <- if (nrow(background)) 
        column_means(background)
    else rep(NA_real_, ncol(relevant))
    data.frame(gene = effect$mapping$gene, depmap_column = effect$mapping$column, 
        relevant_mean_gene_effect = relevant_mean, relevant_median_gene_effect = relevant_median, 
        relevant_min_gene_effect = relevant_min, relevant_fraction_below_cutoff = fraction_dependent, 
        background_mean_gene_effect = background_mean, dependency_selectivity = background_mean - 
            relevant_mean, relevant_models_n = colSums(is.finite(relevant)), 
        stringsAsFactors = FALSE)
}
summarize_specificity <-
function (data, comparators, fdr = 0.050000000000000003, min_abs_log2fc = 1) 
{
    if (!is.finite(fdr) || fdr <= 0 || fdr > 1 || !is.finite(min_abs_log2fc) || 
        min_abs_log2fc < 0) 
        stop("Require 0 < fdr <= 1 and min_abs_log2fc >= 0", 
            call. = FALSE)
    rows <- lapply(split(data, data$gene), function(d) {
        tested <- is.finite(d$adjusted_p) & is.finite(d$log2fc)
        significant <- tested & d$adjusted_p <= fdr & abs(d$log2fc) >= 
            min_abs_log2fc & d$log2fc != 0
        complete <- setequal(d$comparator[tested], comparators)
        same_direction <- all(d$log2fc[tested] > 0) || all(d$log2fc[tested] < 
            0)
        specific <- complete && all(significant) && same_direction
        data.frame(gene = d$gene[1], comparisons_tested = sum(tested), 
            comparisons_required = length(comparators), comparisons_passed = sum(significant), 
            consistent_direction = same_direction && any(tested), 
            subgroup_specific = specific && length(comparators) > 
                1L, passes_all_comparators = specific, direction = if (!any(tested) || 
                !same_direction) 
                "mixed_or_unknown"
            else if (d$log2fc[which(tested)[1]] > 0) 
                "up"
            else "down", min_abs_log2fc = if (any(tested)) 
                min(abs(d$log2fc[tested]))
            else NA_real_, max_adjusted_p = if (any(tested)) 
                max(d$adjusted_p[tested])
            else NA_real_, stringsAsFactors = FALSE)
    })
    if (!length(rows)) 
        stop("No differential-expression rows available", call. = FALSE)
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result[order(!result$passes_all_comparators, result$max_adjusted_p, 
        -result$min_abs_log2fc), ]
}
utc_now <-
function () 
format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
verify_r2_survival_payload <-
function (html, js, retrieval) 
{
    if (!grepl("proportion *= (nrRemaining - status) / nrRemaining;", 
        js, fixed = TRUE)) 
        stop("R2 event semantics changed; re-audit kaplan.js")
    p <- r2_script_payload(html, "/d3/plots/kaplan.js")$plotData
    if (!grepl("month", p$xLabel, ignore.case = TRUE) || !grepl("overall survival", 
        p$yLabel, ignore.case = TRUE)) 
        stop("Survival endpoint/time unit not confirmed")
    d <- p$data
    require_columns(d, c("id", "status", "xValue"))
    base <- retrieval$data
    keys <- toupper(d$id)
    reference <- toupper(base$sample_id)
    if (anyDuplicated(keys) || anyDuplicated(reference) || any(!keys %in% 
        reference)) 
        stop("Survival IDs are duplicate or unmatched")
    i <- match(keys, reference)
    if (any(!d$status %in% c(0, 1)) || any(d$status != r2_numeric(base$survival_status_raw[i], 
        "dead"))) 
        stop("Survival event flags disagree with source dead metadata")
    complete <- is.finite(base$survival_time_years) & !r2_missing(base$survival_status_raw)
    if (!setequal(reference[complete], keys)) 
        stop("Survival endpoint samples differ from complete annotated cases")
    if (any(!is.finite(d$xValue)) || any(d$xValue < 0) || any(abs(d$xValue/12 - 
        base$survival_time_years[i]) > 0.11)) 
        stop("Survival times disagree beyond annotation rounding")
    list(data = data.frame(sample_id = base$sample_id[i], time = d$xValue/12, 
        event = as.integer(d$status), source_months = d$xValue, 
        source_years_rounded = base$survival_time_years[i]), 
        provenance = list(endpoint = "Overall survival", time = "R2 Kaplan follow-up months / 12, checked against rounded os_(years)", 
            event = "status 1 decreases survival in R2 kaplan.js; 0 censored; per-patient status agrees with dead metadata", 
            identifier_join = "Unique case-normalised GEO accessions", 
            r2_statistics = p$statistics))
}
write_gene_page <-
function (retrieval, validation, publications = NULL, publication_error = NULL, 
    output_path = NULL, extra_html = "") 
{
    gene <- retrieval$provenance$queried_gene
    output_path <- output_path %||% file.path(retrieval$provenance$output_dir, 
        paste0(gene, "_gene_page.html"))
    stats <- validation$statistics
    figure <- if (!is.null(validation$plot)) 
        paste0("<div class='figure' role='img' aria-label='", 
            html_escape(paste(gene, "expression by subgroup")), 
            "'>", plot_inline_svg(validation$plot), "</div>")
    else "<p>The subgroup plot is unavailable.</p>"
    meta <- retrieval$audit
    fields <- do.call(rbind, lapply(meta$fields, function(x) data.frame(field = x$field, 
        samples = x$n, missing = x$missing_n)))
    clinical_note <- paste0("Metastasis coding: ", meta$metastasis$evidence, 
        ". Survival: ", meta$survival$note)
    pub_header <- if (!is.null(publications)) 
        paste0("<p>", length(publications$articles), " latest results shown from ", 
            publications$total_matches, " matches. Scope: <strong>", 
            html_escape(publications$disease_filter %||% "all diseases and biological contexts"), 
            "</strong>.</p>", "<p class='small'>Fetched ", html_escape(publications$provenance$retrieved_at), 
            ". ", html_escape(publications$provenance$sort), 
            ". Journal citation dates and electronic dates can differ.</p><details><summary>Search terms, aliases and limitations</summary><p><code>", 
            html_escape(publications$provenance$query), "</code></p><p>", 
            html_escape(publications$provenance$alias_resolution$status), 
            "</p><p>", html_escape(publications$provenance$limitations), 
            "</p><p>Candidate screening: ", html_escape(json_text(publications$provenance$screening)), 
            "</p></details><label for='paper-search'>Filter these results</label><input id='paper-search' type='search' placeholder='Title, journal or description'><p id='paper-count' class='small' aria-live='polite'></p>")
    else paste0("<p class='notice'>Publications unavailable: ", 
        html_escape(publication_error %||% "Not requested"), 
        ". Refresh the page through the MCP tool to retry.</p>")
    css <- "*{box-sizing:border-box}body{margin:0;background:#f3f5f8;color:#162435;font:16px/1.6 system-ui,-apple-system,Segoe UI,sans-serif}main{max-width:1120px;margin:auto;padding:30px 24px 70px}header{padding:24px 0}h1{font-size:2.7rem;line-height:1.1;margin:8px 0}h2{font-size:1.5rem;margin:0 0 16px}h3{font-size:1.1rem;line-height:1.45;margin:0 0 8px}a{color:#095f75;text-decoration-thickness:1px;text-underline-offset:3px}.eyebrow{font-size:.8rem;letter-spacing:.15em;text-transform:uppercase;color:#466576}nav{display:flex;gap:20px;flex-wrap:wrap;margin:22px 0}.card{background:white;border:1px solid #dee5eb;border-radius:14px;padding:26px;margin:22px 0;box-shadow:0 3px 12px #153b5006}.metrics{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}.metric{background:#eaf2f3;border-radius:10px;padding:15px}.metric strong{display:block;font-size:1.7rem}.small,.paper-meta{font-size:.88rem;color:#526474}.figure svg{display:block;width:100%;height:auto}.table-wrap{overflow:auto;max-height:560px}table{border-collapse:collapse;width:100%;font-size:.9rem}th,td{text-align:left;padding:9px 12px;border-bottom:1px solid #e7edf1;white-space:nowrap}thead{position:sticky;top:0;background:#eef3f6}th{font-weight:600}.paper{padding:24px 0;border-bottom:1px solid #dee5eb}.paper:last-child{border:0}.description{border-left:3px solid #8fbec4;padding-left:16px}.paper-links{font-size:.88rem;overflow-wrap:anywhere}.notice{padding:12px 16px;background:#fff4df;border-left:3px solid #b8811b}details{margin:14px 0}summary{cursor:pointer;color:#095f75;font-weight:600}code{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.86rem}input{width:100%;padding:12px;border:1px solid #bfcdd5;border-radius:7px;margin-top:6px;font:inherit}label{font-size:.9rem;font-weight:600}footer{color:#526474;font-size:.86rem}section{scroll-margin-top:20px}@media(max-width:650px){main{padding:18px 12px}h1{font-size:2rem}.card{padding:18px}.metrics{grid-template-columns:1fr}}"
    html <- paste0("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>", 
        "<meta http-equiv='Content-Security-Policy' content=\"default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'; img-src data:; base-uri 'none'; form-action 'none'\">", 
        "<title>", html_escape(gene), " &mdash; medulloblastoma gene page</title><style>", 
        css, "</style></head><body><main>", "<header><div class='eyebrow'>Medulloblastoma &middot; R2 patient evidence + PubMed</div><h1>", 
        html_escape(gene), " gene overview</h1>", "<p>Patient expression and recent literature, with source data and analysis details.</p><nav><a href='#expression'>Expression</a><a href='#publications'>Publications</a><a href='#metadata'>Metadata</a><a href='#patients'>Patient data</a></nav></header>", 
        "<div class='metrics'><div class='metric'><strong>", 
        nrow(retrieval$data), "</strong>R2 samples</div><div class='metric'><strong>", 
        length(unique(retrieval$data$subgroup[!r2_missing(retrieval$data$subgroup)])), 
        "</strong>annotated subgroups</div><div class='metric'><strong>", 
        if (!is.null(publications)) 
            length(publications$articles)
        else "&mdash;", "</strong>recent publications</div></div>", 
        "<section class='card' id='expression'><h2>Expression across broad subgroups</h2>", 
        figure, "<p>One point per patient; original R2 subgroup labels. The overall test compares distributions and does not establish subgroup specificity.</p>", 
        html_table(validation$summaries), "<details><summary>Statistical test and limitations</summary>", 
        html_table(stats), "<p>Predefined unpaired Kruskal&ndash;Wallis test. No pairwise testing, clinical analysis or adjustment for batch/age was performed; unequal distribution shapes prevent a pure median interpretation.</p></details></section>", 
        extra_html, "<section class='card' id='publications'><h2>Recent PubMed publications</h2>", 
        pub_header, if (!is.null(publications)) 
            paste0("<div id='papers'>", publication_cards(publications), 
                "</div>"), "</section>", "<section class='card' id='metadata'><h2>Metadata and source provenance</h2><p>", 
        html_escape(retrieval$provenance$dataset_label), "<br>Reporter ", 
        html_escape(retrieval$provenance$reporter), " &middot; Retrieved ", 
        html_escape(retrieval$provenance$retrieved_at), "</p><p>", 
        html_escape(clinical_note), "</p>", "<details><summary>Metadata availability and original subtype labels</summary>", 
        html_table(fields), html_table(meta$subtype_by_subgroup), 
        "</details>", "<details><summary>Identifiers and expression transformation</summary><p>", 
        html_escape(meta$identifier_note), "</p><p>", html_escape(retrieval$provenance$transformation), 
        "</p><p>", html_escape(retrieval$provenance$sample_join), 
        "</p></details></section>", "<section class='card' id='patients'><h2>Patient-level data</h2><details><summary>Inspect all ", 
        nrow(retrieval$data), " samples</summary>", "<label for='patient-search'>Filter samples or subgroup labels</label><input id='patient-search' type='search' placeholder='Sample ID, subgroup or subtype'>", 
        html_table(retrieval$data[, intersect(c("sample_id", 
            "expression", "subgroup", "subtype", "metastasis_raw", 
            "metastasis_label", "survival_time_years", "survival_status_raw"), 
            names(retrieval$data)), drop = FALSE], "patient-table"), 
        "</details></section><footer>This local HTML file includes its plot and data; no web server or external scripts are required. PubMed and DOI links open online. Original source files and R code remain alongside the analysis outputs.</footer></main>", 
        "<script>const p=document.getElementById('paper-search');if(p){const rows=[...document.querySelectorAll('.paper')];const f=()=>{let n=0;for(const r of rows){r.hidden=!r.textContent.toLowerCase().includes(p.value.toLowerCase());if(!r.hidden)n++;}document.getElementById('paper-count').textContent=n+' publications visible';};p.addEventListener('input',f);f();}const s=document.getElementById('patient-search');if(s)s.addEventListener('input',()=>{for(const r of document.querySelectorAll('#patient-table tbody tr'))r.hidden=!r.textContent.toLowerCase().includes(s.value.toLowerCase());});</script></body></html>")
    writeLines(enc2utf8(html), output_path, useBytes = TRUE)
    normalizePath(output_path, winslash = "/", mustWork = TRUE)
}
xml_value <-
function (node, xpath) 
{
    n <- xml2::xml_find_first(node, xpath)
    if (inherits(n, "xml_missing")) 
        return(NA_character_)
    value <- trimws(gsub("[\\s\\x{00A0}]+", " ", xml2::xml_text(n), 
        perl = TRUE))
    if (nzchar(value)) 
        value
    else NA_character_
}
