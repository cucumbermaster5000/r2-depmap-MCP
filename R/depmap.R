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
