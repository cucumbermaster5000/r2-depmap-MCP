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
clean_table_names <-
function (x) 
{
    x <- trimws(gsub("[\r\n]+", " ", x))
    x <- gsub("[^A-Za-z0-9]+", "", x)
    tolower(x)
}
encode_form <-
function (fields) 
{
    encode_one <- function(x) utils::URLencode(as.character(x %||% 
        ""), reserved = TRUE)
    paste0(vapply(names(fields), encode_one, character(1)), "=", 
        vapply(fields, encode_one, character(1)), collapse = "&")
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
r2_session <-
function (timeout_seconds = 240) 
{
    h <- curl::new_handle()
    curl::handle_setopt(h, followlocation = TRUE, timeout = timeout_seconds, 
        connecttimeout = 20, cookiefile = "", useragent = "r2-depmap-mcp/0.1.0 (+read-only research connector)")
    h
}
