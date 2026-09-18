# NCBI E-utilities, public read-only requests. No expression/clinical data sent.
ncbi_get <- local({
  last_request <- as.POSIXct("1970-01-01", tz = "UTC")
  function(endpoint, params) {
    if (!endpoint %in% c("esearch.fcgi", "esummary.fcgi", "efetch.fcgi")) stop("Unsupported NCBI endpoint")
    params$tool <- "r2_mb_gene_profile"
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", endpoint, "?", encode_form(params))
    for (attempt in 1:3) {
      delay <- 0.36 - as.numeric(difftime(Sys.time(), last_request, units = "secs"))
      if (delay > 0) Sys.sleep(delay)
      last_request <<- Sys.time()
      response <- curl::curl_fetch_memory(url, curl::new_handle(timeout = 60, connecttimeout = 15, followlocation = TRUE,
        useragent = "r2-mb-gene-profile/0.3 (+NCBI E-utilities)"))
      if (response$status_code == 200L) return(rawToChar(response$content))
      if (!response$status_code %in% c(429L, 500L, 502L, 503L, 504L) || attempt == 3L)
        stop("NCBI returned HTTP ", response$status_code, call. = FALSE)
      Sys.sleep(attempt)
    }
  }
})

pubmed_gene_terms <- function(gene, aliases = character()) {
  if (length(gene) != 1L || !grepl("^[A-Za-z][A-Za-z0-9._-]*$", gene)) stop("Provide one gene symbol", call. = FALSE)
  aliases <- unlist(aliases, use.names = FALSE)
  terms <- unique(c(toupper(gene), aliases))
  if (anyNA(terms) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9 ._()-]{0,100}$", terms))) stop("Invalid gene alias", call. = FALSE)
  terms
}

pubmed_query <- function(terms, disease = NULL) {
  q <- paste0("(", paste0('"', terms, '"[Title/Abstract]', collapse = " OR "), ")")
  if (!is.null(disease)) {
    if (length(disease) != 1 || !grepl("^[A-Za-z0-9][A-Za-z0-9 -]{0,100}$", disease)) stop("Invalid disease text", call. = FALSE)
    q <- paste0(q, ' AND "', disease, '"[Title/Abstract]')
  }
  q
}

resolve_pubmed_aliases <- function(gene, fetch = ncbi_get) {
  raw_search <- fetch("esearch.fcgi", list(db = "gene", term = paste0(gene, "[Gene Name] AND Homo sapiens[Organism]"), retmode = "json", retmax = 20))
  search <- r2_json(raw_search, "NCBI Gene search")
  ids <- search$esearchresult$idlist
  if (!length(ids)) return(list(aliases = character(), status = "No matching NCBI human Gene record; using supplied symbol", raw_search = raw_search))
  raw_summary <- fetch("esummary.fcgi", list(db = "gene", id = paste(ids, collapse = ","), retmode = "json"))
  summary <- r2_json(raw_summary, "NCBI Gene summary")$result
  hits <- lapply(as.character(ids), function(id) summary[[id]])
  hits <- Filter(function(x) !is.null(x$name) && toupper(x$name) == gene && as.character(x$organism$taxid) == "9606", hits)
  if (length(hits) != 1L) return(list(aliases = character(), status = "Official symbol could not be uniquely verified; using supplied symbol", raw_search = raw_search, raw_summary = raw_summary))
  hit <- hits[[1]]
  aliases <- trimws(unlist(strsplit(hit$otheraliases %||% "", ",", fixed = TRUE)))
  aliases <- aliases[nzchar(aliases) & grepl("^[A-Za-z0-9][A-Za-z0-9 ._()-]{0,100}$", aliases)]
  list(aliases = unique(aliases), gene_id = hit$uid, official_symbol = hit$name, description = hit$description,
    status = "Human gene symbol verified against NCBI Gene; documented aliases included", raw_search = raw_search, raw_summary = raw_summary)
}

xml_value <- function(node, xpath) {
  n <- xml2::xml_find_first(node, xpath)
  if (inherits(n, "xml_missing")) return(NA_character_)
  value <- trimws(gsub("[\\s\\x{00A0}]+", " ", xml2::xml_text(n), perl = TRUE))
  if (nzchar(value)) value else NA_character_
}

pubmed_date <- function(node) {
  if (inherits(node, "xml_missing")) return(NA_character_)
  medline <- xml_value(node, "./MedlineDate")
  if (!is.na(medline)) return(medline)
  parts <- vapply(c("Year", "Month", "Day", "Season"), function(x) xml_value(node, paste0("./", x)), character(1))
  if (!any(!is.na(parts))) return(NA_character_)
  paste(parts[!is.na(parts)], collapse = " ")
}

pubmed_mentions <- function(text, terms) {
  text <- gsub("(*UTF)[\\x{2010}-\\x{2015}\\x{2212}]", "-", enc2utf8(text), perl = TRUE)
  any(vapply(terms, function(term) grepl(paste0("(?i)(?<![A-Za-z0-9_-])\\Q", term, "\\E(?![A-Za-z0-9_-])"), text, perl = TRUE), logical(1)))
}

pubmed_description <- function(article, terms, title) {
  blocks <- xml2::xml_find_all(article, "./MedlineCitation/Article/Abstract/AbstractText")
  candidates <- character(); scores <- numeric()
  for (b in blocks) {
    text <- trimws(gsub("[\\s\\x{00A0}]+", " ", xml2::xml_text(b), perl = TRUE))
    # An extractive summary: never synthesizes a new biological claim.
    sentences <- unlist(strsplit(text, "(?<=[.!?])\\s+(?=[A-Z0-9])", perl = TRUE))
    sentences <- sentences[nzchar(sentences)]
    label <- paste(xml2::xml_attr(b, "Label"), xml2::xml_attr(b, "NlmCategory"))
    matches <- vapply(sentences, pubmed_mentions, logical(1), terms = terms)
    # Prefer explicit gene findings, then conclusions/results, preserving sentence text.
    score <- 10 * matches + 3 * grepl("conclu", label, ignore.case = TRUE) + 2 * grepl("result", label, ignore.case = TRUE)
    candidates <- c(candidates, sentences); scores <- c(scores, score)
  }
  if (!length(candidates)) return(list(text = paste0("No abstract is available in PubMed; the publication is titled \u201c", sub("[.]$", "", title), "\u201d."), source = "title_only_no_abstract"))
  list(text = candidates[which.max(scores)], source = "abstract_sentence_excerpt")
}

parse_pubmed_articles <- function(xml, terms) {
  doc <- xml2::read_xml(charToRaw(xml), options = "NONET")
  errors <- xml2::xml_find_all(doc, "//ERROR")
  if (length(errors)) stop("NCBI EFetch error: ", paste(xml2::xml_text(errors), collapse = "; "))
  if (!identical(xml2::xml_name(doc), "PubmedArticleSet")) stop("Unexpected PubMed XML root", call. = FALSE)
  articles <- xml2::xml_find_all(doc, "/PubmedArticleSet/PubmedArticle")
  lapply(articles, function(a) {
    pmid <- xml_value(a, "./MedlineCitation/PMID")
    title <- xml_value(a, "./MedlineCitation/Article/ArticleTitle")
    if (is.na(pmid) || !grepl("^[0-9]+$", pmid) || is.na(title)) stop("Incomplete PubMed record")
    doi <- xml_value(a, './PubmedData/ArticleIdList/ArticleId[@IdType="doi"]')
    if (is.na(doi)) doi <- xml_value(a, './MedlineCitation/Article/ELocationID[@EIdType="doi"]')
    if (!is.na(doi) && !grepl("^10\\.[0-9]{4,9}/[^[:space:]]+$", doi)) doi <- NA_character_
    description <- pubmed_description(a, terms, title)
    types <- xml2::xml_text(xml2::xml_find_all(a, "./MedlineCitation/Article/PublicationTypeList/PublicationType"))
    notices <- xml2::xml_find_all(a, './MedlineCitation/CommentsCorrectionsList/CommentsCorrections[@RefType="RetractionIn" or @RefType="ErratumIn" or @RefType="ExpressionOfConcernIn"]')
    list(pmid = pmid, title = title, journal = xml_value(a, "./MedlineCitation/Article/Journal/Title"),
      publication_date = pubmed_date(xml2::xml_find_first(a, "./MedlineCitation/Article/Journal/JournalIssue/PubDate")),
      electronic_date = pubmed_date(xml2::xml_find_first(a, './MedlineCitation/Article/ArticleDate[@DateType="Electronic"]')),
      doi = doi, doi_url = if (!is.na(doi)) paste0("https://doi.org/", utils::URLencode(doi, reserved = TRUE)) else NULL,
      pubmed_url = paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/"), description = description$text, description_source = description$source,
      abstract = paste(xml2::xml_text(xml2::xml_find_all(a, "./MedlineCitation/Article/Abstract/AbstractText")), collapse = "\n"),
      publication_types = as.list(types), notices = lapply(notices, function(n) list(type = xml2::xml_attr(n, "RefType"), pmid = xml_value(n, "./PMID"))),
      retracted = "Retracted Publication" %in% types || any(xml2::xml_attr(notices, "RefType") == "RetractionIn"))
  })
}

get_gene_publications <- function(cfg, gene, limit = 10L, disease = NULL, aliases = character(), include_aliases = TRUE, refresh = FALSE) {
  terms <- pubmed_gene_terms(gene, aliases); gene <- terms[1]
  if (length(limit) != 1L || !is.numeric(limit) || !is.finite(limit) || limit != floor(limit) || limit < 1L || limit > 30L) stop("limit must be 1-30")
  pubmed_query(terms, disease) # validate before network activity
  key <- cache_key("pubmed_v2", list(gene = gene, disease = disease %||% "all", aliases = paste(terms, collapse = "|"), expand = include_aliases, limit = limit))
  cached <- if (!refresh) cache_get(cfg, key) else NULL
  if (!is.null(cached)) {
    age <- as.numeric(difftime(Sys.time(), as.POSIXct(cached$provenance$retrieved_at, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), units = "hours"))
    if (is.finite(age) && age >= 0 && age < 24) { cached$provenance$cache_hit <- TRUE; return(cached) }
  }
  resolved <- if (include_aliases) tryCatch(resolve_pubmed_aliases(gene), error = function(e) list(aliases = character(), status = paste("Alias lookup unavailable:", conditionMessage(e)))) else list(aliases = character(), status = "Automatic aliases disabled")
  terms <- pubmed_gene_terms(gene, unique(c(aliases, resolved$aliases)))
  query <- pubmed_query(terms, disease)
  screening_limit <- min(200L, max(50L, limit * 5L))
  raw_search <- ncbi_get("esearch.fcgi", list(db = "pubmed", term = query, retmode = "json", retmax = screening_limit, sort = "pub_date"))
  search <- r2_json(raw_search, "PubMed search")
  if (!is.null(search$error) || !is.null(search$esearchresult$errorlist)) stop("PubMed rejected the query: ", json_text(search$error %||% search$esearchresult$errorlist))
  if (is.null(search$esearchresult$count) || is.null(search$esearchresult$idlist)) stop("Unexpected PubMed search result")
  ids <- as.character(search$esearchresult$idlist)
  raw_xml <- NULL; articles <- list()
  if (length(ids)) {
    raw_xml <- ncbi_get("efetch.fcgi", list(db = "pubmed", id = paste(ids, collapse = ","), retmode = "xml"))
    articles <- parse_pubmed_articles(raw_xml, terms)
    returned <- vapply(articles, `[[`, character(1), "pmid")
    if (anyDuplicated(returned) || !setequal(ids, returned)) stop("PubMed search/fetch record mismatch; refusing partial results")
    articles <- articles[match(ids, returned)]
  }
  relevant <- vapply(articles, function(a) pubmed_mentions(paste(a$title, a$abstract), terms), logical(1))
  excluded_ids <- vapply(articles[!relevant], `[[`, character(1), "pmid")
  articles <- head(articles[relevant], limit)
  raw_dir <- file.path(cfg$cache_dir, "pubmed_sources", basename(tempfile("retrieval-")))
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  writeBin(charToRaw(raw_search), file.path(raw_dir, "search.json"))
  if (!is.null(raw_xml)) writeBin(charToRaw(raw_xml), file.path(raw_dir, "articles.xml"))
  if (!is.null(resolved$raw_search)) writeBin(charToRaw(resolved$raw_search), file.path(raw_dir, "gene-search.json"))
  if (!is.null(resolved$raw_summary)) writeBin(charToRaw(resolved$raw_summary), file.path(raw_dir, "gene-summary.json"))
  resolved$raw_search <- resolved$raw_summary <- NULL
  result <- list(gene = gene, disease_filter = disease, total_matches = as.integer(search$esearchresult$count), articles = articles,
    provenance = list(source = "NCBI PubMed E-utilities", query = query, query_translation = search$esearchresult$querytranslation,
      searched_terms = as.list(terms), alias_resolution = resolved, sort = "PubMed publication date descending (pub_date), not indexing date",
      date_definition = "Journal citation date displayed verbatim, with electronic publication date separately when provided; partial dates are not invented",
      summary_method = "One verbatim abstract sentence, preferring gene mentions and conclusions/results; title-only fallback when no abstract is supplied",
      screening = list(records_checked = length(ids), excluded_pmids = as.list(excluded_ids), rule = "Require a standalone gene symbol or alias in title/abstract; hyphenated drug codes such as HLX-02 do not qualify", maximum_records = screening_limit),
      limitations = "Keyword/alias matches can include unrelated uses or non-human studies; no implication that a paper studies medulloblastoma unless separately checked",
      search_warnings = search$esearchresult$warninglist, retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), cache_hit = FALSE,
      raw_source_dir = normalizePath(raw_dir, winslash = "/")))
  cache_put(cfg, key, result)
  result
}
