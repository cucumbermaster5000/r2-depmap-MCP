html_escape <- function(x) {
  x <- as.character(x); x[is.na(x)] <- "Not reported"
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

html_table <- function(d, id = NULL) {
  if (is.null(d) || !nrow(d)) return("<p>No rows available.</p>")
  display <- d
  for (key in names(display)) if (is.numeric(display[[key]])) display[[key]] <- format(signif(display[[key]], 6), trim = TRUE)
  header <- paste0("<th scope='col'>", html_escape(names(display)), "</th>", collapse = "")
  rows <- apply(display, 1, function(row) paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>"))
  paste0("<div class='table-wrap'><table", if (!is.null(id)) paste0(" id='", html_escape(id), "'"), "><thead><tr>", header,
    "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>")
}

plot_inline_svg <- function(plot) {
  path <- tempfile(fileext = ".svg")
  on.exit(unlink(path))
  grDevices::svg(path, width = 9, height = 6.3, onefile = TRUE)
  tryCatch(print(plot), finally = grDevices::dev.off())
  svg <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr("<svg", svg, fixed = TRUE)[1]
  if (start < 1) stop("Plot rendering did not produce SVG")
  svg <- substring(svg,start)
  prefix <- paste0('svg',basename(tempfile()),'-')
  svg <- gsub('id="',paste0('id="',prefix),svg,fixed=TRUE)
  svg <- gsub('href="#',paste0('href="#',prefix),svg,fixed=TRUE)
  gsub('url(#',paste0('url(#',prefix),svg,fixed=TRUE)
}

publication_cards <- function(publications) {
  if (!length(publications$articles)) return("<p>No publications matched this search. The disease filter, alias choices and search date are shown above.</p>")
  paste(vapply(publications$articles, function(p) {
    notices <- c(if (isTRUE(p$retracted)) "Retracted publication", vapply(p$notices, function(n) n$type, character(1)))
    paste0("<article class='paper'>",
      if (length(notices)) paste0("<p class='notice'>", html_escape(paste(unique(notices), collapse = "; ")), " &mdash; inspect the PubMed record.</p>"),
      "<h3><a href='", html_escape(p$pubmed_url), "' target='_blank' rel='noopener noreferrer'>", html_escape(p$title), "</a></h3>",
      "<p class='paper-meta'><strong>", html_escape(p$journal), "</strong><br>Publication date: ", html_escape(p$publication_date),
      if (!is.na(p$electronic_date)) paste0(" &middot; Electronic publication: ", html_escape(p$electronic_date)), "</p>",
      "<p class='description'>", html_escape(p$description), "</p>",
      "<p class='small'>", if (identical(p$description_source, "abstract_sentence_excerpt")) "One-sentence abstract excerpt; the authors' wording." else "No abstract supplied; title-only description.", "</p>",
      "<p class='paper-links'>PMID ", html_escape(p$pmid), " &middot; DOI: ",
      if (is.null(p$doi_url)) "Not reported in PubMed" else paste0("<a href='", html_escape(p$doi_url), "' target='_blank' rel='noopener noreferrer'>", html_escape(p$doi), "</a>"),
      "</p></article>")
  }, character(1)), collapse = "\n")
}

write_gene_page <- function(retrieval, validation, publications = NULL, publication_error = NULL, output_path = NULL, extra_html = '') {
  gene <- retrieval$provenance$queried_gene
  output_path <- output_path %||% file.path(retrieval$provenance$output_dir, paste0(gene, "_gene_page.html"))
  stats <- validation$statistics
  figure <- if (!is.null(validation$plot)) paste0("<div class='figure' role='img' aria-label='", html_escape(paste(gene, "expression by subgroup")), "'>", plot_inline_svg(validation$plot), "</div>") else "<p>The subgroup plot is unavailable.</p>"
  meta <- retrieval$audit
  fields <- do.call(rbind, lapply(meta$fields, function(x) data.frame(field = x$field, samples = x$n, missing = x$missing_n)))
  clinical_note <- paste0("Metastasis coding: ", meta$metastasis$evidence, ". Survival: ", meta$survival$note)
  pub_header <- if (!is.null(publications)) paste0("<p>", length(publications$articles), " latest results shown from ", publications$total_matches,
    " matches. Scope: <strong>", html_escape(publications$disease_filter %||% "all diseases and biological contexts"), "</strong>.</p>",
    "<p class='small'>Fetched ", html_escape(publications$provenance$retrieved_at), ". ", html_escape(publications$provenance$sort),
    ". Journal citation dates and electronic dates can differ.</p><details><summary>Search terms, aliases and limitations</summary><p><code>",
    html_escape(publications$provenance$query), "</code></p><p>", html_escape(publications$provenance$alias_resolution$status), "</p><p>",
    html_escape(publications$provenance$limitations), "</p><p>Candidate screening: ", html_escape(json_text(publications$provenance$screening)), "</p></details><label for='paper-search'>Filter these results</label><input id='paper-search' type='search' placeholder='Title, journal or description'><p id='paper-count' class='small' aria-live='polite'></p>") else
    paste0("<p class='notice'>Publications unavailable: ", html_escape(publication_error %||% "Not requested"), ". Refresh the page through the MCP tool to retry.</p>")
  css <- "*{box-sizing:border-box}body{margin:0;background:#f3f5f8;color:#162435;font:16px/1.6 system-ui,-apple-system,Segoe UI,sans-serif}main{max-width:1120px;margin:auto;padding:30px 24px 70px}header{padding:24px 0}h1{font-size:2.7rem;line-height:1.1;margin:8px 0}h2{font-size:1.5rem;margin:0 0 16px}h3{font-size:1.1rem;line-height:1.45;margin:0 0 8px}a{color:#095f75;text-decoration-thickness:1px;text-underline-offset:3px}.eyebrow{font-size:.8rem;letter-spacing:.15em;text-transform:uppercase;color:#466576}nav{display:flex;gap:20px;flex-wrap:wrap;margin:22px 0}.card{background:white;border:1px solid #dee5eb;border-radius:14px;padding:26px;margin:22px 0;box-shadow:0 3px 12px #153b5006}.metrics{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}.metric{background:#eaf2f3;border-radius:10px;padding:15px}.metric strong{display:block;font-size:1.7rem}.small,.paper-meta{font-size:.88rem;color:#526474}.figure svg{display:block;width:100%;height:auto}.table-wrap{overflow:auto;max-height:560px}table{border-collapse:collapse;width:100%;font-size:.9rem}th,td{text-align:left;padding:9px 12px;border-bottom:1px solid #e7edf1;white-space:nowrap}thead{position:sticky;top:0;background:#eef3f6}th{font-weight:600}.paper{padding:24px 0;border-bottom:1px solid #dee5eb}.paper:last-child{border:0}.description{border-left:3px solid #8fbec4;padding-left:16px}.paper-links{font-size:.88rem;overflow-wrap:anywhere}.notice{padding:12px 16px;background:#fff4df;border-left:3px solid #b8811b}details{margin:14px 0}summary{cursor:pointer;color:#095f75;font-weight:600}code{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.86rem}input{width:100%;padding:12px;border:1px solid #bfcdd5;border-radius:7px;margin-top:6px;font:inherit}label{font-size:.9rem;font-weight:600}footer{color:#526474;font-size:.86rem}section{scroll-margin-top:20px}@media(max-width:650px){main{padding:18px 12px}h1{font-size:2rem}.card{padding:18px}.metrics{grid-template-columns:1fr}}"
  html <- paste0("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<meta http-equiv='Content-Security-Policy' content=\"default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'; img-src data:; base-uri 'none'; form-action 'none'\">",
    "<title>", html_escape(gene), " &mdash; medulloblastoma gene page</title><style>", css, "</style></head><body><main>",
    "<header><div class='eyebrow'>Medulloblastoma &middot; R2 patient evidence + PubMed</div><h1>", html_escape(gene), " gene overview</h1>",
    "<p>Patient expression and recent literature, with source data and analysis details.</p><nav><a href='#expression'>Expression</a><a href='#publications'>Publications</a><a href='#metadata'>Metadata</a><a href='#patients'>Patient data</a></nav></header>",
    "<div class='metrics'><div class='metric'><strong>", nrow(retrieval$data), "</strong>R2 samples</div><div class='metric'><strong>", length(unique(retrieval$data$subgroup[!r2_missing(retrieval$data$subgroup)])), "</strong>annotated subgroups</div><div class='metric'><strong>", if (!is.null(publications)) length(publications$articles) else "&mdash;", "</strong>recent publications</div></div>",
    "<section class='card' id='expression'><h2>Expression across broad subgroups</h2>", figure,
    "<p>One point per patient; original R2 subgroup labels. The overall test compares distributions and does not establish subgroup specificity.</p>",
    html_table(validation$summaries), "<details><summary>Statistical test and limitations</summary>", html_table(stats),
    "<p>Predefined unpaired Kruskal&ndash;Wallis test. No pairwise testing, clinical analysis or adjustment for batch/age was performed; unequal distribution shapes prevent a pure median interpretation.</p></details></section>",
    extra_html, "<section class='card' id='publications'><h2>Recent PubMed publications</h2>", pub_header,
    if (!is.null(publications)) paste0("<div id='papers'>", publication_cards(publications), "</div>"), "</section>",
    "<section class='card' id='metadata'><h2>Metadata and source provenance</h2><p>", html_escape(retrieval$provenance$dataset_label), "<br>Reporter ",
    html_escape(retrieval$provenance$reporter), " &middot; Retrieved ", html_escape(retrieval$provenance$retrieved_at), "</p><p>", html_escape(clinical_note), "</p>",
    "<details><summary>Metadata availability and original subtype labels</summary>", html_table(fields), html_table(meta$subtype_by_subgroup), "</details>",
    "<details><summary>Identifiers and expression transformation</summary><p>", html_escape(meta$identifier_note), "</p><p>",
    html_escape(retrieval$provenance$transformation), "</p><p>", html_escape(retrieval$provenance$sample_join), "</p></details></section>",
    "<section class='card' id='patients'><h2>Patient-level data</h2><details><summary>Inspect all ", nrow(retrieval$data), " samples</summary>",
    "<label for='patient-search'>Filter samples or subgroup labels</label><input id='patient-search' type='search' placeholder='Sample ID, subgroup or subtype'>",
    html_table(retrieval$data[, intersect(c("sample_id", "expression", "subgroup", "subtype", "metastasis_raw", "metastasis_label", "survival_time_years", "survival_status_raw"), names(retrieval$data)), drop = FALSE], "patient-table"),
    "</details></section><footer>This local HTML file includes its plot and data; no web server or external scripts are required. PubMed and DOI links open online. Original source files and R code remain alongside the analysis outputs.</footer></main>",
    "<script>const p=document.getElementById('paper-search');if(p){const rows=[...document.querySelectorAll('.paper')];const f=()=>{let n=0;for(const r of rows){r.hidden=!r.textContent.toLowerCase().includes(p.value.toLowerCase());if(!r.hidden)n++;}document.getElementById('paper-count').textContent=n+' publications visible';};p.addEventListener('input',f);f();}const s=document.getElementById('patient-search');if(s)s.addEventListener('input',()=>{for(const r of document.querySelectorAll('#patient-table tbody tr'))r.hidden=!r.textContent.toLowerCase().includes(s.value.toLowerCase());});</script></body></html>")
  writeLines(enc2utf8(html), output_path, useBytes = TRUE)
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}

