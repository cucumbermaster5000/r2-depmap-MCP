#!/usr/bin/env Rscript
# Read-only source audit. Preserve the exact public HTML and its observed controls.
arg <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- normalizePath(file.path(dirname(arg[1]), ".."), winslash = "/")
for (f in c("common", "r2_client")) source(file.path(root, "R", paste0(f, ".R")))
out <- file.path(root, "outputs", "audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cfg <- jsonlite::fromJSON(file.path(root, "config", "defaults.json"))
requests_path <- file.path(out, "requests.json")
requests <- if (file.exists(requests_path)) jsonlite::fromJSON(requests_path, simplifyVector = FALSE) else list(
  list(name = "dataset", url = paste0(cfg$r2_base_url, "?table=", cfg$r2_dataset_table)))
h <- r2_session(120)
for (request in requests) {
  html <- r2_request(h, request$url %||% cfg$r2_base_url, request$fields)
  writeLines(html, file.path(out, paste0(request$name, ".html")), useBytes = TRUE)
  doc <- xml2::read_html(file.path(out, paste0(request$name, ".html")))
  controls <- lapply(xml2::xml_find_all(doc, ".//select|.//input|.//textarea"), function(node) {
    options <- xml2::xml_find_all(node, ".//option")
    list(tag = xml2::xml_name(node), name = xml2::xml_attr(node, "name"), type = xml2::xml_attr(node, "type"),
      value = xml2::xml_attr(node, "value"), options = lapply(options, function(n) list(value = xml2::xml_attr(n, "value"), label = trimws(xml2::xml_text(n)))))
  })
  links <- xml2::xml_find_all(doc, ".//a[@href]")
  link_table <- data.frame(label = trimws(xml2::xml_text(links)), href = xml2::xml_attr(links, "href"))
  writeLines(json_text(list(controls = controls, links = link_table)), file.path(out, paste0(request$name, "-controls.json")))
  cat(request$name, ": ", nchar(html), " chars; controls: ", paste(vapply(controls, function(x) x$name %||% "", character(1)), collapse = ", "), "\n", sep = "")
}
