# Acceptance artifacts using real R2 data (existing cache is reused).
source("R/load_core.R")
core <- load_r2_core(".")
source("shiny/interface.R", encoding = "UTF-8")
cfg <- core$read_server_config("config/defaults.json")
cfg$output_dir <- normalizePath("outputs", winslash = "/")
folder <- "outputs/shiny-v1-acceptance"
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
rows <- lapply(c("HLX", "SLC2A1"), function(gene) {
  r <- core$get_r2_expression(cfg, gene)
  a <- core$analyze_r2_subgroups(r)
  x <- list(gene = gene, retrieval = r, analysis = a)
  for (kind in c("patients", "statistics", "pdf", "png")) {
    extension <- if (kind %in% c("pdf", "png")) kind else "csv"
    path <- file.path(folder, paste0(gene, "_", kind, ".", extension))
    write_explorer_download(x, kind, path)
    stopifnot(file.info(path)$size > 100)
    if (kind == "patients") stopifnot(nrow(data.table::fread(path)) == nrow(r$data))
    if (kind == "statistics") stopifnot(isTRUE(all.equal(data.table::fread(path)$p_value, a$statistics$p_value)))
  }
  data.frame(gene = gene, retrieved = nrow(r$data), analyzed = a$statistics$n,
    H = a$statistics$statistic, p = a$statistics$p_value,
    cache_hit = isTRUE(r$provenance$cache_hit))
})
results <- do.call(rbind, rows)
data.table::fwrite(results, file.path(folder, "acceptance.csv"))
print(results)
