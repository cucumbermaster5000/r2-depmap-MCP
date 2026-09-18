#!/usr/bin/env Rscript
file_arg <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- normalizePath(file.path(dirname(file_arg[1]), ".."), winslash = "/")
.libPaths(c(file.path(root, ".R-library"), .libPaths()))
for (f in c("common", "cache", "r2_client", "study", "enrichment")) source(file.path(root, "R", paste0(f, ".R")))
cfg <- read_server_config(file.path(root, "config", "defaults.json"))
catalog <- list_enrichr_libraries(cfg, pattern = "Reactome", refresh = TRUE)
if (!nrow(catalog$libraries)) stop("No Reactome library in Enrichr catalogue")
library <- tail(sort(catalog$libraries$libraryName), 1)
sets <- load_gene_sets(cfg, library, refresh = TRUE)
cat("Enrichr OK: ", library, ", ", length(sets$sets), " terms\n", sep = "")
if ("--r2" %in% commandArgs(trailingOnly = TRUE)) {
  result <- r2_compare_groups(cfg, top_n = 10L, p_threshold = 1, refresh = TRUE)
  cat("R2 OK: ", nrow(result$data), " differential-expression rows\n", sep = "")
}
