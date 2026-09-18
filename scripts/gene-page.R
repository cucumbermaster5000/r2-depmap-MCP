#!/usr/bin/env Rscript
# In RStudio, source this from the project root; the default gene is HLX.
file_args <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- if (length(file_args)) normalizePath(file.path(dirname(file_args[1]), ".."), winslash = "/") else normalizePath(".", winslash = "/")
source(file.path(root, "R", "load_core.R"))
core <- load_r2_core(root)
cfg <- core$read_server_config(file.path(root, "config", "defaults.json"))
cfg$output_dir <- file.path(root, "outputs")
args <- commandArgs(trailingOnly = TRUE)
gene <- if (length(args) && !startsWith(args[1], "--")) args[1] else "HLX"
page <- core$build_gene_page(cfg, gene, publication_disease = if ("--medulloblastoma" %in% args) "medulloblastoma" else NULL,
  refresh_publications = "--refresh-publications" %in% args)
cat(core$json_text(page), "\n")
