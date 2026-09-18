#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
file_arg <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- normalizePath(file.path(dirname(file_arg[1]), ".."), winslash = "/")
lib <- file.path(root, ".R-library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(timeout = 300)
required <- c("jsonlite", "curl", "xml2", "rvest", "data.table", "survival", "ggplot2")
if ("--tests" %in% args) required <- c(required, "testthat")
if ("--shiny" %in% args) required <- c(required, "shiny", "bslib")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
}

if ("--analysis" %in% args) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
  analysis <- c("limma", "edgeR")
  missing_analysis <- analysis[!vapply(analysis, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_analysis)) BiocManager::install(missing_analysis, lib = lib, ask = FALSE, update = FALSE)
  required <- c(required, analysis)
}
if (any(!vapply(required, requireNamespace, logical(1), quietly = TRUE))) stop("Some dependencies could not be installed")

message("Installed MCP runtime dependencies: ", paste(required, collapse = ", "))
message("Project library: ", lib)
message("Use --analysis for raw matrix analysis (limma/edgeR); --tests for testthat; --shiny for the web app.")
