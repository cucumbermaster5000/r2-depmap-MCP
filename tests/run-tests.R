#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
test_dir <- dirname(normalizePath(file_arg[[1]], mustWork = TRUE))
plugin_dir <- normalizePath(file.path(test_dir, ".."), mustWork = TRUE)

.libPaths(c(file.path(plugin_dir, ".R-library"), .libPaths()))
if (!requireNamespace("testthat", quietly = TRUE)) stop("Run scripts/install-dependencies.R --tests")
source(file.path(plugin_dir, "R", "load_core.R"))
load_r2_core(plugin_dir, globalenv())
for (file in c("tools.R", "study_tools.R", "v1_tools.R")) {
  source(file.path(plugin_dir, "R", file), local = globalenv())
}

require_runtime_packages()
Sys.setenv(R2_DEPMAP_TEST_PLUGIN_DIR = plugin_dir)
testthat::test_dir(file.path(plugin_dir, "tests", "testthat"), reporter = "summary")
