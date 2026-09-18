#!/usr/bin/env Rscript
# Also usable from RStudio: source('scripts/version1.R'). Defaults to HLX.
file_args <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- if (length(file_args)) normalizePath(file.path(dirname(file_args[1]), ".."), winslash = "/") else normalizePath(".", winslash = "/")
.libPaths(c(file.path(root, ".R-library"), .libPaths()))
for (f in c("common", "cache", "r2_client", "r2_expression")) source(file.path(root, "R", paste0(f, ".R")))
cfg <- read_server_config(file.path(root, "config", "defaults.json"))
cfg$output_dir <- file.path(root, "outputs")
arguments <- commandArgs(trailingOnly = TRUE)
gene <- if (length(arguments) && !startsWith(arguments[1], "--")) arguments[1] else "HLX"
retrieval <- get_r2_expression(cfg, gene, refresh = "--refresh" %in% arguments)
validation <- plot_r2_subgroup_test(retrieval)
print(validation$statistics)
print(validation$summaries)
if (interactive()) print(validation$plot)
cat("Version 1 complete; no downstream modules executed.\n", retrieval$provenance$output_dir, "\n")
