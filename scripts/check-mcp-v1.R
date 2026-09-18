#!/usr/bin/env Rscript
# Exercise the actual stdio server with the real cached HLX retrieval.
file_arg <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
root <- normalizePath(file.path(dirname(file_arg[1]), ".."), winslash = "/")
source(file.path(root, "R", "common.R"))
requests <- list(
  list(jsonrpc = "2.0", id = 1, method = "initialize", params = list(protocolVersion = "2025-06-18")),
  list(jsonrpc = "2.0", method = "notifications/initialized"),
  list(jsonrpc = "2.0", id = 2, method = "tools/list"),
  list(jsonrpc = "2.0", id = 3, method = "tools/call", params = list(name = "get_r2_expression", arguments = list(gene = "HLX", limit = 2))),
  list(jsonrpc = "2.0", id = 4, method = "tools/call", params = list(name = "plot_r2_subgroup_test", arguments = list(gene = "HLX"))),
  list(jsonrpc = "2.0", id = 5, method = "tools/call", params = list(name = "build_gene_page", arguments = list(gene = "HLX"))))
out <- file.path(root, "outputs", "audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
input <- file.path(out, "mcp-v1-requests.jsonl")
output <- file.path(out, "mcp-v1-responses.jsonl")
writeLines(vapply(requests, json_text, character(1), pretty = FALSE), input)
executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- system2(executable, c("--vanilla", shQuote(file.path(root, "R", "server.R"))), stdin = input, stdout = output, stderr = file.path(out, "mcp-v1-stderr.txt"))
if (status != 0) stop("Server exited with code ", status)
responses <- lapply(readLines(output), jsonlite::fromJSON, simplifyVector = FALSE)
stopifnot(length(responses) == 5L, length(responses[[2]]$result$tools) == 7L,
  !responses[[3]]$result$isError, responses[[3]]$result$structuredContent$total_samples == 763,
  !responses[[4]]$result$isError, !responses[[5]]$result$isError,
  responses[[5]]$result$structuredContent$publication_status == "ok",
  all(vapply(responses[[5]]$result$structuredContent$modules,function(m)m$status=='ok',logical(1))),
  file.exists(responses[[5]]$result$structuredContent$html_path))
cat("MCP OK: seven tools; 763 real HLX samples; R plot/test and PubMed browser page generated.\n")
