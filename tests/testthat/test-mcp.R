testthat::test_that("stdio MCP negotiates, exposes only Version 1 and retrieves cached patient data", {
  root <- Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR")
  cfg <- test_cfg()
  withr::local_envvar(c(R2_DEPMAP_CACHE_DIR = cfg$cache_dir))
  cached <- v1_assembled_fixture()
  key <- cache_key("r2_patient_v1", list(dataset = "ps_avgpres_gse85217geo763_hugene11t", gene = "HLX", reporter = "unique", transformation = "log2"))
  cache_put(cfg, key, cached)
  requests <- list(
    list(jsonrpc = "2.0", id = 1, method = "initialize", params = list(protocolVersion = "future-version")),
    list(jsonrpc = "2.0", method = "notifications/initialized"),
    list(jsonrpc = "2.0", id = 2, method = "tools/list"),
    list(jsonrpc = "2.0", id = 3, method = "tools/call", params = list(name = "get_r2_expression", arguments = list(gene = "HLX", limit = 2))),
    list(jsonrpc = "2.0", id = 4, method = "tools/call", params = list(name = "prioritize_candidates", arguments = list(analysis_id = "bad"))),
    list(jsonrpc = "2.0", id = 5, method = "ping"))
  input <- tempfile(); output <- tempfile(); errors <- tempfile()
  writeLines(c(vapply(requests, json_text, character(1), pretty = FALSE), "{bad json", "[]"), input)
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(executable, c("--vanilla", shQuote(file.path(root, "R", "server.R"))), stdin = input, stdout = output, stderr = errors)
  testthat::expect_equal(status, 0)
  responses <- lapply(readLines(output), jsonlite::fromJSON, simplifyVector = FALSE)
  testthat::expect_length(responses, 7)
  testthat::expect_equal(responses[[1]]$result$protocolVersion, "2025-06-18")
  testthat::expect_length(responses[[2]]$result$tools, 7)
  testthat::expect_false(responses[[3]]$result$isError)
  testthat::expect_equal(responses[[3]]$result$structuredContent$total_samples, 12)
  testthat::expect_length(responses[[3]]$result$structuredContent$rows, 2)
  testthat::expect_true(responses[[4]]$result$isError)
  testthat::expect_equal(responses[[6]]$error$code, -32700)
  testthat::expect_equal(responses[[7]]$error$code, -32600)
  testthat::expect_match(readLines(output)[5], '"result":\\{\\}')
})
