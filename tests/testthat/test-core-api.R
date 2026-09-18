isolated_core <- function() {
  root <- Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR")
  bootstrap <- new.env(parent = globalenv())
  sys.source(file.path(root, "R", "load_core.R"), bootstrap)
  # Base plus standard attached R packages, with no test/global function fallback.
  bootstrap$load_r2_core(root, new.env(parent = as.environment("package:stats")))
}

testthat::test_that("core analysis runs without MCP, network or artifact writes", {
  core <- isolated_core()
  testthat::expect_false(exists("dispatch_v1_tool", core, inherits = FALSE))
  testthat::expect_false(exists("handle_request", core, inherits = FALSE))
  r <- v1_assembled_fixture()
  r$provenance$output_dir <- tempfile("must-not-write-")
  core$get_r2_expression <- function(...) stop("Unexpected network access")
  x <- core$analyze_gene_sources(list(retrieval = r))
  testthat::expect_false(dir.exists(r$provenance$output_dir))
  expected <- stats::kruskal.test(expression ~ subgroup, r$data)
  testthat::expect_equal(x$validation$statistics$p_value, expected$p.value)
  testthat::expect_s3_class(x$validation$plot, "ggplot")
  testthat::expect_equal(x$profile$subtypes$statistics$n, 12)
  testthat::expect_true(all(x$profile$status$status[-1] == "unavailable"))
  testthat::expect_true(all(vapply(x$profile$clinical,
    function(panel) identical(panel$survival$status, "unavailable"), logical(1))))
  testthat::expect_error(core$analyze_gene_sources(list(retrieval = r), "best"), "arg")
})

testthat::test_that("export and offline reproduction work from an isolated core", {
  core <- isolated_core()
  r <- v1_assembled_fixture()
  x <- core$analyze_gene_sources(list(retrieval = r))
  cfg <- list(output_dir = tempfile("core-report-"))
  page <- core$export_gene_report(x, cfg)
  testthat::expect_true(file.exists(page$html_path))
  saved <- readRDS(page$report_data)
  testthat::expect_equal(saved$validation$statistics, x$validation$statistics)
  offline <- new.env(parent = as.environment("package:stats"))
  sys.source(file.path(dirname(page$html_path), "analysis_functions.R"), offline)
  testthat::expect_false(exists("dispatch_v1_tool", offline, inherits = FALSE))
  regenerated <- file.path(dirname(page$html_path), "regenerated.html")
  offline$write_gene_page(saved$retrieval, saved$validation, saved$publications,
    saved$publication_error, regenerated, offline$profile_extra_html(saved$profile, saved$cutoff))
  testthat::expect_true(file.exists(regenerated))
})

testthat::test_that("MCP and direct R share report service and statistical results", {
  core <- isolated_core()
  r <- v1_assembled_fixture()
  # Inject source fixtures at the I/O boundary, retaining the real analysis/export.
  core$collect_gene_sources <- function(...) list(retrieval = r)
  adapter <- new.env(parent = core)
  sys.source(file.path(Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR"), "R", "v1_tools.R"), adapter)
  cfg <- list(output_dir = tempfile("adapter-report-"))
  direct <- core$build_gene_page(cfg, "HLX")
  via_mcp <- adapter$dispatch_v1_tool("build_gene_page", list(gene = "HLX"), cfg)
  testthat::expect_equal(via_mcp$modules, direct$modules)
  testthat::expect_equal(readRDS(via_mcp$report_data)$validation$statistics,
    readRDS(direct$report_data)$validation$statistics)
  testthat::expect_equal(via_mcp$available_cutoffs, direct$available_cutoffs)
})
