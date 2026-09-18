testthat::test_that("Shiny only runs on click and delegates analysis to the core", {
  testthat::skip_if_not_installed("shiny")
  source(file.path(Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR"), "shiny", "interface.R"), local = TRUE)
  calls <- 0L
  core <- new.env()
  core$r2_missing <- r2_missing
  core$get_r2_expression <- function(cfg, gene) {
    calls <<- calls + 1L
    if (gene == "ABSENT") stop("Expected one exact reporter; found 0")
    if (gene == "OFFLINE") stop("Connection timed out: technical detail")
    r <- v1_assembled_fixture()
    r$provenance$queried_gene <- gene
    if (gene == "NOANNOTATION") r$data$subgroup <- NA_character_
    r
  }
  core$analyze_r2_subgroups <- analyze_r2_subgroups
  core$analyze_r2_subtypes <- analyze_r2_subtypes
  shiny::testServer(explorer_server(core, list()), {
    session$setInputs(gene = "HLX")
    testthat::expect_equal(calls, 0L)
    session$setInputs(run = 1)
    testthat::expect_equal(calls, 1L)
    testthat::expect_equal(result()$gene, "HLX")
    testthat::expect_equal(result()$analysis$statistics$n, 12)
    session$setInputs(gene = "SLC2A1")
    testthat::expect_equal(calls, 1L)
    testthat::expect_equal(result()$gene, "HLX")
    session$setInputs(run = 2)
    testthat::expect_equal(result()$gene, "SLC2A1")
    for (kind in c("patients", "statistics", "pdf", "png")) {
      path <- tempfile(fileext = paste0(".", if (kind %in% c("pdf", "png")) kind else "csv"))
      write_explorer_download(result(), kind, path)
      testthat::expect_gt(file.info(path)$size, 100)
      if (kind == "patients") testthat::expect_equal(nrow(data.table::fread(path)), 12)
      if (kind == "statistics") testthat::expect_equal(data.table::fread(path)$p_value, result()$analysis$statistics$p_value)
    }
    session$setInputs(gene = "ABSENT", run = 3)
    testthat::expect_null(result())
    testthat::expect_match(failure(), "Gene not found")
    session$setInputs(gene = "OFFLINE", run = 4)
    testthat::expect_identical(failure(), "R2 data retrieval failed. Please try again later.")
    session$setInputs(gene = "NOANNOTATION", run = 5)
    testthat::expect_match(failure(), "MB subgroup annotation is unavailable")
    session$setInputs(gene = "../bad", run = 6)
    testthat::expect_match(failure(), "valid gene symbol")
    session$setInputs(gene = "HLX", run = 7)
    testthat::expect_null(failure())
    testthat::expect_equal(result()$gene, "HLX")
    testthat::expect_false(busy())
  })
})
