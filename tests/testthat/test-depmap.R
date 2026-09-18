testthat::test_that("DepMap gene dependencies are filtered to medulloblastoma", {
  plugin_dir <- Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR")
  fixture <- function(name) file.path(plugin_dir, "tests", "fixtures", name)
  result <- depmap_gene_dependency(
    gene = "HLX",
    model_csv = fixture("Model.csv"),
    gene_effect_csv = fixture("CRISPRGeneEffect.csv"),
    model_regex = "(?i)medulloblastoma",
    release = "fixture"
  )

  testthat::expect_equal(result$summary$matched_models, 2L)
  testthat::expect_equal(result$summary$mean_gene_effect, -0.7)
  testthat::expect_equal(sort(result$models$model_name), c("D283", "D341"))
})

testthat::test_that("DepMap summaries calculate selectivity", {
  plugin_dir <- Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR")
  fixture <- function(name) file.path(plugin_dir, "tests", "fixtures", name)
  result <- summarize_depmap_for_genes(
    genes = c("HLX", "SLC16A1"),
    model_csv = fixture("Model.csv"),
    gene_effect_csv = fixture("CRISPRGeneEffect.csv"),
    model_regex = "(?i)medulloblastoma"
  )

  hlx <- result[result$gene == "HLX", ]
  testthat::expect_equal(hlx$relevant_mean_gene_effect, -0.7)
  testthat::expect_equal(hlx$background_mean_gene_effect, -0.1)
  testthat::expect_equal(hlx$dependency_selectivity, 0.6)
})
