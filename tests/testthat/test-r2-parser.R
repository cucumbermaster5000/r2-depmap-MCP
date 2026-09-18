testthat::test_that("R2 differential-expression direction is normalized", {
  plugin_dir <- Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR")
  html <- paste(readLines(file.path(plugin_dir, "tests", "fixtures", "r2_diff.html")), collapse = "\n")
  parsed <- parse_r2_diff_html(html, "group3", "group4", "subgroup")

  testthat::expect_equal(nrow(parsed$data), 2L)
  testthat::expect_equal(parsed$data$patient_effect_group1_minus_group2[parsed$data$gene == "HLX"], 1.5)
  testthat::expect_true(parsed$data$group1_higher[parsed$data$gene == "HLX"])
  testthat::expect_equal(parsed$data$patient_effect_group1_minus_group2[parsed$data$gene == "RBM24"], -2.9)
  testthat::expect_equal(parsed$details$group_1_n, 144L)
  testthat::expect_equal(parsed$details$group_2_n, 326L)
})
