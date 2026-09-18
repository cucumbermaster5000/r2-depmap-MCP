testthat::test_that("live R2 connector returns the expected columns", {
  testthat::skip_if_not(tolower(Sys.getenv("RUN_LIVE_R2_TESTS", "false")) == "true")
  plugin_dir <- Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR")
  cfg <- read_server_config(file.path(plugin_dir, "config", "defaults.json"))
  result <- r2_compare_groups(
    cfg,
    track = "subgroup",
    group_1 = "group3",
    group_2 = "group4",
    top_n = 10L,
    p_threshold = 1,
    refresh = TRUE
  )
  expected <- c("gene", "adjusted_p", "r2_log2fc", "patient_effect_group1_minus_group2")
  testthat::expect_true(all(expected %in% names(result$data)))
  testthat::expect_equal(result$details$group_1_n, 144L)
  testthat::expect_equal(result$details$group_2_n, 326L)
})
