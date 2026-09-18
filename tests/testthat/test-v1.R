testthat::test_that("R2 patient parser validates dataset, reporter and transformation", {
  f <- v1_fixture()
  e <- parse_r2_patient_expression(f$html, "HLX", "123", "ps_avgpres_gse85217geo763_hugene11t")
  testthat::expect_equal(e$expression, as.numeric(1:12))
  testthat::expect_equal(e$sample_id, paste0("S", 1:12))
  testthat::expect_error(parse_r2_patient_expression(f$html, "MYC", "123", "ps_avgpres_gse85217geo763_hugene11t"), "different")
  testthat::expect_error(parse_r2_patient_expression(f$html, "HLX", "999", "ps_avgpres_gse85217geo763_hugene11t"), "different")
  testthat::expect_error(parse_r2_patient_expression(f$html, "HLX", "123", "other"), "different")
  testthat::expect_error(parse_r2_patient_expression(f$html, "HLX", "123", "ps_avgpres_gse85217geo763_hugene11t", "transform_none"), "transformation")
  testthat::expect_error(r2_script_payload("<html>Login required</html>", "/loadDataTableModal.js"), "exactly one")
  testthat::expect_error(r2_numeric(c("1", "bad"), "test"), "Invalid numeric")
})

testthat::test_that("patient join uses exact R2 identifiers and retains raw annotations", {
  f <- v1_fixture()
  e <- parse_r2_patient_expression(f$html, "HLX", "123", "ps_avgpres_gse85217geo763_hugene11t")
  f$annotations$data <- f$annotations$data[12:1, ]
  r <- assemble_r2_patient_data(e, f$annotations, f$info, f$tracks)
  testthat::expect_identical(r$metadata, f$annotations$data)
  testthat::expect_equal(r$data$subgroup[1], "wnt")
  testthat::expect_true(r$audit$metastasis$coding_verified)
  testthat::expect_equal(r$data$metastasis_raw[3], "na")
  testthat::expect_true(is.na(r$data$metastasis_label[3]))
  testthat::expect_equal(r$data$metastasis_label[1:2], c("M0 (no metastasis)", "Metastatic"))
  testthat::expect_true(all(is.na(r$data$survival_event)))
  testthat::expect_false(r$audit$survival$event_coding_verified)
  testthat::expect_equal(nrow(r$audit$subtype_by_subgroup), 4)
  f$annotations$data$samplenames[1] <- "not-a-match"
  testthat::expect_error(assemble_r2_patient_data(e, f$annotations, f$info, f$tracks), "sample ID sets differ")
})

testthat::test_that("unverified or invalid metastasis coding never produces labels", {
  f <- v1_fixture()
  e <- parse_r2_patient_expression(f$html, "HLX", "123", "ps_avgpres_gse85217geo763_hugene11t")
  r <- assemble_r2_patient_data(e, f$annotations, f$info, list(trackDescriptions = "metastasis: 0, 1"))
  testthat::expect_false(r$audit$metastasis$coding_verified)
  testthat::expect_true(all(is.na(r$data$metastasis_label)))
  f$annotations$data[["met_status_(1_met__0_m0)"]][1] <- "2"
  testthat::expect_error(assemble_r2_patient_data(e, f$annotations, f$info, f$tracks), "Unexpected metastasis")
})

testthat::test_that("Version 1 plot exports statistics, figures and exclusions", {
  testthat::skip_if_not_installed("ggplot2")
  r <- v1_assembled_fixture()
  r$data$expression[1] <- NA_real_
  # After one missing sample, retain enough observations per subgroup for this test.
  r$data <- rbind(r$data, transform(r$data[2, ], sample_id = "additional_sample"))
  out <- tempfile(); dir.create(out)
  x <- plot_r2_subgroup_test(r, out)
  testthat::expect_equal(x$statistics$excluded_n, 1)
  testthat::expect_equal(x$statistics$n, 12)
  testthat::expect_true(all(file.exists(unlist(x$files))))
  testthat::expect_equal(x$statistics$p_value, stats::kruskal.test(expression ~ subgroup, r$data)$p.value)
})

testthat::test_that("Version 1 rejects deferred tools and unsafe gene/path arguments", {
  cfg <- test_cfg()
  testthat::expect_length(v1_tool_definitions(), 7)
  testthat::expect_error(dispatch_v1_tool("prioritize_candidates", list(), cfg), "not enabled")
  testthat::expect_error(dispatch_v1_tool("get_r2_expression", list(gene = "HLX", limit = -1), cfg), "range")
  testthat::expect_error(get_r2_expression(cfg, "../HLX"), "gene symbol")
})
