subtype_fixture <- function() {
  r <- v1_assembled_fixture()
  r$data <- r$data[rep(seq_len(12),each=3),]
  r$metadata <- r$metadata[rep(seq_len(12),each=3),]
  r$data$sample_id <- r$metadata$samplenames <- paste0("P",seq_len(36))
  r
}

testthat::test_that("subtypes preserve metadata and reproduce rank statistics and BH family", {
  r <- subtype_fixture()
  a <- analyze_r2_subtypes(r)
  testthat::expect_equal(a$status,"complete")
  testthat::expect_equal(nrow(a$data),36)
  testthat::expect_true(all(names(r$metadata) %in% names(a$data)))
  testthat::expect_equal(a$data$samplenames,r$data$sample_id)
  testthat::expect_equal(a$overall$p_value,kruskal.test(r$data$expression,r$metadata$subtype)$p.value)
  testthat::expect_equal(nrow(a$pairwise),6)
  testthat::expect_equal(a$pairwise$adjusted_p_value,p.adjust(a$pairwise$p_value,"BH"))
  testthat::expect_equal(a$pairwise$rank_biserial,rep(-1,6))
  testthat::expect_match(a$metadata$classification,"Unidentified")
  testthat::expect_s3_class(a$plot,"ggplot")
  # Joining is by ID, not input row position.
  r$metadata <- r$metadata[36:1,]
  testthat::expect_equal(analyze_r2_subtypes(r)$summaries,a$summaries)
})

testthat::test_that("missing, tiny, single and constant subtypes are explicit", {
  r <- subtype_fixture()
  r$metadata$subtype[1] <- "na"; r$data$expression[2] <- NA_real_
  a <- analyze_r2_subtypes(r)
  testthat::expect_equal(a$missing$excluded_n,2)
  testthat::expect_equal(a$metadata$missing_n,1)
  testthat::expect_equal(nrow(a$data),36)
  testthat::expect_equal(sum(a$summaries$n),34)
  r$metadata$subtype <- NULL
  a <- analyze_r2_subtypes(r)
  testthat::expect_false(a$metadata$available)
  testthat::expect_equal(a$overall$n,0)
  r <- subtype_fixture(); r$metadata$subtype <- "one"
  testthat::expect_equal(analyze_r2_subtypes(r)$status,"descriptive_only")
  r <- subtype_fixture(); r$data$expression <- 1
  testthat::expect_match(analyze_r2_subtypes(r)$overall$reason,"constant")
  testthat::expect_match(analyze_r2_subtypes(v1_assembled_fixture())$overall$reason,"fewer than five")
  r$data$sample_id[1] <- r$data$sample_id[2]
  testthat::expect_error(analyze_r2_subtypes(r),"identifiers")
})

testthat::test_that("subtype failures cannot invalidate the Version 1 result", {
  source(file.path(Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR"),"shiny","interface.R"),local=TRUE)
  core <- new.env(); core$r2_missing <- r2_missing
  core$get_r2_expression <- function(...) subtype_fixture()
  core$analyze_r2_subgroups <- analyze_r2_subgroups
  core$analyze_r2_subtypes <- function(...) stop("secret technical traceback")
  shiny::testServer(explorer_server(core,list()), {
    session$setInputs(gene="HLX"); session$setInputs(run=1)
    testthat::expect_null(failure())
    testthat::expect_equal(result()$analysis$statistics$n,36)
    testthat::expect_null(result()$subtypes)
    testthat::expect_false(grepl("secret",result()$subtype_failure))
  })
})
