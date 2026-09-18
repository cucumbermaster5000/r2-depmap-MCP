test_cfg <- function() {
  path <- tempfile("r2-study-")
  dir.create(path)
  list(cache_dir = path, r2_disease = "medulloblastoma")
}

de_fixture <- function() {
  d <- expand.grid(gene = c("HLX", "SLC16A1", "RBM24", "UNKNOWN", "NEUTRAL"), comparator = c("group4", "shh", "wnt"), stringsAsFactors = FALSE)
  d$log2fc <- rep(c(2, 1.5, -2, 1.2, 0), 3)
  d$adjusted_p <- rep(c(0.001, 0.002, 0.003, 0.01, 0.9), 3)
  d
}

write_fixture <- function(d) { p <- tempfile(fileext = ".csv"); data.table::fwrite(d, p); p }

testthat::test_that("DE import preserves statistics and requires every comparator", {
  cfg <- test_cfg()
  d <- de_fixture()
  d$adjusted_p[d$gene == "SLC16A1" & d$comparator == "wnt"] <- 0.5
  d <- d[!(d$gene == "UNKNOWN" & d$comparator == "shh"), ]
  r <- import_de_results(cfg, write_fixture(d), "medulloblastoma", "group3")
  a <- load_study_object(cfg, r$id, "analysis")
  testthat::expect_setequal(a$summary$gene[a$summary$subgroup_specific], c("HLX", "RBM24"))
  testthat::expect_equal(a$data$log2fc, d$log2fc)
  testthat::expect_false("UNKNOWN" %in% a$universe)
  testthat::expect_equal(nrow(analysis_results(cfg, r$id, TRUE)$rows), 2)
  testthat::expect_equal(nrow(analysis_results(cfg, r$id, offset = 999)$rows), 0)
  testthat::expect_equal(length(list_studies(cfg)$studies), 1)
  one <- import_de_results(cfg, write_fixture(d[d$comparator == "group4", ]), "medulloblastoma", "group3")
  testthat::expect_false(any(load_study_object(cfg, one$id, "analysis")$summary$subgroup_specific))
})

testthat::test_that("invalid and ambiguous DE input fails", {
  cfg <- test_cfg(); d <- de_fixture()
  testthat::expect_error(import_de_results(cfg, write_fixture(rbind(d, d[1, ])), "MB", "g3"), "Duplicate")
  d$adjusted_p[1] <- 2
  testthat::expect_error(import_de_results(cfg, write_fixture(d), "MB", "g3"), "between")
  d$adjusted_p[1] <- 0.01; d$gene[1] <- "ENSG000001.2"
  testthat::expect_error(import_de_results(cfg, write_fixture(d), "MB", "g3"), "Map Ensembl")
  testthat::expect_error(load_study_object(cfg, "../anything", "analysis"), "Invalid")
})

testthat::test_that("ranking retains missing dependency evidence and joins case-normalized symbols", {
  cfg <- test_cfg(); d <- de_fixture(); d$gene <- tolower(d$gene)
  a <- import_de_results(cfg, write_fixture(d), "medulloblastoma", "group3")
  fixture <- function(p) file.path(Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR"), "tests", "fixtures", p)
  ranked <- prioritize_candidates(cfg, a$id, fixture("Model.csv"), fixture("CRISPRGeneEffect.csv"), depmap_release = "fixture")
  testthat::expect_equal(ranked$depmap_status, "ok")
  testthat::expect_true(all(c("HLX", "SLC16A1", "UNKNOWN") %in% ranked$candidates$gene))
  testthat::expect_true(is.na(ranked$candidates$priority_score[ranked$candidates$gene == "UNKNOWN"]))
  testthat::expect_equal(ranked$candidates$dependency_selectivity[ranked$candidates$gene == "HLX"], 0.6)
  gene <- query_gene_context(cfg, "hlx", "medulloblastoma", a$id, model_csv = fixture("Model.csv"), gene_effect_csv = fixture("CRISPRGeneEffect.csv"))
  testthat::expect_equal(gene$depmap$result$summary$mean_gene_effect, -0.7)
  testthat::expect_equal(gene$clinical$status, "unavailable")
  testthat::expect_error(query_gene_context(cfg, "HLX", "glioma", a$id), "does not match")
})

testthat::test_that("local enrichment uses measured background and adjusts zero-overlap terms", {
  sets <- list(hit = c("A", "B", "OUTSIDE"), zero = c("C", "D"))
  r <- overrepresentation(c("A", "B"), LETTERS[1:10], sets, 1, 10)
  expected <- stats::phyper(1, 2, 8, 2, lower.tail = FALSE)
  testthat::expect_equal(r$p_value[1], expected)
  testthat::expect_equal(r$adjusted_p[1], expected * 2)
  testthat::expect_equal(r$set_n[1], 2)
  testthat::expect_equal(nrow(r), 2)
  cfg <- test_cfg()
  a <- import_de_results(cfg, write_fixture(de_fixture()), "medulloblastoma", "group3")
  p <- tempfile(fileext = ".gmt")
  writeLines(c("up\tdescription\tHLX\tSLC16A1", "down\tdescription\tRBM24"), p)
  out <- enrich_pathways(cfg, a$id, gmt_path = p, min_set_size = 1, fdr = 1)
  testthat::expect_equal(length(out$results), 2)
  testthat::expect_equal(out$results[[1]]$input_genes_n, 3)
  testthat::expect_equal(out$results[[2]]$input_genes_n, 1)
})

make_cohort <- function(cfg, scale = "log2") {
  set.seed(402)
  groups <- rep(c("group3", "group4", "shh"), each = 20)
  mat <- matrix(rnorm(120 * 60, 7, 0.7), 120, 60)
  mat[1, groups == "group3"] <- mat[1, groups == "group3"] + 3
  if (scale == "counts") mat <- round(2^mat)
  colnames(mat) <- paste0("s", 1:60)
  e <- data.frame(gene = c("HLX", paste0("GENE", 2:120)), mat, check.names = FALSE)
  m <- data.frame(sample_id = colnames(mat), subgroup = groups, age = runif(60, 3, 18),
    survival_time = rexp(60, 0.03), survival_event = rbinom(60, 1, 0.65), metastasis = rbinom(60, 1, 0.5))
  register_cohort(cfg, write_fixture(e), write_fixture(m), "medulloblastoma", scale = scale)
}

testthat::test_that("cohort analysis works for log2 matrices and raw counts", {
  testthat::skip_if_not_installed("limma"); testthat::skip_if_not_installed("edgeR")
  for (scale in c("log2", "counts")) {
    cfg <- test_cfg(); c <- make_cohort(cfg, scale)
    a <- analyze_subgroup(cfg, c$id, "group3", covariates = "age")
    result <- load_study_object(cfg, a$id, "analysis")
    testthat::expect_true(result$summary$subgroup_specific[result$summary$gene == "HLX"])
    testthat::expect_true(all(result$data$log2fc[result$data$gene == "HLX"] > 2))
    testthat::expect_equal(length(result$comparators), 2)
  }
})

testthat::test_that("clinical models return uncertainty and explicit unavailable endpoints", {
  cfg <- test_cfg(); c <- make_cohort(cfg)
  r <- clinical_associations(cfg, c$id, c("HLX", "GENE2"), covariates = c("age", "subgroup"))
  testthat::expect_equal(nrow(r$rows), 4)
  testthat::expect_true(all(r$rows$status == "ok"))
  testthat::expect_true(all(r$rows$ci_lower < r$rows$effect_ratio & r$rows$ci_upper > r$rows$effect_ratio))
  testthat::expect_true(all(r$rows$adjusted_p >= r$rows$p_value))
  absent <- clinical_associations(cfg, c$id, "MISSING")
  testthat::expect_true(all(absent$rows$status == "unavailable"))
  testthat::expect_true(all(is.na(absent$rows$p_value)))
  missing_outcome <- clinical_associations(cfg, c$id, "HLX", metastasis_column = "unrecorded")
  testthat::expect_equal(missing_outcome$rows$status[2], "unavailable")
})

testthat::test_that("MCP schemas validate required fields, ranges and types", {
  schema <- jsonlite::fromJSON(json_text(tool_definitions(), FALSE), simplifyVector = FALSE)
  testthat::expect_true(is.list(schema[[3]]$inputSchema$required))
  testthat::expect_error(validate_tool_arguments("import_de_results", list(path = "x")), "Missing")
  testthat::expect_error(validate_tool_arguments("analysis_results", list(analysis_id = "x", limit = -1)), "range")
  testthat::expect_error(validate_tool_arguments("clinical_associations", list(cohort_id = "x", genes = list(2))), "type")
  testthat::expect_error(validate_tool_arguments("analysis_results", list(analysis_id = "x", unwanted = TRUE)), "Unknown")
  testthat::expect_silent(validate_tool_arguments("list_studies", list()))
})
