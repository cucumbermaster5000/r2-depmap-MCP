pub_fixture <- function() paste(readLines(file.path(Sys.getenv("R2_DEPMAP_TEST_PLUGIN_DIR"), "tests", "fixtures", "pubmed.xml")), collapse = "\n")

testthat::test_that("PubMed parsing preserves dates, absent DOI, sourced descriptions and notices", {
  p <- parse_pubmed_articles(pub_fixture(), "HLX")
  testthat::expect_length(p, 2)
  testthat::expect_equal(p[[1]]$publication_date, "2026 Sep")
  testthat::expect_equal(p[[1]]$electronic_date, "2026 08 02")
  testthat::expect_equal(p[[1]]$description, "HLX expression was associated with the experimental phenotype.")
  testthat::expect_equal(p[[1]]$doi, "10.1234/example")
  testthat::expect_true(p[[1]]$retracted)
  testthat::expect_equal(p[[2]]$publication_date, "2025 Winter")
  testthat::expect_true(is.na(p[[2]]$doi))
  testthat::expect_null(p[[2]]$doi_url)
  testthat::expect_equal(p[[2]]$description_source, "title_only_no_abstract")
  testthat::expect_error(parse_pubmed_articles("<ERROR>invalid</ERROR>", "HLX"), "error")
})

testthat::test_that("literature queries are explicit and aliases require matching human gene identity", {
  testthat::expect_false(pubmed_mentions("HLX-02 and HLX\u201102", "HLX"))
  testthat::expect_true(pubmed_mentions("HLX, a gene", "HLX"))
  testthat::expect_true(pubmed_mentions("HLX\u00a0expression", "HLX"))
  testthat::expect_equal(pubmed_query(c("SLC16A1", "MCT1"), "medulloblastoma"), '("SLC16A1"[Title/Abstract] OR "MCT1"[Title/Abstract]) AND "medulloblastoma"[Title/Abstract]')
  testthat::expect_error(pubmed_gene_terms("../x"), "gene symbol")
  testthat::expect_error(pubmed_query("HLX", 'x" OR cancer'), "disease")
  mock <- function(endpoint, params) {
    if (endpoint == "esearch.fcgi") return('{"esearchresult":{"idlist":["1"]}}')
    '{"result":{"1":{"uid":"1","name":"SLC16A1","organism":{"taxid":9606},"otheraliases":"MCT1, MCT","description":"transporter"}}}'
  }
  testthat::expect_equal(resolve_pubmed_aliases("SLC16A1", mock)$aliases, c("MCT1", "MCT"))
  testthat::expect_length(resolve_pubmed_aliases("HLX", mock)$aliases, 0)
})

testthat::test_that("HTML embeds plot and patient data while escaping external paper text", {
  testthat::skip_if_not_installed("ggplot2")
  r <- v1_assembled_fixture(); out <- tempfile(); dir.create(out)
  r$provenance$output_dir <- out
  v <- plot_r2_subgroup_test(r, out)
  testthat::expect_false(any(grepl("pdf|png", names(v$files))))
  pubs <- list(articles = parse_pubmed_articles(pub_fixture(), "HLX"), total_matches = 2,
    provenance = list(query = '"HLX"[Title/Abstract]', retrieved_at = "2026-09-18", sort = "pub_date", alias_resolution = list(status = "fixture"), limitations = "fixture"))
  path <- write_gene_page(r, v, pubs)
  html <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  testthat::expect_match(html, "<svg", fixed = TRUE)
  testthat::expect_match(html, "&lt;script&gt;test&lt;/script&gt;", fixed = TRUE)
  testthat::expect_false(grepl("<script>test</script>", html, fixed = TRUE))
  testthat::expect_match(html, "https://doi.org/", fixed = TRUE)
  testthat::expect_match(html, "Not reported in PubMed", fixed = TRUE)
  testthat::expect_match(html, "S12", fixed = TRUE)
  testthat::expect_match(html, "Retracted publication", fixed = TRUE)
  testthat::expect_false(grepl("<script src=", html, fixed = TRUE))
  testthat::expect_match(paste(readLines(write_gene_page(r, v, publication_error = "Offline")), collapse = ""), "Publications unavailable: Offline", fixed = TRUE)
})
