v1_fixture <- function() {
  ids <- paste0("S", 1:12)
  table <- list(tableData = list(data = data.frame(samplenames = ids, x = as.character(1:12))))
  plot <- list(table = "ps_avgpres_gse85217geo763_hugene11t", corType = "transform_log2", plotData = list(reporter = "123", reporterSymbol = "HLX"))
  html <- paste0("<html><body><script>import('/loadDataTableModal.js').then(m => m['default'](", json_text(table, FALSE), "));</script>",
    "<script>import('/d3/plots/plot.js').then(m => m['default'](", json_text(plot, FALSE), "));</script></body></html>")
  metadata <- data.frame(samplenames = ids, subgroup = rep(c("wnt", "shh", "group3", "group4"), each = 3),
    subtype = rep(c("original_A", "original_B", "original_C", "original_D"), each = 3), check.names = FALSE)
  metadata[["met_status_(1_met__0_m0)"]] <- rep(c("0", "1", "na"), 4)
  metadata[["os_(years)"]] <- rep(c("1", "2", "na"), 4)
  metadata$dead <- rep(c("0", "1", "na"), 4)
  metadata[["geo-gse85217"]] <- paste0("DIFFERENT_GEO_ID_", 1:12)
  list(html = html, annotations = list(data = metadata), info = list(sampleSize = "12", pubmedId = "fixture", design = "Synthetic fixture"),
    tracks = list(trackDescriptions = "met_status_(1_met__0_m0): 0, 1, na"))
}

v1_assembled_fixture <- function() {
  f <- v1_fixture()
  e <- parse_r2_patient_expression(f$html, "HLX", "123", "ps_avgpres_gse85217geo763_hugene11t")
  r <- assemble_r2_patient_data(e, f$annotations, f$info, f$tracks)
  r$provenance <- list(queried_gene = "HLX", reporter = "123", dataset_table = "ps_avgpres_gse85217geo763_hugene11t", source = "SYNTHETIC TEST FIXTURE", output_dir = tempdir())
  r$files <- list(patient_data = "synthetic-test-fixture.csv")
  r
}
