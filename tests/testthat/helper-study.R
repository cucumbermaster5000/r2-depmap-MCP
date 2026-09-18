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
