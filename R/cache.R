cache_key <- function(prefix, parts) {
  raw <- paste(c(prefix, unlist(parts, use.names = TRUE)), collapse = "__")
  safe <- gsub("[^A-Za-z0-9._-]+", "-", raw)
  if (nchar(safe) > 180L) {
    tmp <- tempfile()
    writeBin(charToRaw(raw), tmp)
    on.exit(unlink(tmp), add = TRUE)
    safe <- paste0(substr(safe, 1L, 120L), "-", unname(tools::md5sum(tmp)))
  }
  safe
}

cache_path <- function(cfg, key) {
  file.path(cfg$cache_dir, paste0(key, ".rds"))
}

cache_get <- function(cfg, key) {
  path <- cache_path(cfg, key)
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

cache_put <- function(cfg, key, value) {
  path <- cache_path(cfg, key)
  tmp <- paste0(path, ".tmp")
  saveRDS(value, tmp, version = 3)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not write cache file: ", path, call. = FALSE)
  }
  invisible(path)
}

cache_inventory <- function(cfg) {
  files <- list.files(cfg$cache_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) return(data.frame())
  rows <- lapply(files, function(path) {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    data.frame(
      file = basename(path),
      created_at = obj$provenance$retrieved_at %||% format(file.info(path)$mtime, tz = "UTC"),
      dataset = obj$provenance$dataset_table %||% NA_character_,
      contrast = obj$provenance$contrast %||% NA_character_,
      rows = if (is.data.frame(obj$data)) nrow(obj$data) else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
