`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

log_stderr <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), ..., "\n", file = stderr())
}

require_runtime_packages <- function() {
  required <- c("jsonlite", "curl", "xml2", "rvest", "data.table")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Run scripts/install-dependencies.R.",
      call. = FALSE
    )
  }
}

read_server_config <- function(path) {
  if (!file.exists(path)) stop("Configuration file not found: ", path, call. = FALSE)
  cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  cfg$cache_dir <- path.expand(Sys.getenv("R2_DEPMAP_CACHE_DIR", cfg$cache_dir))
  if (!grepl("^(/|[A-Za-z]:|\\\\)", cfg$cache_dir)) cfg$cache_dir <- file.path(dirname(path), "..", cfg$cache_dir)
  cfg$r2_base_url <- Sys.getenv("R2_BASE_URL", cfg$r2_base_url)
  cfg$r2_dataset_table <- Sys.getenv("R2_DATASET_TABLE", cfg$r2_dataset_table)
  dir.create(cfg$cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cfg$cache_dir)) stop("Cannot create cache directory: ", cfg$cache_dir, call. = FALSE)
  cfg$cache_dir <- normalizePath(cfg$cache_dir, winslash = "/", mustWork = TRUE)
  cfg
}

as_number <- function(x, default = NA_real_) {
  out <- suppressWarnings(as.numeric(x))
  if (!length(out) || is.na(out[[1]])) default else out[[1]]
}

as_integer_scalar <- function(x, default, minimum = NULL, maximum = NULL) {
  out <- suppressWarnings(as.integer(x %||% default))
  if (is.na(out)) out <- as.integer(default)
  if (!is.null(minimum)) out <- max(out, minimum)
  if (!is.null(maximum)) out <- min(out, maximum)
  out
}

as_logical_scalar <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(default)
  if (is.logical(x)) return(isTRUE(x[[1]]))
  tolower(as.character(x[[1]])) %in% c("1", "true", "yes", "y")
}

safe_zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

clean_symbol <- function(x) {
  trimws(sub("\\s*\\([^)]*\\)\\s*$", "", as.character(x)))
}

json_text <- function(x, pretty = TRUE) {
  jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null",
    na = "null",
    digits = 8,
    pretty = pretty
  )
}

normalize_local_path <- function(path, env_name = NULL) {
  value <- path %||% if (!is.null(env_name)) Sys.getenv(env_name, "") else ""
  if (!nzchar(value)) {
    stop("A file path is required", if (!is.null(env_name)) paste0(" (or set ", env_name, ")") else "", call. = FALSE)
  }
  value <- path.expand(value)
  if (!file.exists(value)) stop("File not found: ", value, call. = FALSE)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}
