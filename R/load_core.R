# Shared bootstrap for R scripts, Shiny and MCP. Never starts a server.
load_r2_core <- function(root = ".", envir = new.env(parent = globalenv())) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  .libPaths(c(file.path(root, ".R-library"), .libPaths()))
  modules <- c("common", "cache", "r2_client", "r2_expression", "subtype_analysis", "pubmed",
    "study", "clinical", "enrichment", "depmap", "ranking", "integration",
    "profile_statistics", "profile_sources", "profile_enrichr", "profile_depmap",
    "gene_analysis", "browser_report", "profile_report", "clinical_v3")
  for (module in modules) sys.source(file.path(root, "R", paste0(module, ".R")), envir)
  envir$.r2_core_functions <- unique(unlist(lapply(modules, function(module) {
    expressions <- parse(file.path(root, "R", paste0(module, ".R")))
    vapply(expressions, function(expr) {
      if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
          is.symbol(expr[[2]]) && is.call(expr[[3]]) &&
          identical(expr[[3]][[1]], as.name("function"))) as.character(expr[[2]]) else ""
    }, character(1))
  })))
  envir$.r2_core_functions <- setdiff(envir$.r2_core_functions, "")
  invisible(envir)
}
