#!/usr/bin/env Rscript

plugin_dir <- Sys.getenv("R2_DEPMAP_PLUGIN_DIR", "")
if (!nzchar(plugin_dir)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  plugin_dir <- normalizePath(file.path(dirname(file_arg[[1]]), ".."), mustWork = TRUE)
}

source(file.path(plugin_dir, "R", "load_core.R"))
core <- load_r2_core(plugin_dir)
adapter <- new.env(parent = core)
sys.source(file.path(plugin_dir, "R", "v1_tools.R"), adapter)
tool_definitions <- adapter$v1_tool_definitions
dispatch_tool <- adapter$dispatch_v1_tool
json_text <- core$json_text
log_stderr <- core$log_stderr
`%||%` <- core$`%||%`
core$require_runtime_packages()
config_path <- Sys.getenv("R2_DEPMAP_CONFIG", file.path(plugin_dir, "config", "defaults.json"))
cfg <- core$read_server_config(config_path)

write_message <- function(message) {
  cat(json_text(message, pretty = FALSE), "\n", sep = "")
  flush.console()
}

rpc_error <- function(id, code, message, data = NULL) {
  error <- list(code = code, message = message)
  if (!is.null(data)) error$data <- data
  list(jsonrpc = "2.0", id = id, error = error)
}

tool_response <- function(value) {
  list(
    content = list(list(type = "text", text = json_text(value, pretty = TRUE))),
    structuredContent = value,
    isError = FALSE
  )
}

handle_request <- function(request) {
  id <- request$id %||% NULL
  method <- request$method %||% ""

  if (method == "initialize") {
    return(list(
      jsonrpc = "2.0",
      id = id,
      result = list(
        protocolVersion = if ((request$params$protocolVersion %||% "") %in% c("2024-11-05", "2025-03-26", "2025-06-18")) request$params$protocolVersion else "2025-06-18",
        capabilities = list(tools = list(listChanged = FALSE)),
        serverInfo = list(name = "r2-depmap", version = "0.4.0"),
        instructions = "Use build_gene_page for the integrated browser gene report: Cavalli subgroups/subtypes, Pfister CNS tumours, verified survival/metastasis with selectable cutoffs, Enrichr membership/TF target enrichment, pinned public DepMap dependencies and PubMed. Report unavailable sections explicitly. Optimised cutoffs are exploratory; annotations do not establish pathway deregulation."
      )
    ))
  }

  if (method == "ping") return(list(jsonrpc = "2.0", id = id, result = structure(list(), names = character())))
  if (method == "tools/list") {
    return(list(jsonrpc = "2.0", id = id, result = list(tools = tool_definitions())))
  }
  if (method == "tools/call") {
    params <- request$params %||% list()
    value <- tryCatch(tool_response(dispatch_tool(params$name %||% "", params$arguments %||% list(), cfg)),
      error = function(e) list(content = list(list(type = "text", text = conditionMessage(e))), isError = TRUE))
    return(list(jsonrpc = "2.0", id = id, result = value))
  }

  if (startsWith(method, "notifications/")) return(NULL)
  rpc_error(id, -32601L, paste0("Method not found: ", method))
}

log_stderr("Starting r2-depmap MCP server")
input <- file("stdin", open = "r")

repeat {
  line <- readLines(input, n = 1L, warn = FALSE)
  if (!length(line)) break
  if (!nzchar(trimws(line))) next

  request <- tryCatch(
    jsonlite::fromJSON(line, simplifyVector = FALSE),
    error = function(e) e
  )
  if (inherits(request, "error")) {
    write_message(rpc_error(NULL, -32700L, "Parse error", conditionMessage(request)))
    next
  }

  if (!is.list(request) || is.null(names(request)) || !identical(request$jsonrpc, "2.0") ||
      !is.character(request$method) || length(request$method) != 1L) {
    write_message(rpc_error(NULL, -32600L, "Invalid Request"))
    next
  }

  response <- tryCatch(
    handle_request(request),
    error = function(e) {
      log_stderr("Request failed:", conditionMessage(e))
      rpc_error(request$id %||% NULL, -32603L, conditionMessage(e))
    }
  )
  if (!is.null(response) && !is.null(request$id)) write_message(response)
}

log_stderr("Stopping r2-depmap MCP server")
close(input)
