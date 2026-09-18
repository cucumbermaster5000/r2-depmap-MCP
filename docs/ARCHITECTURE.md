# Reusable R analysis architecture

The biological and statistical functions are ordinary R functions. MCP is an
optional adapter; no AI agent, MCP process or Shiny session is required to run
an analysis. `R/load_core.R` is the shared bootstrap for scripts and interfaces.
It loads definitions into an environment without starting a server or retrieving
data. It adds the project's installed R library to `.libPaths()`.

## Boundaries

| Layer | Files / entry points | Responsibility |
| --- | --- | --- |
| Sources and cache | `r2_client.R`, `r2_expression.R`, `profile_sources.R`, `pubmed.R`, `profile_depmap.R` | Retrieve, verify, normalize and cache source data with provenance |
| Reusable analyses | `analyze_r2_subgroups()`, `profile_statistics.R`, `clinical.R`, `study.R`, `enrichment.R`, `profile_enrichr.R`, `depmap.R` | Group tests, clinical associations, enrichment and dependency summaries |
| Application services | `gene_analysis.R` | Collect sources and assemble structured gene analysis results |
| Report presentation | `browser_report.R`, `profile_report.R` | Render and export computed results to HTML, RDS and CSV |
| MCP adapter | `v1_tools.R`, `server.R` | Tool schemas, request validation, dispatch, response formatting and JSON-RPC transport |

Legacy `tools.R` and `study_tools.R` are adapters for the broader study API;
they are not loaded by the core or enabled by the current MCP server.
Core functions must never call a dispatcher, inspect MCP requests, or depend
on an agent. New analyses belong in reusable R modules with direct-R tests.
Adapters may paginate or omit large objects for transport but must not implement
biological calculations or change scientific interpretation.

## Direct R API

```r
source("R/load_core.R")
core <- load_r2_core(".")
cfg <- core$read_server_config("config/defaults.json")
cfg$output_dir <- "outputs"

# Network/cache access; returns verified source objects and optional failures.
sources <- core$collect_gene_sources(cfg, "HMGCR",
  publication_disease = "medulloblastoma")
saveRDS(sources, "HMGCR_sources.rds")

# In-memory, no network access or file writes. Can be called again offline.
result <- core$analyze_gene_sources(sources, cutoff = "median")
result$validation$statistics
result$validation$plot
result$profile$clinical[["all:median"]]$survival

# Optional export; no source retrieval or re-fitting.
page <- core$export_gene_report(result, cfg)

# Convenience services using those same stages:
result <- core$analyze_gene(cfg, "HMGCR")
page <- core$build_gene_page(cfg, "HMGCR")
```

`collect_gene_sources()` returns `retrieval`, `survival`, `pfister`, `enrichr`,
`depmap`, and `publications`. Cavalli retrieval is required; optional modules
carry `status = "unavailable"` and a reason on failure. Source services may
perform normalization and annotation enrichment and save their own caches.
`analyze_gene_sources()` accepts these objects, including omitted optional
sources, and returns `retrieval`, `validation`, `publications`,
`publication_error`, `profile`, and `cutoff`. Results retain exclusions,
provenance, statistical tables, clinical panels and ggplot objects. All three
clinical cutoff methods are computed; `cutoff` selects the initial HTML panel.
The statistical methods and correction families are unchanged by this refactor.

For an individual comparison, `analyze_r2_subgroups(retrieval)` returns tables,
exclusions and a plot without writing files. The existing
`plot_r2_subgroup_test()` wraps it with the historical artifact exports.

## Future Shiny integration

Load the core once during application startup. Keep each user's sources and
results in session-local reactives; pass configuration explicitly. A server
function can call the API directly:

```r
server <- function(input, output, session) {
  result <- shiny::eventReactive(input$run, {
    shiny::req(input$gene)
    core$analyze_gene(cfg, input$gene)
  })
  output$subgroups <- shiny::renderPlot(result()$validation$plot)
  output$statistics <- shiny::renderTable(result()$validation$statistics)
  output$coverage <- shiny::renderTable(result()$profile$status)
}
```

This is an integration example, not a deployed Shiny app. Slow source downloads
can later run in a background worker without moving analysis into the UI. A
GUI should show optional-module failures and retain the report's exploratory
interpretation and cutoff caveats.

## Verification

Run `Rscript --vanilla tests/run-tests.R`. Core API tests load an isolated
environment with no MCP functions, analyze synthetic data without I/O, export
and regenerate a report offline, and compare direct-R and MCP-dispatched results.
Existing statistical tests verify the underlying estimators and corrections.
