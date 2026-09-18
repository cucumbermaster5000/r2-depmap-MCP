# Open this file in RStudio and click Run App, or shiny::runApp("shiny").
root <- if (file.exists(file.path("..", "R", "load_core.R"))) ".." else "."
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "load_core.R"), local = TRUE)
core <- load_r2_core(root)
for (package in c("shiny", "bslib")) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install web dependencies with: Rscript scripts/install-dependencies.R --shiny")
}
cfg <- core$read_server_config(file.path(root, "config", "defaults.json"))
cfg$output_dir <- file.path(root, "outputs")
source(file.path(root, "shiny", "interface.R"), local = TRUE, encoding = "UTF-8")
shiny::shinyApp(explorer_ui(), explorer_server(core, cfg))
