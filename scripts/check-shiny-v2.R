# Real-data acceptance; only expression, V1 subgroup and V2 subtype analyses.
source("R/load_core.R")
core <- load_r2_core(".")
source("shiny/interface.R", encoding="UTF-8")
cfg <- core$read_server_config("config/defaults.json")
cfg$output_dir <- normalizePath("outputs",winslash="/")
folder <- "outputs/shiny-v2-acceptance"
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
baseline <- as.data.frame(data.table::fread("outputs/shiny-v1-acceptance/acceptance.csv"))
rows <- lapply(c("HLX","MYC","SLC2A1"),function(gene) {
  message("Checking ",gene)
  r <- core$get_r2_expression(cfg,gene)
  a <- core$analyze_r2_subtypes(r)
  v1 <- core$analyze_r2_subgroups(r)
  expected <- baseline[baseline$gene==gene,]
  if(nrow(expected)) stopifnot(isTRUE(all.equal(v1$statistics$p_value,expected$p)),
    isTRUE(all.equal(v1$statistics$statistic,expected$H)))
  stopifnot(nrow(a$data)==763,nrow(a$summaries)==12,nrow(a$pairwise)==66,
    a$metadata$missing_n==0,sum(a$summaries$n)==763,a$status=="complete",
    isTRUE(all.equal(a$pairwise$adjusted_p_value,p.adjust(a$pairwise$p_value,"BH"))))
  # Drive the actual Shiny server with real retrievals and exercise renderers.
  local_core <- list2env(as.list(core), parent=globalenv())
  local_core$get_r2_expression <- function(...) r
  shiny::testServer(explorer_server(local_core,cfg), {
    session$setInputs(gene=gene); session$setInputs(run=1)
    stopifnot(is.null(failure()),result()$gene==gene,nrow(result()$subtypes$pairwise)==66)
    stopifnot(nchar(output$subtype_pairs)>100,nchar(output$subtype_descriptives)>100,
      nchar(output$subtype_overall)>100,!is.null(output$subtype_plot))
    x <- result()
    for(kind in c("patients","statistics","pdf","png","subtype_patients","subtype_summaries",
      "subtype_pairwise","subtype_overall","subtype_pdf","subtype_png")) {
      ext <- if(grepl("pdf$",kind)) "pdf" else if(grepl("png$",kind)) "png" else "csv"
      path <- file.path(folder,paste0(gene,"_",kind,".",ext))
      write_explorer_download(x,kind,path)
      stopifnot(file.info(path)$size>100)
      if(kind=="subtype_pairwise") stopifnot(isTRUE(all.equal(data.table::fread(path)$p_value,a$pairwise$p_value)))
      if(kind=="subtype_patients") stopifnot(nrow(data.table::fread(path))==763)
    }
  })
  saveRDS(a,file.path(folder,paste0(gene,"_analysis.rds")))
  data.frame(gene=gene,n=a$overall$n,subtypes=nrow(a$summaries),pairs=nrow(a$pairwise),
    H=a$overall$statistic,p=a$overall$p_value,epsilon_squared=a$overall$epsilon_squared,
    v1_H=v1$statistics$statistic,v1_p=v1$statistics$p_value,cache_hit=isTRUE(r$provenance$cache_hit))
})
results <- do.call(rbind,rows)
data.table::fwrite(results,file.path(folder,"acceptance.csv"))
print(results)
