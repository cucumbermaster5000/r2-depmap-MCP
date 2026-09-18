source("R/load_core.R")
core <- load_r2_core(".")
source("shiny/interface.R",encoding="UTF-8")
cfg <- core$read_server_config("config/defaults.json")
cfg$output_dir <- normalizePath("outputs",winslash="/")
folder <- "outputs/shiny-v3-acceptance"
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
baseline <- as.data.frame(data.table::fread("outputs/shiny-v2-acceptance/acceptance.csv"))
rows <- lapply(c("HLX","MYC","SLC2A1"),function(gene) {
  message("Checking V3 ",gene)
  r <- core$get_r2_expression(cfg,gene)
  v1 <- core$analyze_r2_subgroups(r); v2 <- core$analyze_r2_subtypes(r)
  old <- baseline[baseline$gene==gene,]
  stopifnot(isTRUE(all.equal(v1$statistics$statistic,old$v1_H)),isTRUE(all.equal(v1$statistics$p_value,old$v1_p)),
    isTRUE(all.equal(v2$overall$statistic,old$H)),isTRUE(all.equal(v2$overall$p_value,old$p)))
  old_v2 <- readRDS(file.path("outputs/shiny-v2-acceptance",paste0(gene,"_analysis.rds")))
  stopifnot(isTRUE(all.equal(v2$pairwise,old_v2$pairwise)))
  a <- core$analyze_r2_clinical_v3(cfg,r,v1)
  stopifnot(a$verification$metastasis$verified,a$verification$survival$verified,
    a$metastasis$all$status=="complete",a$survival$all$status=="complete")
  stopifnot(setequal(names(a$survival),c("all",a$selection$selected)))
  for(ep in c("metastasis","survival")) for(co in names(a[[ep]])) {
    obj <- a[[ep]][[co]]
    stopifnot(obj$missing$cohort_n==obj$missing$used_n+obj$missing$excluded_n)
    if(ep=="survival") stopifnot(isTRUE(all.equal(obj$cutoff$value,mean(obj$data$expression[is.finite(obj$data$expression)]))))
    for(kind in c("patient_data","statistics",if(!is.null(obj$plot)) c("plot_pdf","plot_png"))) {
      ext <- if(kind=="plot_pdf") "pdf" else if(kind=="plot_png") "png" else "csv"
      path <- file.path(folder,paste(gene,co,ep,paste0(kind,".",ext),sep="_"))
      core$write_r2_clinical_download(obj,kind,path)
      stopifnot(file.info(path)$size>50)
    }
  }
  mock <- list2env(as.list(core),parent=globalenv())
  mock$get_r2_expression <- function(...) r
  mock$analyze_r2_clinical_v3 <- function(...) a
  shiny::testServer(explorer_server(mock,cfg),{
    session$setInputs(gene=gene);session$setInputs(run=1)
    stopifnot(is.null(failure()),!is.null(result()$clinical),!is.null(output$overview),
      !is.null(output$metastasis_content),!is.null(output$survival_content),
      !is.null(output$metastasis_all_plot),!is.null(output$survival_all_plot),
      !is.null(output$subtype_plot),!is.null(output$subgroup_plot))
  })
  saveRDS(a,file.path(folder,paste0(gene,"_clinical.rds")))
  data.frame(gene=gene,selected=paste(a$selection$selected,collapse=";"),metastasis_n=a$metastasis$all$missing$used_n,
    metastasis_p=a$metastasis$all$statistics$p_value,survival_n=a$survival$all$statistics$n,
    events=a$survival$all$statistics$events,mean_cutoff=a$survival$all$cutoff$value,
    logrank_p=a$survival$all$statistics$logrank_p,HR=a$survival$all$statistics$hazard_ratio_high_vs_low,
    cox_p=a$survival$all$statistics$cox_p,PH_p=a$survival$all$statistics$ph_test_p)
})
results <- do.call(rbind,rows)
data.table::fwrite(results,file.path(folder,"acceptance.csv"));print(results)
