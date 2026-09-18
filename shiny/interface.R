# Presentation and orchestration only. Scientific computations live in R/.
explorer_ui <- function() {
  future <- function(title) bslib::nav_panel(title,
    bslib::card(bslib::card_header(title), shiny::p("Not implemented yet")))
  bslib::page_sidebar(
    title = "Medulloblastoma Gene Explorer",
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#216b75"),
    fillable = FALSE,
    sidebar = bslib::sidebar(width = 270,
      shiny::p(class = "text-uppercase text-muted small", "Cavalli cohort / GSE85217"),
      shiny::textInput("gene", "Gene symbol", value = "HLX", placeholder = "e.g. HLX"),
      bslib::input_task_button("run", "Run analysis", label_busy = "Analyzing-", width = "100%"),
      shiny::p(class = "small text-muted", "Enter an official gene symbol. Results update only when you run an analysis."),
      shiny::hr(), shiny::p(class = "small", "Version 3 - Expression, metastasis and survival")),
    shiny::uiOutput("status"),
    bslib::navset_card_tab(id = "section",
      bslib::nav_panel("Overview", shiny::uiOutput("overview")),
      bslib::nav_panel("MB Subgroups",
        shiny::uiOutput("plot_heading"),
        shiny::div(style = "overflow-x:auto", shiny::div(style = "min-width:720px",
          shiny::plotOutput("subgroup_plot", height = "540px"))),
        shiny::h4("Statistical comparison"),
        shiny::div(style = "overflow-x:auto", shiny::tableOutput("statistics")),
        shiny::p(class = "small text-muted", "The original omnibus result is unchanged. Additional pairwise evidence for clinical cohort selection is shown below."),
        shiny::uiOutput("selection_summary"),
        shiny::div(style="overflow-x:auto",shiny::tableOutput("selection_pairs")),
        shiny::h4("Subgroup summaries"), shiny::tableOutput("summaries")),
      bslib::nav_panel("MB Subtypes",
        shiny::uiOutput("subtype_summary"),
        shiny::div(style="overflow-x:auto", shiny::div(style="min-width:1100px",
          shiny::plotOutput("subtype_plot", height="670px"))),
        shiny::h4("Descriptive statistics"),
        shiny::div(style="overflow-x:auto", shiny::tableOutput("subtype_descriptives")),
        shiny::h4("Overall comparison"),
        shiny::div(style="overflow-x:auto", shiny::tableOutput("subtype_overall")),
        shiny::h4("All pairwise subtype comparisons"),
        shiny::p("Positive rank-biserial effect means subtype A tends higher than B; negative means lower. Median difference is A minus B on the source expression scale. BH adjustment covers all subtype pairs for this gene. These results do not establish subtype specificity."),
        shiny::div(style="overflow:auto;max-height:600px", shiny::tableOutput("subtype_pairs"))),
      bslib::nav_panel("Survival",shiny::uiOutput("survival_content")),
      bslib::nav_panel("Metastasis",shiny::uiOutput("metastasis_content")),
      future("Functional Biology"), future("TF Targets"), future("DepMap"),
      bslib::nav_panel("Downloads", shiny::uiOutput("downloads"))
    )
  )
}

explorer_error <- function(error, stage) {
  detail <- conditionMessage(error)
  if (grepl("gene symbol|Gene absent|found 0", detail))
    return("Gene not found in the selected R2 dataset.")
  if (grepl("Specify reporter|exact reporter", detail))
    return("This gene has multiple R2 reporters. Reporter selection is not available in Version 1.")
  if (grepl("subgroup|annotated subgroups", detail, ignore.case = TRUE))
    return("MB subgroup annotation is unavailable for this dataset.")
  if (stage == "retrieval") return("R2 data retrieval failed. Please try again later.")
  "The subgroup comparison could not be completed. Please try another gene."
}

# Uses the stored data/plot; these handlers never run statistical analysis.
write_explorer_download <- function(result, kind, file) {
  switch(kind,
    subtype_patients = data.table::fwrite(result$subtypes$data, file, na="NA"),
    subtype_summaries = data.table::fwrite(result$subtypes$summaries, file, na="NA"),
    subtype_pairwise = data.table::fwrite(result$subtypes$pairwise, file, na="NA"),
    subtype_overall = data.table::fwrite(result$subtypes$overall, file, na="NA"),
    subtype_pdf = ggplot2::ggsave(file, result$subtypes$plot, device=grDevices::pdf, width=13, height=7.5),
    subtype_png = ggplot2::ggsave(file, result$subtypes$plot, device="png", width=13, height=7.5, dpi=300),
    patients = data.table::fwrite(result$retrieval$data, file, na = "NA"),
    statistics = data.table::fwrite(result$analysis$statistics, file, na = "NA"),
    pdf = ggplot2::ggsave(file, result$analysis$plot, device = grDevices::pdf,
      width = 8, height = 5.8),
    png = ggplot2::ggsave(file, result$analysis$plot, device = "png",
      width = 8, height = 5.8, dpi = 300),
    stop("Unknown download type"))
  invisible(file)
}

explorer_server <- function(core, cfg) {
  force(core); force(cfg)
  function(input, output, session) {
    result <- shiny::reactiveVal(NULL)
    failure <- shiny::reactiveVal(NULL)
    busy <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(input$run, {
      if (busy()) return()
      busy(TRUE)
      on.exit(busy(FALSE), add = TRUE)
      result(NULL); failure(NULL)
      gene <- toupper(trimws(input$gene))
      if (length(gene) != 1L || is.na(gene) || !grepl("^[A-Z][A-Z0-9._-]*$", gene)) {
        failure("Please enter one valid gene symbol, for example HLX or SLC2A1.")
        return()
      }
      stage <- "retrieval"
      tryCatch(shiny::withProgress(message = paste("Analyzing", gene), value = 0, {
        shiny::incProgress(.15, detail = "Retrieving R2 patient expression")
        retrieval <- core$get_r2_expression(cfg, gene)
        stage <- "analysis"
        shiny::incProgress(.55, detail = "Comparing broad MB subgroups")
        if (!"subgroup" %in% names(retrieval$data) || all(core$r2_missing(retrieval$data$subgroup)))
          stop("MB subgroup annotation unavailable")
        analysis <- core$analyze_r2_subgroups(retrieval)
        shiny::incProgress(.15, detail="Comparing molecular subtypes")
        subtype_failure <- NULL
        subtypes <- tryCatch(core$analyze_r2_subtypes(retrieval), error=function(e) {
          message(sprintf("[Gene Explorer] gene=%s stage=subtypes: %s", gene, conditionMessage(e)))
          subtype_failure <<- "Molecular subtype analysis could not be completed. Broad subgroup results remain available."
          NULL
        })
        clinical <- NULL
        if(is.function(core$analyze_r2_clinical_v3)) clinical <- tryCatch(
          core$analyze_r2_clinical_v3(cfg,retrieval,analysis),error=function(e) {
            message("[Gene Explorer clinical] ",conditionMessage(e)); NULL
          })
        result(list(gene = retrieval$provenance$queried_gene, retrieval = retrieval, analysis = analysis,
          subtypes=subtypes, subtype_failure=subtype_failure,clinical=clinical))
        shiny::incProgress(.3, detail = "Preparing figure and tables")
      }), error = function(e) {
        # Technical details go to the R console/server log, never the web page.
        message(sprintf("[Gene Explorer] gene=%s stage=%s: %s", gene, stage, conditionMessage(e)))
        failure(explorer_error(e, stage))
      })
    }, ignoreInit = TRUE)

    output$status <- shiny::renderUI({
      if (!is.null(failure())) return(shiny::div(class = "alert alert-danger", role = "alert", failure()))
      if (is.null(result())) return(shiny::div(class = "alert alert-info", "Enter a gene and click Run analysis to begin."))
      shiny::div(class = "alert alert-success", role = "status", paste("Analysis complete:", result()$gene))
    })
    output$overview <- shiny::renderUI({
      x <- result()
      if (is.null(x)) return(shiny::p("No analysis yet."))
      s <- x$analysis$statistics
      metadata <- names(x$retrieval$metadata)
      shiny::tagList(
        bslib::layout_columns(
          bslib::value_box("Queried gene", x$gene),
          bslib::value_box("Patient samples retrieved", nrow(x$retrieval$data)),
          bslib::value_box("Samples analyzed", s$n)),
        shiny::h4("Dataset"), shiny::p(x$retrieval$provenance$dataset_label),
        shiny::p("GSE85217 - R2 log2 expression - WNT, SHH, Group 3 and Group 4"),
        shiny::h4("Broad subgroup result"),
        shiny::p(sprintf("Kruskal-Wallis comparison across %d subgroups: H(%d) = %.3f, raw p = %s. %d samples excluded for missing expression or subgroup annotation.",
          s$groups, s$df, s$statistic, format(s$p_value, digits = 4, scientific = TRUE), s$excluded_n)),
        shiny::div(class = "alert alert-warning", "This is an unadjusted comparison of expression distributions. It does not establish subgroup specificity or therapeutic benefit; age, batch and other confounders are not adjusted."),
        shiny::h4("Molecular subtype result"),
        shiny::p(if(is.null(x$subtypes)) "Unavailable" else sprintf("%d subtypes; overall p = %s",nrow(x$subtypes$summaries),format(x$subtypes$overall$p_value,digits=4))),
        shiny::h4("Clinical results"),
        if(is.null(x$clinical)) shiny::p("Clinical analyses unavailable; expression results remain available.") else shiny::tagList(
          shiny::p(paste("Selected subgroups:",if(length(x$clinical$selection$selected)) paste(x$clinical$selection$selected,collapse=", ") else "None; complete MB cohort only")),
          shiny::p(x$clinical$selection$criteria),
          lapply(c("metastasis","survival"),function(endpoint) shiny::tagList(
            shiny::h5(endpoint),lapply(x$clinical[[endpoint]],function(a) {
              s <- a$statistics
              if(a$status!="complete") return(shiny::p(paste(a$cohort,paste(a$warnings,collapse=" "),sep=": ")))
              if(endpoint=="metastasis") shiny::p(sprintf("%s: metastatic N=%d; M0 N=%d; raw p=%.4g; rank-biserial effect=%.3f.",a$cohort,s$n_metastatic,s$n_m0,s$p_value,s$rank_biserial_metastatic_vs_m0)) else
                shiny::tagList(shiny::p(sprintf("%s: survival N=%d; events=%d; mean cutoff=%.10g; log-rank p=%.4g; HR High/Low=%.3g (95%% CI %.3g-%.3g).",a$cohort,s$n,s$events,s$cutoff,s$logrank_p,s$hazard_ratio_high_vs_low,s$ci_lower,s$ci_upper)),
                  lapply(a$warnings[-1],shiny::p))
            }))),shiny::p(class="alert alert-warning",x$clinical$warnings)),
        shiny::p("Functional Biology: Not yet implemented | DepMap: Not yet implemented"),
        shiny::h4("Available metadata fields"), shiny::p(paste(metadata, collapse = ", ")),
        shiny::p(class = "small text-muted", paste("Retrieved:", x$retrieval$provenance$retrieved_at,
          if (isTRUE(x$retrieval$provenance$cache_hit)) "- Existing R2 cache reused" else "- Retrieved from R2"))
      )
    })
    output$plot_heading <- shiny::renderUI({
      shiny::req(result()); shiny::h3(paste(result()$gene, "expression across MB subgroups"))
    })
    output$subgroup_plot <- shiny::renderPlot({
      shiny::req(result()); result()$analysis$plot
    }, res = 120)
    output$statistics <- shiny::renderTable({
      shiny::req(result()); s <- result()$analysis$statistics
      data.frame(Comparison = s$test, N = s$n, `Effect size` = signif(s$epsilon_squared, 4),
        Statistic = signif(s$statistic, 5), `Raw p-value` = format(s$p_value, scientific = TRUE, digits = 4),
        `Adjusted p-value` = "Not applicable (one test)", check.names = FALSE)
    }, striped = TRUE, bordered = FALSE, spacing = "m", digits = 4)
    output$summaries <- shiny::renderTable({
      shiny::req(result()); result()$analysis$summaries
    }, striped = TRUE, digits = 3)
    output$subtype_summary <- shiny::renderUI({
      x <- result()
      if (is.null(x)) return(shiny::p("Run an analysis to view molecular subtypes."))
      if (is.null(x$subtypes)) return(shiny::div(class="alert alert-warning", x$subtype_failure))
      a <- x$subtypes
      shiny::tagList(shiny::h3(paste(x$gene,"molecular subtype expression")),
        shiny::p(a$metadata$classification),
        shiny::p(sprintf("Metadata field: %s | Annotated: %d | Missing subtype: %d | Subtypes: %d | Analyzed: %d | Excluded: %d",
          a$metadata$field,a$metadata$annotated_n,a$metadata$missing_n,nrow(a$summaries),a$overall$n,a$missing$excluded_n)),
        shiny::p(a$diagnostics$rationale),
        lapply(a$warnings, function(w) shiny::div(class="alert alert-warning",w)))
    })
    output$subtype_plot <- shiny::renderPlot({
      shiny::req(result()$subtypes); result()$subtypes$plot
    }, res=120)
    output$subtype_descriptives <- shiny::renderTable({
      shiny::req(result()$subtypes); result()$subtypes$summaries
    }, striped=TRUE, digits=3)
    # Format only at presentation time; CSVs retain full numeric precision.
    subtype_table <- function(x) {
      for (field in intersect(c("p_value","adjusted_p_value"), names(x)))
        x[[field]] <- ifelse(is.na(x[[field]]), "Not applicable", format(x[[field]],scientific=TRUE,digits=4))
      x
    }
    output$subtype_overall <- shiny::renderTable({
      shiny::req(result()$subtypes); subtype_table(result()$subtypes$overall)
    }, striped=TRUE, digits=4)
    output$subtype_pairs <- shiny::renderTable({
      shiny::req(result()$subtypes); subtype_table(result()$subtypes$pairwise)
    }, striped=TRUE, digits=4)
    output$downloads <- shiny::renderUI({
      if (is.null(result())) return(shiny::p("Run an analysis to enable downloads."))
      shiny::tagList(shiny::h3(paste("Download", result()$gene, "results")),
        shiny::p("Exports contain the same patient data, statistics and figure shown in this analysis."),
        shiny::div(class = "d-flex flex-wrap gap-3",
          shiny::downloadButton("patients", "Patient expression CSV"),
          shiny::downloadButton("stats_csv", "Subgroup statistics CSV"),
          shiny::downloadButton("figure_pdf", "Figure PDF"),
          shiny::downloadButton("figure_png", "Figure PNG")),
        if (!is.null(result()$subtypes)) shiny::tagList(shiny::h4("Molecular subtypes"),
          shiny::div(class="d-flex flex-wrap gap-3",
            shiny::downloadButton("subtype_patients", "Subtype patient CSV"),
            shiny::downloadButton("subtype_summaries", "Subtype descriptive CSV"),
            shiny::downloadButton("subtype_pairwise", "Subtype pairwise CSV"),
            shiny::downloadButton("subtype_overall_csv", "Subtype overall test CSV"),
            shiny::downloadButton("subtype_pdf", "Subtype PDF"),
            shiny::downloadButton("subtype_png", "Subtype PNG (300 dpi)"))),
        shiny::uiOutput("clinical_downloads"))
    })
    downloads <- list(patients = c("patients", "csv", "text/csv"),
      stats_csv = c("statistics", "csv", "text/csv"),
      figure_pdf = c("pdf", "pdf", "application/pdf"),
      figure_png = c("png", "png", "image/png"),
      subtype_patients=c("subtype_patients","csv","text/csv"),
      subtype_summaries=c("subtype_summaries","csv","text/csv"),
      subtype_pairwise=c("subtype_pairwise","csv","text/csv"),
      subtype_overall_csv=c("subtype_overall","csv","text/csv"),
      subtype_pdf=c("subtype_pdf","pdf","application/pdf"),
      subtype_png=c("subtype_png","png","image/png"))
    for (id in names(downloads)) local({
      spec <- downloads[[id]]
      output[[id]] <- shiny::downloadHandler(
        filename = function() { shiny::req(result()); paste0(result()$gene, "_", spec[1], ".", spec[2]) },
        contentType = spec[3],
        content = function(file) { shiny::req(result()); write_explorer_download(result(), spec[1], file) })
    })
    output$selection_summary <- shiny::renderUI({
      shiny::req(result()$clinical); s <- result()$clinical$selection
      shiny::tagList(shiny::h4("Clinical cohort selection"),shiny::p(s$criteria),
        shiny::p(paste("Selected:",if(length(s$selected)) paste(s$selected,collapse=", ") else "None")))
    })
    output$selection_pairs <- shiny::renderTable({
      shiny::req(result()$clinical); result()$clinical$selection$pairwise
    },digits=6,striped=TRUE)
    for(endpoint in c("metastasis","survival")) local({
      ep <- endpoint
      output[[paste0(ep,"_content")]] <- shiny::renderUI({
        x <- result()$clinical
        if(is.null(x)) return(shiny::p("Clinical results are not available. Run an analysis or retry retrieval."))
        v <- x$verification[[ep]]
        shiny::tagList(shiny::h3(paste(result()$gene,ep)),shiny::h4("Metadata verification"),
          shiny::p(v$evidence),
          if(ep=="metastasis") shiny::p(sprintf("Field: %s | Type: %s | Missing: %d | Annotated: %d",v$field,v$source_type,v$missing_n,v$usable_n)) else
            shiny::p(sprintf("Time: %s | Event: %s | Units: %s | Missing event: %d | Missing time: %d | Verified usable: %d",v$time_field,v$event_field,v$units,v$missing_event_n,v$missing_time_n,v$usable_n)),
          shiny::tableOutput(paste0(ep,"_codes")),
          shiny::p(class="alert alert-warning",x$warnings),
          lapply(names(x[[ep]]),function(cohort) {
            a <- x[[ep]][[cohort]]; id <- paste(ep,cohort,sep="_")
            shiny::tagList(shiny::h4(if(cohort=="all") "Complete MB cohort" else paste("Subgroup:",cohort)),
              lapply(a$warnings,function(w)shiny::p(w)),
              shiny::div(style="overflow-x:auto",shiny::tableOutput(paste0(id,"_counts"))),
              shiny::div(style="overflow-x:auto",shiny::tableOutput(paste0(id,"_stats"))),
              if(!is.null(a$plot)) shiny::plotOutput(paste0(id,"_plot"),height="530px"),
              shiny::div(style="overflow-x:auto",shiny::tableOutput(paste0(id,"_details"))))
          }))
      })
      output[[paste0(ep,"_codes")]] <- shiny::renderTable({
        shiny::req(result()$clinical)
        v <- result()$clinical$verification[[ep]]
        if(ep=="metastasis") v$values else v$event_values
      })
      for(cohort in c("all","wnt","shh","group3","group4")) local({
        co <- cohort; id <- paste(ep,co,sep="_")
        get_analysis <- function() {shiny::req(result()$clinical[[ep]][[co]]);result()$clinical[[ep]][[co]]}
        output[[paste0(id,"_plot")]] <- shiny::renderPlot({get_analysis()$plot},res=120)
        output[[paste0(id,"_counts")]] <- shiny::renderTable({get_analysis()$missing})
        output[[paste0(id,"_stats")]] <- shiny::renderTable({
          d <- get_analysis()$statistics
          d <- d[,setdiff(names(d),c("warnings","coding_evidence")),drop=FALSE]
          for(f in names(d)[vapply(d,is.numeric,logical(1))]) d[[f]] <- format(d[[f]],digits=16,trim=TRUE)
          d
        },striped=TRUE)
        output[[paste0(id,"_details")]] <- shiny::renderTable({
          a <- get_analysis(); if(ep=="metastasis") a$summaries else a$risk_table
        },digits=4,striped=TRUE)
        for(kind in c("patient_data","statistics","plot_pdf","plot_png")) local({
          k <- kind; ext <- if(k=="plot_pdf") "pdf" else if(k=="plot_png") "png" else "csv"
          output[[paste0(id,"_download_",k)]] <- shiny::downloadHandler(
            filename=function()paste0(result()$gene,if(co=="all") "" else paste0("_",co),"_",ep,"_",sub("_pdf$|_png$","",k),".",ext),
            content=function(file)core$write_r2_clinical_download(get_analysis(),k,file))
        })
      })
    })
    output$clinical_downloads <- shiny::renderUI({
      shiny::req(result()$clinical)
      shiny::tagList(shiny::h4("Clinical analyses"),lapply(c("metastasis","survival"),function(ep)
        lapply(names(result()$clinical[[ep]]),function(co) {
          a <- result()$clinical[[ep]][[co]]
          if(is.null(a$data)) return(NULL)
          kinds <- c("patient_data","statistics",if(!is.null(a$plot)) c("plot_pdf","plot_png"))
          shiny::tagList(shiny::h5(paste(ep,co)),shiny::div(class="d-flex flex-wrap gap-3",
            lapply(kinds,function(k)shiny::downloadButton(paste(ep,co,"download",k,sep="_"),gsub("_"," ",k)))))
        })))
    })
  }
}
