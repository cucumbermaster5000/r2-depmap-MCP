profile_plot_html <- function(x) {
  if(is.null(x$plot))return(paste0("<p class='notice'>Unavailable: ",html_escape(x$reason %||% 'No eligible observations'),'</p>'))
  paste0("<div class='figure'>",plot_inline_svg(x$plot),'</div>',html_table(x$statistics),
    if(!is.null(x$risk_table))paste0('<details><summary>Number at risk</summary>',html_table(x$risk_table),'</details>'),
    if(!is.null(x$pairwise))paste0('<details><summary>Pairwise comparisons (BH-adjusted)</summary>',html_table(x$pairwise),'</details>'),
    if(!is.null(x$summaries))paste0('<details><summary>Group summaries</summary>',html_table(x$summaries),'</details>'))
}

profile_extra_html <- function(profile,default_cutoff='median') {
  clinical<-profile$clinical
  panels<-paste(vapply(clinical,function(x) {
    paste0("<div class='clinical-panel' data-cohort='",html_escape(x$cohort),"' data-cutoff='",html_escape(x$method),"'",
      if(x$cohort!='all'||x$method!=default_cutoff)' hidden' else '','>',
      '<h3>Overall survival</h3>',profile_plot_html(x$survival),
      if(x$method=='optimal')"<p class='notice'>Optimised in these same patients: the cutoff-adjusted p-value uses Bonferroni correction over all eligible cutoffs. Hazard ratios and confidence intervals are descriptive after selection and may be optimistic. This is not independent validation.</p>",
      '<h3>Metastasis status</h3>',profile_plot_html(x$metastasis),
      '<details><summary>Group counts and cutoff search</summary>',html_table(x$metastasis$counts),html_table(x$cutoff$scan),'</details></div>')
  },character(1)),collapse='')
  cohorts<-unique(vapply(clinical,`[[`,character(1),'cohort'))
  controls<-paste0("<label for='clinical-cohort'>Cohort</label> <select id='clinical-cohort'>",paste0('<option value="',html_escape(cohorts),'">',html_escape(cohorts),'</option>',collapse=''),'</select> ',
    "<label for='clinical-cutoff'>Expression cutoff</label> <select id='clinical-cutoff'>",paste(vapply(c('median','mean','optimal'),function(m)paste0('<option value="',m,'"',if(m==default_cutoff)' selected' else '','>',m,'</option>'),character(1)),collapse=''),'</select>')
  pf<-profile$pfister;en<-profile$enrichr;dep<-profile$depmap
  en_html<-if(identical(en$status,'unavailable'))paste0('<p class="notice">',html_escape(en$reason),'</p>') else {
    sig<-en$enrichment;if(nrow(sig))sig<-sig[sig$adjusted_p<=.05,,drop=FALSE]
    if(nrow(sig))sig<-sig[order(sig$adjusted_p),,drop=FALSE]
    paste0('<p>Known gene-set memberships, not evidence that a pathway is activated, inhibited or deregulated.</p>',html_table(en$membership),
      '<h3>Transcription factor annotation</h3><p>',html_escape(en$tf$classification),' (',html_escape(en$tf$version %||% 'source unavailable'),')</p>',
      '<h3>Binding-associated TF targets</h3>',if(!nrow(en$targets))'<p>No exact human TF-target experiment for this gene was found in ChEA 2022. This does not mean the gene has no targets.</p>' else paste0('<details><summary>Inspect target genes and source experiments</summary>',html_table(en$targets),'</details>'),
      '<h3>TF-target enrichment</h3>',if(!nrow(sig))'<p>No eligible/significant TF-target enrichment results at BH FDR 0.05.</p>' else html_table(head(sig,100)),
      '<p>',html_escape(en$provenance$method),'</p><p>Top 100 significant rows are shown; all tested terms are retained in the R object and CSV files.</p>')
  }
  dep_html<-if(identical(dep$status,'unavailable'))paste0('<p class="notice">',html_escape(dep$reason),'</p>') else paste0('<p>',html_escape(dep$provenance$note),'</p>',
    profile_plot_html(dep),"<div class='figure'>",plot_inline_svg(dep$mb_plot),'</div>',html_table(dep$models),'<p>',html_escape(dep$provenance$model_selection),'. ',html_escape(dep$provenance$interpretation),'</p>')
  paste0("<section class='card' id='profile-summary'><h2>Report coverage</h2>",html_table(profile$status),
    '<p>Separate evidence layers: expression distributions, unadjusted clinical associations, functional annotations, in-vitro dependency and literature. No composite therapeutic score.</p>',
    '<nav><a href="#subtypes">Subtypes</a><a href="#pfister">CNS tumours</a><a href="#clinical">Survival and metastasis</a><a href="#enrichr">Enrichr / TF targets</a><a href="#depmap">DepMap</a></nav></section>',
    "<section class='card' id='subtypes'><h2>Cavalli molecular subtypes</h2>",profile_plot_html(profile$subtypes),
    '<p>All original subtype labels are preserved. Pairwise rank tests compare distributions; they do not prove subtype specificity.</p></section>',
    "<section class='card' id='pfister'><h2>Expression across CNS tumours: Pfister 272</h2>",profile_plot_html(pf),
    '<p>FPKM is transformed locally as log2(1 + FPKM). This panel is separate from Cavalli microarray expression. CNS labels are explicitly listed in provenance; only primary samples enter the comparison.</p>',
    '<details><summary>Source group counts and excluded samples</summary>',html_table(pf$all_counts),html_table(pf$excluded_samples),'</details></section>',
    "<section class='card' id='clinical'><h2>Survival and metastasis</h2>",controls,
    '<p>Mean and median are calculated on all finite log2 expression values within the selected cohort, before clinical missing-data exclusions; the same threshold is used for survival and metastasis. High &gt; cutoff; Low &lt;= cutoff. The optimal threshold uses survival-complete cases. Subgroup cutoffs are recalculated within each subgroup.</p>',
    '<p>Family-adjusted p-values use BH across all five cohorts and all three cutoff methods, separately for survival and metastasis. Survival first corrects the optimised cutoff search. Metastasis at the survival-selected threshold is descriptive only: inferential p-values and confidence intervals are omitted because the endpoints may be correlated. Pooled results are not adjusted for subgroup, age or treatment. Check cohort-specific estimates and the Cox proportional-hazards diagnostic.</p>',panels,
    '<details><summary>Survival source verification</summary><pre>',html_escape(json_text(profile$survival_provenance)),'</pre></details></section>',
    "<section class='card' id='enrichr'><h2>Enrichr functional evidence</h2>",en_html,'</section>',
    "<section class='card' id='depmap'><h2>DepMap dependencies</h2>",dep_html,'</section>',
    '<script>function updateClinical(){const c=document.getElementById("clinical-cohort").value;const m=document.getElementById("clinical-cutoff").value;document.querySelectorAll(".clinical-panel").forEach(x=>{x.hidden=x.dataset.cohort!==c||x.dataset.cutoff!==m;});}document.getElementById("clinical-cohort").addEventListener("change",updateClinical);document.getElementById("clinical-cutoff").addEventListener("change",updateClinical);updateClinical();</script>')
}

save_profile_tables <- function(x,path) {
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  if(is.data.frame(x)) {if(ncol(x)) data.table::fwrite(x,file.path(path,'data.csv'),na='NA');return(invisible(NULL))}
  if(!is.list(x)||inherits(x,'ggplot'))return(invisible(NULL))
  for(n in names(x)) if(n!='plot'&&n!='mb_plot') {
    child<-file.path(path,gsub('[^a-zA-Z0-9_-]','_',n))
    if(is.data.frame(x[[n]])||is.list(x[[n]]))save_profile_tables(x[[n]],child)
  }
}

# Export an already computed result; no retrieval or statistical calculations.
export_gene_report <- function(analysis, cfg) {
  retrieval <- analysis$retrieval
  validation <- analysis$validation
  pubs <- analysis$publications
  pub_error <- analysis$publication_error
  profile <- analysis$profile
  cutoff <- analysis$cutoff
  gene <- retrieval$provenance$queried_gene
  status <- profile$status
  source <- list(provenance = profile$survival_provenance)
  pfister <- profile$pfister
  en <- profile$enrichr
  dep <- profile$depmap
  root<-cfg$output_dir %||% file.path(dirname(cfg$cache_dir),'..','outputs')
  out<-file.path(root,toupper(gene),basename(tempfile(paste0('report-',format(Sys.time(),'%Y%m%dT%H%M%S',tz='UTC'),'-'))))
  dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,winslash='/')
  saveRDS(list(retrieval=retrieval,validation=validation,publications=pubs,publication_error=pub_error,profile=profile,cutoff=cutoff),file.path(out,'report_data.rds'))
  save_profile_tables(profile,file.path(out,'tables'))
  data.table::fwrite(retrieval$data,file.path(out,'Cavalli_patient_data.csv'),na='NA')
  writeLines(json_text(list(generated_at=utc_now(),cavalli=retrieval$provenance,survival=source$provenance,pfister=pfister$provenance,enrichr=en$provenance,depmap=dep$provenance,pubmed=pubs$provenance)),file.path(out,'provenance.json'))
  if(!is.null(pubs))writeLines(json_text(pubs),file.path(out,'pubmed_results.json'))
  path<-write_gene_page(retrieval,validation,pubs,pub_error,file.path(out,paste0(toupper(gene),'_gene_page.html')),profile_extra_html(profile,cutoff))
  # Snapshot every report function and a network-free rendering entrypoint.
  core_env <- environment(export_gene_report)
  functions <- get('.r2_core_functions', envir=core_env, inherits=FALSE)
  dump(functions,file.path(out,'analysis_functions.R'),envir=core_env)
  cat('\n.r2_core_functions <- ', file=file.path(out,'analysis_functions.R'), append=TRUE)
  cat(paste(capture.output(dput(functions)), collapse='\n'), '\n', file=file.path(out,'analysis_functions.R'), append=TRUE)
  writeLines(c("source('analysis_functions.R')","x <- readRDS('report_data.rds')","write_gene_page(x$retrieval, x$validation, x$publications, x$publication_error, 'regenerated_report.html', profile_extra_html(x$profile, x$cutoff))"),file.path(out,'reproduce_report.R'))
  writeLines(capture.output(sessionInfo()),file.path(out,'session-info.txt'))
  list(gene=toupper(gene),html_path=path,samples=nrow(retrieval$data),publication_count=length(pubs$articles),publication_status=if(is.null(pubs))'unavailable' else 'ok',
    modules=status,clinical_cutoff_default=cutoff,available_cutoffs=as.list(c('median','mean','optimal')),report_data=file.path(out,'report_data.rds'),
    instructions='Open html_path in a browser; use cohort and cutoff selectors. Source data, tables, R functions and offline reproduction script are saved alongside the report.')
}

# Stable entry point shared by MCP and the command-line script.
build_gene_page <- function(cfg, gene, reporter = NULL, publication_limit = 10L,
    publication_disease = NULL, aliases = character(), include_aliases = TRUE,
    refresh_publications = FALSE, cutoff = "median", pfister_reporter = NULL) {
  analysis <- analyze_gene(cfg, gene, reporter, publication_limit,
    publication_disease, aliases, include_aliases, refresh_publications,
    cutoff, pfister_reporter)
  export_gene_report(analysis, cfg)
}
