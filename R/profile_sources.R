get_pfister_expression <- function(cfg,gene,reporter=NULL,refresh=FALSE) {
  pcfg <- cfg;pcfg$r2_dataset_table <- 'ps_avgpres_pfisterb272_informp3'
  pcfg$r2_transformation <- 'transform_none';pcfg$r2_grouping_track <- 'cancer_type'
  pcfg$output_dir <- file.path(cfg$output_dir %||% file.path(dirname(cfg$cache_dir),'..','outputs'),'pfister')
  r <- get_r2_expression(pcfg,gene,reporter,refresh)
  if(nrow(r$data)!=272L || !grepl('Pfister',r$provenance$dataset_label,fixed=TRUE)) stop('Pfister dataset identity/count mismatch')
  d <- r$data; a <- r$metadata[match(d$sample_id,r$metadata$samplenames),]
  require_columns(a,c('cancer_type','type'))
  if(any(d$expression<0,na.rm=TRUE)) stop('Negative source FPKM values')
  d$fpkm <- d$expression;d$expression<-log2(1+d$fpkm);d$cancer_type<-a$cancer_type;d$sample_type<-a$type
  cns <- c('atrt','epd_it','etmr','hgg_k27m','hggother','mb_group3','mb_group4','mb_shh','mb_wnt','pa')
  known <- c(cns,'ews','nb','os','rms','t-all')
  if(any(!a$cancer_type %in% known)) stop('Unrecognised Pfister cancer label; update documented CNS mapping before filtering')
  d$include_cns <- a$cancer_type %in% cns & a$type=='primary'
  d$exclusion_reason <- ifelse(!a$cancer_type %in% cns,'Non-CNS cancer label',ifelse(a$type!='primary','Non-primary sample; avoid paired relapse inference',''))
  r$data <- d;r$provenance$analysis_transform <- 'log2(1 + FPKM), applied locally to R2 transform_none values'
  r$provenance$cns_labels <- as.list(cns);r$provenance$selection <- 'Explicit cancer_type CNS label whitelist and type=primary; full 272 samples retained in source files'
  r
}

verify_r2_survival_payload <- function(html,js,retrieval) {
  # R2's own renderer decreases the curve by status / nrRemaining; status=1 is an event.
  if(!grepl('proportion *= (nrRemaining - status) / nrRemaining;',js,fixed=TRUE)) stop('R2 event semantics changed; re-audit kaplan.js')
  p <- r2_script_payload(html,'/d3/plots/kaplan.js')$plotData
  if(!grepl('month',p$xLabel,ignore.case=TRUE)||!grepl('overall survival',p$yLabel,ignore.case=TRUE)) stop('Survival endpoint/time unit not confirmed')
  d <- p$data;require_columns(d,c('id','status','xValue'))
  base <- retrieval$data
  # R2 Kaplan IDs lowercase GEO accessions; case normalisation only, never positional joining.
  keys <- toupper(d$id);reference <- toupper(base$sample_id)
  if(anyDuplicated(keys)||anyDuplicated(reference)||any(!keys %in% reference)) stop('Survival IDs are duplicate or unmatched')
  i <- match(keys,reference)
  if(any(!d$status %in% c(0,1))||any(d$status!=r2_numeric(base$survival_status_raw[i],'dead'))) stop('Survival event flags disagree with source dead metadata')
  complete <- is.finite(base$survival_time_years)&!r2_missing(base$survival_status_raw)
  if(!setequal(reference[complete],keys)) stop('Survival endpoint samples differ from complete annotated cases')
  # The annotation rounds years, while R2 Kaplan carries more precise follow-up in months.
  if(any(!is.finite(d$xValue))||any(d$xValue<0)||any(abs(d$xValue/12-base$survival_time_years[i])>.11)) stop('Survival times disagree beyond annotation rounding')
  list(data=data.frame(sample_id=base$sample_id[i],time=d$xValue/12,event=as.integer(d$status),source_months=d$xValue,source_years_rounded=base$survival_time_years[i]),
    provenance=list(endpoint='Overall survival',time='R2 Kaplan follow-up months / 12, checked against rounded os_(years)',event='status 1 decreases survival in R2 kaplan.js; 0 censored; per-patient status agrees with dead metadata',
      identifier_join='Unique case-normalised GEO accessions',r2_statistics=p$statistics))
}

get_cavalli_survival <- function(cfg,retrieval,refresh=FALSE) {
  key <- cache_key('cavalli_survival_v2',list(dataset=cfg$r2_dataset_table,reporter=retrieval$provenance$reporter))
  cached<-if(!refresh)cache_get(cfg,key) else NULL
  if(!is.null(cached))return(cached)
  h<-r2_session(120);landing<-r2_request(h,cfg$r2_base_url)
  if(grepl('Use R2 without an account',landing,fixed=TRUE))r2_request(h,cfg$r2_base_url,list(method='guest',open_page='auth'))
  html<-r2_request(h,cfg$r2_base_url,list(table=cfg$r2_dataset_table,option='kaplanscan',inputprobeset=retrieval$provenance$reporter,scanmodus='median',survival='overall',mingrpsize='8',subset=''))
  doc<-xml2::read_html(charToRaw(html),encoding='UTF-8');scripts<-xml2::xml_text(xml2::xml_find_all(doc,'.//script'))
  line<-scripts[grepl('/d3/plots/kaplan.js',scripts,fixed=TRUE)]
  if(length(line)!=1)stop('R2 Kaplan source module unavailable')
  url<-sub(".*import\\('([^']+)'.*",'\\1',line)
  if(!grepl('^https://hgserver1[.]amc[.]nl/r2/assets/[^/]+/js/src/d3/plots/kaplan[.]js$',url))stop('Unexpected Kaplan renderer URL')
  js<-r2_request(h,url); result<-verify_r2_survival_payload(html,js,retrieval)
  folder<-file.path(retrieval$provenance$output_dir,'survival_source');dir.create(folder,showWarnings=FALSE)
  writeBin(charToRaw(html),file.path(folder,'r2-kaplan.html'));writeBin(charToRaw(js),file.path(folder,'kaplan.js'))
  data.table::fwrite(result$data,file.path(folder,'verified_survival.csv'))
  result$provenance$source_dir<-folder;result$provenance$retrieved_at<-utc_now();result$provenance$renderer_url<-url
  cache_put(cfg,key,result);result
}
