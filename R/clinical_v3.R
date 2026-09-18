# V3 orchestration reuses V1 data, rank comparisons, and verified R2 survival.
# No clinical outcome is used to select subgroups or expression thresholds.
select_r2_clinical_subgroups <- function(retrieval, subgroup_analysis) {
  comparisons <- profile_group_comparison(retrieval$data, "subgroup", "Broad subgroup associations")$pairwise
  groups <- sort(unique(comparisons$group1))
  groups <- sort(unique(c(groups, comparisons$group2)))
  evidence <- do.call(rbind,lapply(groups,function(g) {
    p <- comparisons[comparisons$group1==g | comparisons$group2==g,,drop=FALSE]
    forward <- p$group1==g
    effect <- ifelse(forward,p$rank_biserial,-p$rank_biserial)
    delta <- ifelse(forward,p$median_difference,-p$median_difference)
    supported <- is.finite(p$adjusted_p) & p$adjusted_p<.05 & abs(effect)>=.30 &
      p$n1>=10 & p$n2>=10 & sign(effect)==sign(delta)
    consistent <- all(effect>0) || all(effect<0)
    selected <- length(groups)==4 && nrow(p)==3 && all(supported) && consistent &&
      is.finite(subgroup_analysis$statistics$p_value) && subgroup_analysis$statistics$p_value<.05
    data.frame(subgroup=g,selected=selected,direction=if(consistent && all(effect>0)) "Higher" else if(consistent) "Lower" else "Mixed",
      supported_comparisons=sum(supported),comparisons=nrow(p),minimum_absolute_effect=min(abs(effect)),
      maximum_adjusted_p=max(p$adjusted_p))
  }))
  list(selected=evidence$subgroup[evidence$selected], evidence=evidence, pairwise=comparisons,
    criteria="Overall subgroup p < 0.05; all three pairwise BH p < 0.05, |rank-biserial| >= 0.30, n >= 10 per group, consistent rank and median direction versus every other subgroup. Exploratory cohort selection, not evidence of a clinical association.")
}

clinical_v3_data <- function(retrieval, cohort="all") {
  d <- as.data.frame(retrieval$data); m <- retrieval$metadata
  if(anyDuplicated(d$sample_id) || any(r2_missing(d$sample_id)) || anyDuplicated(m$samplenames) ||
      !setequal(d$sample_id,m$samplenames)) stop("Invalid clinical metadata identifiers")
  m <- m[match(d$sample_id,m$samplenames),,drop=FALSE]
  for(f in setdiff(names(m),names(d))) d[[f]] <- m[[f]]
  if(cohort!="all") d <- d[!r2_missing(d$subgroup) & d$subgroup==cohort,,drop=FALSE]
  d$cohort <- cohort
  d
}

verify_r2_clinical_metadata <- function(retrieval, survival_source=NULL) {
  m <- retrieval$metadata
  field <- "met_status_(1_met__0_m0)"
  raw <- if(field %in% names(m)) as.character(m[[field]]) else rep(NA_character_,nrow(m))
  verified <- field %in% names(m) && isTRUE(retrieval$audit$metastasis$coding_verified) &&
    all(raw[!r2_missing(raw)] %in% c("0","1"))
  met <- list(verified=verified, field=field, values=as.data.frame(table(raw,useNA="ifany")),
    source_type=if(field %in% names(m)) class(m[[field]])[1] else "absent",
    missing_n=sum(r2_missing(raw)),usable_n=sum(!r2_missing(raw)),
    evidence=if(!field %in% names(m)) "Metastasis metadata are absent." else if(verified) "R2 annotation and track field explicitly states met_status_(1_met__0_m0): 1=metastatic, 0=M0. No additional metastatic-stage field is present." else
      "Metastasis metadata were found, but the biological meaning of the coding could not be reliably verified.")
  surv <- list(verified=FALSE,time_field="os_(years)",event_field="dead",units="Unresolved",
    evidence="Survival event coding or units could not be verified; no survival analysis performed.")
  if(!all(c("os_(years)","dead") %in% names(m))) surv$evidence <- "Required survival metadata are absent; no survival analysis performed."
  if(!is.null(survival_source)) {
    folder <- survival_source$provenance$source_dir
    # Revalidate saved authoritative source, rather than trusting a cached boolean.
    checked <- tryCatch(verify_r2_survival_payload(
      paste(readLines(file.path(folder,"r2-kaplan.html"),warn=FALSE),collapse="\n"),
      paste(readLines(file.path(folder,"kaplan.js"),warn=FALSE),collapse="\n"),retrieval),error=function(e) {
        message("[Clinical source verification] ",conditionMessage(e)); NULL
      })
    if(!is.null(checked)) {
      surv$verified <- TRUE; surv$data <- checked$data; surv$units <- "years"
      surv$evidence <- paste(checked$provenance$time,checked$provenance$event,sep="; ")
      surv$source_dir <- folder; surv$renderer_url <- survival_source$provenance$renderer_url
      surv$source_md5 <- tools::md5sum(file.path(folder,c("r2-kaplan.html","kaplan.js")))
    }
  }
  surv$event_values <- if("dead" %in% names(m)) as.data.frame(table(m$dead,useNA="ifany")) else data.frame()
  surv$missing_event_n <- if("dead" %in% names(m)) sum(r2_missing(m$dead)) else nrow(m)
  surv$missing_time_n <- if("os_(years)" %in% names(m)) sum(r2_missing(m[["os_(years)"]])) else nrow(m)
  surv$usable_n <- if(surv$verified) nrow(surv$data) else 0L
  list(metastasis=met,survival=surv)
}

clinical_v3_counts <- function(d, retrieval, clinical_ok, expression_ok) {
  data.frame(total_retrieved=nrow(retrieval$data),cohort_n=nrow(d),usable_expression=sum(expression_ok),
    usable_clinical=sum(clinical_ok),used_n=sum(clinical_ok & expression_ok),excluded_n=sum(!clinical_ok | !expression_ok))
}

analyze_r2_metastasis <- function(retrieval, verification, cohort="all") {
  d <- clinical_v3_data(retrieval,cohort)
  raw <- if(verification$field %in% names(d)) as.character(d[[verification$field]]) else rep(NA_character_,nrow(d))
  d$metastasis_raw <- raw; d$metastasis_label <- NA_character_
  ok <- !r2_missing(raw) & raw %in% c("0","1") & isTRUE(verification$verified)
  expr <- is.finite(d$expression)
  if(verification$verified) d$metastasis_label[ok] <- ifelse(raw[ok]=="1","Metastatic","M0 (non-metastatic)")
  d$included <- ok & expr
  d$exclusion_reason <- ifelse(!ok & !expr,"Missing/invalid clinical data and expression",ifelse(!ok,"Missing or unverified metastasis",ifelse(!expr,"Non-finite expression","")))
  counts <- clinical_v3_counts(d,retrieval,ok,expr)
  out <- list(gene=retrieval$provenance$queried_gene,cohort=cohort,status="unavailable",data=d,
    verification=verification,missing=counts,statistics=data.frame(cohort=cohort,status="unavailable",p_value=NA_real_),
    summaries=data.frame(),plot=NULL,warnings=character())
  if(!verification$verified) {out$warnings <- verification$evidence; return(out)}
  used <- d[d$included,,drop=FALSE]
  out$summaries <- do.call(rbind,lapply(c("M0 (non-metastatic)","Metastatic"),function(g) {
    x <- used$expression[used$metastasis_label==g]
    data.frame(group=g,n=length(x),mean=if(length(x)) mean(x) else NA_real_,median=if(length(x)) median(x) else NA_real_,
      sd=if(length(x)>1) sd(x) else NA_real_,iqr=if(length(x)) IQR(x) else NA_real_,
      variance=if(length(x)>1) var(x) else NA_real_)
  }))
  out$plot <- ggplot2::ggplot(used,ggplot2::aes(metastasis_label,expression,fill=metastasis_label))+
    ggplot2::geom_boxplot(outlier.shape=NA,alpha=.4)+ggplot2::geom_point(position=ggplot2::position_jitter(width=.15,seed=85217),alpha=.5,size=1)+
    ggplot2::scale_x_discrete(labels=setNames(paste0(out$summaries$group,"\n(n=",out$summaries$n,")"),out$summaries$group))+
    ggplot2::theme_classic(base_size=12)+ggplot2::guides(fill="none")+
    ggplot2::labs(title=paste(out$gene,"metastasis:",cohort),x="Verified metastasis status",y="Expression (R2 log2)")
  if(any(out$summaries$n<5) || length(unique(used$expression))<2) {
    out$warnings <- "Insufficient variable expression or fewer than five patients per metastasis group; descriptive results only."
    return(out)
  }
  # Reuse the established pairwise rank framework, following size/spread checks.
  a <- profile_group_comparison(used,"metastasis_label","Metastasis")$pairwise
  x <- used$expression[used$metastasis_label=="Metastatic"]; y <- used$expression[used$metastasis_label=="M0 (non-metastatic)"]
  # Orient effects explicitly as metastatic minus M0 regardless of alphabetical order.
  direction <- if(a$group1=="Metastatic") 1 else -1
  out$statistics <- data.frame(cohort=cohort,status="complete",test="Two-sided Wilcoxon rank-sum; tie-adjusted normal approximation",
    statistic=(direction*a$rank_biserial+1)*length(x)*length(y)/2,p_value=a$p_value,
    rank_biserial_metastatic_vs_m0=direction*a$rank_biserial,median_difference=median(x)-median(y),
    n_metastatic=length(x),n_m0=length(y))
  out$status <- "complete"
  out$warnings <- "Rank-based distribution comparison after inspecting group sizes, variance and missingness; not necessarily a median-only effect. No age/batch adjustment."
  out
}

analyze_r2_survival <- function(retrieval, verification, cohort="all", cutoff="mean") {
  if(!identical(cutoff,"mean")) stop("Version 3 supports only the predefined mean cutoff")
  d <- clinical_v3_data(retrieval,cohort)
  d$time <- d$event <- NA_real_
  if(isTRUE(verification$verified)) {
    i <- match(d$sample_id,verification$data$sample_id)
    d$time <- verification$data$time[i]; d$event <- verification$data$event[i]
  }
  expr <- is.finite(d$expression)
  ok <- is.finite(d$time) & d$time>=0 & !is.na(d$event) & d$event %in% c(0,1)
  value <- if(any(expr)) mean(d$expression[expr]) else NA_real_
  d$expression_group <- ifelse(!expr,NA_character_,ifelse(d$expression>value,"High","Low"))
  d$included <- ok & expr
  d$exclusion_reason <- ifelse(!ok & !expr,"Missing/invalid clinical data and expression",ifelse(!ok,"Missing or unverified survival time/event",ifelse(!expr,"Non-finite expression","")))
  counts <- clinical_v3_counts(d,retrieval,ok,expr)
  used <- d[d$included,,drop=FALSE]
  out <- list(gene=retrieval$provenance$queried_gene,cohort=cohort,status="unavailable",data=d,missing=counts,
    verification=verification,cutoff=list(method="mean",value=value,reference_n=sum(expr),candidates=1L),
    statistics=data.frame(cohort=cohort,status="unavailable",cutoff_method="mean",cutoff=value,
      cutoff_reference_n=sum(expr),n=nrow(used),high_n=sum(used$expression_group=="High"),low_n=sum(used$expression_group=="Low"),
      events=sum(used$event),censored=sum(used$event==0),logrank_p=NA_real_,cox_p=NA_real_),
    plot=NULL,km_fit=NULL,cox_fit=NULL,logrank_test=NULL,ph_diagnostic=NULL,risk_table=data.frame(),warnings=character())
  if(!isTRUE(verification$verified)) {out$warnings <- verification$evidence; return(out)}
  sizes <- table(factor(used$expression_group,levels=c("Low","High")))
  if(nrow(used)<20 || any(sizes<10) || sum(used$event)<5 || length(unique(used$time))<2) {
    out$warnings <- "Insufficient data for reliable survival analysis: require 20 patients, 10 per expression group, five events and variable follow-up."
    return(out)
  }
  events <- vapply(c("Low","High"),function(g)sum(used$event[used$expression_group==g]),numeric(1))
  fit_cox <- sum(events)>=10 && all(events>=3)
  fit <- profile_survival(d,out$cutoff,cohort,fit_cox=fit_cox)
  out$status <- "complete"; out$plot <- fit$plot; out$km_fit <- fit$km_fit
  out$logrank_test <- fit$logrank_test; out$cox_fit <- fit$cox_fit; out$risk_table <- fit$risk_table
  out$statistics <- cbind(out$statistics[,!names(out$statistics) %in% names(fit$statistics),drop=FALSE],fit$statistics)
  out$statistics$status <- "complete"; out$statistics$logrank_statistic <- fit$logrank_test$chisq
  out$statistics$logrank_df <- 1L
  out$warnings <- "High > mean; Low <= mean. Mean uses all finite expression values in this cohort before clinical exclusions. Unadjusted exploratory association; independent/non-informative censoring assumed."
  if(!fit_cox) out$warnings <- c(out$warnings,"Cox model withheld: require ten events overall and three events per expression group.")
  if(nzchar(fit$statistics$cox_note)) out$warnings <- c(out$warnings,"Cox model was unstable; hazard ratio is withheld.")
  if(!is.null(fit$cox_fit) && !nzchar(fit$statistics$cox_note)) {
    cf <- summary(fit$cox_fit)$coefficients
    out$statistics$cox_p <- cf[1,"Pr(>|z|)"]
    out$ph_diagnostic <- tryCatch(survival::cox.zph(fit$cox_fit),error=function(e) NULL)
    if(is.null(out$ph_diagnostic)) out$warnings <- c(out$warnings,"Proportional-hazards diagnostic could not be computed.")
    else if(out$ph_diagnostic$table[1,"p"]<.05) out$warnings <- c(out$warnings,"Proportional-hazards assumption may be violated (Schoenfeld p < 0.05); a single HR may not describe the time-varying association.")
    if(any(!is.finite(unlist(out$statistics[c("hazard_ratio_high_vs_low","ci_lower","ci_upper","cox_p")])))) {
      out$statistics[c("hazard_ratio_high_vs_low","ci_lower","ci_upper","cox_p")] <- NA_real_
      out$warnings <- c(out$warnings,"Non-finite Cox estimate or interval: hazard ratio withheld.")
    }
  }
  out$plot <- out$plot + ggplot2::labs(title=paste(out$gene,"overall survival:",cohort),
    subtitle=sprintf("Mean cutoff %.10g; N=%d; log-rank p=%.3g",value,nrow(used),out$statistics$logrank_p),
    caption=paste0("High > mean; Low <= mean. HR High/Low = ",format(out$statistics$hazard_ratio_high_vs_low,digits=4),
      " (95% CI ",format(out$statistics$ci_lower,digits=4)," - ",format(out$statistics$ci_upper,digits=4),").\n",
      "Dotted: pointwise 95% CI; + censoring. Unadjusted exploratory association; see PH diagnostic."))
  out
}

analyze_r2_clinical_v3 <- function(cfg,retrieval,subgroup_analysis) {
  selection <- select_r2_clinical_subgroups(retrieval,subgroup_analysis)
  source <- tryCatch(get_cavalli_survival(cfg,retrieval),error=function(e) {
    message("[Clinical survival retrieval] ",conditionMessage(e)); NULL
  })
  verification <- verify_r2_clinical_metadata(retrieval,source)
  cohorts <- c("all",selection$selected)
  safe <- function(fun,ver,cohort) tryCatch(fun(retrieval,ver,cohort),error=function(e) {
    message("[Clinical analysis] ",cohort,": ",conditionMessage(e))
    list(status="unavailable",cohort=cohort,warnings="Clinical analysis could not be completed. Earlier expression results remain available.",plot=NULL)
  })
  metastasis <- setNames(lapply(cohorts,function(g)safe(analyze_r2_metastasis,verification$metastasis,g)),cohorts)
  survival <- setNames(lapply(cohorts,function(g)safe(analyze_r2_survival,verification$survival,g)),cohorts)
  # Descriptive multiplicity control within each endpoint; does not fix selection bias.
  adjust <- function(results,column) {
    p <- vapply(results,function(x) if(is.null(x$statistics[[column]])) NA_real_ else x$statistics[[column]],numeric(1))
    adj <- p.adjust(p,"BH",n=length(cohorts))
    for(i in seq_along(results)) if(!is.null(results[[i]]$statistics)) {
      results[[i]]$statistics[[paste0(column,"_BH")]] <- adj[i]
      for(f in names(results[[i]]$missing)) results[[i]]$statistics[[f]] <- results[[i]]$missing[[f]]
      results[[i]]$statistics$warnings <- paste(results[[i]]$warnings,collapse="; ")
      results[[i]]$statistics$coding_evidence <- results[[i]]$verification$evidence
    }
    results
  }
  list(selection=selection,verification=verification,metastasis=adjust(metastasis,"p_value"),
    survival=adjust(adjust(survival,"logrank_p"),"cox_p"),
    warnings="Subgroup selection and clinical testing use overlapping patients. Clinical p-values are exploratory; BH across tested cohorts does not correct selection bias or multiple queried genes. Complete MB results may reflect subgroup composition.")
}

write_r2_clinical_download <- function(analysis, kind, file) {
  switch(kind,patient_data=data.table::fwrite(analysis$data,file,na="NA"),
    statistics=data.table::fwrite(analysis$statistics,file,na="NA"),
    plot_pdf=ggplot2::ggsave(file,analysis$plot,device=grDevices::pdf,width=9,height=6),
    plot_png=ggplot2::ggsave(file,analysis$plot,device="png",width=9,height=6,dpi=300),
    stop("Unknown clinical download type"))
  invisible(file)
}
