# Transparent exploratory analyses; original labels and explicit exclusions.
profile_group_comparison <- function(data, field, title, y_label = 'Expression (R2 log2)', pairwise = TRUE) {
  d <- data.frame(sample_id=data$sample_id, expression=data$expression, group=as.character(data[[field]]))
  ok <- is.finite(d$expression) & !r2_missing(d$group)
  excluded <- d[!ok,,drop=FALSE]; d <- d[ok,,drop=FALSE]
  sizes <- table(d$group)
  if (length(sizes)<2 || anyDuplicated(d$sample_id)) stop('Need independent samples in at least two groups')
  summaries <- do.call(rbind,lapply(names(sizes),function(g) { x <- d$expression[d$group==g]; data.frame(group=g,n=length(x),mean=mean(x),median=median(x),iqr=IQR(x)) }))
  kw <- kruskal.test(expression~group,d)
  comparisons <- data.frame()
  if(pairwise) {
    pairs <- combn(names(sizes),2,simplify=FALSE)
    comparisons <- do.call(rbind,lapply(pairs,function(g) {
      x <- d$expression[d$group==g[1]]; y <- d$expression[d$group==g[2]]
      test <- suppressWarnings(wilcox.test(x,y,exact=FALSE))
      data.frame(group1=g[1],group2=g[2],n1=length(x),n2=length(y),median_difference=median(x)-median(y),
        rank_biserial=2*as.numeric(test$statistic)/(length(x)*length(y))-1,p_value=test$p.value)
    }))
    comparisons$adjusted_p <- p.adjust(comparisons$p_value,'BH')
  }
  statistics <- data.frame(test='Kruskal-Wallis',n=nrow(d),excluded_n=nrow(excluded),groups=length(sizes),statistic=unname(kw$statistic),p_value=kw$p.value,
    epsilon_squared=max(0,(unname(kw$statistic)-length(sizes)+1)/(nrow(d)-length(sizes))))
  p <- ggplot2::ggplot(d,ggplot2::aes(x=group,y=expression,fill=group)) + ggplot2::geom_boxplot(outlier.shape=NA,alpha=.4) +
    ggplot2::geom_point(position=ggplot2::position_jitter(width=.17,height=0,seed=85217),size=.8,alpha=.5) +
    ggplot2::scale_x_discrete(labels=setNames(paste0(names(sizes),' (n=',sizes,')'),names(sizes))) +
    ggplot2::coord_flip() + ggplot2::guides(fill='none') + ggplot2::theme_classic(base_size=11) +
    ggplot2::labs(title=title,subtitle=sprintf('Kruskal-Wallis p = %.3g; %d samples',kw$p.value,nrow(d)),x='Original source labels',y=y_label,
      caption='Exploratory unadjusted distribution comparison; small groups remain visible. Pairwise tests use BH correction.')
  list(plot=p,statistics=statistics,summaries=summaries,pairwise=comparisons,exclusions=excluded,data=d)
}

profile_cutoff <- function(expression, method, survival_data=NULL) {
  method <- match.arg(method,c('median','mean','optimal'))
  x <- expression[is.finite(expression)]
  if(length(x)<10 || length(unique(x))<2) stop('Insufficient variable expression for a high/low split')
  if(method!='optimal') return(list(value=if(method=='median') median(x) else mean(x),method=method,candidates=1L,scan=data.frame(),reference_n=length(x)))
  if(is.null(survival_data)) stop('Optimised cutoff requires verified survival data')
  d <- survival_data
  values <- sort(unique(d$expression)); cuts <- head(values,-1) # high > cutoff, low <= cutoff; ties never split
  minimum <- max(10L,ceiling(.1*nrow(d)))
  rows <- lapply(cuts,function(cut) {
    high <- d$expression>cut
    if(min(sum(high),sum(!high))<minimum) return(NULL)
    test <- tryCatch(survival::survdiff(survival::Surv(time,event)~high,data=d),error=function(e) NULL)
    if(is.null(test)||!is.finite(test$chisq)) return(NULL)
    data.frame(cutoff=cut,low_n=sum(!high),high_n=sum(high),p_value=pchisq(test$chisq,1,lower.tail=FALSE))
  })
  scan <- do.call(rbind,rows)
  if(is.null(scan)||!nrow(scan)) stop('No eligible optimal cutoffs (minimum 10 patients and 10% per group)')
  scan$adjusted_p <- p.adjust(scan$p_value,'bonferroni')
  best <- which.min(scan$p_value)
  list(value=scan$cutoff[best],method=method,candidates=nrow(scan),scan=scan,reference_n=nrow(d))
}

profile_survival <- function(data, cutoff, cohort, fit_cox = TRUE) {
  d <- data
  eligible <- is.finite(d$expression)&is.finite(d$time)&d$time>=0&!is.na(d$event)&d$event %in% c(0,1)
  excluded <- d[!eligible,,drop=FALSE]; d <- d[eligible,,drop=FALSE]
  if(nrow(d)<10 || sum(d$event)<5) stop('Survival needs >=10 complete samples and >=5 deaths')
  d$group <- factor(ifelse(d$expression>cutoff$value,'High','Low'),levels=c('Low','High'))
  if(any(table(d$group)<5)) stop('Survival split needs >=5 samples in each expression group')
  fit <- survival::survfit(survival::Surv(time,event)~group,data=d,conf.type='log-log')
  test <- survival::survdiff(survival::Surv(time,event)~group,data=d)
  rawp <- pchisq(test$chisq,1,lower.tail=FALSE)
  notes <- character()
  cox <- if (!fit_cox) NULL else withCallingHandlers(tryCatch(survival::coxph(survival::Surv(time,event)~group,data=d,x=TRUE),error=function(e){notes<<-c(notes,conditionMessage(e));NULL}),
    warning=function(w){notes<<-c(notes,conditionMessage(w));invokeRestart('muffleWarning')})
  hr <- lo <- hi <- ph <- NA_real_
  if(!is.null(cox)&&!length(notes)) {
    cf <- summary(cox)$conf.int; hr<-cf[1,'exp(coef)'];lo<-cf[1,'lower .95'];hi<-cf[1,'upper .95']
    ph <- tryCatch(survival::cox.zph(cox)$table[1,'p'],error=function(e) NA_real_)
  }
  curves <- data.frame(time=fit$time,survival=fit$surv,lower=fit$lower,upper=fit$upper,censored=fit$n.censor,
    group=rep(sub('group=','',names(fit$strata),fixed=TRUE),as.integer(fit$strata)))
  curves <- rbind(data.frame(time=0,survival=1,lower=1,upper=1,censored=0,group=c('Low','High')),curves)
  curves <- curves[order(curves$group,curves$time),]
  times <- pretty(c(0,max(d$time)),n=6);times<-times[times>=0&times<=max(d$time)]
  risk <- do.call(rbind,lapply(c('Low','High'),function(g) data.frame(group=g,time_years=times,n_at_risk=vapply(times,function(t)sum(d$group==g&d$time>=t),integer(1)))))
  statistics <- data.frame(cohort=cohort,cutoff_method=cutoff$method,cutoff=cutoff$value,n=nrow(d),deaths=sum(d$event),low_n=sum(d$group=='Low'),high_n=sum(d$group=='High'),
    excluded_n=nrow(excluded),logrank_p=rawp,cutoff_adjusted_p=min(1,rawp*cutoff$candidates),candidate_cutoffs=cutoff$candidates,
    hazard_ratio_high_vs_low=hr,ci_lower=lo,ci_upper=hi,ph_test_p=ph,cox_note=paste(notes,collapse='; '))
  p <- ggplot2::ggplot(curves,ggplot2::aes(time,survival,color=group,group=group)) + ggplot2::geom_step(linewidth=.8)+
    ggplot2::geom_step(ggplot2::aes(y=lower),linetype=3,alpha=.4,na.rm=TRUE)+ggplot2::geom_step(ggplot2::aes(y=upper),linetype=3,alpha=.4,na.rm=TRUE)+
    ggplot2::geom_point(data=curves[curves$censored>0,],shape=3,size=1.8)+ggplot2::scale_color_manual(values=c(High='#b43b47',Low='#236fa1'))+
    ggplot2::coord_cartesian(ylim=c(0,1))+ggplot2::theme_classic(base_size=12)+
    ggplot2::labs(title=paste('Overall survival:',cohort),subtitle=sprintf('%s cutoff %.4g; low %d / high %d; log-rank p %.3g',cutoff$method,cutoff$value,statistics$low_n,statistics$high_n,rawp),
      x='Follow-up (years)',y='Survival probability',color='Expression',caption=paste0('High > cutoff; Low <= cutoff. Dashed: pointwise 95% CI where estimable; + censoring.\n',if(cutoff$method=='optimal')sprintf('Cutoff-search Bonferroni p = %.3g; selected estimates are exploratory.',statistics$cutoff_adjusted_p) else 'Unadjusted exploratory association.'))
  list(plot=p,statistics=statistics,risk_table=risk,curves=curves,data=d,exclusions=excluded,cutoff=cutoff,
    km_fit=fit,logrank_test=test,cox_fit=cox)
}

profile_metastasis <- function(data,cutoff,cohort) {
  if(!all(na.omit(data$metastasis_raw[!r2_missing(data$metastasis_raw)]) %in% c('0','1'))) stop('Unexpected metastasis codes')
  eligible <- is.finite(data$expression)&!r2_missing(data$metastasis_raw)
  d <- data[eligible,,drop=FALSE]
  d$group <- factor(ifelse(d$expression>cutoff$value,'High','Low'),levels=c('Low','High'))
  d$metastasis <- factor(ifelse(d$metastasis_raw=='1','Metastatic','M0'),levels=c('M0','Metastatic'))
  tab <- table(d$group,d$metastasis)
  if(any(rowSums(tab)<5)||any(colSums(tab)==0)) stop('Metastasis comparison needs >=5 per expression group and both metastasis states')
  # Matrix orientation gives odds(Metastatic | High) / odds(Metastatic | Low).
  fisher <- fisher.test(tab)
  counts <- as.data.frame(tab);names(counts)<-c('expression_group','metastasis','n')
  counts$fraction <- counts$n / as.numeric(rowSums(tab)[as.character(counts$expression_group)])
  statistics <- data.frame(cohort=cohort,cutoff_method=cutoff$method,cutoff=cutoff$value,n=nrow(d),excluded_n=sum(!eligible),
    odds_ratio_high_vs_low=unname(fisher$estimate),ci_lower=fisher$conf.int[1],ci_upper=fisher$conf.int[2],p_value=fisher$p.value,
    interpretation=if(cutoff$method=='optimal') 'Exploratory: threshold selected using survival in overlapping patients; not an independent validation' else 'Fisher exact test; unadjusted for clinical confounding')
  if(cutoff$method=='optimal') {statistics$p_value<-NA_real_;statistics$ci_lower<-NA_real_;statistics$ci_upper<-NA_real_}
  p <- ggplot2::ggplot(counts,ggplot2::aes(expression_group,fraction,fill=metastasis))+ggplot2::geom_col(width=.6)+
    ggplot2::geom_text(ggplot2::aes(label=paste0('n=',n)),position=ggplot2::position_stack(vjust=.5),color='white')+
    ggplot2::scale_fill_manual(values=c(M0='#477b9e',Metastatic='#b43b47'))+ggplot2::theme_classic(base_size=12)+
    ggplot2::labs(title=paste('Metastasis by expression:',cohort),subtitle=if(cutoff$method=='optimal')sprintf('Survival-selected cutoff %.4g; descriptive only (no inferential p-value)',cutoff$value) else sprintf('%s cutoff %.4g; Fisher p %.3g',cutoff$method,cutoff$value,fisher$p.value),x='Expression group',y='Fraction within expression group',fill='Verified status')
  list(plot=p,statistics=statistics,counts=counts,data=d,exclusions=data[!eligible,,drop=FALSE])
}

profile_clinical <- function(retrieval,survival_source,methods=c('median','mean','optimal')) {
  results <- list(); d <- retrieval$data
  d$time <- d$event <- NA_real_
  if(!is.null(survival_source)) {
    i <- match(d$sample_id,survival_source$data$sample_id)
    d$time <- survival_source$data$time[i];d$event<-survival_source$data$event[i]
  }
  for(cohort in c('all',sort(unique(d$subgroup[!r2_missing(d$subgroup)])))) for(method in methods) {
    subset <- if(cohort=='all') d else d[which(d$subgroup==cohort),,drop=FALSE]
    complete <- subset[is.finite(subset$expression)&is.finite(subset$time)&!is.na(subset$event),,drop=FALSE]
    attempt <- function(expr) tryCatch(expr,error=function(e)list(status='unavailable',reason=conditionMessage(e)))
    cut <- attempt(profile_cutoff(subset$expression,method,complete))
    survival <- if(!is.null(cut$value)) attempt(profile_survival(subset,cut,cohort)) else cut
    met <- if(!isTRUE(retrieval$audit$metastasis$coding_verified)) list(status='unavailable',reason='Metastasis coding is not verified') else if(!is.null(cut$value)) attempt(profile_metastasis(subset,cut,cohort)) else cut
    results[[paste(cohort,method,sep=':')]] <- list(cohort=cohort,method=method,cutoff=cut,survival=survival,metastasis=met)
  }
  # One exploratory family across all attempted cohorts and all three cutoff modes.
  for(endpoint in c('survival','metastasis')) {
    ps <- vapply(results,function(r) {s<-r[[endpoint]]$statistics;if(is.null(s))return(NA_real_);if(endpoint=='survival')s$cutoff_adjusted_p else s$p_value},numeric(1))
    adjusted <- p.adjust(ps,'BH',n=length(results))
    for(i in seq_along(results)) if(!is.null(results[[i]][[endpoint]]$statistics)) results[[i]][[endpoint]]$statistics$family_adjusted_p <- adjusted[i]
  }
  results
}
