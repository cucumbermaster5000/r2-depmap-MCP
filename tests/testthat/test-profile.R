testthat::test_that('cutoffs preserve ties and correct the complete cutoff search', {
  x<-rep(1:20,each=2)
  d<-data.frame(expression=x,time=rev(seq_along(x))/10,event=rep(c(1,1,0,1),10))
  testthat::expect_equal(profile_cutoff(x,'median')$value,median(x))
  testthat::expect_equal(profile_cutoff(x,'mean')$value,mean(x))
  cut<-profile_cutoff(x,'optimal',d)
  testthat::expect_true(all(cut$scan$low_n>=10 & cut$scan$high_n>=10))
  testthat::expect_equal(cut$scan$adjusted_p,p.adjust(cut$scan$p_value,'bonferroni'))
  testthat::expect_equal(cut$value,cut$scan$cutoff[which.min(cut$scan$p_value)])
  testthat::expect_error(profile_cutoff(rep(1,20),'median'),'variable')
})

testthat::test_that('metastasis odds ratio is high vs low with correct denominators', {
  d<-data.frame(sample_id=paste0('s',1:40),expression=rep(c(1,2),each=20),metastasis_raw=c(rep('0',18),rep('1',2),rep('0',5),rep('1',15)))
  m<-profile_metastasis(d,list(value=1.5,method='median'),'test')
  expected<-fisher.test(matrix(c(18,5,2,15),nrow=2))
  testthat::expect_equal(m$statistics$odds_ratio_high_vs_low,unname(expected$estimate))
  testthat::expect_gt(m$statistics$odds_ratio_high_vs_low,1)
  testthat::expect_equal(as.numeric(tapply(m$counts$fraction,m$counts$expression_group,sum)),c(1,1))
  selected<-profile_metastasis(d,list(value=1.5,method='optimal'),'test')
  testthat::expect_true(is.na(selected$statistics$p_value))
  testthat::expect_true(is.na(selected$statistics$ci_lower))
  d$metastasis_raw[1]<-'na'
  testthat::expect_equal(profile_metastasis(d,list(value=1.5,method='mean'),'test')$statistics$excluded_n,1)
})

testthat::test_that('survival source verification rejects event and identifier mismatches', {
  retrieval<-list(data=data.frame(sample_id=c('GSM1','GSM2'),survival_time_years=c(1,2),survival_status_raw=c('0','1')))
  payload<-list(plotData=list(xLabel='Follow up in months',yLabel='overall survival probability',data=data.frame(id=c('gsm1','gsm2'),status=c(0,1),xValue=c(12,24))))
  html<-function(p)paste0("<script>import('/d3/plots/kaplan.js').then(m => m['default'](",json_text(p,pretty=FALSE),"));</script>")
  js<-'proportion *= (nrRemaining - status) / nrRemaining;'
  testthat::expect_equal(verify_r2_survival_payload(html(payload),js,retrieval)$data$time,c(1,2))
  payload$plotData$data$status<-c(1,0)
  testthat::expect_error(verify_r2_survival_payload(html(payload),js,retrieval),'event flags')
  payload$plotData$data$id<-c('other','gsm2')
  testthat::expect_error(verify_r2_survival_payload(html(payload),js,retrieval),'unmatched')
})

testthat::test_that('survival KM, censoring and risk tables agree with survival package', {
  d<-data.frame(sample_id=paste0('s',1:40),expression=rep(c(1,2),each=20),time=c(1:20,1:20/2),event=rep(c(1,0,1,1),10))
  cut<-list(value=1.5,method='median',candidates=1)
  s<-profile_survival(d,cut,'test')
  d$group<-factor(ifelse(d$expression>1.5,'High','Low'),levels=c('Low','High'))
  lr<-survival::survdiff(survival::Surv(time,event)~group,d)
  testthat::expect_equal(s$statistics$logrank_p,pchisq(lr$chisq,1,lower.tail=FALSE))
  testthat::expect_equal(s$statistics$deaths,sum(d$event))
  testthat::expect_equal(s$risk_table$n_at_risk[s$risk_table$time_years==0],c(20L,20L))
  testthat::expect_equal(sum(s$curves$censored),sum(d$event==0))
  d$time[1]<-NA
  testthat::expect_equal(profile_survival(d,cut,'test')$statistics$excluded_n,1)
})

testthat::test_that('all subtype labels and BH pairwise family are retained', {
  d<-data.frame(sample_id=paste0('s',1:12),expression=1:12,subtype=rep(c('a','b','c'),each=4))
  r<-profile_group_comparison(d,'subtype','fixture')
  testthat::expect_equal(nrow(r$pairwise),3)
  testthat::expect_equal(r$pairwise$adjusted_p,p.adjust(r$pairwise$p_value,'BH'))
  testthat::expect_equal(sort(r$summaries$group),c('a','b','c'))
})

testthat::test_that('multiple inline SVGs have distinct identifier namespaces', {
  p<-ggplot2::ggplot(data.frame(x=1:3,y=1:3),ggplot2::aes(x,y))+ggplot2::geom_point()
  a<-xml2::read_xml(plot_inline_svg(p));b<-xml2::read_xml(plot_inline_svg(p))
  ids<-function(x)xml2::xml_attr(xml2::xml_find_all(x,'//*[@id]'),'id')
  testthat::expect_length(intersect(ids(a),ids(b)),0)
  testthat::expect_error(dispatch_v1_tool('build_gene_page',list(gene='HLX',cutoff='best'),test_cfg()),'choice')
})
