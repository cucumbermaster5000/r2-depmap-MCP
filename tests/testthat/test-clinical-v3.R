clinical_fixture <- function() {
  r <- v1_assembled_fixture()
  r$data <- r$data[rep(seq_len(12),each=3),]
  r$metadata <- r$metadata[rep(seq_len(12),each=3),]
  r$data$sample_id <- r$metadata$samplenames <- paste0("P",seq_len(36))
  r$data$expression <- rep(1:12,3)
  r$metadata[["met_status_(1_met__0_m0)"]] <- rep(c("0","1"),18)
  r$audit$metastasis$coding_verified <- TRUE
  r
}

testthat::test_that("metastasis requires explicit coding and retains missing rows", {
  r <- clinical_fixture(); v <- verify_r2_clinical_metadata(r)$metastasis
  a <- analyze_r2_metastasis(r,v)
  testthat::expect_equal(a$status,"complete")
  x <- a$data$expression[a$data$metastasis_raw=="1"]; y <- a$data$expression[a$data$metastasis_raw=="0"]
  w <- wilcox.test(x,y,exact=FALSE)
  testthat::expect_equal(a$statistics$p_value,w$p.value)
  testthat::expect_equal(a$statistics$rank_biserial_metastatic_vs_m0,2*as.numeric(w$statistic)/(length(x)*length(y))-1)
  r$metadata[[v$field]][1] <- "na";r$data$expression[2] <- NA_real_
  a <- analyze_r2_metastasis(r,v)
  testthat::expect_equal(a$missing$excluded_n,2)
  testthat::expect_equal(nrow(a$data),36)
  v$verified <- FALSE;a <- analyze_r2_metastasis(r,v)
  testthat::expect_null(a$plot)
  testthat::expect_true(all(is.na(a$data$metastasis_label)))
  r$metadata[[v$field]][3] <- "unexpected"
  testthat::expect_false(verify_r2_clinical_metadata(r)$metastasis$verified)
  r$metadata[[v$field]] <- NULL
  testthat::expect_match(verify_r2_clinical_metadata(r)$metastasis$evidence,"absent")
})

testthat::test_that("survival uses full-cohort mean and returns valid independent fits", {
  r <- clinical_fixture()
  v <- list(verified=TRUE,evidence="Synthetic verified fixture",data=data.frame(sample_id=r$data$sample_id,
    time=rep(1:18,2)/2,event=rep(c(1,0,1),12)))
  v$data$time[1] <- NA_real_
  a <- analyze_r2_survival(r,v)
  testthat::expect_equal(a$status,"complete")
  testthat::expect_equal(a$cutoff$value,mean(r$data$expression))
  testthat::expect_equal(a$statistics$n,35)
  used <- a$data[a$data$included,];used$group <- factor(used$expression_group,levels=c("Low","High"))
  fit <- survival::coxph(survival::Surv(time,event)~group,used)
  testthat::expect_equal(a$statistics$hazard_ratio_high_vs_low,unname(exp(coef(fit))))
  testthat::expect_equal(a$statistics$cox_p,summary(fit)$coefficients[1,"Pr(>|z|)"])
  testthat::expect_equal(a$statistics$logrank_p,pchisq(survival::survdiff(survival::Surv(time,event)~group,used)$chisq,1,lower.tail=FALSE))
  testthat::expect_s3_class(a$km_fit,"survfit")
  testthat::expect_s3_class(a$ph_diagnostic,"cox.zph")
  testthat::expect_error(analyze_r2_survival(r,v,cutoff="optimal"),"only")
  r$data$expression <- rep(c(1,2,3),12)
  a <- analyze_r2_survival(r,v)
  testthat::expect_true(any(a$data$expression==a$cutoff$value))
  testthat::expect_true(all(a$data$expression_group[a$data$expression==a$cutoff$value]=="Low"))
  v$verified <- FALSE
  testthat::expect_null(analyze_r2_survival(r,v)$km_fit)
})

testthat::test_that("low-event cohorts do not force Cox and absent outcomes do not fit KM", {
  r <- clinical_fixture()
  v <- list(verified=TRUE,data=data.frame(sample_id=r$data$sample_id,time=1:36,event=rep(0,36)))
  testthat::expect_null(analyze_r2_survival(r,v)$km_fit)
  v$data$event[c(1,5,10,15,20,25)] <- 1
  a <- analyze_r2_survival(r,v)
  testthat::expect_equal(a$status,"complete")
  testthat::expect_null(a$cox_fit)
  testthat::expect_true(is.na(a$statistics$cox_p))
})

testthat::test_that("subgroup selection requires consistent effect, FDR, sizes and omnibus", {
  r <- clinical_fixture()
  # Larger fixture with WNT uniquely high; no assumption favouring Group 3.
  i <- rep(seq_len(nrow(r$data)),2)
  r$data <- r$data[i,];r$data$sample_id <- paste0("x",seq_len(nrow(r$data)))
  r$data$expression <- ave(seq_len(nrow(r$data)),r$data$subgroup,FUN=seq_along)/100 + ifelse(r$data$subgroup=="wnt",10,0)
  s <- select_r2_clinical_subgroups(r,analyze_r2_subgroups(r))
  testthat::expect_identical(s$selected,"wnt")
  testthat::expect_equal(nrow(s$pairwise),6)
  testthat::expect_equal(s$pairwise$adjusted_p,p.adjust(s$pairwise$p_value,"BH"))
  a <- analyze_r2_subgroups(r);a$statistics$p_value <- 1
  testthat::expect_length(select_r2_clinical_subgroups(r,a)$selected,0)
})
