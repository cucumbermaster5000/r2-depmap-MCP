# Independent of Shiny; consumes the already validated retrieval, without network I/O.
analyze_r2_subtypes <- function(retrieval, subtype_field = "subtype") {
  d <- as.data.frame(retrieval$data)
  m <- as.data.frame(retrieval$metadata)
  if (!all(c("sample_id", "expression") %in% names(d))) stop("Missing patient identifiers/expression")
  if (anyDuplicated(d$sample_id) || any(r2_missing(d$sample_id))) stop("Invalid patient identifiers")
  gene <- retrieval$provenance$queried_gene
  warnings <- character()
  if (!"samplenames" %in% names(m) || anyDuplicated(m$samplenames) ||
      !setequal(d$sample_id, m$samplenames)) stop("Subtype metadata identifiers do not match expression")
  m <- m[match(d$sample_id, m$samplenames), , drop = FALSE]
  # Preserve every original field under its exact name, including collisions.
  for (field in names(m)) {
    if (field %in% names(d)) d[[paste0("retrieved_", field)]] <- d[[field]]
    d[[field]] <- m[[field]]
  }
  available <- subtype_field %in% names(m)
  labels <- if (available) as.character(m[[subtype_field]]) else rep(NA_character_, nrow(d))
  d$subtype <- labels
  if (!"subgroup" %in% names(d)) d$subgroup <- NA_character_
  d$queried_gene <- gene
  expected <- c("wnt_alpha", "wnt_beta", "shh_alpha", "shh_beta", "shh_gamma", "shh_delta",
    "group3_alpha", "group3_beta", "group3_gamma", "group4_alpha", "group4_beta", "group4_gamma")
  observed <- unique(labels[!r2_missing(labels)])
  documented <- identical(as.character(retrieval$audit$classification$source_publication), "28609654") &&
    identical(retrieval$provenance$dataset_table, "ps_avgpres_gse85217geo763_hugene11t") &&
    identical(subtype_field, "subtype") && length(observed) > 0 && all(observed %in% expected)
  classification <- if (documented) "Cavalli 2017: original R2 molecular subtypes (PMID 28609654)" else
    "Unidentified source classification; original labels only"
  if (!available) warnings <- c(warnings, "No molecular subtype metadata is available.")
  if (available && !documented) warnings <- c(warnings, "Subtype classification could not be identified. Labels are shown without assigning a published classification.")
  d$subtype_classification <- classification
  bad_subtype <- r2_missing(labels)
  bad_expression <- !is.finite(d$expression)
  d$included_in_subtype_analysis <- !bad_subtype & !bad_expression
  d$exclusion_reason <- ifelse(bad_subtype & bad_expression, "Missing subtype and non-finite expression",
    ifelse(bad_subtype, "Missing subtype", ifelse(bad_expression, "Non-finite expression", "")))
  used <- d[d$included_in_subtype_analysis, , drop = FALSE]
  cross <- unique(d[!bad_subtype, c("subtype", "subgroup"), drop = FALSE])
  # Ordering and colour use observed metadata membership, never subtype string parsing.
  group_order <- c("wnt", "shh", "group3", "group4")
  order_labels <- observed[order(match(vapply(observed, function(s) {
    g <- unique(as.character(cross$subgroup[cross$subtype == s]))
    if (length(g) == 1 && !r2_missing(g)) g else "Unresolved"
  }, character(1)), group_order), match(observed, expected), observed, na.last = TRUE)]
  summaries <- do.call(rbind, lapply(order_labels, function(s) {
    x <- used$expression[used$subtype == s]
    g <- unique(as.character(cross$subgroup[cross$subtype == s]))
    data.frame(subtype = s, subgroup = if (length(g) == 1 && !r2_missing(g)) g else "Unresolved",
      annotated_n = sum(labels == s, na.rm = TRUE), n = length(x),
      mean = if(length(x)) mean(x) else NA_real_, median = if(length(x)) median(x) else NA_real_,
      sd = if(length(x)>1) sd(x) else NA_real_, variance = if(length(x)>1) var(x) else NA_real_,
      iqr = if(length(x)) IQR(x) else NA_real_, min = if(length(x)) min(x) else NA_real_,
      max = if(length(x)) max(x) else NA_real_, unique_values = length(unique(x)))
  }))
  if (is.null(summaries)) summaries <- data.frame(subtype=character(), subgroup=character(), annotated_n=integer(), n=integer(), mean=numeric(), median=numeric(), sd=numeric(), variance=numeric(), iqr=numeric(), min=numeric(), max=numeric(), unique_values=integer())
  sizes <- summaries$n[summaries$n > 0]
  k <- length(sizes); n <- nrow(used)
  rationale <- "Predefined rank-based comparison of independent tumour expression distributions; no normality-test-driven selection. Group sizes, spread, ties and missingness are inspected first. Unequal shapes/spreads preclude interpreting this solely as a median test."
  overall <- data.frame(test = "Not performed", statistic = NA_real_, df = NA_real_, p_value = NA_real_,
    adjusted_p_value = NA_real_, epsilon_squared = NA_real_, n = n, groups = k,
    reason = "", multiple_testing = "One omnibus test per gene; no across-gene correction")
  reason <- if (k < 2) "At least two subtypes with finite expression are required." else if (any(sizes < 5))
    "At least one subtype has fewer than five observations: descriptive results only; asymptotic inference withheld." else if (length(unique(used$expression)) < 2)
    "Expression is constant: inferential tests are undefined." else ""
  pairwise <- data.frame(subtype_a=character(), subtype_b=character(), n_a=integer(), n_b=integer(),
    test=character(), statistic=numeric(), df=numeric(), median_difference=numeric(), rank_biserial=numeric(),
    direction=character(), p_value=numeric(), adjusted_p_value=numeric())
  if (nzchar(reason)) { warnings <- c(warnings, reason); overall$reason <- reason } else {
    kw <- stats::kruskal.test(used$expression, used$subtype)
    overall$test <- kw$method; overall$statistic <- unname(kw$statistic)
    overall$df <- unname(kw$parameter); overall$p_value <- kw$p.value
    overall$epsilon_squared <- max(0, (overall$statistic-k+1)/(n-k))
    overall$reason <- rationale
    pairs <- combn(summaries$subtype[summaries$n > 0], 2, simplify = FALSE)
    pairwise <- do.call(rbind, lapply(pairs, function(pair) {
      x <- used$expression[used$subtype == pair[1]]; y <- used$expression[used$subtype == pair[2]]
      u <- sum(rank(c(x,y))[seq_along(x)]) - length(x)*(length(x)+1)/2
      effect <- 2*u/(length(x)*length(y))-1
      constant <- length(unique(c(x,y))) == 1
      exact <- length(x)<50 && length(y)<50 && !anyDuplicated(c(x,y))
      w <- if (!constant) stats::wilcox.test(x,y, exact=exact, correct=TRUE) else NULL
      data.frame(subtype_a=pair[1], subtype_b=pair[2], n_a=length(x), n_b=length(y),
        test=if(constant) "Identical constant values; p=1" else w$method,
        statistic=u, df=NA_real_, median_difference=median(x)-median(y), rank_biserial=effect,
        direction=if(effect>0) "A higher rank tendency" else if(effect<0) "A lower rank tendency" else "No rank tendency",
        p_value=if(constant) 1 else w$p.value, adjusted_p_value=NA_real_)
    }))
    pairwise$adjusted_p_value <- stats::p.adjust(pairwise$p_value, method="BH")
  }
  associations <- pairwise[, c("subtype_a","subtype_b","direction","n_a","n_b","rank_biserial","median_difference","p_value","adjusted_p_value")]
  names(associations)[1:2] <- c("subtype", "comparison")
  used$plot_subtype <- factor(used$subtype, levels=order_labels)
  used$plot_subgroup <- summaries$subgroup[match(used$subtype, summaries$subtype)]
  p <- ggplot2::ggplot(used, ggplot2::aes(plot_subtype, expression, fill=plot_subgroup)) +
    ggplot2::geom_boxplot(width=.6, outlier.shape=NA, alpha=.45) +
    ggplot2::geom_point(position=ggplot2::position_jitter(width=.18, height=0, seed=85217), size=.9, alpha=.5) +
    ggplot2::scale_x_discrete(drop=FALSE, labels=setNames(paste0(summaries$subtype,"\n(n = ",summaries$n,")"), summaries$subtype)) +
    ggplot2::scale_fill_manual(values=c(wnt="#377EB8",shh="#E69F00",group3="#984EA3",group4="#009E73",Unresolved="#777777"), na.value="#777777") +
    ggplot2::labs(title=paste(gene,"expression across molecular subtypes"), subtitle=classification,
      x="Original R2 subtype labels", y=if(grepl("log2",retrieval$provenance$transformation %||% "")) "Expression (R2 log2)" else "Expression (source scale)",
      fill="Observed subgroup", caption=paste0("Each point is one sample; excluded: ",sum(!d$included_in_subtype_analysis),
        ". Full pairwise results are reported separately.\nUnadjusted observational comparisons do not establish subtype specificity.")) +
    ggplot2::theme_classic(base_size=12) + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1),
      legend.position="top", plot.caption=ggplot2::element_text(hjust=0,size=9))
  list(gene=gene, status=if(nzchar(reason)) "descriptive_only" else "complete", data=d,
    metadata=list(field=subtype_field, available=available, classification=classification,
      source_publication=retrieval$audit$classification$source_publication, original_fields=names(m),
      original_labels=order_labels, subtype_by_subgroup=cross, annotated_n=sum(!bad_subtype), missing_n=sum(bad_subtype)),
    missing=list(total_n=nrow(d), missing_subtype=sum(bad_subtype), nonfinite_expression=sum(bad_expression),
      excluded_n=sum(!d$included_in_subtype_analysis)), exclusions=d[!d$included_in_subtype_analysis,,drop=FALSE],
    summaries=summaries, diagnostics=list(group_sizes=sizes, tied_observations=n-length(unique(used$expression)),
      rationale=rationale), overall=overall, pairwise=pairwise, associations=associations, plot=p,
    warnings=unique(c(warnings,"Comparisons are unadjusted for age, batch and other confounders. BH correction covers subtype pairs for this gene only; repeated gene queries need separate multiplicity control.")))
}
