r2_missing <-
function (x) 
is.na(x) | trimws(as.character(x)) %in% c("", "na", "NA")
plot_r2_subgroup_test <-
function (result, output_dir = result$provenance$output_dir) 
{
    if (!requireNamespace("ggplot2", quietly = TRUE)) 
        stop("Install ggplot2 for the Version 1 plot", call. = FALSE)
    d <- result$data
    eligible <- is.finite(d$expression) & !r2_missing(d$subgroup)
    exclusions <- d[!eligible, c("sample_id", "expression", "subgroup"), 
        drop = FALSE]
    used <- d[eligible, , drop = FALSE]
    sizes <- table(used$subgroup)
    if (length(sizes) < 2 || any(sizes < 3)) 
        stop("Need at least two annotated subgroups with three samples each", 
            call. = FALSE)
    if (anyDuplicated(used$sample_id)) 
        stop("Repeated sample IDs; unpaired comparison is invalid", 
            call. = FALSE)
    test <- stats::kruskal.test(expression ~ subgroup, used)
    summaries <- do.call(rbind, lapply(split(used$expression, 
        used$subgroup), function(x) data.frame(n = length(x), 
        mean = mean(x), median = stats::median(x), sd = stats::sd(x), 
        iqr = stats::IQR(x), min = min(x), max = max(x))))
    summaries$subgroup <- rownames(summaries)
    rownames(summaries) <- NULL
    statistics <- data.frame(test = test$method, statistic = unname(test$statistic), 
        df = unname(test$parameter), p_value = test$p.value, 
        n = nrow(used), groups = length(sizes), excluded_n = nrow(exclusions), 
        epsilon_squared = max(0, (unname(test$statistic) - length(sizes) + 
            1)/(nrow(used) - length(sizes))), multiple_testing = "One predefined omnibus test for one gene; no pairwise tests or multiplicity adjustment in Version 1")
    labels <- setNames(paste0(names(sizes), "\n(n = ", as.integer(sizes), 
        ")"), names(sizes))
    p <- ggplot2::ggplot(used, ggplot2::aes(x = subgroup, y = expression, 
        fill = subgroup)) + ggplot2::geom_boxplot(width = 0.5, 
        outlier.shape = NA, alpha = 0.34999999999999998, linewidth = 0.45000000000000001) + 
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.16, 
            height = 0, seed = 85217), size = 1.1000000000000001, 
            alpha = 0.5, shape = 16) + ggplot2::scale_x_discrete(labels = labels) + 
        ggplot2::guides(fill = "none") + ggplot2::labs(title = paste(unique(used$queried_gene), 
        "expression across medulloblastoma subgroups"), subtitle = sprintf("Kruskal-Wallis H(%d) = %.2f; p = %.3g", 
        test$parameter, test$statistic, test$p.value), x = "Original R2 subgroup labels", 
        y = "Expression (R2 log2 transformation)", caption = paste0("GSE85217 / Cavalli; reporter ", 
            result$provenance$reporter, ". Each point is one R2 sample.\n", 
            "Overall distribution comparison; this test does not establish subgroup specificity. Excluded: ", 
            nrow(exclusions), ".")) + ggplot2::theme_classic(base_size = 12) + 
        ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, 
            size = 9))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    gene <- result$provenance$queried_gene
    files <- list(pdf = file.path(output_dir, paste0(gene, "_subgroup_test.pdf")), 
        png = file.path(output_dir, paste0(gene, "_subgroup_test.png")), 
        statistics = file.path(output_dir, paste0(gene, "_subgroup_test.csv")), 
        summaries = file.path(output_dir, paste0(gene, "_subgroup_summaries.csv")), 
        exclusions = file.path(output_dir, "subgroup_test_exclusions.csv"))
    ggplot2::ggsave(files$pdf, p, width = 8, height = 5.7999999999999998, 
        device = grDevices::pdf)
    ggplot2::ggsave(files$png, p, width = 8, height = 5.7999999999999998, 
        dpi = 300)
    data.table::fwrite(statistics, files$statistics)
    data.table::fwrite(summaries, files$summaries)
    data.table::fwrite(exclusions, files$exclusions)
    writeLines("Predefined unpaired Kruskal-Wallis comparison of distributions. Samples are treated as independent primary tumours as described by the source. Batch, age and other confounding are not adjusted in this Version 1 smoke test. No pairwise or subtype/clinical analyses were performed. Missing annotations/expression are reported in the exclusions file. Epsilon-squared is the non-negative rank-based omnibus effect estimate (H-k+1)/(n-k).", 
        file.path(output_dir, "subgroup_test_method.txt"))
    dump(c("r2_missing", "plot_r2_subgroup_test"), file = file.path(output_dir, 
        "analysis_functions.R"))
    writeLines(c("# Open this output directory as the working directory in RStudio.", 
        "# Uses the saved patient data; no network request and no new clinical analysis.", 
        "source('analysis_functions.R')", "retrieval <- readRDS('retrieval.rds')", 
        "validation <- plot_r2_subgroup_test(retrieval, output_dir = 'rerun')", 
        "print(validation$plot)", "validation$statistics"), file.path(output_dir, 
        "reproduce_plot.R"))
    writeLines(capture.output(sessionInfo()), file.path(output_dir, 
        "analysis_session-info.txt"))
    list(statistics = statistics, summaries = summaries, files = files, 
        plot = p)
}
