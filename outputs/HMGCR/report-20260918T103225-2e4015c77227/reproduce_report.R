source('analysis_functions.R')
x <- readRDS('report_data.rds')
write_gene_page(x$retrieval, x$validation, x$publications, x$publication_error, 'regenerated_report.html', profile_extra_html(x$profile, x$cutoff))
