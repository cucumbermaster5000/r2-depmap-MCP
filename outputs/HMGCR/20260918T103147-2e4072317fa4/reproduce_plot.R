# Open this output directory as the working directory in RStudio.
# Uses the saved patient data; no network request and no new clinical analysis.
source('analysis_functions.R')
retrieval <- readRDS('retrieval.rds')
validation <- plot_r2_subgroup_test(retrieval, output_dir = 'rerun')
print(validation$plot)
validation$statistics
