source('R/common.R')
responses<-lapply(readLines('outputs/audit/mcp-v1-responses.jsonl'),jsonlite::fromJSON,simplifyVector=FALSE)
result<-responses[[5]]$result$structuredContent
path<-result$html_path;folder<-dirname(path);x<-readRDS(file.path(folder,'report_data.rds'))
stopifnot(all(x$profile$status$status=='ok'),nrow(x$profile$pfister$all_samples)==272,nrow(x$profile$pfister$data)==246,
  nrow(x$profile$depmap$models)==10,all(is.na(x$profile$clinical[['all:optimal']]$metastasis$statistics$p_value)))
doc<-xml2::read_html(path)
ids<-xml2::xml_attr(xml2::xml_find_all(doc,'//*[@id]'),'id')
stopifnot(!anyDuplicated(ids),length(xml2::xml_find_all(doc,"//*[@class='clinical-panel']"))==15,
  length(xml2::xml_find_all(doc,"//*[@class='clinical-panel' and not(@hidden)]"))==1)
cat('Artifact OK:',path,'\n')
cat('DepMap medulloblastoma models with gene effect:',sum(is.finite(x$profile$depmap$models$gene_effect)),'of',nrow(x$profile$depmap$models),'\n')
cat('TF annotation:',x$profile$enrichr$tf$classification,'\n')
cat('Enrichr memberships:',nrow(x$profile$enrichr$membership),'; TF target rows:',nrow(x$profile$enrichr$targets),'\n')
old<-setwd(folder)
source('reproduce_report.R')
setwd(old)
stopifnot(file.exists(file.path(folder,'regenerated_report.html')))
cat('Offline HTML reproduction OK\n')
