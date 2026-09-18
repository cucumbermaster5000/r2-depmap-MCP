source('R/common.R');source('R/r2_client.R')
out <- 'outputs/audit/profile-sources';dir.create(out,recursive=TRUE,showWarnings=FALSE)
h<-r2_session(90)
urls <- list(figshare='https://api.figshare.com/v2/articles/27993248/versions/1',
  enrichr_catalog='https://maayanlab.cloud/Enrichr/datasetStatistics')
for(n in names(urls)) tryCatch({txt<-r2_request(h,urls[[n]]);writeBin(charToRaw(txt),file.path(out,paste0('depmap-',n,'.txt')));cat(n,nchar(txt),'characters\n')},error=function(e)cat(n,conditionMessage(e),'\n'))
