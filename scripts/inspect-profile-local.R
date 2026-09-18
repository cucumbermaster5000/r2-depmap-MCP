for(f in c('common','r2_client','r2_expression')) source(file.path('R',paste0(f,'.R')))
x <- paste(readLines('outputs/audit/profile-sources/survival-result.html',warn=FALSE),collapse='\n')
p <- r2_script_payload(x,'/d3/plots/kaplan.js')
print(names(p$plotData)); print(names(p$plotData$data)); print(head(p$plotData$data,3)); print(p$plotData$legend)
saveRDS(p,'outputs/audit/profile-sources/kaplan-payload.rds')
a <- jsonlite::fromJSON('outputs/audit/profile-sources/pfister-json_cg_sampleannotation_v1.html')$data
print(table(a$entity,a$cancer_type)); print(table(a$entity,a$subgroup))
