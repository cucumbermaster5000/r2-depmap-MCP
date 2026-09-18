#!/usr/bin/env Rscript
for (f in c('common', 'r2_client', 'r2_expression')) source(file.path('R', paste0(f, '.R')))
cfg <- read_server_config('config/defaults.json')
out <- 'outputs/audit/profile-sources'; dir.create(out, recursive=TRUE, showWarnings=FALSE)
h <- r2_session(120)
fetch <- function(name, params, post=FALSE) {
  txt <- r2_request(h, if (post) cfg$r2_base_url else paste0(cfg$r2_base_url, '?', encode_form(params)), if (post) params else NULL)
  writeBin(charToRaw(txt), file.path(out,paste0(name,'.html')))
  doc <- xml2::read_html(charToRaw(txt), encoding='UTF-8')
  controls <- lapply(xml2::xml_find_all(doc,'.//select|.//input'), function(n) list(name=xml2::xml_attr(n,'name'),value=xml2::xml_attr(n,'value'),options=lapply(xml2::xml_find_all(n,'.//option'),function(o) list(value=xml2::xml_attr(o,'value'),label=xml2::xml_text(o)))))
  scripts <- xml2::xml_text(xml2::xml_find_all(doc,'.//script'))
  writeLines(json_text(list(controls=controls,scripts=as.list(scripts))),file.path(out,paste0(name,'-controls.json')))
  cat(name,':',nchar(txt),'characters\n'); invisible(txt)
}
fetch('landing',list(table=cfg$r2_dataset_table))
fetch('guest',list(method='guest',open_page='auth'),TRUE)
fetch('main',list(table=cfg$r2_dataset_table))
fetch('survival',list(table=cfg$r2_dataset_table,option='kaplan_main'))
fetch('survival-gene',list(table=cfg$r2_dataset_table,option='kaplanscan_gene',gene='HLX'))
js <- r2_request(h,'https://hgserver1.amc.nl/r2/assets/v393/js/src/datasetSelectionModal.js')
writeLines(js,file.path(out,'datasetSelectionModal.js'),useBytes=TRUE)
catalog_raw <- fetch('catalog',list(json_option='json_dataset',selected='',filter_datasets=''))
catalog <- r2_json(gsub('[[:cntrl:]]',' ',catalog_raw),'catalog (control characters replaced for catalogue discovery only)')
print(catalog[grepl('Pfister',catalog$showname),c('dataset','showname')])
fetch('survival-result',list(table=cfg$r2_dataset_table,option='kaplanscan',inputprobeset='7909890',scanmodus='median',survival='overall',mingrpsize='8',subset=''),TRUE)
hit <- catalog$dataset[grepl('Pfister',catalog$showname) & catalog$dataset_samples == 272 & catalog$platform == 'informp3']
stopifnot(length(hit)==1)
for (kind in c('json_dataset_info','json_dataset_info_tracks','json_cg_sampleannotation_v1')) fetch(paste0('pfister-',kind),list(json_option=kind,table=hit,ctable=hit,subset=''))
writeLines(r2_request(h,'https://hgserver1.amc.nl/r2/assets/v393/js/src/d3/plots/kaplan.js'),file.path(out,'kaplan.js'),useBytes=TRUE)
for (entry in list(c('depmap-openapi','https://depmap.org/breadbox/openapi.json'),c('depmap-files','https://depmap.org/portal/api/download/files'))) {
  txt <- r2_request(h,entry[2]); writeBin(charToRaw(txt),file.path(out,paste0(entry[1],'.txt'))); cat(entry[1],nchar(txt),'characters\n')
}
