for(f in c('common','cache','r2_client','r2_expression','study','enrichment','depmap','profile_sources','profile_depmap'))source(file.path('R',paste0(f,'.R')))
cfg<-read_server_config('config/defaults.json');cfg$output_dir<-normalizePath('outputs',winslash='/')
r<-get_r2_expression(cfg,'HLX');s<-get_cavalli_survival(cfg,r);cat('Verified survival:',nrow(s$data),'\n')
p<-get_pfister_expression(cfg,'HLX');cat('Pfister:',nrow(p$data),'CNS primary:',sum(p$data$include_cns),'\n')
for(l in c('GO_Biological_Process_2026','Reactome_Pathways_2024','KEGG_2021_Human','WikiPathways_2024_Human','ChEA_2022')) {
  gs<-load_gene_sets(cfg,l);cat(l,length(gs$sets),'terms\n')
}
d<-get_public_depmap_gene(cfg,'HLX');cat('DepMap:',nrow(d$data),'models;',sum(d$data$is_medulloblastoma),'MB models\n')
