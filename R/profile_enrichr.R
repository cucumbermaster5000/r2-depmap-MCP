get_tf_annotation <- function(cfg,gene) {
  key<-'human_tfs_lambert_v101';x<-cache_get(cfg,key)
  if(is.null(x)) {
    url<-'https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt'
    raw<-enrichr_get(url);genes<-trimws(strsplit(raw,'\n',fixed=TRUE)[[1]])
    if(length(genes)<1000||any(grepl('<html',genes,ignore.case=TRUE)))stop('Unexpected Human TFs response')
    path<-file.path(cfg$cache_dir,'TF_names_v_1.01.txt');writeBin(charToRaw(raw),path)
    x<-list(genes=genes,source=url,version='Lambert human TF list v1.01',retrieved_at=utc_now(),md5=unname(tools::md5sum(path)),raw_file=path);cache_put(cfg,key,x)
  }
  list(gene=gene,classification=if(toupper(gene) %in% toupper(x$genes))'Listed human transcription factor' else 'Not listed in this version; absence is not proof of non-TF function',source=x$source,version=x$version,raw_file=x$raw_file,md5=x$md5)
}

profile_enrichr <- function(cfg,gene,libraries=c('GO_Biological_Process_2026','Reactome_Pathways_2024','KEGG_2021_Human','WikiPathways_2024_Human')) {
  gene<-toupper(gene)
  inputs<-setNames(lapply(libraries,function(l)load_gene_sets(cfg,l)),libraries)
  rows<-lapply(names(inputs),function(l) {
    sets<-inputs[[l]]$sets;terms<-names(sets)[vapply(sets,function(s)gene %in% s,logical(1))]
    data.frame(library=rep(l,length(terms)),term=terms,set_size=as.integer(lengths(sets[terms])),evidence=rep('Annotated gene membership; not enrichment or pathway activity',length(terms)))
  })
  membership<-do.call(rbind,rows)
  chea<-load_gene_sets(cfg,'ChEA_2022')
  labels<-names(chea$sets)
  selected<-labels[startsWith(labels,paste0(gene,'_')) & grepl('HUMAN',labels,ignore.case=TRUE)]
  target_sets<-chea$sets[selected]
  targets<-if(length(target_sets)) do.call(rbind,lapply(names(target_sets),function(term)data.frame(source_term=term,gene=target_sets[[term]],evidence='ChEA ChIP binding-associated target; not proof of regulation in MB'))) else data.frame()
  enrichments<-list()
  # Per-experiment target sets are kept separate; background is the human ChEA annotation universe.
  human<-chea$sets[grepl('HUMAN',names(chea$sets),ignore.case=TRUE)]
  universe<-unique(unlist(human,use.names=FALSE))
  for(term in names(target_sets)) for(l in names(inputs)) {
    tab<-overrepresentation(target_sets[[term]],universe,inputs[[l]]$sets)
    if(nrow(tab)) {tab$source_term<-term;tab$library<-l;enrichments[[paste(term,l,sep=':')]]<-tab}
  }
  tf<-tryCatch(get_tf_annotation(cfg,gene),error=function(e)list(classification=paste('Unavailable:',conditionMessage(e))))
  list(membership=membership,targets=targets,tf=tf,enrichment=if(length(enrichments))do.call(rbind,enrichments) else data.frame(),
    provenance=list(libraries=lapply(inputs,`[[`,'provenance'),tf_source=chea$provenance,tf_target_sets=as.list(selected),background_n=length(universe),
      method='TF target over-representation: hypergeometric; BH over all eligible terms separately per target experiment and library; human ChEA annotated genes as background. This is an annotation background, not a measured ChIP assay universe.',
      interpretation='Single-gene membership does not establish pathway deregulation. TF ChIP-associated targets are kept separate from perturbation, predictions and correlations; no MB-specific regulation is inferred. No user transcriptomics data were supplied.'))
}
