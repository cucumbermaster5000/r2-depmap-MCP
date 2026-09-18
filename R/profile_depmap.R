depmap_public_files <- function(cfg) {
  folder<-file.path(cfg$cache_dir,'depmap-24Q4-v1');dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  manifest_path<-file.path(folder,'manifest.json')
  if(!file.exists(manifest_path)) {
    txt<-enrichr_get('https://api.figshare.com/v2/articles/27993248/versions/1')
    manifest<-jsonlite::fromJSON(txt)
    if(!identical(manifest$title,'DepMap 24Q4 Public')||manifest$version!=1L)stop('DepMap archive identity mismatch')
    writeBin(charToRaw(txt),manifest_path)
  }
  manifest<-jsonlite::fromJSON(manifest_path)
  wanted<-c('Model.csv','CRISPRGeneEffect.csv')
  paths<-setNames(character(length(wanted)),wanted)
  for(name in wanted) {
    f<-manifest$files[manifest$files$name==name,,drop=FALSE]
    if(nrow(f)!=1)stop('Ambiguous DepMap archive file')
    path<-file.path(folder,name)
    if(!file.exists(path)) {
      message('Downloading pinned DepMap ',name,' (',round(f$size/1e6),' MB), once per installation')
      part<-paste0(path,'.part')
      response<-curl::curl_fetch_disk(f$download_url,part,curl::new_handle(timeout=900,connecttimeout=30,followlocation=TRUE))
      if(response$status_code!=200 || unname(file.info(part)$size)!=f$size || unname(tools::md5sum(part))!=f$computed_md5)stop('DepMap download failed integrity validation: ',name)
      if(!file.rename(part,path))stop('Could not finalize DepMap download')
    }
    if(unname(file.info(path)$size)!=f$size||unname(tools::md5sum(path))!=f$computed_md5)stop('Cached DepMap file integrity mismatch')
    paths[name]<-normalizePath(path,winslash='/')
  }
  list(paths=paths,provenance=list(release=manifest$title,doi=manifest$doi,manifest=normalizePath(manifest_path,winslash='/'),
    note='Pinned official archived 24Q4 v1 release, not the latest release. Portal automated access currently requires human verification.',files=manifest$files[manifest$files$name %in% wanted,c('name','download_url','computed_md5','size')]))
}

get_public_depmap_gene <- function(cfg,gene) {
  key<-cache_key('public_depmap_24q4_v1',list(gene=toupper(gene)))
  cached<-cache_get(cfg,key);if(!is.null(cached))return(cached)
  files<-depmap_public_files(cfg)
  models<-read_model_metadata(files$paths[['Model.csv']]); require_columns(models,c('ModelID','CellLineName','OncotreePrimaryDisease'))
  effects<-read_gene_effect_columns(files$paths[['CRISPRGeneEffect.csv']],gene)
  effect_col<-effects$mapping$column[1];ids<-as.character(effects$data[[effects$id_col]])
  if(any(!ids %in% models$ModelID))stop('DepMap effect IDs absent from same-release metadata')
  values<-as.numeric(effects$data[[effect_col]]);idx<-match(models$ModelID,ids)
  data<-models[,intersect(c('ModelID','CellLineName','OncotreeLineage','OncotreePrimaryDisease','OncotreeSubtype','OncotreeCode'),names(models)),drop=FALSE]
  data$gene_effect<-values[idx]
  disease_fields<-intersect(c('OncotreePrimaryDisease','OncotreeSubtype'),names(data))
  data$is_medulloblastoma<-Reduce(`|`,lapply(disease_fields,function(f)!is.na(data[[f]])&grepl('^medulloblastoma(,|$)',data[[f]],ignore.case=TRUE)))
  if(!any(data$is_medulloblastoma))stop('No explicit medulloblastoma disease annotation found')
  data$gene<-toupper(gene)
  result<-list(data=data,provenance=c(files$provenance,list(gene_column=effect_col,model_selection='Explicit OncotreePrimaryDisease or OncotreeSubtype label Medulloblastoma (including comma-qualified subtypes); no subgroup inferred from model names',retrieved_at=utc_now(),interpretation='Chronos gene effect: more negative indicates greater in-vitro dependency; 0 and -1 are reference values, not clinical thresholds')))
  cache_put(cfg,key,result);result
}

plot_public_depmap <- function(result) {
  d<-result$data;d<-d[is.finite(d$gene_effect),,drop=FALSE]
  d<-d[order(d$gene_effect,d$ModelID),];d$rank<-seq_len(nrow(d))
  d$context<-ifelse(d$is_medulloblastoma,'Medulloblastoma','Other models')
  p<-ggplot2::ggplot(d,ggplot2::aes(rank,gene_effect,color=context))+ggplot2::geom_point(alpha=.65,size=1.5)+
    ggplot2::geom_hline(yintercept=c(0,-1),linetype=3,color='grey50')+
    ggplot2::geom_text(data=d[d$is_medulloblastoma,],ggplot2::aes(label=CellLineName),size=2.6,angle=45,hjust=0,vjust=-.4,check_overlap=FALSE)+
    ggplot2::scale_color_manual(values=c(Medulloblastoma='#b43b47','Other models'='#b6c0c7'))+ggplot2::theme_classic()+
    ggplot2::labs(title=paste(unique(d$gene),'dependency across DepMap models'),subtitle=result$provenance$release,x='Models ordered by gene effect',y='Chronos gene effect',color='Disease context')
  mb<-result$data[result$data$is_medulloblastoma,,drop=FALSE]
  mp<-ggplot2::ggplot(mb[is.finite(mb$gene_effect),],ggplot2::aes(reorder(CellLineName,gene_effect),gene_effect))+ggplot2::geom_point(color='#b43b47',size=3)+
    ggplot2::geom_hline(yintercept=c(0,-1),linetype=3,color='grey50')+ggplot2::coord_flip()+ggplot2::theme_classic()+
    ggplot2::labs(title='Medulloblastoma models',x='DepMap cell-line name',y='Chronos gene effect')
  list(plot=p,mb_plot=mp,models=mb,all_models=result$data,provenance=result$provenance)
}
