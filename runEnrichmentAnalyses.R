library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)

runEnrichmentAnalyses <- function(dds, gene_id_format, custom_label, counts, entrez_df, org_db, outdir, GOBP = TRUE, KEGG = TRUE, WIKIPATHWAYS = TRUE, REACTOME = TRUE){
  
  res_dirname <- custom_label
  res_dir <- file.path(outdir, res_dirname)
  
  # If results dir doesn't exist, create it
  if(!file.exists(res_dir)){
    message("Creating directory: ", res_dir)
    dir.create(res_dir, recursive = T)
  }else{
    message("Outputting files to directory: ", res_dir)
  }
  
  message("Creating log file...")
  sink(file(paste0(res_dir, "/log.txt"), open = "wt"), 
       append=FALSE, 
       type = "message")
  
  res <- DESeq2::results(dds) %>% 
    as.data.frame() %>% 
    rownames_to_column(gene_id_format)
  
  message("~~~ Running enrichment on DEGs from:", custom_label, "~~~")
  
  # Vector of significant genes
  message("Retrieving significant DEGs...")
  
  sigDEGs_pos <- res %>%
    dplyr::filter(padj < 0.05 & log2FoldChange > 0) %>%
    dplyr::pull(!!sym(gene_id_format))
  message(length(sigDEGs_pos), " DEGs found enriched in treatment group")
  
  sigDEGs_neg <- res %>%
    dplyr::filter(padj < 0.05 & log2FoldChange < 0) %>%
    dplyr::pull(!!sym(gene_id_format))
  message(length(sigDEGs_neg), " DEGs found enriched in control group")
  
  # Retrieve samples involved in the comparison
  
  metadata <- colData(dds) %>% as.data.frame()
  target_var <- dds@design[[2]] %>% 
    as.character() %>% 
    str_split("\\+") %>% 
    unlist() %>% 
    last()
  
  reference_level <- levels(metadata[,target_var])[1]
  treatment_level <- levels(metadata[,target_var])[2]
  
  samples_pos <- metadata %>% 
    dplyr::filter(!!sym(target_var) == treatment_level) %>% 
    rownames()
  
  message("Using treatment level: ", treatment_level)
  message("Using treatment group samples: ", paste(samples_pos, collapse = ", "))
  
  samples_neg <- metadata %>% 
    dplyr::filter(!!sym(target_var) == reference_level) %>% 
    rownames()
  message("Using reference level: ", reference_level)
  message("Using reference group samples: ", paste(samples_neg, collapse = ", "))
  
  # Filter counts to only keep those samples
  counts_filtered <- counts[,colnames(counts) %in% c(samples_pos, samples_neg)]
  
  # Get the background universes
  message("Generating background universe...")
  background_universe <- counts_filtered %>% 
    as.data.frame() %>% 
    dplyr::filter(rowSums(.) > 0) %>% 
    rownames()
  message(length(background_universe), " genes used as background universe")
  
  if(GOBP){
    
    message("Running hypergeometric tests against GO:BP...")
    
    gobp_res_pos <- clusterProfiler::enrichGO(gene = sigDEGs_pos,
                                              OrgDb = org_db,
                                              keyType = gene_id_format,
                                              ont = "BP",
                                              pAdjustMethod = "BH",
                                              qvalueCutoff = 0.05,
                                              universe = background_universe)
    
    if(is.null(gobp_res_pos)){gobp_res_pos <- NA}
    gobp_res_pos_df <- as.data.frame(gobp_res_pos)
    write.csv(gobp_res_pos_df, paste0(res_dir, "/GOBP_pos_", custom_label, ".csv"), row.names = F)
    
    gobp_res_neg <- clusterProfiler::enrichGO(gene = sigDEGs_neg,
                                              OrgDb = org_db,
                                              keyType = gene_id_format,
                                              ont = "BP",
                                              pAdjustMethod = "BH",
                                              qvalueCutoff = 0.05,
                                              universe = background_universe)
    
    if(is.null(gobp_res_neg)){gobp_res_neg <- NA}
    gobp_res_neg_df <- as.data.frame(gobp_res_neg)
    write.csv(gobp_res_neg_df, paste0(res_dir, "/GOBP_neg_", custom_label, ".csv"), row.names = F)
  }
  
  # KEGG, WikiPathways and reactome require entrez ids
  if(KEGG | WIKIPATHWAYS | REACTOME){
    
    message("Translating universe from ", gene_id_format, " to entrez...")
    background_universe_entrez <- entrez_df %>% 
      dplyr::filter(!!sym(gene_id_format) %in% background_universe) %>% 
      dplyr::filter(!is.na(ENTREZID)) %>% 
      dplyr::pull(ENTREZID) %>% 
      as.character()
    message(length(background_universe_entrez), " genes used as entrez id background universe")
    
    message("Translating DEGs from ", gene_id_format, " to entrez...")
    sigDEGs_pos_entrez <- entrez_df %>% 
      dplyr::filter(!!sym(gene_id_format) %in% sigDEGs_pos) %>% 
      dplyr::filter(!is.na(ENTREZID)) %>% 
      dplyr::pull(ENTREZID)
    message(length(sigDEGs_pos_entrez), " entrez DEGs enriched in treatment group mapped")
    
    sigDEGs_neg_entrez <- entrez_df %>% 
      dplyr::filter(!!sym(gene_id_format) %in% sigDEGs_neg) %>% 
      dplyr::filter(!is.na(ENTREZID)) %>% 
      dplyr::pull(ENTREZID)
    message(length(sigDEGs_neg_entrez), " entrez DEGs enriched in reference group mapped")
    
    if(org_db$packageName == "org.Hs.eg.db"){
      organism_kegg <- "hsa"
      organism_wp <- "Homo sapiens"
      organism_reactome <- "human"
    }else if(org_db$packageName == "org.Mm.eg.db"){
      organism_kegg <- "mmu"
      organism_wp <- "Mus musculus"
      organism_reactome <- "mouse"
    }
  }
  
  
  if(KEGG){
    
    message("Running hypergeometric tests against KEGG db...")
    
    kegg_res_pos <- clusterProfiler::enrichKEGG(gene = sigDEGs_pos_entrez,
                                                organism = organism_kegg,
                                                keyType = "ncbi-geneid",
                                                pAdjustMethod = "BH",
                                                qvalueCutoff = 0.05,
                                                universe = background_universe_entrez)
    
    if(is.null(kegg_res_pos)){kegg_res_pos <- NA}
    kegg_res_pos_df <- as.data.frame(kegg_res_pos)
    write.csv(kegg_res_pos_df, paste0(res_dir, "/KEGG_pos_", custom_label, ".csv"), row.names = F)
    
    kegg_res_neg <- clusterProfiler::enrichKEGG(gene = sigDEGs_neg_entrez,
                                                organism = organism_kegg,
                                                keyType = "ncbi-geneid",
                                                pAdjustMethod = "BH",
                                                qvalueCutoff = 0.05,
                                                universe = background_universe_entrez)
    
    if(is.null(kegg_res_neg)){kegg_res_neg <- NA}
    kegg_res_neg_df <- as.data.frame(kegg_res_neg)
    write.csv(kegg_res_neg_df, paste0(res_dir, "/KEGG_neg_", custom_label, ".csv"), row.names = F)
    
  }
  
  
  if(WIKIPATHWAYS){
    
    message("Running hypergeometric tests against WikiPathways...")
    
    wiki_res_pos <- clusterProfiler::enrichWP(gene = sigDEGs_pos_entrez,
                                              organism = organism_wp,
                                              pAdjustMethod = "BH",
                                              qvalueCutoff = 0.05,
                                              universe = background_universe_entrez)
    
    if(is.null(wiki_res_pos)){wiki_res_pos <- NA}
    wiki_res_pos_df <- as.data.frame(wiki_res_pos)
    write.csv(wiki_res_pos_df, paste0(res_dir, "/WIKIPATHWAYS_pos_", custom_label, ".csv"), row.names = F)
    
    
    wiki_res_neg <- clusterProfiler::enrichWP(gene = sigDEGs_neg_entrez,
                                              organism = organism_wp,
                                              pAdjustMethod = "BH",
                                              qvalueCutoff = 0.05,
                                              universe = background_universe_entrez)
    
    if(is.null(wiki_res_neg)){wiki_res_neg <- NA}
    wiki_res_neg_df <- as.data.frame(wiki_res_neg)
    write.csv(wiki_res_neg_df, paste0(res_dir, "/WIKIPATHWAYS_neg_", custom_label, ".csv"), row.names = F)
    
  }
  
  
  if(REACTOME){
    
    message("Running hypergeometric test against REACTOMEdb...")
    
    reactome_res_pos <- ReactomePA::enrichPathway(gene = sigDEGs_pos_entrez,
                                                  organism = organism_reactome, 
                                                  pAdjustMethod = "BH",
                                                  qvalueCutoff = 0.05,
                                                  universe = background_universe_entrez)
    
    if(is.null(reactome_res_pos)){reactome_res_pos <- NA}
    reactome_res_pos_df <- as.data.frame(reactome_res_pos)
    write.csv(reactome_res_pos_df, paste0(res_dir, "/REACTOME_pos_", custom_label, ".csv"), row.names = F)
    
    reactome_res_neg <- ReactomePA::enrichPathway(gene = sigDEGs_neg_entrez,
                                                  organism = organism_reactome, 
                                                  pAdjustMethod = "BH",
                                                  qvalueCutoff = 0.05,
                                                  universe = background_universe_entrez)
    
    if(is.null(reactome_res_neg)){reactome_res_neg <- NA}
    reactome_res_neg_df <- as.data.frame(reactome_res_neg)
    write.csv(reactome_res_neg_df, paste0(res_dir, "/REACTOME_neg_", custom_label, ".csv"), row.names = F)
    
  }
  
  closeAllConnections()
  
  res <- list()
  
  if(exists("gobp_res_pos")){res[[1]] <- gobp_res_pos}else{res[[1]] <- "GO_not_run"}
  if(exists("gobp_res_neg")){res[[2]] <- gobp_res_neg}else{res[[2]] <- "GO_not_run"}
  if(exists("gobp_res_pos")){res[[3]] <- kegg_res_pos}else{res[[3]] <- "KEGG_not_run"}
  if(exists("gobp_res_neg")){res[[4]] <- kegg_res_neg}else{res[[4]] <- "KEGG_not_run"}
  if(exists("wiki_res_pos")){res[[5]] <- wiki_res_pos}else{res[[5]] <- "wikipathways_not_run"}
  if(exists("wiki_res_neg")){res[[6]] <- wiki_res_neg}else{res[[6]] <- "wikipathways_not_run"}
  if(exists("reactome_res_pos")){res[[7]] <- reactome_res_pos}else{res[[7]] <- "reactome_not_run"}
  if(exists("reactome_res_neg")){res[[8]] <- reactome_res_neg}else{res[[8]] <- "reactome_not_run"}
  
  names(res) <- c("GOBP_pos", "GOBP_neg", 
                  "KEGG_pos", "KEGG_neg", 
                  "WIKIPATHWAYS_pos", "WIKIPATHWAYS_neg", 
                  "REACTOME_pos", "REACTOME_neg")
  return(res)
}


batchRunEnrichmentAnalyses <- function(dds_list, gene_id_format, counts, entrez_df, org_db, outdir, GOBP = TRUE, KEGG = TRUE, WIKIPATHWAYS = TRUE, REACTOME = TRUE){
  
  # dds_list = subc_dds
  
  res_list <- list() 
  
  for(i in 1:length(dds_list)){
    
    custom_label <- names(dds_list[i])
    
    res_list[[i]] <- runEnrichmentAnalyses(dds = dds_list[[i]],
                                           gene_id_format = gene_id_format,
                                           custom_label = custom_label,
                                           counts = counts,
                                           entrez_df = symbol2ensembl2entrez, 
                                           org_db = org.Mm.eg.db,
                                           outdir = outdir)
  }
}
