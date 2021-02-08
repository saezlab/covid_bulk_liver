

run_normalize <- function(tdf, dir){
  # Generate targets data-frame
  gene_names <- rownames(tdf)
  # Trim violin plots
  tdf[log2(tdf) < 4 ] <- NA
  # Remove rows (genes) which aren't well-measured in enough samples
  msk <- rowSums(is.na(tdf)) < 2
  tdf <- tdf[msk,]
  gene_names <- gene_names[msk]
  #now we can normalise the cleaned dataframe using vsn
  fit <- vsnMatrix(as.matrix(tdf)) #train vsn parameters
  count_df_vsn <- as.data.frame(vsn::predict(fit,as.matrix(tdf)))
  to_write <- count_df_vsn
  to_write$gene <- gene_names
  to_write <- to_write[,c(length(to_write[1,]),1:(length(to_write[1,])-1))]
  # Create dir where to store results
  path <- file.path("..", "results", dir)
  dir.create(path, showWarnings = FALSE, recursive = T)
  write_csv(to_write, file.path(path, "count_df_vsn.csv"))
}

run_tf <- function(dir){
  path <- file.path("..", "results", dir)
  ## We read the normalised counts and the experimental design 
  Normalised_counts <- read_csv(file.path(path, "count_df_vsn.csv"))
  Normalised_counts_matrix <- Normalised_counts %>% 
    dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
    tibble::column_to_rownames(var = "gene") %>% 
    as.matrix()
  ## We load Dorothea Regulons
  data(dorothea_hs, package = "dorothea")
  regulons <- dorothea_hs %>%
    dplyr::filter(confidence %in% c("A", "B","C"))
  tf_activities_counts <- 
    dorothea::run_viper(Normalised_counts_matrix, regulons,
                        options =  list(minsize = 5, eset.filter = FALSE, 
                                        cores = 1, verbose = FALSE, method = c("scale")))
  tf_out <- as_tibble(tf_activities_counts, rownames=NA) %>% 
    rownames_to_column('tf')
  write_csv(tf_out, 
            file.path(path, "tf_act.csv"))
}

run_pathway <- function(dir){
  path <- file.path("..", "results", dir)
  ## We read the normalised counts and the experimental design 
  Normalised_counts <- read_csv(file.path(path, "count_df_vsn.csv"))
  Normalised_counts_matrix <- Normalised_counts %>% 
    dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
    tibble::column_to_rownames(var = "gene") %>% 
    as.matrix()
  PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                    organism="Human", top = 100)
  pathway_out <- as_tibble(t(PathwayActivity_counts), rownames=NA) %>% 
    rownames_to_column('pathway')
  write_csv(pathway_out, 
            file.path(path, "pw_act.csv"))
}

get_plot_df <- function(model, type_name){
  coef <- model$coefficients[,type_name]
  pvals <- model$p.value[,type_name]
  names <- names(coef)
  df_plot <- tibble(
    func=names,
    coef=coef,
    pvals=pvals,
    thr=ifelse(pvals < 0.05, 'sign', 'non-sign')
  ) %>% arrange(pvals) %>% rowid_to_column('order')
}

limma_volcano <- function(model){
  t_names <- model$coefficients %>% colnames()
  p <- lapply(t_names, function(t_name){
    df_plot <- get_plot_df(model, t_name)
    p <- ggplot(df_plot) +
      geom_point(aes(x=coef, y=-log10(pvals), colour=thr)) +
      geom_text_repel(aes(x = coef, y = -log10(pvals), 
                          label = ifelse(order < 15 & thr == 'sign', 
                                         func, ""))) +
      ggtitle(t_name) +
      xlab("coeff") + 
      ylab("-log10(p-value)") +
      theme_classic() +
      theme(legend.position = "none",
            text = element_text(size=15),
            plot.title = element_text(size = rel(0.9), hjust = 0.5),
            axis.title = element_text(size = rel(1)))
  })
  max_y <- max(unlist(lapply(p, function(x){-log10(x$data$pvals)})))
  max_y <- max_y + (max_y * 0.05)
  max_x <- max(abs(unlist(lapply(p, function(x){x$data$coef}))))
  max_x <- max_x + (max_x * 0.05)
  p <- lapply(p, function(x){
    x + coord_cartesian(xlim = c(-max_x, max_x), 
                        ylim = c(0, max_y))
  })
  #p[[1]] <- p[[1]] + ggtitle('COVID+/PCR- vs COVID-/PCR-')
  #p[[2]] <- p[[2]] + ggtitle('COVID+/PCR+ vs COVID-/PCR-')
  #p[[3]] <- p[[3]] + ggtitle('Difference')
  return(p)
}

GSE_analysis = function(geneList,Annotation_DB){
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]
  
  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")
  
  DB_genecontent = length(unique(unlist(Annotation_DB)))
  
  GenesDB = DB_genecontent 
  SelectedGenes = length(geneList)
  
  for(gset in rownames(ResultsDF)){
    GP = length(Annotation_DB[[gset]])
    GL = length(intersect(Annotation_DB[[gset]],geneList))
    
    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }
  
  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]
  
  ResultsDF = ResultsDF %>% 
    rownames_to_column("gset") %>% 
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"), 
              as.numeric) %>% 
    dplyr::arrange(corr_p_value,GenesInList)
  
  return(ResultsDF)
  
}
