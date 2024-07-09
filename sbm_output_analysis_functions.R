#helper functions hsbm/ nsbm output analysis
#ack 17.08.23
#anne-claire.kroger@curie.fr

#------------------------------------------------------------------------------#

get_nsbmclust_df <- function(nsbm_clust, anno, anno_samplecol, anno_groupcol) {
  clust_list <- apply(nsbm_clust, 2, function(clust){anno[match(clust, anno[,anno_samplecol]),anno_groupcol]})
  clust_df <- as.data.frame(matrix(nrow = ncol(clust_list), ncol = length(table(anno[,anno_groupcol]))))
  colnames(clust_df) <- names(table(anno[,anno_groupcol]))
  rownames(clust_df) <- colnames(clust_list)
  for (i in 1:ncol(clust_list)) {
    c <- match(names(table(clust_list[,i])), colnames(clust_df))
    clust_df[i,c] <- as.vector(table(clust_list[,i]))
  }
  return(clust_df)
}

#------------------------------------------------------------------------------#

clean_todist <- function(todist) {
  rownames(todist) <- todist$doc
  todist <- todist[,-c(1,2)]
  colnames(todist) <- paste0("Topic.", 1:ncol(todist))
  return(todist)
}

#------------------------------------------------------------------------------#

clean_worddist <- function(worddist, tri = FALSE) {
  if (tri == TRUE) {
    rownames(worddist) <- str_replace_all(worddist$X, "#", "")
    rownames(worddist) <- str_remove(rownames(worddist), "_P$")
    worddist <- worddist[,-1]
    colnames(worddist) <- paste0("Topic.", 1:ncol(worddist))
  } else{
    rownames(worddist) <- worddist$word_names
    worddist <- worddist[,-c(ncol(worddist))]
    colnames(worddist) <- paste0("Topic.", 1:ncol(worddist))
  }
  return(worddist)
}

#------------------------------------------------------------------------------#


get_sankey_input <- function(ft_custering) {
  #change cluster names
  for (i in 1:ncol(ft_custering)) {
    ft_custering[,i] <- paste0(paste0(rep("_", i-1), collapse = ""), ft_custering[,i]) 
  }
  
  #turn into river plot table
  ft_custering <- as.matrix(ft_custering)
  links <- ft_custering[,1:2]
  for (i in 3:ncol(ft_custering)) { 
    links <- rbind(links, ft_custering[,(i-1):i]) }
  links <- as.data.frame(links)%>%group_by_all()%>%dplyr::count()
  
  colnames(links) <- c("source", "target", "value")
  
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  return(list("links"=links, "nodes"=nodes))
}

#------------------------------------------------------------------------------#

read_sbm_output <- function(folder_path, level, type = "", omics_layers = 1, force = FALSE){
  
  if (omics_layers == 3) {
    
    clust <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_clusters.csv"))
    if(dim(clust)[2] > 20 & force == FALSE){
      stop(paste("There are more then 20 clusters,", dim(clust)[2], "in total, are you sure you want to load all this? If yes set force=TRUE"))
      }

    
    todist_rna <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_topic-dist.csv"))
    todist_prot <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_2_metadatum-dist.csv"))
    todist_phos <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_3_metadatum-dist.csv"))
    todist_rna <- clean_todist(todist_rna)
    todist_prot <- clean_todist(todist_prot)
    todist_phos <- clean_todist(todist_phos) 

    worddist_rna <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_word-dist.csv"))
    worddist_prot <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_2_keyword-dist.csv"))
    worddist_phos <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_3_keyword-dist.csv"))
    worddist_rna <- clean_worddist(worddist_rna, tri = TRUE)
    worddist_prot <- clean_worddist(worddist_prot, tri = TRUE)
    worddist_phos <- clean_worddist(worddist_phos, tri = TRUE) 
        
    topics_rna <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_topics.csv"))
    topics_prot <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_2_metadata.csv"))
    topics_phos <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_3_metadata.csv"))
    
    summary <- list("n_clust"= dim(clust)[2], 
                    "n_topics_rna"= dim(topics_rna)[2], 
                    "max_topic_length_rna"=dim(topics_rna)[1], 
                    "average_topic_length_rna" = mean(apply(topics_rna, 2, function(x) {length(which(x!=""))})),
                    "n_topics_prot"= dim(topics_prot)[2], 
                    "max_topic_length_prot"=dim(topics_prot)[1],
                    "average_topic_length_prot" = mean(apply(topics_prot, 2, function(x) {length(which(x!=""))})),
                    "n_topics_phos"= dim(topics_phos)[2], 
                    "max_topic_length_phos"=dim(topics_phos)[1],
                    "average_topic_length_phos" = mean(apply(topics_phos, 2, function(x) {length(which(x!=""))})),
                    "n_samples"=dim(todist_rna)[1])
    outlist <- list("summary" = summary, 
                    "clust" = clust, 
                    "todist_rna" = todist_rna, 
                    "topics_rna"=topics_rna, 
                    "worddist_rna"=worddist_rna,
                    "todist_prot" = todist_prot, 
                    "topics_prot"=topics_prot, 
                    "worddist_prot"=worddist_prot,
                    "todist_phos" = todist_phos, 
                    "topics_phos"=topics_phos, 
                    "worddist_phos"=worddist_phos)
  } 
  
  if (omics_layers == 2) {
    
    clust <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_clusters.csv"))
    if(dim(clust)[2] > 20 & force == FALSE){
      stop(paste("There are more then 20 clusters,", dim(clust)[2], "in total, are you sure you want to load all this? If yes set force=TRUE"))
    }
    
    
    todist_rna <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_topic-dist.csv"))
    todist_prot <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_2_metadatum-dist.csv"))
    todist_rna <- clean_todist(todist_rna)
    todist_prot <- clean_todist(todist_prot)
    
    worddist_rna <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_word-dist.csv"))
    worddist_prot <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_2_keyword-dist.csv"))
    worddist_rna <- clean_worddist(worddist_rna, tri = TRUE)
    worddist_prot <- clean_worddist(worddist_prot, tri = TRUE)
    
    topics_rna <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_topics.csv"))
    topics_prot <- read.csv(paste0(folder_path, "/trisbm_level_", level,"_kind_2_metadata.csv"))
    
    summary <- list("n_clust"= dim(clust)[2], 
                    "n_topics_rna"= dim(topics_rna)[2], 
                    "max_topic_length_rna"=dim(topics_rna)[1],
                    "average_topic_length_rna" = mean(apply(topics_rna, 2, function(x) {length(which(x!=""))})),
                    "n_topics_prot"= dim(topics_prot)[2], 
                    "max_topic_length_prot"=dim(topics_prot)[1],
                    "average_topic_length_prot" = mean(apply(topics_prot, 2, function(x) {length(which(x!=""))})),
                    "n_samples"=dim(todist_rna)[1])
    outlist <- list("summary" = summary, 
                    "clust" = clust, 
                    "todist_rna" = todist_rna, 
                    "topics_rna"=topics_rna, 
                    "worddist_rna"=worddist_rna,
                    "todist_prot" = todist_prot, 
                    "topics_prot"=topics_prot, 
                    "worddist_prot"=worddist_prot)
  } 
  
  if (omics_layers == 1) {
  clust <- read.csv(paste0(folder_path, "/topsbm_level_", level,"_clusters.csv"))
  todist <- read.csv(paste0(folder_path, "/topsbm_level_", level,"_topic-dist.csv"))
  todist <- clean_todist(todist)
  topics <- read.csv(paste0(folder_path, "/topsbm_level_", level,"_topics.csv"))
  
  summary <- list("n_clust"= dim(clust)[2], 
                  "n_topics"= dim(topics)[2], 
                  "max_topic_length"=dim(topics)[1],
                  "average_topic_length" = mean(apply(topics, 2, function(x) {length(which(x!=""))})),
                  "n_samples"=dim(todist)[1])
  
  f <- paste0(folder_path, "/", type, "_l", level, "_worddist.csv")
  if(type != "" & file.exists(f)) {
    worddist <- read.csv(f)
    worddist <- clean_worddist(worddist)
    summary$n_features <- dim(worddist)[1]
    outlist <- list("summary" = summary, "clust" = clust, "todist" = todist, "topics"=topics, "worddist"=worddist)
  } else {
    outlist <- list("summary" = summary, "clust" = clust, "todist" = todist, "topics"=topics)
  }
  }
  
  
  return(outlist)
}

#------------------------------------------------------------------------------#

analyse_topics <- function(sbm_output_list, omics_layers = 1){
  
  if (omics_layers == 1) {
    topic_imp <- colSums(sbm_output_list$todist)
    mean_topic_imp <- colSums(sbm_output_list$todist)/ nrow(sbm_output_list$todist)
    topic_length <- apply(sbm_output_list$worddist, 2, function(x) {length(which(x != 0))})
    
    topic_info <- data.frame(topic_length, topic_imp, mean_topic_imp, row.names = names(topic_length))
    
    clust_l <- lapply(as.list(sbm_output_list$clust), function(x){x[nzchar(x)]})
    clust_mean_topic_imp_l <- lapply(clust_l, function(clust_samp){
      return(colSums(sbm_output_list$todist[clust_samp,])/ length(clust_samp))
    })
    clust_mean_topic_imp <- as.data.frame(clust_mean_topic_imp_l)
    
    clust_ratio_topic_imp <- clust_mean_topic_imp/mean_topic_imp
    
    clust_centered_topic_imp <- clust_mean_topic_imp - mean_topic_imp
    
    colnames(clust_mean_topic_imp) <- paste("mean_imp_",colnames(clust_mean_topic_imp), sep = "")
    colnames(clust_ratio_topic_imp) <- paste("ratio_imp_",colnames(clust_ratio_topic_imp), sep = "")
    colnames(clust_centered_topic_imp) <- paste("centered_imp_",colnames(clust_centered_topic_imp), sep = "")
    
    topic_info <- cbind(topic_info, clust_mean_topic_imp, clust_ratio_topic_imp, clust_centered_topic_imp)
    
    clust_cent <- clust_centered_topic_imp
    clust_cent$topic <- rownames(clust_cent)
    clust_cent$topic_length <- topic_info$topic_length
    clust_cent$topic_mean_imp <- topic_info$mean_topic_imp
    clust_cent <- melt(clust_cent, id.vars = c("topic", "topic_length", "topic_mean_imp"))
    clust_cent$cluster <- unlist(lapply(str_split(clust_cent$variable, "_imp_"), "[[", 2))
    
    return(list("info_all"=topic_info, "info_centered"=clust_cent))
  }
  
  if (omics_layers == 2) {
    #overall topic info
    topic_imp_list <- lapply(sbm_output_list[c(3,6)], colSums)
    mean_topic_imp_list <- list()
    for (i in c(3,6)) {
      mean_topic_imp_list[[i/3]] <- colSums(sbm_output_list[[i]])/ nrow(sbm_output_list[[i]])
    }
    topic_length_list <- lapply(sbm_output_list[c(5,8)], function(x){apply(x, 2, function(y) {length(which(y != 0))})})
    
    #topic info for each cluster
    clust_l <- lapply(as.list(sbm_output_list$clust), function(x){x[nzchar(x)]})
    
    topic_info_l <- list()
    clust_cent_l <- list()
    
    for (i in c(1:2)) {
      
      topic_info_l[[i]] <- data.frame(topic_length_list[[i]], 
                                      topic_imp_list[[i]], 
                                      mean_topic_imp_list[[i]], 
                                      row.names = names(topic_length_list[[i]]))
      colnames(topic_info_l[[i]]) <- c("topic_length", "topic_imp", "mean_topic_imp")
      
      
      clust_mean_topic_imp_l <- lapply(clust_l, function(clust_samp){
        return(colSums(sbm_output_list[[i*3]][clust_samp,])/ length(clust_samp))
      })
      clust_mean_topic_imp <- as.data.frame(clust_mean_topic_imp_l)
      clust_ratio_topic_imp <- clust_mean_topic_imp/mean_topic_imp_list[[i]]
      clust_centered_topic_imp <- clust_mean_topic_imp - mean_topic_imp_list[[i]]
      colnames(clust_mean_topic_imp) <- paste("mean_imp_",colnames(clust_mean_topic_imp), sep = "")
      colnames(clust_ratio_topic_imp) <- paste("ratio_imp_",colnames(clust_ratio_topic_imp), sep = "")
      colnames(clust_centered_topic_imp) <- paste("centered_imp_",colnames(clust_centered_topic_imp), sep = "")
      
      topic_info_l[[i]] <- cbind(topic_info_l[[i]] , 
                                 clust_mean_topic_imp,
                                 clust_ratio_topic_imp,
                                 clust_centered_topic_imp)
      
      
      clust_cent <- clust_centered_topic_imp
      clust_cent$topic <- rownames(clust_cent)
      clust_cent$topic_length <- topic_info_l[[i]]$topic_length
      clust_cent$topic_mean_imp <- topic_info_l[[i]]$mean_topic_imp
      clust_cent <- melt(clust_cent, id.vars = c("topic", "topic_length", "topic_mean_imp"))
      clust_cent$cluster <- unlist(lapply(str_split(clust_cent$variable, "_imp_"), "[[", 2))
      clust_cent_l[[i]] <- clust_cent
    }
    
    return(list("rna"=list("info_all"=topic_info_l[[1]], "info_centered"=clust_cent_l[[1]]), 
                "prot"=list("info_all"=topic_info_l[[2]], "info_centered"=clust_cent_l[[2]])))
  }
  
  if (omics_layers == 3) {
    #overall topic info
    topic_imp_list <- lapply(sbm_output_list[c(3,6,9)], colSums)
    mean_topic_imp_list <- list()
    for (i in c(3,6,9)) {
      mean_topic_imp_list[[i/3]] <- colSums(sbm_output_list[[i]])/ nrow(sbm_output_list[[i]])
    }
    topic_length_list <- lapply(sbm_output_list[c(5,8,11)], function(x){apply(x, 2, function(y) {length(which(y != 0))})})
    
    #topic info for each cluster
    clust_l <- lapply(as.list(sbm_output_list$clust), function(x){x[nzchar(x)]})
    
    topic_info_l <- list()
    clust_cent_l <- list()
    
    for (i in c(1:3)) {
      
      topic_info_l[[i]] <- data.frame(topic_length_list[[i]], 
                                      topic_imp_list[[i]], 
                                      mean_topic_imp_list[[i]], 
                                      row.names = names(topic_length_list[[i]]))
      colnames(topic_info_l[[i]]) <- c("topic_length", "topic_imp", "mean_topic_imp")
      
      
      clust_mean_topic_imp_l <- lapply(clust_l, function(clust_samp){
        return(colSums(sbm_output_list[[i*3]][clust_samp,])/ length(clust_samp))
      })
      clust_mean_topic_imp <- as.data.frame(clust_mean_topic_imp_l)
      clust_ratio_topic_imp <- clust_mean_topic_imp/mean_topic_imp_list[[i]]
      clust_centered_topic_imp <- clust_mean_topic_imp - mean_topic_imp_list[[i]]
      colnames(clust_mean_topic_imp) <- paste("mean_imp_",colnames(clust_mean_topic_imp), sep = "")
      colnames(clust_ratio_topic_imp) <- paste("ratio_imp_",colnames(clust_ratio_topic_imp), sep = "")
      colnames(clust_centered_topic_imp) <- paste("centered_imp_",colnames(clust_centered_topic_imp), sep = "")
      
      topic_info_l[[i]] <- cbind(topic_info_l[[i]] , 
                                 clust_mean_topic_imp,
                                 clust_ratio_topic_imp,
                                 clust_centered_topic_imp)

      
      clust_cent <- clust_centered_topic_imp
      clust_cent$topic <- rownames(clust_cent)
      clust_cent$topic_length <- topic_info_l[[i]]$topic_length
      clust_cent$topic_mean_imp <- topic_info_l[[i]]$mean_topic_imp
      clust_cent <- melt(clust_cent, id.vars = c("topic", "topic_length", "topic_mean_imp"))
      clust_cent$cluster <- unlist(lapply(str_split(clust_cent$variable, "_imp_"), "[[", 2))
      clust_cent_l[[i]] <- clust_cent
    }
    
    return(list("rna"=list("info_all"=topic_info_l[[1]], "info_centered"=clust_cent_l[[1]]), 
                "prot"=list("info_all"=topic_info_l[[2]], "info_centered"=clust_cent_l[[2]]),
                "phos"=list("info_all"=topic_info_l[[3]], "info_centered"=clust_cent_l[[3]])))
    }
}

#------------------------------------------------------------------------------#

topic_length_corr <- function(stats_topics){
  n <- ncol(stats_topics$info_all)
  for (i in (n-((n-3)/3)+1):n) {
    print(colnames(stats_topics$info_all)[i])
    print(cor(stats_topics$info_all$topic_length, abs(stats_topics$info_all[,i])))
  }
}

#------------------------------------------------------------------------------#

rand_index_sbmcluster <- function(nsbmclust_df, adjusted_rand=FALSE){
  colnames(nsbmclust_df) <- as.numeric(as.factor(colnames(nsbmclust_df)))
  nsbmclust_df$cluster <- colnames(nsbmclust_df)[apply(nsbmclust_df, 1,which.max)]
  nsbmclust_df <- melt(nsbmclust_df, by = "cluster")
  nsbmclust_df <- nsbmclust_df[!is.na(nsbmclust_df$value),]
  
  g1 <- unlist(apply(nsbmclust_df, 1, function(x){return(rep(x[1], x[3]))}))
  g2 <- unlist(apply(nsbmclust_df, 1, function(x){return(rep(x[2], x[3]))}))
  
  if(adjusted_rand==TRUE){ return(adj.rand.index(as.numeric(g1),as.numeric(g2)))
  } else {
    return(round(rand.index(as.numeric(g1),as.numeric(g2)), 3))
    }
  
}


