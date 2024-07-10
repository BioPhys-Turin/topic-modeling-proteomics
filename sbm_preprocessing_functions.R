#helper functions preprocessing CPTAC data for SBM analysis
#ack 17.08.23
#anne-claire.kroger@curie.fr


#------------------------------------------------------------------------------#

annotate_cptac <- function(cptac_df, phospho = FALSE){
  
  if (phospho) {
    #clean ensembl ID
    localisation <- unlist(lapply(str_split(cptac_df$idx, "\\|"), "[[", 3))
    gene_ensid <- unlist(lapply(str_split(cptac_df$idx, "\\|"), "[[", 1))
    gene_ensid <- str_replace(gene_ensid, "\\..*$", "")
    
    #annotate with symbol and transcript type
    annots <- select(org.Hs.eg.db, keys=gene_ensid, 
                     columns=c("SYMBOL", "GENETYPE"), keytype="ENSEMBL")
    annots <- unique(annots)
    annots <- annots[match(gene_ensid, annots$ENSEMBL),]
    annots$new_names <- paste0(annots$SYMBOL, "_", localisation)
    
  }else{
    #clean ensembl ID
    gene_ensid <- str_replace(cptac_df$idx, "\\..*$", "")
    
    #annotate with symbol and transcript type
    annots <- AnnotationDbi::select(org.Hs.eg.db, keys=gene_ensid, 
                                    columns=c("SYMBOL", "GENETYPE"), keytype="ENSEMBL")
    annots <- unique(annots)
    annots <- annots[match(gene_ensid, annots$ENSEMBL),]
    
  }
  cptac_df <- cbind(annots, cptac_df)
  return(cptac_df)
  
}

#------------------------------------------------------------------------------#

clean_annotated_cptac <- function(cptac_df, phospho = FALSE) {
  
  #throw out NAs that don't have a symbol assigned
  cptac_df <- cptac_df[!is.na(cptac_df$SYMBOL),]
  
  
  
  if (phospho) {
    cptac_df <- cptac_df[rowSums(cptac_df[,6:ncol(cptac_df)], na.rm = TRUE) !=0,]
    
    #get doubles still left
    dubs <- names(which(table(cptac_df$new_names)>1))
    dubs <- cptac_df[cptac_df$new_names %in% dubs,]
    #remove everything that does not have the ending 1
    dubs <- dubs$idx[unlist(lapply(str_split(dubs$idx, "\\|"), "[[", 5)) != 1]
    cptac_df <- cptac_df[!cptac_df$idx %in% dubs,]
    
    rownames(cptac_df) <- cptac_df$new_names
    cptac_df <- cptac_df[,6:ncol(cptac_df)]
    
  } else {
    cptac_df <- cptac_df[rowSums(cptac_df[,5:ncol(cptac_df)], na.rm = TRUE) > 0,]
    cptac_df <- cptac_df[!str_detect(cptac_df$idx, "PAR"),]
    cptac_df <- cptac_df[cptac_df$GENETYPE == "protein-coding",]
    
    
    dubs <- cptac_df[cptac_df$SYMBOL %in% names(which(table(cptac_df$SYMBOL) != 1)), ]
    sums <- cbind(dubs[,1:4], "sum" = rowSums(dubs[,5:ncol(dubs)], na.rm = TRUE))
    out <-  sapply(unique(sums$SYMBOL), aggregate_max, sums = sums)
    out <- dubs$idx[!dubs$idx %in% array(out)]
    cptac_df <- cptac_df[!cptac_df$idx %in% out,]
    
    rownames(cptac_df) <- cptac_df$SYMBOL
    cptac_df <- cptac_df[,5:ncol(cptac_df)]
  }
  return(cptac_df)
}

#------------------------------------------------------------------------------#

aggregate_max <- function(symbol, sums) {
  d <- sums[sums$SYMBOL== symbol,]
  return(d$idx[which.max(d$sum)])}

#------------------------------------------------------------------------------#

prep_nsbm <- function(cptac_df, samples, log2trans = TRUE, shift = FALSE){
  if(shift){
    cptac_df <- cptac_df - 10
    cptac_df[cptac_df < 0] <- 0
  }
  if (log2trans){
    cptac_df <- round(2^cptac_df)}
  cptac_df <- cptac_df[rowSums(cptac_df, na.rm = TRUE) != 0,]
  cptac_df <- data.matrix(cptac_df,  rownames.force = NULL)
  cptac_df[is.infinite(cptac_df)] <- NA
  cptac_df <- as.data.frame(cptac_df[, samples])
  return(cptac_df)
  
}

#------------------------------------------------------------------------------#
