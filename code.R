##########source: https://www.biostars.org/p/449688/
#Istall phyloseq, DESeq2, and curatedMetagenomicData
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("DESeq2")
#BiocManager::install("curatedMetagenomicData")

#Loading the library
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(curatedMetagenomicData)

#Getting the source file
mphlanin <- read.csv("merged_abundance_table_reformatted.txt", sep = "\t", strip.white = T, stringsAsFactors = F, row.names = 1)
metadata <- read.delim("metadata.txt", header=TRUE, sep = "\t")
metadatadf <- data.frame(metadata)
row.names(metadatadf) <- metadatadf$sample
samples_df <- metadatadf %>% select (-sample)
sample <- sample_data(samples_df)

phyloseqin= metaphlanToPhyloseq(mphlanin, metadat = sample)
phyloseqin
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 229 taxa and 22 samples ]
#sample_data() Sample Data:       [ 22 samples by 1 sample variables ]
#tax_table()   Taxonomy Table:    [ 229 taxa by 7 taxonomic ranks ]

#' convert the output of a metaphlan2_taxonomic_table_joined.tsv object to a otu_table + tax_table object
#' 
#' 
#' @param phyloseq object 
#' @param 
#' @export
#' @examples
#' mtph_tbru_phy <- metaphlanToPhyloseq(tax = mtph_tbru, split = "|")
metaphlanToPhyloseq <- function(
  tax,
  metadat=NULL,
  simplenames=TRUE,
  roundtointeger=FALSE,
  split="|"){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ## if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}