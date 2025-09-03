# ----------------------
# GSEA tutorial
# ----------------------

# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(123456)

# Set project library
.libPaths('C:/Users/laura/Documents/Biostatsquid/Scripts/R4.2.3')
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)

# Set relevant paths
list.files()
in_path <- "Datasets/"
out_path <- "PEA/Results/"

# Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Get all the genes in your dataset and assign them to my_genes 
my_genes <- df$gene_symbol

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# Copy the .gmt file to your folder, in my case, its 'PEA/Background_genes/'
# Then read in the .gmt file
gmt_file <- "PEA/Background_genes/c2.cp.kegg.v7.5.1.symbols.gmt"
hidden <- unique(unlist(gmt))

# Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
              nrow = length(hidden), ncol = length(gmt))
for (i in 1:dim(mat)[2]){
  mat[,i] <- as.numeric(hidden %in% gmt[[i]])
}
#Subset to the genes that are present in our data to avoid bias
hidden1 <- intersect(genes_in_data, hidden)
mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
# And get the list again using the function we previously defined
final_list <- matrix_to_list(mat)

# Set relevant paths
list.files()
in_path <- "Datasets/"
out_path <- "PEA/Results/"
bg_path <- "PEA/Background_genes/"

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

# Analysis ====================================================

## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(paste0(in_path, 'severevshealthy_degresults.csv'), row.names = 1)

## 2. Prepare background genes -----------------------------------------------

# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$gene_symbol
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
bg_genes <- prepare_gmt(gmt_files[2], genes_in_data, savefile = FALSE)

# You can use a list of filenames, and just change [1], [2], [3]...
bg_genes <- prepare_gmt(gmt_files[1], my_genes, savefile = FALSE)
            
# Or write it out
bg_genes <- prepare_gmt("PEA/Background_genes/c5.go.bp.v7.5.1.symbols.gmt", my_genes, savefile = FALSE)

rankings <- sign(df$log2fc)*(-log10(df$pval)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df$gene_symbol # genes as names#

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation


## 6. Check results ------------------------------------------------------
# Top 6 enriched pathways (ordered by p-val)
head(GSEAres[order(pval), ])

plotEnrichment(bg_genes[['REACTOME_FORMATION_OF_FIBRIN_CLOT_CLOTTING_CASCADE']],
               rankings) + 
  labs(title = 'Reactome pathway: Formation of fibrin clot / clotting cascade') + 
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) +
  scale_y_continuous('Enrichment score (ES)') +
  geom_line(col = 'purple', linewidth = 2) 









