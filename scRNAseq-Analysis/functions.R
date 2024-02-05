# ~~~ Functions File ~~~
# Desc: Contains all functions called in downstream analysis
# 
# Best Practices:
# * Carefully define the purpose of each function, any arguments it takes, and
# *   and anything it returns
# * Single Responsibility: Each function should perform one job and do it well
# * Manage Scope Carefully
# * Parametrize your functions to make them configurable and reusable

LoadPackages <- function() {
  # Desc: Loads all necessary packages cached in compiled CTHC transferred file
  #       Args: None
  #       Returns: Void
  library(jsonlite)
  library(tidyverse)
  library(DropletUtils)
  library(Seurat)
  library(ShinyCell)
  library(Matrix)
  library(scales)
  library(rjson)
  library(R2HTML)
  library(DT)
  library(SeuratDisk)
  # to install SeuratDisk, which is not on CRAN, use the following:
  # if (!requireNamespace("remotes", quietly = TRUE)) {
  #   install.packages("remotes")
  # }
  # remotes::install_github("mojaveazure/seurat-disk")
}

ImportCountData <- function(sample_prefix_array) {
  # Desc: Iterates through each sample filename prefix and creates a raw read mtx
  # Args: Array of sample prefixes specified in JSON Params
  # Returns: List of raw matrices for later processing
  
  raw_mtx_list <- list()
  
  for(sample in sample_prefix_array){
    # load raw data matrix using the readMM function from the Matrix package
    raw_mtx <- readMM(paste0(sample,'matrix.mtx.gz'))
    
    # load genes
    genes <- read.csv(paste0(sample,'features.tsv.gz'), sep = '\t', header = F)
    
    # add ensemble gene_ids to the data matrix as rownames
    rownames(raw_mtx) <- genes[,1] 
    
    # add cell barcodes as column names
    barcodes <- read.csv(paste0(sample,'barcodes.tsv.gz'), sep = '\t', header = F) 
    colnames(raw_mtx) <- barcodes[,1]
    
    # Add the read and processed mtx to the list
    raw_mtx_list[[paste0(sample,".raw_mtx")]] <- raw_mtx
  }
  
  #return the list of matrices
  return(raw_mtx_list)
}
