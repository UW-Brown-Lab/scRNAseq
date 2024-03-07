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
  #library(DT)
  library(SingleCellExperiment)
  library(scDblFinder)
  #library(SeuratDisk)
  # to install SeuratDisk, which is not on CRAN, use the following:
  # if (!requireNamespace("remotes", quietly = TRUE)) {
  #   install.packages("remotes")
  # }
  # remotes::install_github("mojaveazure/seurat-disk")
}

ImportCountData <- function(samples,dir_info) {
    # Desc: Iterates through each sample file path and reads 10X counts for
    #       the raw and filtered data, then combined them into a list
    # Args: samples: The array of samples specified in the params
    #       dir_info: Info about the directory setup from params
    # Returns: A list of items, each containing the raw droplets and filtered
    #          counts
  
  # Get output directory from the alignment step from params (based on STAR)
  output_path <- paste0("./",dir_info$alignment_output_dir,"/star/")
  
  # Create list of lists to output
  output_list <- list()
  
  for(i in 1:nrow(samples)){

        # Get Sample
        sample <- samples[i,]

        # Get Sample Name
        sample_name <- sample$name

        # Get Sample File Path (based on STAR output)
        sample_path <- paste0(output_path,
                              sample_name,"/",sample_name,".Solo.out/Gene/")
        
        # Read raw counts
        droplets.raw <- Read10X(paste0(sample_path,"raw/"))
        
        # Read filtered counts
        counts.filt <- Read10X(paste0(sample_path,"filtered/"))
        
        assign(sample_name,
               list(
                    filtered = counts.filt,
                    raw = droplets.raw
                    )
               )
        
        output_list[[sample_name]] <- get(sample_name)
  }
  
  return(output_list)
}

DecontaminateSoup <- function(data_list, dim_reduct_config) {
  # Desc: Estimates cell free mRNA contamination with SoupX and performs a
  #       correction
  # Args: data_list:  a list of lists, one item for each sample, each containing
  #                   raw and filtered counts
  #       dim_reduct_config:  settings imported from params
  # Returns: a list of adjusted count matrices, one per sample
  
  # Create list to return
  output_list <- list()
  
  # Iterate through each sample
  for(i in seq_along(data_list)){
    # Get this sample item from the list
    sample <- data_list[[i]]
    
    # Get sample name
    sample_name <- names(data_list)[i]
    
    print(paste0("Executing on sample ", sample_name))
    
    # Get counts
    counts.raw <- sample$raw
    counts.filt <- sample$filtered
    
    # Create seurat object from filtered counts for clustering
    filt_seurat <- CreateSeuratObject(counts.filt)
    
    # Perform Normalization, Scaling, and Variable Feature selection with SCTransform
    filt_seurat <- SCTransform(filt_seurat, verbose = F)
    
    # Perform PCA and UMAP
    filt_seurat <- RunPCA(filt_seurat, verbose = F)
    filt_seurat <- RunUMAP(filt_seurat, dims = 1:dim_reduct_config$PCA_component_count)
    
    # Find Neighbors, then clusters for SoupX
    filt_seurat <- FindNeighbors(filt_seurat, dims=1:dim_reduct_config$PCA_component_count)
    filt_seurat <- FindClusters(filt_seurat, verbose = T)
    
    # Get clustering and dimensional reduction data from seurat object
    meta.data <- filt_seurat@meta.data
    umap.data <- filt_seurat@reductions$umap@cell.embeddings
    
    # Create SoupX SoupChannel object
    SoupX_channel <- SoupChannel(counts.raw, counts.filt)
    
    # Apply seurat clusters
    SoupX_channel <- setClusters(
      SoupX_channel, 
      setNames(meta.data$seurat_clusters,rownames(meta.data))
      )
    
    # Apply dimensional reduction for visualizations
    SoupX_channel <- setDR(SoupX_channel, umap.data)
    
    # Estimate Contamination and Remove it
    SoupX_channel <- autoEstCont(SoupX_channel, verbose = T)
    corrected_matrix <- adjustCounts(SoupX_channel, verbose = T)
    
    # Add corrected count matrix to output
    output_list[[sample_name]] <- corrected_matrix
  }
  
  return(output_list)
}

CountMitochondrialFeatures <- function(merged_seurat) {
  # Desc: Counts percent of HUMAN mitochondrial features for each cell
  #       and appends this to seurat object
  # Args: merged_seurat: a merged seurat object generated by 
  #       CreateMergedSeurat()
  # Returns: A seurat object annotated with percent.mt (human only)
  
  
    # Get percentage of Mitochondrial features and append to object
    # NOTE: Currently only works with HUMAN features!
    merged_seurat[["percent.mt"]] <- PercentageFeatureSet(object = merged_seurat,
                                                          pattern = "^chrM")
    
  # Return annotated seurat list
  return(merged_seurat)
}

FilterSeuratObject <- function(merged_seurat, thresholds) {
  # Desc: Subsets seurat object by specified parameters
  # Args: merged_seurat: A merged suerat object.
  #       thresholds: the params$filter_options$thresholds
  # Returns: A subsetted seurat object and before/after vln plots
  
  # Count Mitochondrial Features
  mt_seurat <- CountMitochondrialFeatures(merged_seurat)
  
  # Get number of cells before filtering
  total_cells <- ncol(mt_seurat)
    
  # Get number of features before filtering
  total_features <- nrow(mt_seurat)
    
  # Filter out cells
  filt_seurat <- subset(mt_seurat, 
                            subset = nCount_RNA < thresholds$nCount_RNA_upper & 
                              nCount_RNA > thresholds$nCount_RNA_lower &
                              nFeature_RNA > thresholds$nFeature_RNA_lower &
                              nFeature_RNA < thresholds$nFeature_RNA_upper & 
                              percent.mt < thresholds$percent_mitochondrial
                            )
    
    # Get number of cells after filtering
    remaining_cells <- ncol(filt_seurat)
    
    # Get number of features after filtering
    remaining_features <- nrow(filt_seurat)
    
    print(paste0(remaining_cells, " cells remaining from ", total_cells))
    print(paste0(remaining_features, " features remaining from ", total_features))
    
    
    before_plot <- VlnPlot(mt_seurat, 
                           features = c("nFeature_RNA",
                                        "nCount_RNA",
                                        "percent.mt"),
                           ncol = 3
                           )
    after_plot <- VlnPlot(filt_seurat, 
                           features = c("nFeature_RNA",
                                        "nCount_RNA",
                                        "percent.mt"),
                           ncol = 3
                          )
  
    output_list <- list(plot_b = before_plot, plot_a = after_plot, seurat = filt_seurat)
    
  # Return annotated seurat list
  return(output_list)
}

AnalyzeDims <- function(seurat_list, dim_red_params) {
  # Desc: Performs scaling, PCA, and dim reduction, and clustering on seurat objects
  # Args: seurat_list: a list of seurat objects to modify
  #       dim_red_params: dimensional reduction parameters from params
  # Returns: a list of modified seurat objects
  
  # Create empty seurat list to return:
  seurat_objects <- list()
  
  # Iterate through each seurat object
  for(i in seq_along(seurat_list)) {
    
    # Get Seurat object
    seurat_object <- seurat_list[[i]]
    
    # Get Sample name
    sample_name <- names(seurat_list)[i]
    
    # Perform Scaling
    # To quote Daniel Beiting of DIYTranscriptomics :
    # it is standard practice to apply a linear transformation ('scaling') before PCA. For single cell data this includes:
    # 1. Shifting the expression of each gene, so that the mean expression across cells is 0
    # 2. Scaling the expression of each gene, so that the variance across cells is 1
    # This gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
    seurat_object <- ScaleData(seurat_object, verbose = FALSE)
    
    # Perform PCA
    seurat_object <- RunPCA(seurat_object, 
                            npcs = dim_red_params$PCA_component_count, 
                            verbose = FALSE)
    
    # Perform UMAP if requested
    if(dim_red_params$run_UMAP){
      seurat_object <- RunUMAP(seurat_object,
                               reduction = "pca", 
                               dims = 1:dim_red_params$PCA_component_count)
      
    }
    
    # Perform tSNE if requested
    if(dim_red_params$run_tSNE){
      seurat_object <- RunTSNE(seurat_object,
                               reduction = "pca",
                               dims = 1:dim_red_params$PCA_component_count)
      
    }
    
    # Perform Nearest Neighbor and Clustering analysis
    seurat_object <- FindNeighbors(seurat_object,
                                   reduction = "pca",
                                   dims = 1:dim_red_params$PCA_component_count)
    seurat_object <- FindClusters(seurat_object, resolution = 0.5)
    
    # Add filtered seurat object to ouput list
    seurat_objects[[sample_name]] <- seurat_object
  }
  
  # Return annotated seurat list
  return(seurat_objects)
}

CreateWebApp <- function(seurat_list, parameters) {
  # Desc: Creates interactive ShinyCell webapp for data exploration
  # Args: seurat_list: A list of seurat objects
  #       parameters: Imported params
  # Returns: Void
  
  # Create directory for ShinyCell Files
  dir.create("./ShinyCellFiles")
  
  # Get ShinyCell specific params (for readability)
  sc_params <- parameters$shinycell_settings
  
  # Create empty vector of Sample prefixes for combining them into 1 app
  sc_prefixes <- c()
  
  # Iterate through each seurat object
  for(i in seq_along(seurat_list)) {
    
    # Get Seurat object
    seurat_object <- seurat_list[[i]]
    
    # Get Sample name
    sample_name <- names(seurat_list)[i]
    
    # Create ShinyCell Config
    shiny_cell_config <- createConfig(seurat_object)
    
    # INSERT CODE HERE FOR FUTURE SHINYCELL CONFIG MODIFICATION
    
    # Check if seurat object is integrated
    is_integrated <- "integrated" %in% Assays(seurat_object)
    
    # Set assay to RNA or Integrated as necessary
    if(is_integrated){
      shiny_assay <- "integrated"
    }else{
        shiny_assay <- "RNA"
        }
    
    makeShinyFiles(seurat_object, 
                   shiny_cell_config,
                   gex.assay = shiny_assay,
                   gex.slot = "data",
                   gene.mapping = TRUE,
                   shiny.prefix = sample_name,
                   shiny.dir = "./ShinyCellFiles/")
    
    sc_prefixes <- c(sc_prefixes, sample_name)
  }
  
  # Assemble Shiny Cell App
  makeShinyCodesMulti(
    shiny.title = "test",
    shiny.headers = sc_prefixes,
    shiny.prefix = sc_prefixes,
    shiny.dir = "./ShinyCellFiles/",
    shiny.footnotes = sc_params$footnote
  )
  
}

IntegrateSeuratData <- function(seurat_list, parameters) {
  # Desc: Integrates Seurat Data Sets according to params
  # Args: seurat_list: A list of all seurat objects
  #       parameters: imported parameters=
  # Returns: Modified list of seurat objects
  
  # Create empty list to return
  seurat_objects <- list()
  
  # Create empty list of integrated Seurats that need Dim analysis
  int_seurat_objects <- list()

  # Iterate through all samples and add them to output
  for(i in seq_along(seurat_list)) {
    
    # Re-add the each unintegrated seurat object
    seurat_objects[[names(seurat_list)[i]]] <- seurat_list[[i]]
    
  }
  
  # Create a dataframe of samples with integration groups
  samples_df <- parameters$study_design$samples %>%
    filter(!is.na(integration_group))
  
  # Get a list of integration groups
  unique_groups <- unique(samples_df$integration_group)
  
  # Iterate through each group
  for(group in unique_groups) {
    # Subset the dataframe
    subset_df <- samples_df[samples_df$integration_group == group]
    
    # Create empty standalone list
    object_list <- c()
    
    # Refactor later 
    # Create Standalone Seurat Objects
    for(object in subset_df$name) {
      assign(object, seurat_list[[object]])
      object_list <- c(object_list, get(object))
    }
    
    # Find integration features for samples
    integration_features <- SelectIntegrationFeatures(object_list
      )
    
    # Find integration anchors for samples
    integration_anchors <- FindIntegrationAnchors(
      object.list = object_list,
      anchor.features = integration_features
      )
    
    # Create integrated seurat set
    integrated_seurat <- IntegrateData(integration_anchors)
    
    # Add to output list
    int_seurat_objects[[group]] <- integrated_seurat
    
  }
  
  int_seurat_objects <- AnalyzeDims(int_seurat_objects, parameters$dimensional_reduction_settings)
  
  seurat_objects <- c(int_seurat_objects, seurat_objects)
  
  # Return updated seurat list
  return(seurat_objects)
  
}

GenerateSCEobjects <- function(cleaned_mtx_list) {
  # Desc: Takes a list of SoupX cleaned counts and returns 
  # Args: cleaned_mtx_list: A list of cleaned matrices from DecontaminateSoup()
  # Returns: a list of sce objects
  
  # Create empty list to return
  sce_objects <- list()
  
  # Create empty seurat list
  seurat_list <- list()
  
  # Create seurat objects from cleaned matrices
  for(i in seq_along(cleaned_mtx_list)) {
    
    # Get Sample Name
    mtx_name <- names(cleaned_mtx_list)[i]
    
    seurat_list[[mtx_name]] <- 
      CreateSeuratObject(cleaned_mtx_list[[mtx_name]], project = mtx_name)
    
  }
  
  # Merge Seurat and create QC plot
  merged_seurat <- seurat_list[[1]]  # Start with the first Seurat object
  if(length(seurat_list) > 1) {
    for(i in 2:length(seurat_list)) {
      merged_seurat <- merge(merged_seurat, seurat_list[[i]])
    }
  }
  plot <- FeatureScatter(merged_seurat, feature1="nCount_RNA", feature2="nFeature_RNA") + 
    ggtitle("Before Doublet Removal")
  
  # For each seurat object
  for(i in seq_along(seurat_list)) {
    
    # Get Seurat object
    seurat_object <- seurat_list[[i]]
    
    # Get Sample name
    sample_name <- names(seurat_list)[i]
    
    # Convert seurat to sce
    sce_obj <- as.SingleCellExperiment(seurat_object)
    
    # Add sce object to output list
    sce_objects[[sample_name]] <- sce_obj
  }
  
  # Return sce list
  return(list(sce_list=sce_objects, before_plot=plot))
}

IdentifyDoublets <- function(sce_list) {
  # Desc: Takes a list of sce objects and returns a list of annotated sce objects
  # Args: sce_list: A list of sce objects
  # Returns: a list of sce annotated objects
  
  # Create empty list to return
  sce_objects <- list()
  
  # Iterate through each sce object
  for(i in seq_along(sce_list)) {
    
    # Get sce object
    sce_object <- sce_list[[i]]
    
    # Get Sample name
    sample_name <- names(sce_list)[i]
    
    # Identify Doublets (uses random method)
    sce_object <- scDblFinder(sce_object)
    
    # Add sce object to output list
    sce_objects[[sample_name]] <- sce_object
  }
  
  # Return sce list
  return(sce_objects)
}

RemoveDoublets <- function(sce_list) {
  # Desc: Takes a list of sce objects, indentifies doublets, and removes them.
  # Args: sce_list: A list of sce objects
  # Returns: a list of filtered sce objects
  
  # Create empty list to return
  sce_objects <- list()
  
  # Identify Doublets
  sce_list <- IdentifyDoublets(sce_list)
  
  # Iterate through each sce object
  for(i in seq_along(sce_list)) {
    
    # Get sce object
    sce_object <- sce_list[[i]]
    
    # Get Sample name
    sample_name <- names(sce_list)[i]
    
    # Identify Doublets (uses random method)
    sce_object <- sce_object[,sce_object$scDblFinder.class == "singlet"]
    
    
    # Add sce object to output list
    sce_objects[[sample_name]] <- sce_object
  }
  
  # Return sce list
  return(sce_objects)
}

CreateMergedSeurat <- function(sce_list) {
  # Desc: Takes a list of sce objects and creates a merged seurat object from it.
  # Args: sce_list: A list of sce objects
  # Returns: one merged seurat object
  
  # Create temp list of seurat objects
  seurat_list <- list()

  # Iterate through single-cell experiment objects
  for(i in seq_along(sce_list)) {
    
    # Get sce object
    sce_object <- sce_list[[i]]
    
    # Get Sample name
    sample_name <- names(sce_list)[i]
    
    # Convert to seurat
    this_seurat <- as.Seurat(sce_object, data = NULL)
    
    # Add to list of seurat objects
    seurat_list[[sample_name]] <- this_seurat
  }

  # Merge all seurat objects into one (each object will have 2 layers)
  merged_seurat <- seurat_list[[1]]  # Start with the first Seurat object
  if(length(seurat_list) > 1) {

    for(i in 2:length(seurat_list)) {
      merged_seurat <- merge(merged_seurat, seurat_list[[i]])
    }
  }

  # Create after doublet removal QC plot
  plot <- FeatureScatter(merged_seurat, feature1="nCount_RNA", feature2="nFeature_RNA") + 
    ggtitle("After Doublet Removal")
  
  return(list(seurat=merged_seurat,after_plot=plot))
}