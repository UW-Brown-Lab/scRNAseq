# ~~~ Functions File ~~~
# Desc: Contains all functions called in downstream analysis
# 
# Best Practices:
# * Carefully define the purpose of each function, any arguments it takes, and
# *   and anything it returns
# * Single Responsibility: Each function should perform one job and do it well
# * Manage Scope Carefully
# * Parametrize your functions to make them configurable and reusable

ImportCountData <- function(samples,dir_info) {
  # Desc: Iterates through each sample file path and reads 10X counts for
  #       the raw and filtered data, then combined them into a list
  # Args: samples: The array of samples specified in the params
  #       dir_info: Info about the directory setup from params
  # Returns: A list of items, each containing the raw droplets and filtered
  #          counts
  
  # Get output directory from the alignment step from params (based on STAR)
  output_path <- paste0("./",dir_info$alignment_output_dir, "/")
  
  # Create list of lists to output
  output_list <- list()
  
  for(i in 1:nrow(samples)){
    
    # Get Sample
    sample <- samples[i,]
    
    # Get Sample Name
    sample_name <- sample$name
    
    # Get Sample File Path
    sample_path <- paste0(output_path, sample$file_path)
    
    # Read raw counts
    droplets.raw <- Read10X(paste0(sample_path,"raw_feature_bc_matrix.h5"))
    
    # Read filtered counts
    if(sample$antigen_capture & sample$antibody_capture) {
      all.counts <- Read10X(paste0(sample_path,"sample_filtered_feature_bc.h5"))
      counts.filt <- all.counts[["Gene Expression"]]
      antibody.counts <- all.counts[["Antibody Capture"]]
      antigen.counts <- all.counts[["Antigen Capture"]]
    } else if(sample$antigen_capture) {
      all.counts <- Read10X(paste0(sample_path,"sample_filtered_feature_bc.h5"))
      counts.filt <- all.counts[["Gene Expression"]]
      antibody.counts <- NULL
      antigen.counts <- all.counts[["Antigen Capture"]]
    } else if(sample$antibody_capture) {
      all.counts <- Read10X(paste0(sample_path,"sample_filtered_feature_bc.h5"))
      counts.filt <- all.counts[["Gene Expression"]]
      antibody.counts <- all.counts[["Antibody Capture"]]
      antigen.counts <- NULL
    } else {
      counts.filt <- Read10X(paste0(sample_path,"sample_filtered_feature_bc.h5"))
      antibody.counts <- NULL
      antigen.counts <- NULL
    }
    
    # Handle multi output raw/filt discrepancy
    if(sample$chemistry = "FRP") {
      print("Data type is FRP ~~~~~ resolving matrix dimensions")
      exclude <- setdiff(rownames(droplets.raw),rownames(counts.filt))
      keep <- !rownames(droplets.raw) %in% exclude
      droplets.raw <- droplets.raw[keep,]
    }
    
    assign(sample_name,
           list(
             filtered = counts.filt,
             raw = droplets.raw,
             antigen = antigen.counts,
             antibody = antibody.counts
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
    antigen.counts <- sample$antigen
    antibody.counts <- sample$antibody
    
    if(is.null(counts.raw)) {
      stop("Exiting execution -- Raw data does not exist")
    }
    
    # Create seurat object from filtered counts for clustering
    filt_seurat <- CreateSeuratObject(counts.filt)
    
    # Perform Normalization, Scaling, and Variable Feature selection with SCTransform
    filt_seurat <- SCTransform(filt_seurat, verbose = F, return.only.var.genes = F)
    
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
    output_list[[sample_name]] <- list(
      counts = corrected_matrix,
      antigen = antigen.counts,
      antibody = antibody.counts
      
    )
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
                                                        pattern = "^MT-")
  
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

GenerateSCEobjects <- function(cleaned_mtx_list, config) {
  # Desc: Takes a list of SoupX cleaned counts and returns 
  # Args: cleaned_mtx_list: A list of cleaned matrices from DecontaminateSoup()
  #       config: the parameters passed
  # Returns: a list of sce objects
  
  # Create empty list to return
  sce_objects <- list()
  
  # Create empty seurat list
  seurat_list <- list()
  
  # Create seurat objects from cleaned matrices
  for(i in seq_along(cleaned_mtx_list)) {
    
    # Get Sample Name
    mtx_name <- names(cleaned_mtx_list)[i]
    
    # Get sample config info from parameters
    sample <- config$samples[[mtx_name]]
    
    # Create base suerat
    seu <- CreateSeuratObject(cleaned_mtx_list[[mtx_name]], project = mtx_name)
    
    if(sample$antigen_capture) {
      # Add antigen data if present
      seu[["Antigen Capture"]] <- CreateAssay5Object(counts=cleaned_mtx_list[[mtx_name]]$antigen)
    }
    
    if(sample$antibody_capture) {
      # DEMUX antibody capture if needed
      seurat_list[[mtx_name]] <- DemuxSample(seu, cleaned_mtx_list[[mtx_name]]$antibody)
    } else {
      seurat_list[[mtx_name]] <- seu
    }
    
  }
  
  if(config$merge_samples){
    # Merge Seurat and create QC plot
    seu <- seurat_list[[1]]  # Start with the first Seurat object
    if(length(seurat_list) > 1) {
      for(i in 2:length(seurat_list)) {
        seu <- merge(seu, seurat_list[[i]])
      }
    }
    plot <- FeatureScatter(seu, feature1="nCount_RNA", feature2="nFeature_RNA") + 
      ggtitle("Before Doublet Removal")
  } else{
      plot <- list()
  }
  
  
  
  # For each seurat object
  for(i in seq_along(seurat_list)) {
    
    # Get Seurat object
    seurat_object <- seurat_list[[i]]
    
    # Get Sample name
    sample_name <- names(seurat_list)[i]
    
    if(!config$merge_samples) {
      plot[[i]] <- FeatureScatter(seu, feature1="nCount_RNA", feature2="nFeature_RNA") + 
        ggtitle("Before Doublet Removal") + plot_annotation(sample_name)
    }
    
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

CreateMergedSeurat <- function(sce_list, config) {
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
  
  if(config$merge_samples){
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
    
    DefaultAssay(merged_seurat) <- "RNA"
    
    return(list(seurat=merged_seurat,after_plot=plot))
  }
  for(i in seq_along(seurat_list)) {
    
    # Get Seurat object
    seurat_object <- seurat_list[[i]]
    
    # Get Sample name
    sample_name <- names(seurat_list)[i]
    
      plot[[i]] <- FeatureScatter(seu, feature1="nCount_RNA", feature2="nFeature_RNA") + 
        ggtitle("After Doublet Removal") + plot_annotation(sample_name)
    
    DefaultAssay(seurat_object) <- "RNA"
    
    seurat_list[[i]] <- seurat_object
  }
  
  return(list(seurat=seurat_list,after_plot=plot))
}

ClusterAnalysis <- function(merged_seurat,parameters) {
  # Desc: Perform clustering and generate feature plots
  # Args: merged_seurat: Merged seurat object or list of seurat objects
  #       parameters: Clustering parameters
  # Returns: A list containing modified seurat, and outputted plots
  

  # make empty plot list
  output_plots <- list()
  
  if(is.list(merged_seurat)) {
    
    output_seu_list <- list()
    for( i in seq_along(merged_seurat)){
      
      # Get Seurat object
      seu <- merged_seurat[[i]]
      
      # Get Sample name
      sample_name <- names(merged_seurat)[i]
      
      if(parameters$regress_mitochondrial){
        seu <- SCTransform(seu, verbose=FALSE, vars.to.regress="percent.mt", return.only.var.ganes = F) %>% 
          RunPCA() %>% 
          FindNeighbors(dims = 1:30) %>%
          RunUMAP(dims = 1:30) %>%
          FindClusters()
      } else {
        seu <- SCTransform(seu, verbose=FALSE, return.only.var.ganes = F) %>% 
          RunPCA() %>% 
          FindNeighbors(dims = 1:30) %>%
          RunUMAP(dims = 1:30) %>%
          FindClusters()
      }
      
      output_plots[[paste0("UMAP-",sample_name)]] <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters") + ggtitle("UMAP Plot") + plot_annotation(sample_name)
      feat_plots <- MakeFeaturePlots(seu, parameters)
      output_plots <- c(output_plots, feat_plots)
      output_seu_list[[sample_name]] <- seu
    }
    
    output_list <- list(
      seu = output_seu_list,
      plots = output_plots
    )
    return(output_list)
  }
  
  
  # Regularize with SCTransform
  seu <- SCTransform(merged_seurat, verbose=FALSE, vars.to.regress="percent.mt", return.only.var.ganes = F) %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters()
  
  # Plot UMAP
  output_plots[["UMAP"]] <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters") + ggtitle("UMAP Plot")
  
  
  # Generate Feature Plots
  feat_plots <- MakeFeaturePlots(seu, parameters)
  
  # Combine outputted plots
  output_plots <- c(output_plots, feat_plots)
  
  # Combine output list
  output_list <- list(
    seu = seu,
    plots = output_plots
  )
  
  # Return combined output list
  return(output_list)
}

MakeFeaturePlots <- function(merged_seurat,parameters) {
  # Desc: Generate feature plots from merged/dimensionally reduced seurat
  # Args: merged_seurat: Merged seurat object (already transformed)
  #       parameters: Clustering parameters
  # Returns: A list of outputted plots
  
  # Get seurat object
  seu <- merged_seurat
  
  # Make output list
  output_plots <- list()
  
  # Get list of Feature Plots
  feat_groups <- parameters$featureplot_settings$feature_groups
  
  # Iterate through feature plots
  for(group in names(feat_groups)) {
    
    # Get Features to show
    feats <- feat_groups[[group]]  
    
    # Make this feature plot
    output_plots[[group]] <- FeaturePlot(seu, features = feats, order = TRUE) + plot_annotation(group)
    
    # Make Heatmap
    output_plots[[paste0(group,"_heatmap")]] <- DoHeatmap(seu, features = feats) + plot_annotation(group)
    
  }
  
  return(output_plots)
}

DemuxSample <- function(seu, hto.mtx) {
  # Add HTO count matrix as assay to seurat obj
  seu[["HTO"]] <- CreateAssay5Object(counts=hto.mtx)
  # Normalize HTO data with CLR method
  seu <- NormalizeData(seu, assay="HTO", normalization.method = "CLR")
  # Call each cell as +/- for each tag
  seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.99)
  # View table of calls
  table(seu$HTO_classification.global)
  # Assign classification as main ident for subsetting
  Idents(seu) <- "HTO_classification.global"
  # subset to singlets
  seu <- subset(seu, idents = "Singlet")
  
  DefaultAssay(seu) <- "RNA"
  
  return(seu)
}

GenerateLoupeFile <- function(seu, filename) {
  # Creates a 10x Loupe Browser file from a Seurat object
  library(loupeR)
  create_loupe_from_seurat(seu, output_name = filename)
}