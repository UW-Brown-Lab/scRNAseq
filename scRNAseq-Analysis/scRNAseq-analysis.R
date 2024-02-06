#Import all functions from file
source('./functions.R')

#Load all package dependencies
LoadPackages()

#Load parameter file
params <- jsonlite::fromJSON("r-analysis-params.json")

#Load post-alignment samples from json
ImportCountData(params$study_design$samples, 
                params$filter_options$filter_with_DropletUtils,
                params$filter_thresholds$DropletUtils_call_threshold)

# Create list of seurat objects from imported counts
seurat_objects <- CreateSeuratObjects(params$seurat_settings$min_cells,
                    params$seurat_settings$min_features)

# Filter Seurat Objects by thresholds in params
if(params$filter_options$filter_out_mitochondrial) {
  # Determine percentage of MT Features
  seurat_objects <- CountMitochondrialFeatures(seurat_objects)
  # Filter by chosen parameters including MT percent
  seurat_objects <- FilterSeuratObject(seurat_objects, params$filter_thresholds)
} else {
  # Override mitochondrial threshold to not remove any
  params$filter_thresholds$percent_mitochondrial = 100
  
  seurat_objects <- FilterSeuratObject(seurat_objects, params$filter_thresholds)
}

# Perform Dimensional Analysis
seurat_objects <- AnalyzeDims(seurat_objects, params$dimensional_reduction_settings)

seurat_objects <- IntegrateSeuratData(seurat_objects,params)

CreateWebApp(seurat_objects, params$shinycell_settings)

