library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)


# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}

data_path <- "C:/Users/georg/Documents/2025_26Proj/RawData/output_combined/all-sample/DGE_unfiltered"
fig_path <- "/figures"

mat_path <- "C:/Users/georg/Documents/2025_26Proj/RawData/output_combined/all-sample/DGE_unfiltered"
mat <- ReadParseBio(mat_path)

# Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

# Create object
pbmc <- CreateSeuratObject(mat, min.features = 100, min.cells = 100,
                           names.field = 0, meta.data = cell_meta)



# Setting our initial cell class to a single type, this will changer after clustering. 
pbmc@meta.data$orig.ident <- factor(rep("pbmc", nrow(pbmc@meta.data)))
Idents(pbmc) <- pbmc@meta.data$orig.ident

SaveObject(pbmc, "seurat_obj_before_QC")
pbmc <- ReadObject("seurat_obj_before_QC")