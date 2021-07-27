library(Seurat)

geo_data_folder <- "<Directory containing the GEO data>"
input_tables_folder <- "<Location of InputTables>"

# Data, metadata and custom colors reload
rawdata <- Read10X(data.dir = geo_data_folder, gene.column = 1)
rawdata[1:5,1:5]
metadata <-read.table(file = paste0(geo_data_folder,"/metadata.tsv.gz"), header = T, sep = "\t", row.names = 1)
head(metadata)
colors.table <- read.table(file = paste0(input_tables_folder,"/ClusterColors.tsv"), sep="\t", header = T, comment.char = "", as.is = T)
colors.use_celltypes   <- setNames(colors.table$color,colors.table$celltype)
colors.use_stages <- setNames(c('#5e5bf4','#f854ee','#f85554','#f8ea54'),c("Day5","Day6","Day7","Day8"))
colors.use_replicates <- c("#eed239","#7070fc")
colors.use_germlayers <- setNames(c("#ffd100","#474747","#9f9f9f","#6f8def","#f56969"),c("ExE","PGC","Epiblast","Ectoderm","Mesendoderm"))
colors.use_epi <- setNames(c("red","blue"),c("Control","Epithelialized"))
colors.use_proportions <- colorRampPalette(c("Red","#CB9898","Grey","#9797CB","Blue"),interpolate="spline")(100)

##### Create Seurat objects, including batch-corrected PCA coordinates (i.e. harmony) and UMAPs #####
SO.overview <- CreateSeuratObject(counts = rawdata, project = "Overview", meta.data = metadata)
SO.overview <- NormalizeData(object = SO.overview, scale.factor = 10000)
SO.overview[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(metadata[,paste0("overview_harmony",1:15)]),key = "harmony_",assay = "RNA")
SO.overview[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(metadata[,c("ecto_closeup_UMAP_1","ecto_closeup_UMAP_2")]),key = "UMAP_",assay = "RNA")

SO.ecto <- SO.overview[,!is.na(SO.overview$ecto_closeup_UMAP_1)]
SO.ecto <- NormalizeData(object = SO.ecto, scale.factor = 10000)
SO.ecto[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(SO.ecto@meta.data[,paste0("ecto_harmony",1:15)]),key = "harmony_",assay = "RNA")
SO.ecto[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(SO.ecto@meta.data[,c("ecto_closeup_UMAP_1","ecto_closeup_UMAP_2")]),key = "UMAP_",assay = "RNA")

SO.meso <- SO.overview[,!is.na(SO.overview$ecto_closeup_UMAP_1)]
SO.meso <- NormalizeData(object = SO.meso, scale.factor = 10000)
SO.meso[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(SO.meso@meta.data[,paste0("ecto_harmony",1:15)]),key = "harmony_",assay = "RNA")
SO.meso[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(SO.meso@meta.data[,c("ecto_closeup_UMAP_1","ecto_closeup_UMAP_2")]),key = "UMAP_",assay = "RNA")

SO.exe <- SO.overview[,!is.na(SO.overview$ecto_closeup_UMAP_1)]
SO.exe <- NormalizeData(object = SO.exe, scale.factor = 10000)
SO.exe[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(SO.exe@meta.data[,paste0("ecto_harmony",1:15)]),key = "harmony_",assay = "RNA")
SO.exe[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(SO.exe@meta.data[,c("ecto_closeup_UMAP_1","ecto_closeup_UMAP_2")]),key = "UMAP_",assay = "RNA")

SO.neurons <- SO.overview[,!is.na(SO.overview$ecto_closeup_UMAP_1)]
SO.neurons <- NormalizeData(object = SO.neurons, scale.factor = 10000)
SO.neurons[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(SO.neurons@meta.data[,paste0("ecto_harmony",1:15)]),key = "harmony_",assay = "RNA")
SO.neurons[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(SO.neurons@meta.data[,c("ecto_closeup_UMAP_1","ecto_closeup_UMAP_2")]),key = "UMAP_",assay = "RNA")

SO.overview@meta.data[,grep("harmony|UMAP",colnames(SO.overview@meta.data),value=T)] <- NULL
SO.ecto@meta.data[,grep("harmony|UMAP",colnames(SO.ecto@meta.data),value=T)] <- NULL
SO.meso@meta.data[,grep("harmony|UMAP",colnames(SO.meso@meta.data),value=T)] <- NULL
SO.exe@meta.data[,grep("harmony|UMAP",colnames(SO.exe@meta.data),value=T)] <- NULL
SO.neurons@meta.data[,grep("harmony|UMAP",colnames(SO.neurons@meta.data),value=T)] <- NULL

##### Choose which view: overview (germ layers), ectoderm + epiblasts, mesendoderm, extraembryonic #####
SO <- SO.overview
SO <- SO.ecto
SO <- SO.meso
SO <- SO.exe
SO <- SO.neurons

# Visualize timepoints:
Idents(SO) <- SO$Stage
DimPlot(SO, cols = colors.use_stages)

# Visualize germ layers:
Idents(SO) <- SO$GermLayer
DimPlot(SO, label = T, cols = colors.use_germlayers)+NoLegend()

# Visualize cell types:
Idents(SO) <- SO$celltype
DimPlot(SO, label = T, cols = colors.use_celltypes)+NoLegend()

# Visualize cell types as transferred by scmap from the Pijuan-Sala 2019 atlas:
Idents(SO) <- SO$pijuansala_scmap_celltype
DimPlot(SO, label = T)

# Visualize epithelialized vs control epiblasts (individual cells red vs blue):
Idents(SO) <- SO$Epithelialized
DimPlot(SO, label = T, cols = colors.use_epi)+NoLegend()

# Visualize epithelialized vs control proportions by cluster (continuous scale from 100% control to 100% epi):
FeaturePlot(SO, features = "proportion_epi", label = T, cols = colors.use_proportions[round(min(SO$proportion_epi)):round(max(SO$proportion_epi))], pt.size=1.1)+NoLegend()

# Plot gene expression:
gene_list <- c("Pou5f1","T","Sox1","Meis1")
FeaturePlot(SO, gene_list, sort=T)

# Plot gene expression split by condition:
gene_list <- c("Cdx2")
stage.use <- c("Day5","Day6","Day7","Day8")[1:4] # select a number or a range like 1:4 for all. 
FeaturePlot(SO, gene_list, sort=T, split.by = "Epithelialized",pt.size=1.1, cells = colnames(SO)[SO$Stage==stage.use])

# Plot gene expression in one condition only: 
gene_list <- c("Egr2")
condition <- c("Control","Epithelialized")[2] # Change between 1 and 2 for control or epi
FeaturePlot(SO, gene_list, sort=T,pt.size=1.1, cells = colnames(SO)[SO$Epithelialized==condition])

# Check for markers between two cell types (both conditions combined):
levels(SO$celltype)
Idents(SO) <- SO$celltype
markers <- FindMarkers(SO, ident.1 = "Spinal Cord", ident.2 = "Midbrain/Hindbrain")
max_pval <- 1e-2
min_logFC <- log(1.25)
min_dropout_rate_difference <- 0.3
markers.sign <- markers[markers$p_val_adj < max_pval & abs(markers$avg_logFC) > min_logFC & abs(markers$pct.1-markers$pct.2) > min_dropout_rate_difference,]
markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
markers.sign

# Check for differentially expressed genes between epi and control (for a given cell type):
levels(SO$celltype)
celltype.use <- "Epiblast"
Idents(SO) <- SO$Epithelialized
markers <- FindMarkers(SO[,SO$celltype==celltype.use], ident.1 = "Epithelialized", ident.2 = "Control")
max_pval <- 1e-2
min_logFC <- log(1.25)
min_dropout_rate_difference <- 0.1
markers.sign <- markers[markers$p_val_adj < max_pval & abs(markers$avg_logFC) > min_logFC & abs(markers$pct.1-markers$pct.2) > min_dropout_rate_difference,]
markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
markers.sign

# Check for markers between one cell type and all the rest (both conditions combined): 
levels(SO$celltype)
Idents(SO) <- SO$celltype
markers <- FindMarkers(SO, ident.1 = "Midbrain/Hindbrain")
max_pval <- 1e-2
min_logFC <- log(1.25)
min_dropout_rate_difference <- 0.3
markers.sign <- markers[markers$p_val_adj < max_pval & abs(markers$avg_logFC) > min_logFC & abs(markers$pct.1-markers$pct.2) > min_dropout_rate_difference,]
markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
markers.sign
