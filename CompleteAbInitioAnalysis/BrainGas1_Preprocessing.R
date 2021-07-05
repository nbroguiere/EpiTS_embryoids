library(rPython)
library(dplyr)
library(matrixStats)
library(Seurat)
library(BiocGenerics)
library(ggplot2)
library(edgeR)
library(ggrepel)
library(reticulate)
library(scmap)
library(scater)
library(harmony)
library(plotly)
library(uwot)
library(DoubletFinder)
library(readr)
library(harmony)

##### Data and output folders #####
dataset_folder <- "D:/Raw_data_LSCB/scRNAseqMehmet/RawCounts" # Set to folder containing cellranger outputs for all datasets as subfolders
output_folder <- "C:/Users/nicki/Documents/PostDoc_LSCB/20-05-30_Mehmet_Brain_Gastruloids/20-07-16_Analysis1_mergeAllDatasets" # Set to folder of the analysis, which initially only contains the InputTables subfolder, and then will contain also subfolders with exported graphs
classifier_folder <- "C:/Users/nicki/Documents/PostDoc_LSCB/20-05-10_TrainScMapOnAtlas" # Set to folder containing the scmap classifier trained on in vivo atlas data.
Days_to_analyze <- c(5,6,7,8) # choose among 5,6,7,8
Batches_to_analyze <- c(1,2) # choose among 1,2

##### Data loading and QC - including a preliminary automated annotation transfer to be able to remove doublets #####
setwd(dir = output_folder)
datasets.all <- read.table(file = "InputTables/DatasetsMetadata.tsv", sep = "\t",header = TRUE, stringsAsFactors = F)
datasets.all <- datasets.all[,1:(ncol(datasets.all)-2)]
datasets.all
datasets <- datasets.all[datasets.all$Day %in% Days_to_analyze & datasets.all$Replicate %in% Batches_to_analyze & !datasets.all$Spikes,]
rownames(datasets) <- 1:nrow(datasets)
datasets
raw.data <- list()

setwd(dir = dataset_folder)
for(i in 1:nrow(datasets)){
  raw.data[[i]] <- Read10X(data.dir = datasets$Filename[i])
}

##### Load the cell type classifier trained on the in vivo atlas from Pijuan-Sala Griffiths Gottgens et al. Nature 2019 10.1038/s41586-019-0933-9 #####
scmap_classifier <- readRDS(file = paste(classifier_folder,"scmap_classifier3_1000markers.rds",sep="/"))
ref_stages <- c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")

# Create Seurat objects
setwd(dir = output_folder)
SO.list <- list()
homotypic.proportion <- list()
doublet_pct <- list()
pK <- list()
BCmetric <- list()
nDoublets <- list()
nDoublets_nonhomo <- list()
for(i in 1:nrow(datasets)){
  SO.list[[i]] <- CreateSeuratObject(counts = raw.data[[i]],project = datasets$Filename[i],assay = "RNA",min.cells = 3,min.features = 100)
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$Filename[i], col.name = "Dataset")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$Day[i], col.name = "Day")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$Replicate[i], col.name = "Replicate")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$Stage[i], col.name = "Stage")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$Batch[i], col.name = "Batch")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = colSums(SO.list[[i]][grep(pattern = "^mt-", x = rownames(SO.list[[i]]), value = TRUE), ])/colSums(SO.list[[i]])*100, col.name = "percent.mito")
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = datasets$Epithelialized[i], col.name = "Epithelialized")
  
  #ggplot(SO.list[[i]]@meta.data) + geom_point(aes(x=nFeature_RNA, y=percent.mito, color=nCount_RNA))+ scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')
  SO.list[[i]] <- subset(SO.list[[i]], subset = nFeature_RNA > datasets$Min.Genes[i] & percent.mito > datasets$Min.Pt.Mito[i] & percent.mito < datasets$Max.Pt.Mito[i])
  SO.list[[i]] <- NormalizeData(object = SO.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  #SO.list[[i]] <- FindVariableFeatures(object = SO.list[[i]],assay = "RNA",selection.method = "vst",nfeatures = 2000,verbose = TRUE)
  
  # Mini dimensionality reduction on individual datasets just for the doublet finder.
  npcs_df <- 30 # number of principal components used for doublet finder
  SO.list[[i]] <- FindVariableFeatures(SO.list[[i]],nfeatures = 1000)
  SO.list[[i]] <- ScaleData(SO.list[[i]],features = VariableFeatures(SO.list[[i]]))
  SO.list[[i]] <- RunPCA(SO.list[[i]],npcs = npcs_df)
  SO.list[[i]] <- RunUMAP(SO.list[[i]],dims = 1:npcs_df)
  
  # Apply cell type classifier (needed already now, early on, to estimate homotypic doublet proportion)
  set.seed(1234567)
  sce <- as.SingleCellExperiment(x = SO.list[[i]])
  rowData(sce)$feature_symbol <- rownames(sce)
  #sce <- sce[!duplicated(rownames(sce)), ] # Takes 20 GB of RAM and some time (probably using non-sparse matrix intermediate...)
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  scmapCluster_results <- scmapCluster(projection = sce, index_list = scmap_classifier[ref_stages],threshold = 0)
  #plot(getSankey(colData(sce)$celltype, scmapCluster_results$scmap_cluster_labs[,'atlas'],plot_height = 1200))
  SO.list[[i]] <- AddMetaData(object = SO.list[[i]], metadata = scmapCluster_results$combined_labs, col.name = "celltype")
  Idents(SO.list[[i]]) <- SO.list[[i]]$celltype
  
  # # DoubletFinder: best pK identification (use for early exploration, then fix the value)
  # sweep.res.list <- paramSweep_v3(SO.list[[i]], PCs = 1:npcs_df, sct = FALSE)
  # sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  # bcmvn <- find.pK(sweep.stats)
  # BCmetric[[i]] <- bcmvn$BCmetric
  # ggplot(bcmvn,mapping = aes(x=pK,y=BCmetric)) + geom_point()
  # pK[[i]] <- as.numeric(as.character(bcmvn$pK))[bcmvn$BCmetric==max(bcmvn$BCmetric)] # pK corresponding to the max in BC metric
  
  # Estimate number of doublets
  homotypic.proportion[[i]] <- modelHomotypic(SO.list[[i]]$celltype)
  doublet_pct[[i]] <- ncol(SO.list[[i]])/1000*0.76 # The estimate of 10x genomics as per the RNA profiling manual are 0.76% doublet per thousand recovered cells.
  nDoublets[[i]] <- round(doublet_pct[[i]]/100*length(colnames(SO.list[[i]])))
  nDoublets_nonhomo[[i]] <- round(nDoublets[[i]]*(1-homotypic.proportion[[i]]))
}

# Check the BCmetrics, that can be used to automate the choice of the pK:
metric_avg <- BCmetric[[1]]*0
pK_ <- as.numeric(as.character(bcmvn$pK))
for(i in 1:nrow(datasets)){
  metric_avg <- metric_avg+BCmetric[[i]]/nrow(datasets)
  plot(pK_,BCmetric[[i]],type="l",col=c("red","green","blue","cyan","magenta","grey","yellow","black")[i])
  par(new=TRUE)
}
plot(pK_,metric_avg)

## Run DoubletFinder
# First try using pK = pK[[i]] but doesn't always work -> datasets 5,6,8 have no clear max and end up with a random very high pK removing a whole true single cell cluster instead of isolated doublets
# We can have rather small true clusters, so we need rather low pK. The average BCmetric curve shows three local peaks: 0.02 to 0.04, 0.15-0.17, and ~0.23-0.25. Keep the first peak, choose pK=0.03 for all. 
pK.use <- 0.03
for(i in 1:nrow(datasets)){
  SO.list[[i]]@meta.data[,grep("DF.classifications",colnames(SO.list[[i]]@meta.data),value=T)] <- NULL
  SO.list[[i]]@meta.data[,grep("pANN",colnames(SO.list[[i]]@meta.data),value=T)] <- NULL
  SO.list[[i]] <- doubletFinder_v3(SO.list[[i]], PCs = 1:npcs_df, pN = 0.25, pK = pK.use, nExp = nDoublets_nonhomo[[i]], reuse.pANN = FALSE, sct = FALSE)
  FeaturePlot(SO.list[[i]],pt.size=1.5,feature=grep("pANN",colnames(SO.list[[i]]@meta.data),value=T))
  DimPlot(SO.list[[i]],pt.size = 1.5,group.by = grep("DF.classifications",colnames(SO.list[[i]]@meta.data),value=T),cols=c("red","lightgrey"))
}
rm(sce)
gc()

##### Merge all the datasets #####
names(SO.list) <- datasets$Filename
SO <- merge(x = SO.list[[1]],y = SO.list[2:length(SO.list)],add.cell.ids = datasets$Filename, project = "Brain_Gastruloids", merge.data = T)
#rm(raw.data) # consider releasing the memory if limited RAM is available.
#rm(SO.list)

# Create metadata columns with unified doublet annotations accross datasets
tmp <- character(length=ncol(SO))
names(tmp) <- colnames(SO)
for(i in grep("DF.classifications",colnames(SO@meta.data),value=T)){
  j <- !is.na(SO@meta.data[,i])
  tmp[j] <- SO@meta.data[j,i]
}
SO <- AddMetaData(object = SO,metadata = tmp,col.name = "doublets")
tmp <- numeric(length=ncol(SO))
names(tmp) <- colnames(SO)
for(i in grep("pANN",colnames(SO@meta.data),value=T)){
  j <- !is.na(SO@meta.data[,i])
  tmp[j] <- SO@meta.data[j,i]
}
SO <- AddMetaData(object = SO,metadata = tmp,col.name = "doublet_score")

# Remove doublets
Idents(SO) <- SO$doublets
SO <- subset(x = SO, idents = "Doublet", invert = T)

# Remove spurrious columns from SO: 
for(i in grep("DF.classifications",colnames(SO@meta.data),value=T)){
  SO@meta.data[,i] <- NULL
}
for(i in grep("pANN",colnames(SO@meta.data),value=T)){
  SO@meta.data[,i] <- NULL
}
SO$doublets <- NULL

##### Analysis - All cells #####
# Dimensionality reduction and batch correction
n.var.feats <- 2000
SO <- FindVariableFeatures(SO,nfeatures = n.var.feats)
LabelPoints(plot = VariableFeaturePlot(SO), points = VariableFeatures(SO), repel = FALSE, xnudge = 0, ynudge = 0,size=4)
SO <- ScaleData(SO,do.center = T, do.scale=T)
SO <- RunPCA(SO)
ElbowPlot(SO,ndims = 50)
Idents(SO) <- SO$celltype
DimPlot(SO,reduction = "pca",dims = c(15,16))+NoLegend()
dims.use <- 1:15
SO <- RunHarmony(object = SO, group.by.vars = "Replicate", reduction = "pca", dims.use = dims.use, nclust = 400)

# UMAP initialized on the default (diffusion maps) -> not behaving very well, some important trajectories are crossing each other. 
# SO <- RunUMAP(SO,reduction = "harmony",dims = dims.use,n.neighbors = 300,min.dist = 1.75,spread = 1, n.epochs = 1000)
# SO[["umap_dm"]] <- SO[["umap"]] # Save a backup of this umap in case of computing another refined one
# Idents(SO) <- SO$Day
# Idents(SO) <- SO$Replicate
# Idents(SO) <- SO$Epithelialized
# Idents(SO) <- SO$orig.ident
# DimPlot(SO,pt.size = 1.5,cells=sample(colnames(SO)))
# Idents(SO) <- SO$celltype
# DimPlot(SO,pt.size = 1.5,label=T)+NoLegend()

# Clustering
res <- 3
seed2 <- 57301
SO <- FindNeighbors(object = SO, reduction = "harmony",dims = dims.use)
SO <- FindClusters(object = SO, resolution = res,random.seed = seed2)
SO$louvain_clusters <- SO$seurat_clusters
SO$RNA_snn_res.4.6 <- NULL
SO$seurat_clusters <- NULL
#DimPlot(SO,pt.size = 1.5,label=T)+NoLegend()

# Markers per cluster for the overview Louvain clustering (set to TRUE to recompute)
if(!dir.exists("OutputTables")){dir.create("OutputTables")}
if(FALSE){
  markers <- FindAllMarkers(SO)
  markers.sign <- markers[markers$p_val_adj<1e-2 & markers$avg_logFC > log(1.5),]
  markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
  markers.sign <- markers.sign[order(markers.sign$cluster),]
  write.table(x=markers.sign,file = "./OutputTables/OverviewClusterMarkers_pP2_FC1.5.tsv",sep = "\t",row.names=F,col.names = T)
}

## Low dimensional embedding
# Define a layout based on a tractable small subset of 100 cells per cluster:
if(!dir.exists("Layouts")){dir.create("Layouts")}
#for(seed in c(72,157,324,538,749,953,1254,5278,10285,10285)){ # Test for various seeds initially, try to find one which doesn't result in important trajectories crossing each other, then fix the seed.
for(seed in c(953)){
  spread <- 2
  min_dist <- 1
  nn=150
  set.seed(seed)
  Idents(SO) <- SO$louvain_clusters
  cells.use <- WhichCells(object = SO,downsample = 100)
  local_connectivity=1 # Tried 2 and was not so convincing. Should not be more than the local intrinsic dimension of the manifold.
  fast_sgd <- F # Should set it to false ultimately, to get reproducible results, but faster on true for early exploration.
  umap_init <- "spectral" # "normlaplacian", "spectral" (with noise), "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates.
  reduction.use <- "harmony"
  set.seed(seed)
  tmp <- umap(X = Embeddings(SO[[reduction.use]])[cells.use,dims.use],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,init_sdev = 1e-4,local_connectivity=local_connectivity,ret_model=T,verbose = T)
  tmp2 <- 0*Embeddings(SO[["harmony"]])[,1:2]
  tmp2[cells.use,] <- tmp$embedding
  SO[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
  Idents(SO) <- SO$louvain_clusters
  png(filename = paste("./Layouts/Spread",spread,"_mind",min_dist,"_seed",seed,"Louvain_.png",sep=""),width = 1000,height = 1000)
  print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,cells = sample(cells.use)) + NoLegend())
  dev.off()
  Idents(SO) <- SO$celltype
  png(filename = paste("./Layouts/Spread",spread,"_mind",min_dist,"_seed",seed,"_Scmap.png",sep=""),width = 1000,height = 1000)
  print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,cells = sample(cells.use)) + NoLegend())
  dev.off()
}

# Pick up the average position of each cluster, as a reference layout:
Idents(SO) <- SO$louvain_clusters
x_ref <- data.frame(row.names = unique(SO$louvain_clusters))
for(i in unique(SO$louvain_clusters)){
  x_ref[i,1] <- mean(SO[["umap"]]@cell.embeddings[intersect(cells.use,names(SO$louvain_clusters[SO$louvain_clusters==i])),1])
  x_ref[i,2] <- mean(SO[["umap"]]@cell.embeddings[intersect(cells.use,names(SO$louvain_clusters[SO$louvain_clusters==i])),2])
}

# Generate a random position around the layout defined positions for all cells according to their cell types:
noise <- 5
set.seed(seed)
x_ini <- data.frame(row.names = colnames(SO))
x_ini[,1] <- runif(length(colnames(SO)))*noise
x_ini[,2] <- runif(length(colnames(SO)))*noise
Idents(SO) <- SO$louvain_clusters
for(i in unique(SO$louvain_clusters)){
  print(i)
  if(i %in% rownames(x_ref)){
    x_ini[WhichCells(SO,idents=i),1] <- x_ini[WhichCells(SO,idents=i),1]+x_ref[i,1]
    x_ini[WhichCells(SO,idents=i),2] <- x_ini[WhichCells(SO,idents=i),2]+x_ref[i,2]
  }
}

# Do an integrated umap initialized on these layout+noise positions:
if(!dir.exists("Idents_UMAP")){dir.create("Idents_UMAP")}
for(min_dist in c(3)){ # Use these loops to look for a good umap view # min.dist = 1.75,spread = 1
  for(spread in c(1.3)){
    for(nn in c(100)){
      cells.use <- colnames(SO) #WhichCells(object = tmp,downsample = 300) #colnames(SO) #sample(colnames(SO),1000)
      local_connectivity=1 # Should not be more than the local intrinsic dimension of the manifold. I would have imagined 2-3 could be reasonable, but doesn't give good results. 
      fast_sgd <- T # Should set it to false ultimately, to get exactly reproducible results, but can use T to get faster for early exploration. 
      umap_init <- as.matrix(x_ini[cells.use,]) # "normlaplacian", "spectral" (with noise),  "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates. 
      set.seed(seed)
      reduction.use <- "harmony"
      tmp <- umap(X = Embeddings(SO[[reduction.use]])[cells.use,dims.use],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,init_sdev = 1e-4,local_connectivity=local_connectivity,ret_model=T,verbose = T,n_epochs = 1000)
      tmp2 <- 0*Embeddings(SO[["harmony"]])[,1:2]
      tmp2[cells.use,] <- tmp$embedding
      SO[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
      Idents(SO) <- SO$celltype
      DimPlot(SO,pt.size=1,reduction = "umap",label=T,label.size = 8,repel = T) + NoLegend()
      
      # Plot identities. 
      for(ident_plot in c("celltype","louvain_clusters","Stage","Dataset","Replicate")){
        png(filename = paste("./Idents_UMAP/Spread",spread,"_mind",min_dist,"_nn",nn,"_",ident_plot,".png",sep=""),width = 1000,height = 1000)
        #png(filename = paste(ident_plot,"_",model_plot,".png",sep=""),width = 1000,height = 1000)
        cells.plot <- cells.use
        Idents(SO) <- factor(x = SO@meta.data[,ident_plot],levels = sort(unique(SO@meta.data[,ident_plot])),ordered = T)
        if(ident_plot=="celltype"){
          print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,repel=T,cells = sample(cells.plot),label.size = 8) + NoLegend()) #, cols = colors.use_transferred[levels(Idents(SO.align))]
        }else if(ident_plot=="Stage"){
          print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=F,repel=T,cells = sample(cells.plot))+NoLegend()) # , cols = colors.use_stages[levels(Idents(SO.align))]
        }else if(ident_plot=="louvain_clusters"){
          print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=T,repel=T,cells = sample(cells.plot))+NoLegend()) # , cols = colors.use_stages[levels(Idents(SO.align))]
        }else{
          print(DimPlot(SO,pt.size=1.5,reduction = "umap",label=F,cells = sample(cells.plot)))
        }
        dev.off()
      }
    }
  }
}

# # Not used, but for a simple check at a UMAP without batch correction: 
# SO <- RunUMAP(SO,reduction = "pca",reduction.name = "umappc",dims = dims.use,n.neighbors = 300,min.dist = 2,spread = 1, n.epochs = 1000)

# Interactive marker plots within R: 
Idents(SO) <- SO$celltype
FeaturePlot(SO,c("Nanog","Pou5f1","T"),sort.cell = T,pt.size=1.5,label = F,ncol = 3)
FeaturePlot(SO,c("Nkx6-1","Emx2"),sort.cell = T,pt.size=1,label = T,repel=T,split.by = "Epithelialized")

# Export UMAP plots of gene expression directly to the harddrive:
gene_list1 <- c("Pitx1","Bmp4","Tubb3","Hand1","Hand2","Sox2","Eya1","Eya2","Six1","Prrx1","Sox9","Tcf7l1","Tcf7l2","Lhx1","Lhx2","Sox1","Wnt7b","Gbx2","Crabp1","Crabp2","Bex1","Bex4","Rrm2","Hoxb4","Foxg1","Sox10","Six3","Hesx1","En1","Nkx6-1","Nkx6-2","Pax6","Hes3","Fgf15","Pax7","Emx2","Ascl1","Dbx2","Shh","Foxa2","Foxa1","Nes","Otx2","Olig2","Wnt7a","Pax3")
gene_list2 <- c("Six3","Otx2","Cdx2","Meox1","Gata6","Evx1","Fgf8","Cdh2","Cer1","Snai1","Snai2","Tcf15")
gene_list3 <- c("Hesx1","Dmbx1","Tnnt2","Ryr2","Hoxaas3","Bex4","Bex1","Irx2","Meset","Id3","Ube2c","Akr1b3","Alyref","Irx2")
gene_list4 <- c("Pou5f1","Nes","Tubb3","Shh","Olig2","En1","Six3","Sox10","Otx2","Sox17","Spink1","Pax3","Pax6","Pax7","Irx3","Hoxb1","En1","Math1","Dbx2","Ngn1","Ngn2","Mash1","Olig3","Ascl1","Atoh1","Nkx6.1","Dbx1","Nkx2.2","Foxa2","Tcf7l2","Dmbx1","Hesx1","Gbx2")
setwd(output_folder)
gene_list <- unique(c(gene_list1,gene_list2,gene_list3,gene_list4))
if(!dir.exists("GenePlots_All")){dir.create("GenePlots_All")}
for(gene in intersect(gene_list,rownames(SO))){
  for(condition in c("Both","Control","Epithelialized")){
    if(condition=="Both"){
      cells.use <- colnames(SO)
    }else{
      cells.use <- colnames(SO)[SO$Epithelialized==condition]
    }
    png(filename = paste("GenePlots_All/",gene,"_",condition,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO,gene,sort.cell = T,pt.size=1.5,label = T,repel=T,cells = cells.use))
    dev.off()
  }
}
dev.off()

##### Separate germ layers and ExE, PGC #####
tmp <- read.table(file = "InputTables/ClustersNamingGermLayers.tsv",sep="\t",header = T)
Idents(SO) <- SO$louvain_clusters
levels(Idents(SO)) <- tmp$GermLayer
SO$GermLayer <- as.character(Idents(SO))
DimPlot(SO,pt.size=1,label=T,label.size = 4,repel = T) + NoLegend()

##### Write down the preprocessed data in DataExport, use this as a starting point for further analysis ##### 
# Note: need to cd to the DataExport folder and execute "gzip *" in the command line to compress all files, recommended. 
setwd(output_folder)
tmp <- GetAssayData(object = SO,slot = "counts",assay = "RNA")
if(!dir.exists("DataExport")){dir.create("DataExport")}
write(colnames(tmp), file = "./DataExport/barcodes.tsv")
write(rownames(tmp), file = "./DataExport/features.tsv")
writeMM(obj = tmp, file = "./DataExport/matrix.mtx")
SO@meta.data[,grep("RNA_snn_res",colnames(SO@meta.data))] <- NULL
tmp <- cbind(SO@meta.data, Embeddings(SO[["umap"]]), Embeddings(SO[["harmony"]])[,dims.use])
tmp %>% head
write.table(x = tmp, file = "./DataExport/metadata.tsv", sep = "\t", row.names = T, col.names = NA)

##### (Facultative) 3D UMAP interactive visualization #####
if(FALSE){
  # Compute the projection
  SO <- RunUMAP(object = SO, reduction = "harmony", reduction.name = "umap3d",reduction.key = "UMAP3D_",dims = dims.use, umap.method = "uwot", min.dist = 1, n.neighbors = 300, spread = 10, seed.use = 750,n.components = 3L)
  
  # Plot discrete properties (e.g.Dataset, celltype, Day, batch, Replicate... )
  visualize <- "celltype"
  plot.data <- FetchData(object = SO, vars = c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3","celltype",visualize))
  plot.data <- plot.data[sample(1:nrow(plot.data),size = nrow(plot.data)),]
  colnames(plot.data) <- c("UMAP3D_1","UMAP3D_2","UMAP3D_3","celltype","Category")
  plot_ly(data = plot.data, x = ~UMAP3D_1, y = ~UMAP3D_2, z = ~UMAP3D_3, color = ~Category, colors = rainbow(length(unique(plot.data$Category))),type = "scatter3d", mode = "markers", marker = list(size = 2), text=~celltype,hoverinfo="text")
  
  # Plot continuous properties (e.g. gene expression value, nFeature_RNA etc)
  visualize <- "Emx2" 
  plot.data <- FetchData(object = SO, vars = c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3", "celltype",visualize))
  plot.data <- plot.data[sample(1:nrow(plot.data),size = nrow(plot.data)),]
  colnames(plot.data) <- c("UMAP3D_1","UMAP3D_2","UMAP3D_3","celltype","value")
  plot_ly(data = plot.data, x = ~UMAP3D_1, y = ~UMAP3D_2, z = ~UMAP3D_3, color = ~value, colors = c("lightgrey","blue"),type = "scatter3d", mode = "markers", marker = list(size = 2), text=~celltype,hoverinfo="text")
}