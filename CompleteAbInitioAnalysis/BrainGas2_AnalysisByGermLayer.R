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

output_folder <- "C:/Users/nicki/Documents/PostDoc_LSCB/20-05-30_Mehmet_Brain_Gastruloids/20-07-16_Analysis1_mergeAllDatasets" # Same as in the preprocessing step. 
data_folder <- "DataExport"
data_export_folder <- "DataExport_closeups"
setwd(output_folder)

rawdata <- Read10X(data.dir = data_folder, gene.column = 1)
rawdata[1:10,1:4]
metadata <- read.table(file = paste(data_folder,"metadata.tsv.gz",sep="/"), header = T, sep = "\t", row.names = 1)
head(metadata) #celltype is cluster number

SO <- CreateSeuratObject(counts = rawdata, project = "Cardiac Gastruloids", meta.data = metadata)
SO <- NormalizeData(object = SO, scale.factor = 10000)
SO[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(metadata[,grep("harmony",colnames(metadata))]),key = "harmony_",assay = "RNA")
SO[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(metadata[,c("UMAP_1","UMAP_2")]),key = "UMAP_",assay = "RNA")
SO@meta.data[,grep("harmony",colnames(SO@meta.data))] <- NULL
dims.use <- 1:ncol(SO[["harmony"]]@cell.embeddings)
Idents(SO) <- SO$celltype
SO$pijuansala_scmap_celltype <- SO$celltype
SO$celltype <- as.character(SO$celltype)
DimPlot(SO, label = T) + NoLegend()

##### Gene plots on all cells to HD #####
# Custom list
gene_list1 <- c("Pitx1","Bmp4","Tubb3","Hand1","Hand2","Sox2","Eya1","Eya2","Six1","Prrx1","Sox9","Tcf7l1","Tcf7l2","Lhx1","Lhx2","Sox1","Wnt7b","Gbx2","Crabp1","Crabp2","Bex1","Bex4","Rrm2","Hoxb4","Foxg1","Sox10","Six3","Hesx1","En1","Nkx6-1","Nkx6-2","Pax6","Hes3","Fgf15","Pax7","Emx2","Ascl1","Dbx2","Shh","Foxa2","Foxa1","Nes","Otx2","Olig2","Wnt7a","Pax3")
gene_list2 <- c("Six3","Otx2","Cdx2","Meox1","Gata6","Evx1","Fgf8","Cdh2","Cer1","Snai1","Snai2","Tcf15")
gene_list3 <- c("Hesx1","Dmbx1","Tnnt2","Ryr2","Hoxaas3","Bex4","Bex1","Irx2","Meset","Id3","Ube2c","Akr1b3","Alyref","Irx2")
gene_list4 <- c("Cer1","Pou5f1","Nes","Tubb3","Shh","Olig2","En1","Six3","Sox10","Otx2","Sox17","Spink1","Pax3","Pax6","Pax7","Irx3","Hoxb1","En1","Math1","Dbx2","Ngn1","Ngn2","Mash1","Olig3","Ascl1","Atoh1","Nkx6.1","Dbx1","Nkx2.2","Foxa2","Tcf7l2","Dmbx1","Hesx1","Gbx2")
gene_list5 <- c("Egr2","Lamp5","Lhx1","Lhx5","Lbx1","Ptf1a","Pax2","Tlx1","Tlx3","Evx1","Evx2","Lmx1a")
setwd(output_folder)
# Most variable genes
gene_list_vg <- VariableFeatures(SO)
# All Hox, Sox, Tbx, Pax genes
gene_list_hox <- grep("^Hox",rownames(SO),value=T)
gene_list_sox <- grep("^Sox",rownames(SO),value=T)
gene_list_tbx <- grep("^Tbx",rownames(SO),value=T)
gene_list_pax <- grep("^Pax",rownames(SO),value=T)
# Combine the lists:
gene_list <- unique(c(gene_list1,gene_list2,gene_list3,gene_list4,gene_list5,gene_list_hox,gene_list_sox,gene_list_tbx,gene_list_pax))
# Make the plots (set to true if needed)
if(FALSE){
  dir.plots <- "GenePlotsAll"
  if(!dir.exists(dir.plots)){dir.create(dir.plots)}
  for(gene in intersect(gene_list,rownames(SO))){
    for(condition in c("Both","Control","Epithelialized")){
      if(condition=="Both"){
        cells.use <- colnames(SO)
      }else{
        cells.use <- colnames(SO)[SO$Epithelialized==condition]
      }
      png(filename = paste(dir.plots,"/",gene,"_",condition,".png",sep=""),width = 800,height = 800)
      print(FeaturePlot(SO,gene,sort.cell = T,pt.size=1.5,label = T,repel=T,cells = cells.use))
      dev.off()
    }
  }
  dev.off()
}

##### Ectoderm #####
SO.ecto <- SO[,SO$GermLayer %in% c("Ectoderm","Epiblast")]
#SO.ecto <- SO.ecto[,SO.ecto$Stage %in% c("Day6")]
SO.ecto <- FindVariableFeatures(SO.ecto)
SO.ecto <- ScaleData(SO.ecto,features = VariableFeatures(SO.ecto))
SO.ecto <- RunPCA(SO.ecto)
ElbowPlot(SO.ecto)
dims.use.ecto <- 1:15
set.seed(72)
SO.ecto <- RunHarmony(object = SO.ecto, group.by.vars = "Replicate", reduction = "pca", dims.use = dims.use.ecto, nclust = 200)
res <- 3
SO.ecto <- FindNeighbors(object = SO.ecto, reduction = "harmony",dims = dims.use.ecto)
SO.ecto <- FindClusters(object = SO.ecto, resolution = res,random.seed = 72)
SO.ecto$louvain_clusters2 <- SO.ecto$seurat_clusters
SO.ecto@meta.data[,grep("RNA_snn_res",colnames(SO.ecto@meta.data))] <- NULL
SO.ecto$seurat_clusters <- NULL

# Markers for the Louvain reclustering (set to TRUE to recompute)
if(!dir.exists("OutputTables")){dir.create("OutputTables")}
if(F){
  Idents(SO.ecto) <- SO.ecto$louvain_clusters2
  dropout.rate <- rowSums(GetAssayData(SO.ecto,slot = "data")==0)/length(colnames(SO.ecto))*100
  hist(dropout.rate,100)
  tmp <- rownames(SO.ecto)[dropout.rate>33]
  markers.ecto <- FindAllMarkers(SO.ecto[tmp,],only.pos = T,min.pct = 0,pseudocount.use = 1) # Low expressed genes get missed because of low avg_logFC. Reducing the pseudocount would help. 
  markers.sign.ecto <- markers.ecto[markers.ecto$p_val_adj<1e-2 & markers.ecto$avg_logFC > log(1.5),]
  markers.sign.ecto <- markers.sign.ecto[order(-markers.sign.ecto$avg_logFC),]
  markers.sign.ecto <- markers.sign.ecto[order(markers.sign.ecto$cluster),]
  write.table(x=markers.sign.ecto,file = "./OutputTables/EctoClusterMarkers_pP2_FC1.5.tsv",sep = "\t",row.names=F,col.names = T)
}

# Cluster naming in the ectoderm 
tmp <- read.table(file = "InputTables/ClustersNamingEcto.tsv",sep="\t",header = T)
Idents(SO.ecto) <- SO.ecto$louvain_clusters2
levels(Idents(SO.ecto)) <- tmp$celltype
SO.ecto$celltype <- as.character(Idents(SO.ecto))
DimPlot(SO.ecto,pt.size=1,label=T,label.size = 4,repel = F) + NoLegend()

# Find a layout, using 100 cells per cluster
#for(seed in c(72,157,324,538,749,953,1254,5278,10285)){ # Screen for seeds to get a layout without trajectory crossings, then fix the seed. 
#for(seed in c(103,225,432,623,853,1152,2503,7892,12502)){
for(seed in c(1254)){
  spread <- 2
  min_dist <- 1
  nn=150
  set.seed(seed)
  Idents(SO.ecto) <- SO.ecto$louvain_clusters2
  cells.use <- WhichCells(object = SO.ecto,downsample = 100)
  local_connectivity=1 # Tried 2 and was not so convincing. Should not be more than the local intrinsic dimension of the manifold.
  fast_sgd <- F # Should set it to false ultimately, to get reproducible results, but faster for early exploration.
  umap_init <- "spectral" # "normlaplacian", "spectral" (with noise), "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates.
  reduction.use <- "harmony"
  set.seed(seed)
  tmp <- umap(X = Embeddings(SO.ecto[[reduction.use]])[cells.use,dims.use.ecto],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,init_sdev = 1e-4,local_connectivity=local_connectivity,ret_model=T,verbose = T,n_epochs = 300)
  tmp2 <- 0*Embeddings(SO.ecto[["harmony"]])[,1:2]
  tmp2[cells.use,] <- tmp$embedding
  SO.ecto[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
  Idents(SO.ecto) <- SO.ecto$louvain_clusters2
  png(filename = paste("./Layouts/Ecto_Spread",spread,"_mind",min_dist,"_seed",seed,"Louvain_.png",sep=""),width = 1000,height = 1000)
  print(DimPlot(SO.ecto,pt.size=1.5,reduction = "umap",label=T,cells = sample(cells.use)) + NoLegend())
  dev.off()
  Idents(SO.ecto) <- SO.ecto$celltype
  png(filename = paste("./Layouts/Ecto_Spread",spread,"_mind",min_dist,"_seed",seed,"_Scmap.png",sep=""),width = 1000,height = 1000)
  print(DimPlot(SO.ecto,pt.size=1.5,reduction = "umap",label=T,cells = sample(cells.use)) + NoLegend())
  dev.off()
}

# Pick up the average position of each cluster, as a reference layout:
Idents(SO.ecto) <- SO.ecto$louvain_clusters2
x_ref <- data.frame(row.names = unique(SO.ecto$louvain_clusters2))
for(i in unique(SO.ecto$louvain_clusters2)){
  x_ref[i,1] <- mean(SO.ecto[["umap"]]@cell.embeddings[intersect(cells.use,names(SO.ecto$louvain_clusters2[SO.ecto$louvain_clusters2==i])),1])
  x_ref[i,2] <- mean(SO.ecto[["umap"]]@cell.embeddings[intersect(cells.use,names(SO.ecto$louvain_clusters2[SO.ecto$louvain_clusters2==i])),2])
}

# Generate a random position around the layout defined positions for all cells according to their cell types:
noise <- 5
set.seed(seed)
x_ini <- data.frame(row.names = colnames(SO.ecto))
x_ini[,1] <- runif(length(colnames(SO.ecto)))*noise
x_ini[,2] <- runif(length(colnames(SO.ecto)))*noise
Idents(SO.ecto) <- SO.ecto$louvain_clusters2
for(i in unique(SO.ecto$louvain_clusters2)){
  print(i)
  if(i %in% rownames(x_ref)){
    x_ini[WhichCells(SO.ecto,idents=i),1] <- x_ini[WhichCells(SO.ecto,idents=i),1]+x_ref[i,1]
    x_ini[WhichCells(SO.ecto,idents=i),2] <- x_ini[WhichCells(SO.ecto,idents=i),2]+x_ref[i,2]
  }
}

# Do an integrated umap initialized on these layout+noise positions:
if(!dir.exists("Idents_UMAP")){dir.create("Idents_UMAP")}
for(min_dist in c(0.9)){ # Use these loops to look for a good umap view # min.dist = 1.75,spread = 1
  for(spread in c(17)){
    for(nn in c(100)){
      cells.use <- colnames(SO.ecto) #WhichCells(object = tmp,downsample = 300) #colnames(SO.ecto) #sample(colnames(SO.ecto),1000)
      local_connectivity=1 # Should not be more than the local intrinsic dimension of the manifold. I would have imagined 2-3 could be reasonable, but doesn't give good results. 
      fast_sgd <- T # Should set it to false ultimately, to get exactly reproducible results, but can use T to get faster for early exploration. 
      umap_init <- as.matrix(x_ini[cells.use,]) # "normlaplacian", "spectral" (with noise),  "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates. 
      set.seed(seed)
      reduction.use <- "harmony"
      n_epochs <- 500
      tmp <- umap(X = Embeddings(SO.ecto[[reduction.use]])[cells.use,dims.use.ecto],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,init_sdev = 1e-4,local_connectivity=local_connectivity,ret_model=T,verbose = T,n_epochs = n_epochs)
      tmp2 <- 0*Embeddings(SO.ecto[["harmony"]])[,1:2]
      tmp2[cells.use,] <- tmp$embedding
      SO.ecto[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
      Idents(SO.ecto) <- SO.ecto$celltype
      DimPlot(SO.ecto,pt.size=1,reduction = "umap",label=T,label.size = 8,repel = T) + NoLegend()
      
      # Plot identities. For final exports: 16 inches height
      for(ident_plot in c("celltype","louvain_clusters2","Stage","Dataset","Replicate")){
        png(filename = paste("./Idents_UMAP/Ecto_Spread",spread,"_seed",seed,"_mind",min_dist,"_nn",nn,"_",ident_plot,".png",sep=""),width = 1000,height = 1000)
        #png(filename = paste(ident_plot,"_",model_plot,".png",sep=""),width = 1000,height = 1000)
        cells.plot <- cells.use
        Idents(SO.ecto) <- factor(x = SO.ecto@meta.data[,ident_plot],levels = sort(unique(SO.ecto@meta.data[,ident_plot])),ordered = T)
        if(ident_plot=="celltype"){
          print(DimPlot(SO.ecto,pt.size=1.5,reduction = "umap",label=T,repel=T,cells = sample(cells.plot),label.size = 8) + NoLegend()) #, cols = colors.use_transferred[levels(Idents(SO.ecto.align))]
        }else if(ident_plot=="Stage"){
          print(DimPlot(SO.ecto,pt.size=1.5,reduction = "umap",label=F,repel=T,cells = sample(cells.plot))+NoLegend()) # , cols = colors.use_stages[levels(Idents(SO.ecto.align))]
        }else if(ident_plot=="louvain_clusters2"){
          print(DimPlot(SO.ecto,pt.size=1.5,reduction = "umap",label=T,repel=T,cells = sample(cells.plot))+NoLegend()) # , cols = colors.use_stages[levels(Idents(SO.ecto.align))]
        }else{
          print(DimPlot(SO.ecto,pt.size=1.5,reduction = "umap",label=F,cells = sample(cells.plot)))
        }
        dev.off()
      }
    }
  }
}

# Export plots of gene expression on UMAP
# Custom list
gene_list1 <- c("Pitx1","Bmp4","Tubb3","Hand1","Hand2","Sox2","Eya1","Eya2","Six1","Prrx1","Sox9","Tcf7l1","Tcf7l2","Lhx1","Lhx2","Sox1","Wnt7b","Gbx2","Crabp1","Crabp2","Bex1","Bex4","Rrm2","Hoxb4","Foxg1","Sox10","Six3","Hesx1","En1","Nkx6-1","Nkx6-2","Pax6","Hes3","Fgf15","Pax7","Emx2","Ascl1","Dbx2","Shh","Foxa2","Foxa1","Nes","Otx2","Olig2","Wnt7a","Pax3")
gene_list2 <- c("Six3","Otx2","Cdx2","Meox1","Gata6","Evx1","Fgf8","Cdh2","Cer1","Snai1","Snai2","Tcf15")
gene_list3 <- c("Hesx1","Dmbx1","Tnnt2","Ryr2","Hoxaas3","Bex4","Bex1","Irx2","Meset","Id3","Ube2c","Akr1b3","Alyref","Irx2")
gene_list4 <- c("Cer1","Pou5f1","Nes","Tubb3","Shh","Olig2","En1","Six3","Sox10","Otx2","Sox17","Spink1","Pax3","Pax6","Pax7","Irx3","Hoxb1","En1","Math1","Dbx2","Ngn1","Ngn2","Mash1","Olig3","Ascl1","Atoh1","Nkx6.1","Dbx1","Nkx2.2","Foxa2","Tcf7l2","Dmbx1","Hesx1","Gbx2")
gene_list5 <- c("Egr2","Lamp5","Lhx1","Lhx5","Lbx1","Ptf1a","Pax2","Tlx1","Tlx3","Evx1","Evx2","Lmx1a")
setwd(output_folder)
# Top PCA contributors
gene_list_pca <- vector()
for(i in dims.use.ecto){
  gene_list_pca <- unique(c(gene_list_pca,TopFeatures(object = SO.ecto, dim = i, nfeatures = 50)))
}
# Most variable genes
gene_list_vg <- VariableFeatures(SO.ecto)
# All Hox, Sox, Tbx, Pax genes
gene_list_hox <- grep("^Hox",rownames(SO.ecto),value=T)
gene_list_sox <- grep("^Sox",rownames(SO.ecto),value=T)
gene_list_tbx <- grep("^Tbx",rownames(SO.ecto),value=T)
gene_list_pax <- grep("^Pax",rownames(SO.ecto),value=T)
# Combine the lists:
gene_list <- unique(c(gene_list1,gene_list2,gene_list3,gene_list4,gene_list5,gene_list_pca,gene_list_hox,gene_list_sox,gene_list_tbx,gene_list_pax))
# Make the plots
dir.plots <- "GenePlotsEcto"
Idents(SO.ecto) <- SO.ecto$celltype
if(!dir.exists(dir.plots)){dir.create(dir.plots)}
for(gene in intersect(gene_list,rownames(SO.ecto))){
  for(condition in c("Both")){ #,"Control","Epithelialized"
    if(condition=="Both"){
      cells.use <- colnames(SO.ecto)
    }else{
      cells.use <- colnames(SO.ecto)[SO.ecto$Epithelialized==condition]
    }
    png(filename = paste(dir.plots,"/",gene,"_",condition,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.ecto,gene,sort.cell = T,pt.size=1.5,label = T,repel=T,cells = cells.use))
    dev.off()
  }
}
dev.off()

# Interactive plots
if(FALSE){
#SO.ecto <- RunUMAP(SO.ecto,dims = dims.use.ecto, n.neighbors = 100, min.dist = 1, spread = 15)
Idents(SO.ecto) <- SO.ecto$celltype
Idents(SO.ecto) <- SO.ecto$louvain_clusters
Idents(SO.ecto) <- SO.ecto$Stage
Idents(SO.ecto) <- SO.ecto$Dataset
DimPlot(SO.ecto,pt.size=1,label=T,label.size = 4,repel = T) + NoLegend()
FeaturePlot(SO.ecto,c("Pou5f1","Nes","Tubb3","Shh","Olig2","En1","Six3","Sox10","Otx2"),ncol=3,sort=T)
FeaturePlot(SO.ecto,c("Sox17","Spink1"),sort=T)
FeaturePlot(SO.ecto,c("Six3","Otx2","En1"),split.by = "Epithelialized")
FeaturePlot(SO.ecto,c("Pax3","Pax6","Pax7","Irx3","Hoxb1","En1","Math1","Dbx2","Ngn1","Ngn2","Mash1","Olig3","Ascl1","Atoh1","Nkx6.1","Dbx1","Nkx2.2","Foxa2","Tcf7l2"),ncol=3,sort=T,pt.size=1.3)
FeaturePlot(SO.ecto,c("Dmbx1","Hesx1","Gbx2","Otx2"),split.by = "Epithelialized") # Otx2 Gbx2 -> mid/hindbrain border
FeaturePlot(SO.ecto,c("Lamp5","Lhx1","Lhx5","Lbx1","Ptf1a","Pax2"),pt.size=1.3,ncol=3,sort=T) # Gabaergic neurons
FeaturePlot(SO.ecto,c("Egr2","Hoxb1","En1","Irx3","Meis2","Gata6"),pt.size=1.3,ncol=3,sort=T) # Hindbrain/Midbrain
FeaturePlot(SO.ecto,c("Lmx1a","Wnt1","Wnt4","Msx1","Msx2","Msx3"),pt.size=1.3,ncol=3,sort=T) # Roof plate
FeaturePlot(SO.ecto,c("Foxd3","Snai1","Snai2"),pt.size=1.3,sort=T)
FeaturePlot(SO.ecto,c("Pax6","Pax7","Nhlh1","Stmn2","Ascl1","Hes5","Hey1"),ncol = 3,pt.size=1.3,sort=T) # Neural differentiation, see histo for pattern, but mostly starts at the flexure around E8.5 and extends throughout from then on. 
FeaturePlot(SO.ecto,c("Prdm8","Prdm12","Prdm13"),pt.size=1.3,sort=T)
# Cluster naming delineation: 
FeaturePlot(SO.ecto,c("Pou5f1","Fgf5","T","Nanog","Foxa2","Mixl1","Sox1","Fgf8","Hes3")) # Delineate epiblasts
FeaturePlot(SO.ecto,c("Otx2","Sox1","Sox3","Pax2","Nog","Fgf15","Pou3f1","Cdx2","Nkx1-2"),pt.size=1.3,sort=T) # Delineate neuroepithelium
FeaturePlot(SO.ecto,c("Krt5","Bmp2","Krt8","Krt18","Krt19","Tfap2c","Dsp","Msx2")) # Non-neural ectoderm / surface ectoderm (essentially absent)
FeaturePlot(SO.ecto,c("T","Sox2","Cdx4","Nkx1-2","Tbx6","Mesp1","Tbx6","Ncam1")) # Delineate NMP
FeaturePlot(SO.ecto,c("Cyp26a1")) # Caudal lateral epiblast
FeaturePlot(SO.ecto,c("Prom1","Tfap2a","Zic1","Sox1","Sox2","Hoxc9","Pax3","Msx1"),ncol=3) # Delineate neural plate - defined as the clusters during which floor plate and lateral neural plate markers appear. 
FeaturePlot(SO.ecto,c("En1","Gbx2","Fgf15","Crabp1","Zic1","Fgf8","Pax5","En1","Gbx2","Fgf8","Egr2","Wnt1","Foxg1"),ncol=4) # Midbrain/Hindbrain. Foxg1 is even a forebrain marker, and is somehow present. 
FeaturePlot(SO.ecto,c("Top2a","Mki67","Nes","Tubb3","Sox2","Sox9","Fabp7","Wnt7a","Wnt7b"),ncol=3) # Delineate the neurons and neural stem cells
}

# Check custom signatures for AP, DV, and differentiation axis
ecto_signatures <- readLines(con = "InputTables/EctodermCustomSignatures.tsv")
ecto_signatures <- strsplit(x = ecto_signatures,split = "\t")
names(ecto_signatures) <- unlist(lapply(ecto_signatures,"[",1))
ecto_signatures <- lapply(ecto_signatures,"[",-1)
ecto_signatures <- lapply(ecto_signatures,intersect,rownames(SO.ecto))
if(!dir.exists("Signatures_Ecto")){dir.create("Signatures_Ecto")}
Idents(SO.ecto) <- SO.ecto$louvain_clusters2
for(i in names(ecto_signatures)){
  if(!dir.exists(paste("Signatures_Ecto",i,sep="/"))){dir.create(paste("Signatures_Ecto",i,sep="/"))}
  for(j in ecto_signatures[[i]]){
    png(filename = paste("Signatures_Ecto/",i,"/",j,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.ecto,j,sort.cell = T,pt.size=1.5,label = T,repel=T))
    dev.off()
  }
}
for(i in names(ecto_signatures)){
  SO.ecto <- AddMetaData(SO.ecto,metadata = colSums(SO.ecto@assays$RNA@data[ecto_signatures[[i]],]),col.name = i)
  for(condition in c("Both")){ #,"Control","Epithelialized"
    if(condition=="Both"){
      cells.use <- colnames(SO.ecto)
    }else{
      cells.use <- colnames(SO.ecto)[SO.ecto$Epithelialized==condition]
    }
    png(filename = paste("Signatures_Ecto/",i,"_",condition,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.ecto,i,sort.cell = T,pt.size=1.5,label = T,repel=T,cells = cells.use))
    dev.off()
  }
}

if(FALSE){
FeaturePlot(SO.ecto,c("EpiblastAndNeuroepithelium","NeuralStemCellsEarly","NeuralStemCellsLate","Neurons"),ncol=2,pt.size=1.3,sort=T)
FeaturePlot(SO.ecto,c("FloorPlate","RoofPlate"),ncol=2,pt.size=1.3,sort=T)
FeaturePlot(SO.ecto,c("EarlyCaudal","SpinalCord_LatePosterior","MidBrain","Anterior"),ncol=2,pt.size=1.3,sort=T)
FeaturePlot(SO.ecto,c("Neurons1","Neurons2","Neurons3","Neurons4"),ncol=2,pt.size=1.3,sort=T)
FeaturePlot(SO.ecto,c("Excitatory","Inhibitory","Isl1","Lhx3"),ncol=2,pt.size=1.3,sort=T)
# Rework the MidSC signature!!
# Rework the MidBrain signature a bit. 
# Make an early rostral signature
}

# Neurons close-up
#SO.neurons <- SO.ecto[,SO.ecto$celltype %in% c("Neurons","Neural Stem Cells")] #  & SO.ecto$Day==8
SO.neurons <- SO.ecto[,SO.ecto$louvain_clusters2 %in% c(32,33,34,36,38,6,28)] # 6 28 and 7 contain progenitors, the rest are differentiated neurons. 6 might be already neurons. 
SO.neurons <- FindVariableFeatures(SO.neurons)
SO.neurons <- ScaleData(SO.neurons,features = VariableFeatures(SO.neurons)) # Latest projection might be with unscaled data, careful.
SO.neurons <- RunPCA(SO.neurons)
ElbowPlot(SO.neurons)
dims.use.neurons <- 1:30 # 1:15 # DimPlot(SO.neurons,reduction = "pca",dims=c(9,10)) 
set.seed(72)
SO.neurons <- RunHarmony(object = SO.neurons, group.by.vars = "Replicate", reduction = "pca", dims.use = dims.use.neurons, nclust = 200)
res.neurons <- 20 #1
SO.neurons <- FindNeighbors(object = SO.neurons, reduction = "harmony",dims = dims.use.neurons)
SO.neurons <- FindClusters(object = SO.neurons, resolution = res,random.seed = 72)
SO.neurons <- RunUMAP(SO.neurons,dims = dims.use.neurons,reduction = "harmony")
SO.neurons$louvain_clusters3 <- SO.neurons$seurat_clusters
SO.neurons$seurat_clusters <- NULL
Idents(SO.neurons) <- SO.neurons$louvain_clusters3
DimPlot(SO.neurons,label=T,pt.size=1.3)+NoLegend()
if(FALSE){
FeaturePlot(SO.neurons,c("Nes","Tubb3","Sox2"),pt.size = 2,sort=T)+NoLegend()
FeaturePlot(SO.neurons,c("Nes","Tubb3","Lhx1","Lhx4","Lhx5","Tlx1","Tlx3","Evx1","Evx2"),pt.size = 2,sort=T)+NoLegend()
FeaturePlot(SO.neurons,c("Lbx1","Isl1","Isl2","Pax2","Pax5","Nrxn3","Mcc","Nkx6-1","En1"),pt.size = 2,sort=T)+NoLegend()
FeaturePlot(SO.neurons,c("Lmx1a","Lmx1b","Foxa2","Six3","En1","Th","Nr4a2","Slc18a2","Ddc")) # Attempt to locate midbrain dopaminergic neurons
FeaturePlot(SO.neurons,c("Shh","Nkx2-2","Nkx6-1","Th","Htr2a","Pou4f1","Isl1","Fgf8"),pt.size=2,sort=T)+NoLegend() 
FeaturePlot(SO.neurons,c("Pax2","En1","Fgf8"),pt.size=2,sort=T)+NoLegend() 
}

# Signatures from Delile https://dev.biologists.org/content/146/12/dev173807
Delile <- list()
# Roof plate
Delile[["1a_RP"]] <- c("Lmx1a","Msx1","Msx2","Pax3","Wnt1","Olig3")
# Dorsal progenitors (Sox2+) dp1 to dp6
Delile[["2a_dp1-6"]] <- c("Msx1","Msx2","Pax3","Wnt1","Olig3","Irx3","Irx5","Pax6","Pax7","Gsx2","Ascl1","Gbx2","Gsx1","Dbx2","Dbx1","Sp8")
# Ventral progenitors p0 to p2 
Delile[["3a_vp0-2"]] <- c("Dbx2","Dbx1","Sp8","Nkx6-2","Prdm12","Nkx6-1","Foxn4")
# Motor neurons progenitors (Sox2+)
Delile[["4a_pMN"]] <- c("Sp8","Olig2","Nkx6-1")
# Progenitors 3 (next to floor plate) and floor plate
Delile[["5a_FP+p3"]] <- c(c("Nkx6-1","Nkx2-2","Nkx2-9"),c("Arx","Shh","Lmx1b","Nkx6-1","Foxa2"))
# Dorsal interneurons dI1 to dI6
Delile[["2b_dI1-6"]] <- c("Pou4f1","Lhx2","Lhx9","Barhl1","Barhl2","Atoh1","Foxd3","Lhx1","Lhx5","Tlx3","Prrx1","Otp","Pax8","Lbx1","Pax2","Gbx1","Bhlhe22","Pax3","Pax7","Ptf1a","Gsx1","Lmx1b","Ascl1","Dmrt3","Wt1")
# Ventral interneurons V1 to V2b
Delile[["3b_V0-2"]] <- c("Foxd3","Lhx1","Lhx5","Otp","Pax8","Pax2","Bhlhe22","Ascl1","Evx1","Evx2","Pitx2","En1", "Lhx3","Vsx2","Sox14","Sox21","Foxn4","Vsx1","Tal1","Gata2","Gata3","Msx1","Slc18a3")
# Motor neurons
Delile[["4b_MN"]] <- c("Isl1","Lhx3","Mnx1","Slc10a4","Slc18a3","Olig2","Aldh1a2","Arhgap36","Cldn3")
# Ventral interneurons V3
Delile[["5b_V3"]] <- c("Nkx2-2","Sim1")

# Interactive plots for exploration:
if(FALSE){
# To locate hindbrain and rhombomeres:
FeaturePlot(SO.neurons,"Phox2b",pt.size=2,order=T,split.by = "Epithelialized")
FeaturePlot(SO.neurons,"Msx3",pt.size=2,order=T)

# Rhombomeres 1-7 are are Hox1-3 positive and 7-11 are Hox4-7 all except 1 express hox genes -> We don't have rhombomere 1, but we have 2++
FeaturePlot(SO.neurons,c("Hoxb1","Hoxb2","Hoxb3","Hoxb4","Hoxb5","Hoxb6","Hoxb7"),pt.size=2,order=T) 

# Hindbrain Neurons Hernandez_miranda Developmental biology 2017
# Ventral to dorsal:
Idents(SO.neurons) <- SO.neurons$louvain_clusters3
FeaturePlot(SO.neurons,c("Lbx1","Pax2","Lhx1","Lhx5","Wt1","Bhlhe22","Dmrt3"),pt.size=2,order=T,label=T)  # Neurons dI6/dB4 more ventral in spine and hindbrain
FeaturePlot(SO.neurons,c("Lbx1","Tlx3","Lmx1b","Prrxl1","Pou4f1"),pt.size=2,order=T,label=T)  # dI5/dB3
FeaturePlot(SO.neurons,c("Phox2b","Lbx1","Atoh1"),pt.size=2,order=T,label=T) # Simultaneous expression only in hindbrain dB2 domain -> we don't seem to have
FeaturePlot(SO.neurons,c("Lbx1","Pax2","Lhx1","Lhx5"),pt.size=2,order=T,label=T) # Simultaneous expression only in hindbrain dI4/dB1
FeaturePlot(SO.neurons,c("Foxd3","Lhx1","Lhx5","Foxp2"),pt.size=2,order=T,label=T) # Simultaneous expression only in hindbrain dA4 domain
FeaturePlot(SO.neurons,c("Pou4f1","Tlx3","Prrxl1","Isl1","Phox2b","Lmx1b"),ncol=3,label=T,pt.size=2,order=T) # dI3/dA3
FeaturePlot(SO.neurons,c("Pou4f1","Lhx1","Lhx5","Foxd3","Foxp2"),pt.size=2,order=T,label=T) # Simultaneous expression only in hindbrain dI2/dA2 domain
FeaturePlot(SO.neurons,c("Pou4f1","Barhl1","Lhx2","Lhx9","Evx1"),pt.size=2,order=T,label=T) # Simultaneous expression only in hindbrain dI1/dA1 domain

FeaturePlot(SO.neurons,c("Phox2b","Foxa1","Shh","Wnt5a","Th","Ddc","Snca"),pt.size=2,order=T) # DA neuron signatures
FeaturePlot(SO.neurons,c("Nr4a2","Th","Pitx3","Lmx1a","Lmx1b","Wnt1","Msx1"),pt.size=2,order=T) # DA neuron signatures Arenas Villascusa 2015 Development
FeaturePlot(SO.neurons,c("Olig3","Ascl1","Neurog2","Ptf1a"),pt.size=2,order=T)  # Progenitor only in hindbrain pdA4 - Ptf1a is the hindbrain specific marker, present bottom right cluster and enriched in epithelialized organoids

FeaturePlot(SO.neurons,c("Nkx2-2","Phox2b","Ascl1","Tbx20","Nkx6-1"),pt.size=2,order=T)  # Branchial motor neurons Pla BMC neural development 2008
FeaturePlot(SO.neurons,c("Lhx2","Lhx9","Olig3","Atoh1"),pt.size=2,order=T)  # Most dorsal

FeaturePlot(SO.neurons,Delile[["MN"]][1:9],ncol=3) 
}

# Signature plots to hard drive
if(!dir.exists("Signatures_Delile")){dir.create("Signatures_Delile")}
Idents(SO.ecto) <- SO.neurons$seurat_clusters
for(i in names(Delile)){
  SO.neurons <- AddMetaData(SO.neurons,metadata = colSums(SO.neurons@assays$RNA@data[Delile[[i]],]),col.name = i)
  for(condition in c("Both")){ #,"Control","Epithelialized"
    if(condition=="Both"){
      cells.use <- colnames(SO.neurons)
    }else{
      cells.use <- colnames(SO.neurons)[SO.neurons$Epithelialized==condition]
    }
    png(filename = paste("Signatures_Delile/",i,"_",condition,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.neurons,i,sort.cell = T,pt.size=1.5,label = T,repel=T,cells = cells.use))
    dev.off()
  }
}
for(i in names(Delile)){
  if(!dir.exists(paste("Signatures_Delile",i,sep="/"))){dir.create(paste("Signatures_Delile",i,sep="/"))}
  for(j in Delile[[i]]){
    png(filename = paste("Signatures_Delile/",i,"/",j,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.neurons,j,sort.cell = T,pt.size=1.5,label = T,repel=T))
    dev.off()
  }
}

# Signed Delile Signatures Progenitors
DelileSignProg <- read.table(file = "InputTables/SignaturesDelileProgenitors.tsv")
colnames(DelileSignProg) <- gsub("\\.","-",colnames(DelileSignProg))
for(i in rownames(DelileSignProg)){
  SO.neurons <- AddMetaData(SO.neurons,metadata = colSums(sweep(as.matrix(SO.neurons@assays$RNA@data[colnames(DelileSignProg),]),1,as.numeric(DelileSignProg[i,]),'*')),col.name = i)
  for(condition in c("Both")){ #,"Control","Epithelialized"
    if(condition=="Both"){
      cells.use <- colnames(SO.neurons)
    }else{
      cells.use <- colnames(SO.neurons)[SO.neurons$Epithelialized==condition]
    }
    png(filename = paste("Signatures_Delile/SignedProgenitors_",i,"_",condition,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.neurons,i,sort.cell = T,pt.size=2,label = F,repel=T,cells = cells.use,min.cutoff = 'q50',max.cutoff = 'q99')+NoLegend())
    dev.off()
  }
}
# Signed Delile Signatures Neurons
DelileSign <- read.table(file = "InputTables/SignaturesDelileNeurons.tsv")
colnames(DelileSign) <- gsub("\\.","-",colnames(DelileSign))
for(i in rownames(DelileSign)){
  SO.neurons <- AddMetaData(SO.neurons,metadata = colSums(sweep(as.matrix(SO.neurons@assays$RNA@data[colnames(DelileSign),]),1,as.numeric(DelileSign[i,]),'*')),col.name = i)
  for(condition in c("Both")){ #,"Control","Epithelialized"
    if(condition=="Both"){
      cells.use <- colnames(SO.neurons)
    }else{
      cells.use <- colnames(SO.neurons)[SO.neurons$Epithelialized==condition]
    }
    png(filename = paste("Signatures_Delile/SignedNeurons_",i,"_",condition,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.neurons,i,sort.cell = T,pt.size=2,label = F,repel=T,cells = cells.use,min.cutoff = 'q50',max.cutoff = 'q99')+NoLegend())
    dev.off()
  }
}

#### Highlight clusters - see cluster 29 is comprised of two disconnected clusters (and continuous with cluster 23) that are actually V2a and V2b. Split it. 
# 23 and 15 are also actually the ventral progenitors giving rise to V2a/V2b rather than V2a/b themselves. 
Idents(SO.neurons) <- SO.neurons$louvain_clusters3
DimPlot(SO.neurons,pt.size=1.5,label=T)+NoLegend()
SO.neurons$highlight <- 0
SO.neurons$highlight[SO.neurons$louvain_clusters3 %in% c(20,22,39)] <- 1
FeaturePlot(SO.neurons,c("highlight","Lhx8"),pt.size=1.5,label=T)
FeaturePlot(SO.neurons,c("Nkx2-2","Nkx2-9"),pt.size=1.5,label=T,split.by = "Epithelialized")

# V2a and V2b ended up clustered together despite of very distinct umap positions and markers, split it manually:
cluster_29_subcluster <- c("0pc_D8_rep1_GAGCTGCCAGTCCGTG","0pc_D8_rep1_TTAATCCTCACACCCT","0pc_D8_rep2_CACGTTCCATCGGAAG","0pc_D8_rep2_CCGGTGACAGCGTATT","0pc_D8_rep2_CTACGGGGTCTGTCCT","0pc_D8_rep2_CTCAACCTCGAGGCAA","0pc_D8_rep2_GAAGAATTCTTCCTAA","0pc_D8_rep2_GATCATGAGTAAATGC","0pc_D8_rep2_GTCTACCTCGGCTCTT","0pc_D8_rep2_TGATGCAGTAGACAGC","0pc_D8_rep2_TTCGATTCACCTGATA","0pc_D8_rep2_TTGAGTGTCCTTCACG","3pc_D8_rep1_AGCCAATAGCTAGTTC","3pc_D8_rep1_TCAGGGCGTCGTACTA","3pc_D8_rep1_TCAGTGAGTTTACACG","3pc_D8_rep1_TGAGGGAGTGGTCCCA","3pc_D8_rep2_CACTGTCCACGTAACT","3pc_D8_rep2_GTAGCTAGTCGTTATG","3pc_D8_rep2_TAAGCGTCACTACTTT","3pc_D8_rep2_TCGGGCACAGGCGAAT")
SO.neurons$louvain_clusters3 <- as.character(SO.neurons$louvain_clusters3)
SO.neurons$louvain_clusters3[cluster_29_subcluster] <- 41
SO.neurons$louvain_clusters3 <- factor(SO.neurons$louvain_clusters3,levels = 0:41)

# The prdm12+ P1 cells are distinct in the UMAP trajectory from the V1 neurons, but got co-clustered, split them: 
prdm12_cluster <- c("0pc_D7_rep2_GGGACAAGTTCAAGTC","0pc_D8_rep1_AGAAGTAGTTTAGTCG","0pc_D8_rep1_AGTGATCAGGCTAGCA","0pc_D8_rep1_CAGATTGGTTCAAAGA","0pc_D8_rep1_CTCGAGGGTCCAGCCA","0pc_D8_rep1_GAACTGTCATTGTAGC","0pc_D8_rep1_GACCAATAGCGGCTCT","0pc_D8_rep1_GCACTAACAAATTGCC","0pc_D8_rep1_GGGACAAGTGGCTACC","0pc_D8_rep1_GTGGGAAAGACTCTAC","0pc_D8_rep1_TAGTGCAAGGGACCAT","0pc_D8_rep1_TGGGAGAGTGTGTGTT","0pc_D8_rep2_AAGAACAAGAGTCCGA","0pc_D8_rep2_AATGAAGTCAACTCTT","0pc_D8_rep2_AATGACCAGGACACTG","0pc_D8_rep2_ACACTGACACACGGTC","0pc_D8_rep2_AGCTTCCTCACAACCA","0pc_D8_rep2_CAACCTCCAGCCGGTT","0pc_D8_rep2_CATACCCCAAGAGAGA","0pc_D8_rep2_CGAAGTTTCGATACGT","0pc_D8_rep2_CTGTGAAAGAAGGGAT","0pc_D8_rep2_GAGTCTATCCTAGAGT","0pc_D8_rep2_GATAGCTGTACTTGTG","0pc_D8_rep2_TAAGCCATCCATATGG","0pc_D8_rep2_TCAAGTGCAGACCTGC","0pc_D8_rep2_TCGTAGACATAAGCGG","0pc_D8_rep2_TGCCGAGTCCAGTTCC","0pc_D8_rep2_TGCGACGCACGCGGTT","0pc_D8_rep2_TGTCAGATCCATCACC","0pc_D8_rep2_TTACAGGCAGACTGCC","0pc_D8_rep2_TTCTTCCCAGTCACGC","3pc_D7_rep2_TTGGGATCATCTAACG","3pc_D8_rep1_AAAGAACGTTTGAAAG","3pc_D8_rep1_AGGGTCCCAATGAAAC","3pc_D8_rep1_ATATCCTCAGTTCTAG","3pc_D8_rep1_CCGAACGGTTGCCGCA","3pc_D8_rep1_CTCCATGGTCGATTAC","3pc_D8_rep1_GTGGGAAAGGCTAGCA","3pc_D8_rep1_TCAAGCATCTTCCTAA","3pc_D8_rep1_TGTTACTTCCGCTTAC","3pc_D8_rep2_CGTGAATTCCAAGGGA","3pc_D8_rep2_GATCAGTGTTTACGAC","3pc_D8_rep2_GCTTTCGTCTTGAACG","3pc_D8_rep2_TATCCTACATGATAGA","3pc_D8_rep2_TCTATCATCACAATGC","3pc_D8_rep2_TGTTACTCAACGTAAA","3pc_D8_rep2_TTAATCCCATGTTACG","3pc_D8_rep2_TTTCACAAGCTGCGAA")
SO.neurons$louvain_clusters3 <- as.character(SO.neurons$louvain_clusters3)
SO.neurons$louvain_clusters3[prdm12_cluster] <- 42
SO.neurons$louvain_clusters3 <- factor(SO.neurons$louvain_clusters3,levels = 0:42)

# Cluster 16 is in between progenitors and neurons, which complicates identification
cluster16_subcluster <- c("0pc_D8_rep1_AAGCGTTGTAGCTGCC","0pc_D8_rep1_ACGTTCCGTATCGTGT","0pc_D8_rep1_AGACCCGGTATTGGCT","0pc_D8_rep1_ATTCACTGTCTGTCCT","0pc_D8_rep1_CAATGACTCTCCTGAC","0pc_D8_rep1_CAGCACGAGGCTAACG","0pc_D8_rep1_CATTGAGGTGATAGTA","0pc_D8_rep1_CCACTTGGTGAGTAAT","0pc_D8_rep1_GAAGTAAGTTTCGCTC","0pc_D8_rep1_GCATTAGAGACATAAC","0pc_D8_rep1_GTGTCCTCAGCACACC","0pc_D8_rep1_TCAATCTAGAGGTTTA","0pc_D8_rep1_TCCGGGACAGAAATCA","0pc_D8_rep1_TCCTCGAAGGATATGT","0pc_D8_rep1_TGCGACGCAAAGTATG","0pc_D8_rep1_TGCTCCAAGGATTTGA","0pc_D8_rep2_AAGGTAACATGGCGCT","0pc_D8_rep2_AGATGAATCGGCTCTT","0pc_D8_rep2_AGTTCGACATTCATCT","0pc_D8_rep2_ATTACCTAGTTGGAGC","0pc_D8_rep2_CATCGTCAGGCCCGTT","0pc_D8_rep2_CATTGCCCAGGTTCCG","0pc_D8_rep2_CCTCATGTCTTCCTAA","0pc_D8_rep2_CTCGAGGTCCATTGGA","0pc_D8_rep2_GTACAGTTCGTTGCCT","0pc_D8_rep2_GTTAGTGTCGGCATTA","0pc_D8_rep2_TCAGGTAGTCAACATC","3pc_D7_rep2_GGTTGTATCGGATACT","3pc_D7_rep2_GTCTACCGTCAGGCAA","3pc_D7_rep2_GTGCACGGTGCATCTA","3pc_D8_rep1_AACAAAGAGTGGTTGG","3pc_D8_rep1_GGGTCACTCACGATAC","3pc_D8_rep2_ACTGATGCACCAATTG","3pc_D8_rep2_ATCGCCTTCGCAGATT","3pc_D8_rep2_ATGAGGGGTTGTGTAC","3pc_D8_rep2_CGTTCTGTCATCACCC","3pc_D8_rep2_CTACTATGTTCGGGTC","3pc_D8_rep2_CTTTCGGGTTGCTCCT","3pc_D8_rep2_GATTGGTGTCGCACGT","3pc_D8_rep2_GCACGGTCATTCATCT","3pc_D8_rep2_GCGTTTCCACCATTCC","3pc_D8_rep2_GGTTCTCTCATGCTAG","3pc_D8_rep2_GTGGGAACATTAGGCT","3pc_D8_rep2_TGTTCTACAACGTAAA","3pc_D8_rep2_TTTCAGTGTTCAGGTT","3pc_D8_rep2_TTTGTTGAGTCTGTAC")
SO.neurons$louvain_clusters3 <- as.character(SO.neurons$louvain_clusters3)
SO.neurons$louvain_clusters3[cluster16_subcluster] <- 43
SO.neurons$louvain_clusters3 <- factor(SO.neurons$louvain_clusters3,levels = 0:43)

# Dot plot of the hox genes
Idents(SO.neurons) <- factor(SO.neurons$louvain_clusters3,levels = names(sort(AverageExpression(object = SO.neurons,features = "Hoxb8")$RNA)))
DotPlot(SO.neurons,features = sort(setdiff(grep("^Hox",rownames(SO.neurons),value=T),grep("^Hox.*[os].*",rownames(SO.neurons),value=T)),decreasing = T),dot.scale = 8)+theme(axis.text.x = element_text(angle=90))
# Hoxa2 negative clusters 20 and 22 -> r1
# Hoxa4,b4,c4 are up until rhombomere 7 -> this is where the bMN are from. 

# Define brain vs spinal cord
SO.neurons$AllHox <- colSums(SO.neurons[grep("Hox",rownames(SO.neurons),value=T),])/length(grep("Hox",rownames(SO.neurons)))*4
Idents(SO.neurons) <- SO.neurons$louvain_clusters3
FeaturePlot(SO.neurons,c("AllHox","Hoxb7","Hoxb8","Hoxb9"),label=T,pt.size=1.5,sort=T)
SO.neurons$SCvsHB <- "SC"
SO.neurons$SCvsHB[SO.neurons$louvain_clusters3 %in% c(1,2,12,20,22,26,35,37,39,40)] <- "HB"
SO.neurons$SCvsHB[SO.neurons$louvain_clusters3 %in% c(8,21,38)] <- "Cervical SC / HB" # Note that 8,21,38 are hox negative but correspond to bMN which are also generated in the cervical spinal cord
Idents(SO.neurons) <- SO.neurons$SCvsHB
DimPlot(SO.neurons,pt.size=2,label=T)+NoLegend()

# Define Neurons vs neural stem cells
Idents(SO.neurons) <- SO.neurons$louvain_clusters3
FeaturePlot(SO.neurons,c("Tubb3","Sox2","Sox1","Nes"),ncol=2,label=T,pt.size=1.5)
SO.neurons$NSCvsNrn <- "Nrn"
SO.neurons$NSCvsNrn[SO.neurons$louvain_clusters3 %in% c(0,1,2,3,4,6,7,9,10,12,13,15,18,21,23,24,26,31,32,33,34,37,38,42,43)] <- "NSC" # Note some cells are double positive for Sox2 and Tubb3 (clusters 2,23,25,28 and 18,11,10). Decide those based on Nes: positive are NSC even if also Tubb3 is +.  
Idents(SO.neurons) <- SO.neurons$NSCvsNrn
DimPlot(SO.neurons,pt.size=2,label=T)+NoLegend()

SO.neurons$neuron_category <- paste(SO.neurons$SCvsHB,SO.neurons$NSCvsNrn,sep="_")
Idents(SO.neurons) <- SO.neurons$neuron_category
DimPlot(SO.neurons,pt.size=1.5)

# Dot plot neurons
Idents(SO.neurons) <- SO.neurons$louvain_clusters3
neurons_ventral_to_dorsal <- c(8,25,14,29,41,27,36,28,17,16,39,22,20,30,5,11,19,35,40)
Idents(SO.neurons) <- factor(SO.neurons$louvain_clusters3,levels = c(neurons_ventral_to_dorsal,setdiff(levels(SO.neurons$louvain_clusters3),neurons_ventral_to_dorsal)))
g <- DotPlot(SO.neurons[,SO.neurons$NSCvsNrn=="Nrn"],features = unique(rev(c("Hoxb7","Hoxb8","Hoxb9","Zic1","Zic5","Sox1","Sox2","Nes","Elavl3","Tubb3",colnames(DelileSign),"Nkx6-1","Phox2b","Tbx20"))),dot.scale = 8)
scale.max <- 2
g$data$avg.exp[g$data$avg.exp>scale.max] <- scale.max # Set the max on the scale manually
g$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
g + geom_point(aes_string(size = 'pct.exp', color = 'avg.exp')) + guides(color = guide_colorbar(title = 'Average Expression')) + theme(axis.text.x = element_text(angle=90))

# Dot plot progenitors
Idents(SO.neurons) <- SO.neurons$louvain_clusters3
isc_ventral_to_dorsal <- c(38,21,10,18,32,6,23,42,31,15,24,7,0,13,4,3,43,9,34,2,33,12,37,1,26)
Idents(SO.neurons) <- factor(SO.neurons$louvain_clusters3,levels = c(neurons_ventral_to_dorsal,isc_ventral_to_dorsal))
g <- DotPlot(SO.neurons[,SO.neurons$NSCvsNrn=="NSC"],features = unique(rev(c("Hoxb7","Hoxb8","Hoxb9","Ptf1a","Zic1","Zic5","Sox1","Sox2","Nes","Elavl3","Tubb3",colnames(DelileSignProg)[1:(length(colnames(DelileSignProg))-4)],"Nkx6-2","Phox2b","Tbx20"))),dot.scale = 8)
scale.max <- 2
g$data$avg.exp[g$data$avg.exp>scale.max] <- scale.max # Set the max on the scale manually
g$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
g + geom_point(aes_string(size = 'pct.exp', color = 'avg.exp')) + guides(color = guide_colorbar(title = 'Average Expression')) + theme(axis.text.x = element_text(angle=90))

# Name the clusters and redo the plots
tmp <- read.table(file = "InputTables/ClustersNamingNeurons.tsv",sep="\t",header = T)
Idents(SO.neurons) <- SO.neurons$louvain_clusters3
levels(Idents(SO.neurons)) <- tmp$celltype
SO.neurons$celltype <- factor(Idents(SO.neurons),levels=rev(c("HB_dI1","HB_dI2-3","SC_dI2-3","SC_dI4","HB_dI5","SC_dI5","SC_dI6","SC_V0","SC_V1","SC_V2a","SC_V2b","SC_MN","HB_bMN","HB_dp1","HB_dp2","HB_dp3","SC_dp3","HB_dp4","SC_dp4","SC_dp5","SC_dp6","SC_p0","SC_p1","SC_p2","SC_pMN","HB_p3")))
Idents(SO.neurons) <- SO.neurons$celltype
DimPlot(SO.neurons,pt.size=1.5,label=T,label.size = 4,repel = T, cols = c(viridis::cividis(20)[8:20],viridis::magma(20)[8:20])) #+ NoLegend() # ,split.by = "Epithelialized" # Export 7x9

# Neurons, export 5.2x22
g <- DotPlot(SO.neurons[,SO.neurons$NSCvsNrn=="Nrn"],features = unique(rev(c("Hoxb7","Hoxb8","Hoxb9","Zic1","Zic5","Sox2","Nes","Elavl3","Tubb3",colnames(DelileSign),"Nkx6-1","Phox2a","Phox2b","Tbx20"))),dot.scale = 8)
scale.max <- 2
g$data$avg.exp[g$data$avg.exp>scale.max] <- scale.max # Set the max on the scale manually
g$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
g + geom_point(aes_string(size = 'pct.exp', color = 'avg.exp')) + guides(color = guide_colorbar(title = 'Average Expression')) + theme(axis.text.x = element_text(angle=90))

# Progenitors, export 4.6x12.5
g <- DotPlot(SO.neurons[,SO.neurons$NSCvsNrn=="NSC"],features = unique(rev(c("Hoxb7","Hoxb8","Hoxb9","Ptf1a","Zic1","Zic5","Sox2","Nes",colnames(DelileSignProg)[1:(length(colnames(DelileSignProg))-4)],"Phox2a","Phox2b","Tbx20"))),dot.scale = 8)
scale.max <- 2
g$data$avg.exp[g$data$avg.exp>scale.max] <- scale.max # Set the max on the scale manually
g$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
g + geom_point(aes_string(size = 'pct.exp', color = 'avg.exp')) + guides(color = guide_colorbar(title = 'Average Expression')) + theme(axis.text.x = element_text(angle=90))

# Note all neurons are Hoxa2 positive, which means they are caudal from r1:r2 limit. Phox2b is in the ventral r2 to r7 and P
FeaturePlot(SO.neurons,c("Hoxa2","Hoxb1"),pt.size=1.5,label=T)

# Also note DRGs should be Isl1 and Sox10 double positive. There are a few, clustered with neural crest and moving towards neurons, but they are not part of the neuron close-up. 
FeaturePlot(SO.neurons,c("Isl1","Sox10"),pt.size=1.5,label=T)
FeaturePlot(SO.ecto,c("Isl1","Sox10","Tubb3","Ngfr","Pou4f1","Nkx1-2","Neurog2"),pt.size=1.5,label=T)

# Reference for Phox2a, Phox2b, Tbx20, stating these branchial motor neurons (bm) can also be in the cervical SC: Hirsch, Goridis, developmental biology, 2006 Forced expression of Phox2
# "spinal accessory motor (nXI) neurons located in the cervical spinal cord are considered bm (branchial motor) neurons, although this has been debated for some time (Krammer et al.,1987). Like their hindbrain counterparts, they arise in the ventralmost progenitor domain, from which they migrate dorsally to their definite location. Their axons leave the neural tube dorsally to form the SAN that runs first alongside the neural tube and innervates the branchial arch-derived muscles of the neck and shoulder girdle (Krammer et al., 1987; Dillon etal., 2005; Lieberam et al., 2005; Pabst et al., 2003). If the nXI neurons were genuine bm neurons, they should express Phox2a and Phox2b and depend on Phox2b activity for proper development (Brunet and Pattyn, 2002). Indeed, we identified a ventrally located subpopulation of the Islet1,2-positive motoneurons that also expressed Phox2b in the cervical spinal cord of 4d-old (HH 22/23) chicken embryos (Figs. 4A, B) and were absent at brachial levels. They should thus correspond to nXI neurons, while the Phox2b-negative Islet1,2-positive cells are sm neurons present in the cervical spinal cord. 
# The picture in this reference also supports Phox2a/Phox2b in the dI5 neurons (Fig.4F)

## Gene plots from list to HD
# Top PCA contributors
gene_list_pca <- vector()
for(i in dims.use.neurons){
  gene_list_pca <- unique(c(gene_list_pca,TopFeatures(object = SO.neurons, dim = i, nfeatures = 15)))
}
# Most variable genes
gene_list_vg <- VariableFeatures(SO.neurons)[1:100]
gene_list <- unique(c(gene_list_nrn))
# Make the plots
dir.plots <- "GenePlotsNeurons"
if(!dir.exists(dir.plots)){dir.create(dir.plots)}
for(gene in intersect(gene_list,rownames(SO.neurons))){
  for(condition in c("Both")){ # "Control", "Epithelialized"
    if(condition=="Both"){
      cells.use <- colnames(SO.neurons)
    }else{
      cells.use <- colnames(SO.neurons)[SO.neurons$Epithelialized==condition]
    }
    png(filename = paste(dir.plots,"/",gene,"_",condition,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.neurons,gene,sort.cell = T,pt.size=1.5,label = T,repel=T,cells = cells.use))
    dev.off()
  }
}
dev.off()

##### Cell proportions for neurons close-up per celltype (collapse the days) in Epi vs Ctrl shown on UMAP #####
# By cell type
df.all <- DataFrame(matrix(nrow=0,ncol=7))
for(condition in c("Control","Epithelialized")){
  SO.tmp <- SO.neurons[,SO.neurons$Epithelialized==condition]
  print(paste("Condition :",condition))
  df.rep <- list()
  for(rep in unique(SO.tmp$Replicate)){
    print(paste("Rep",rep))
    df <- matrix(data = 0,nrow = length(unique(SO.tmp$Day))*length(unique(SO.tmp$celltype)),ncol = 3)
    colnames(df) <- c("Day","Celltype","Percent")
    df <- as.data.frame(df)
    k=1
    for(day in unique(SO.tmp$Day)){
      print(paste("Day",day))
      if(sum(SO.tmp$Day==day & SO.tmp$Replicate==rep)>0){
        tmp <- SO.tmp[,SO.tmp$Day==day & SO.tmp$Replicate==rep]$celltype
        total_cells <- length(tmp)
        for(celltype in unique(SO.tmp$celltype)){
          df[k,"Day"] <- day
          df[k,"Celltype"] <- celltype
          df[k,"Percent"] <- sum(tmp==celltype)/total_cells*100
          k <- k+1
        }
      }else{
        total_cells <- 0
        for(celltype in unique(SO.tmp$celltype)){
          df[k,"Day"] <- day
          df[k,"Celltype"] <- celltype
          df[k,"Percent"] <- 0
          k <- k+1
        }
      }
    }
    df.rep[[rep]] <- df
  }
  df$Percent <- (df.rep[[1]]$Percent+df.rep[[2]]$Percent)/2 # Recycle the temporary df for the two first identical columns, anyway identical in all dfs.
  df$SE <- abs(df.rep[[1]]$Percent-df.rep[[2]]$Percent)/2
  df$PercentReplicate1 <- df.rep[[1]]$Percent
  df$PercentReplicate2 <- df.rep[[2]]$Percent
  #df$Celltype <- factor(df$Celltype,levels = intersect(names(colors.use_celltypes),unique(df$Celltype)))
  pdf(file = paste("ProportionsVsTime/Neurons_Celltype_",condition,".pdf",sep=""),width = 10,height = 6)
  print(ggplot(as.data.frame(df), aes(x = Day, y = Percent, fill = Celltype)) + geom_area(color = "black") + labs(title = "Cell type proportions",subtitle = paste(condition,"gastruloids"),caption = "LSCB 2020",x = "Day",y = "Proportions of cells",fill = "Cell type") + theme_minimal() + scale_fill_manual(aesthetics = "fill",values = colors.use_celltypes))
  dev.off()
  
  df$condition <- condition
  df.all <- rbind(df.all,df)
}
# Compute the proportions
df.all2 <- matrix(nrow=length(unique(df.all$Celltype)),ncol=length(unique(df.all$condition)),dimnames = list(unique(df.all$Celltype),unique(df.all$condition)))
for(condition in unique(df.all$condition)){
  for(celltype in unique(df.all$Celltype)){
    df.all2[celltype,condition] <- sum(df.all[df.all$Celltype==celltype & df.all$condition==condition,"Percent"])
  }
}
df.all2 <- as.data.frame(df.all2/4)
df.all2$proportion <- df.all2$Epithelialized/(df.all2$Control+df.all2$Epithelialized)*100
# Assign them to a metadata column of the seurat objects
SO.neurons$proportion_epi <- 0
for(celltype in rownames(df.all2)){
  SO.neurons$proportion_epi[SO.neurons$celltype==celltype] <- df.all2[celltype,"proportion"]
}
# Make plots (interactive)
if(FALSE){
  colors.tmp <- colorRampPalette(c("Red","#CB9898","Grey","#9797CB","Blue"),interpolate="spline")(100)
  FeaturePlot(SO.neurons,"proportion_epi",cols = colors.tmp[round(min(SO.neurons$proportion_epi)):round(max(SO.neurons$proportion_epi))],cells = sample(colnames(SO.neurons)),pt.size=1.5)+NoLegend()
  # Just for the purpose of getting the full scale bar, can set two cells to 0 and 100%: 
  SO.neurons$proportion_epi[1] <- 0
  SO.neurons$proportion_epi[2] <- 100
  FeaturePlot(SO.neurons,"proportion_epi",cols = colors.tmp[round(min(SO.neurons$proportion_epi)):round(max(SO.neurons$proportion_epi))],cells = sample(colnames(SO.neurons)),pt.size=3)
}
# Alternatively, I could just show which condition each cell belongs to: 
Idents(SO.neurons) <- SO.neurons$Epithelialized
DimPlot(SO.neurons,cols = c("Blue","Red"),cells = sample(colnames(SO.neurons)),pt.size=1.5)+NoLegend()

# Note we don't seem to have dopaminergic neurons
FeaturePlot(SO.neurons,c("Nkx2-2","Lmx1b","Pet2","Fev","Gata3","Ascl1","Gata2"),pt.size=1.5)
# We don't seem to have V3 SC neurons either
FeaturePlot(SO.neurons,c("Sim1","Uncx","Neurog3"),pt.size=1.5)

##### Mesendoderm #####
SO.meso <- SO[,SO$GermLayer=="Mesendoderm"]
SO.meso <- FindVariableFeatures(SO.meso)
SO.meso <- ScaleData(SO.meso,features = VariableFeatures(SO.meso))
SO.meso <- RunPCA(SO.meso)
ElbowPlot(SO.meso)
dims.use.meso <- 1:15
set.seed(72)
SO.meso <- RunHarmony(object = SO.meso, group.by.vars = "Replicate", reduction = "pca", dims.use = dims.use.meso, nclust = 200)
res <- 3
SO.meso <- FindNeighbors(object = SO.meso, reduction = "harmony",dims = dims.use.meso)
SO.meso <- FindClusters(object = SO.meso, resolution = res,random.seed = 72)
SO.meso$louvain_clusters2 <- SO.meso$seurat_clusters
SO.meso@meta.data[,grep("RNA_snn_res",colnames(SO.meso@meta.data))] <- NULL
SO.meso$seurat_clusters <- NULL

# Find a layout, using 100 cells per cluster
#for(seed in c(72,157,324,538,749,953,1254,5278,10285,10285)){ # Screen for seeds to get a layout without trajectory crossings, then fix the seed.
for(seed in c(72)){
  spread <- 2
  min_dist <- 1
  nn=150
  set.seed(seed)
  Idents(SO.meso) <- SO.meso$louvain_clusters2
  cells.use <- WhichCells(object = SO.meso,downsample = 100)
  local_connectivity=1 # Tried 2 and was not so convincing. Should not be more than the local intrinsic dimension of the manifold.
  fast_sgd <- F # Should set it to false ultimately, to get reproducible results, but faster for early exploration.
  umap_init <- "spectral" # "normlaplacian", "spectral" (with noise), "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates.
  reduction.use <- "harmony"
  set.seed(seed)
  tmp <- umap(X = Embeddings(SO.meso[[reduction.use]])[cells.use,dims.use.meso],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,init_sdev = 1e-4,local_connectivity=local_connectivity,ret_model=T,verbose = T,n_epochs = 300)
  tmp2 <- 0*Embeddings(SO.meso[["harmony"]])[,1:2]
  tmp2[cells.use,] <- tmp$embedding
  SO.meso[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
  Idents(SO.meso) <- SO.meso$louvain_clusters2
  png(filename = paste("./Layouts/Meso_Spread",spread,"_mind",min_dist,"_seed",seed,"Louvain_.png",sep=""),width = 1000,height = 1000)
  print(DimPlot(SO.meso,pt.size=1.5,reduction = "umap",label=T,cells = sample(cells.use)) + NoLegend())
  dev.off()
  Idents(SO.meso) <- SO.meso$celltype
  png(filename = paste("./Layouts/Meso_Spread",spread,"_mind",min_dist,"_seed",seed,"_Scmap.png",sep=""),width = 1000,height = 1000)
  print(DimPlot(SO.meso,pt.size=1.5,reduction = "umap",label=T,cells = sample(cells.use)) + NoLegend())
  dev.off()
}

# Pick up the average position of each cluster, as a reference layout:
Idents(SO.meso) <- SO.meso$louvain_clusters2
x_ref <- data.frame(row.names = unique(SO.meso$louvain_clusters2))
for(i in unique(SO.meso$louvain_clusters2)){
  x_ref[i,1] <- mean(SO.meso[["umap"]]@cell.embeddings[intersect(cells.use,names(SO.meso$louvain_clusters2[SO.meso$louvain_clusters2==i])),1])
  x_ref[i,2] <- mean(SO.meso[["umap"]]@cell.embeddings[intersect(cells.use,names(SO.meso$louvain_clusters2[SO.meso$louvain_clusters2==i])),2])
}

# Generate a random position around the layout defined positions for all cells according to their cell types:
noise <- 5
set.seed(seed)
x_ini <- data.frame(row.names = colnames(SO.meso))
x_ini[,1] <- runif(length(colnames(SO.meso)))*noise
x_ini[,2] <- runif(length(colnames(SO.meso)))*noise
Idents(SO.meso) <- SO.meso$louvain_clusters2
for(i in unique(SO.meso$louvain_clusters2)){
  print(i)
  if(i %in% rownames(x_ref)){
    x_ini[WhichCells(SO.meso,idents=i),1] <- x_ini[WhichCells(SO.meso,idents=i),1]+x_ref[i,1]
    x_ini[WhichCells(SO.meso,idents=i),2] <- x_ini[WhichCells(SO.meso,idents=i),2]+x_ref[i,2]
  }
}

# Do an integrated umap initialized on these layout+noise positions:
if(!dir.exists("Idents_UMAP")){dir.create("Idents_UMAP")}
for(min_dist in c(1.5)){ # Use these loops to look for a good umap view # min.dist = 1.75,spread = 1
  for(spread in c(30)){
    for(nn in c(100)){
      cells.use <- colnames(SO.meso) #WhichCells(object = tmp,downsample = 300) #colnames(SO.meso) #sample(colnames(SO.meso),1000)
      local_connectivity=1 # Should not be more than the local intrinsic dimension of the manifold. I would have imagined 2-3 could be reasonable, but doesn't give good results. 
      fast_sgd <- T # Should set it to false ultimately, to get exactly reproducible results, but can use T to get faster for early exploration. 
      umap_init <- as.matrix(x_ini[cells.use,]) # "normlaplacian", "spectral" (with noise),  "random", "lvrandom" (Gaussian std 1e-4), "laplacian", or a matrix of initial coordinates. 
      set.seed(seed)
      reduction.use <- "harmony"
      n_epochs <- 1000
      tmp <- umap(X = Embeddings(SO.meso[[reduction.use]])[cells.use,dims.use.meso],init = umap_init,n_neighbors = nn,n_components = 2,metric = "cosine",min_dist = min_dist,spread = spread,init_sdev = 1e-4,local_connectivity=local_connectivity,ret_model=T,verbose = T,n_epochs = n_epochs)
      tmp2 <- 0*Embeddings(SO.meso[["harmony"]])[,1:2]
      tmp2[cells.use,] <- tmp$embedding
      SO.meso[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = "RNA")
      Idents(SO.meso) <- SO.meso$celltype
      DimPlot(SO.meso,pt.size=1,reduction = "umap",label=T,label.size = 8,repel = T) + NoLegend()
      
      # Plot identities. For final exports: 16 inches height
      for(ident_plot in c("celltype","louvain_clusters2","Stage","Dataset","Replicate")){
        png(filename = paste("./Idents_UMAP/Meso_Spread",spread,"_mind",min_dist,"_nn",nn,"_",ident_plot,".png",sep=""),width = 1000,height = 1000)
        #png(filename = paste(ident_plot,"_",model_plot,".png",sep=""),width = 1000,height = 1000)
        cells.plot <- cells.use
        Idents(SO.meso) <- factor(x = SO.meso@meta.data[,ident_plot],levels = sort(unique(SO.meso@meta.data[,ident_plot])),ordered = T)
        if(ident_plot=="celltype"){
          print(DimPlot(SO.meso,pt.size=1.5,reduction = "umap",label=T,repel=T,cells = sample(cells.plot),label.size = 8) + NoLegend()) #, cols = colors.use_transferred[levels(Idents(SO.meso.align))]
        }else if(ident_plot=="Stage"){
          print(DimPlot(SO.meso,pt.size=1.5,reduction = "umap",label=F,repel=T,cells = sample(cells.plot))+NoLegend()) # , cols = colors.use_stages[levels(Idents(SO.meso.align))]
        }else if(ident_plot=="louvain_clusters2"){
          print(DimPlot(SO.meso,pt.size=1.5,reduction = "umap",label=T,repel=T,cells = sample(cells.plot))+NoLegend()) # , cols = colors.use_stages[levels(Idents(SO.meso.align))]
        }else{
          print(DimPlot(SO.meso,pt.size=1.5,reduction = "umap",label=F,cells = sample(cells.plot)))
        }
        dev.off()
      }
    }
  }
}

# Interactive marker plots to help identify cell types:
if(FALSE){
#SO.meso <- RunUMAP(SO.meso,dims = dims.use.meso,reduction = "harmony")
Idents(SO.meso) <- SO.meso$louvain_clusters2
DimPlot(SO.meso,pt.size=1,label=T,label.size = 4,repel = T) + NoLegend()
FeaturePlot(SO.meso,c("Shh","Noto","Chrd","T","Sox10","Sox17","Mesp1","Ripply2","Dll1"),ncol=3,sort=T) # Mesendoderm markers
FeaturePlot(SO.meso,c("Gata3","Emx2","Spink1","Epcam","Pax3","Pax8","Pax2","Wfdc2","Pcp4"),ncol=3,sort=T) # Intermediate mesoderm -> kidney
FeaturePlot(SO.meso,c("Top2a","Mki67"),pt.size=1.2,sort=T) # Cell cycle
FeaturePlot(SO.meso,c("Slc2a1","Sox3","Bst2","Pkm","Pou3f1"),pt.size=1.2,sort=T) # Cell cycle
FeaturePlot(SO.meso,c("Ripply2","Pax3"),pt.size=1.2,sort=T) # Differentiation front and Somitic mesoderm
FeaturePlot(SO.meso,c("Pax1","Myog","Pax9","Foxc2","Tbx18","Uncx","Twist1","Pax3","Pax7"),pt.size=1.2,sort=T,ncol=3) # Somite differentiation in sclerotome/dermo-myotome
FeaturePlot(SO.meso,c("Pax1","Pax7","Lum","Myf6"),pt.size=1.2,sort=T) # Somite differentiation in sclerotome/dermo-myotome
FeaturePlot(SO.meso,c("Isl1","Nkx2-5"),pt.size=1.2,sort=T) # Pharyngeal mesoderm
FeaturePlot(SO.meso,c("Wnt10a","Wnt16"),pt.size=1.2,sort=T) # Nascent posterior mesoderm, dermomyotome, also Wnt8a is a PS/caudal epi marker and Wnt6 surface ectoderm marker. 
# From https://doi.org/10.1242/dev.034876
FeaturePlot(SO.meso,c("Prom1","Ncam1","Six2","Epcam"),pt.size=1.2,sort=T) # Kidney stem cells -> Six2 appears later at E10.5
FeaturePlot(SO.meso,c("Osr1","Lhx1"),pt.size=1.2,sort=T) # earlyLPM
FeaturePlot(SO.meso,c("Osr2","Foxc1","Foxc2","Gdnf","Hoxa11"),pt.size=1.2,sort=T) # kidney specification
FeaturePlot(SO.meso,c("Foxp1","Eya1"),pt.size=1.2,sort=T) # kidney specification
FeaturePlot(SO.meso,c("Ret","Gfra1"),pt.size=1.2,sort=T) # kidney specification - GDNF axis, excellent markers!
FeaturePlot(SO.meso,c("Ret","Wnt11"),pt.size=1.2,sort=T) # kidney specification - ureteric bud
FeaturePlot(SO.meso,c("Six2","Cited1","Gdnf"),pt.size=1.2,sort=T) # kidney specification - Cap mesenchyme
FeaturePlot(SO.meso,c("Wnt4","Pax8"),pt.size=1.2,sort=T) # kidney specification - Renal vesicle
FeaturePlot(SO.meso,c("Notch2"),pt.size=1.2,sort=T) # kidney specification - nephron. Note nephron forms by associating with endothelium, explains the tendency to converge in UMAP!
FeaturePlot(SO.meso,c("Wt1","Tle4"),pt.size=1.2,sort=T) # kidney specification - podocytes. 
# From https://doi.org/10.1016/j.ydbio.2013.09.026
FeaturePlot(SO.meso,c("Osr1","Eya1","Pax2"),pt.size=1.2,sort=T) # Mesonephros / Early mesonephric mesenchyme 
FeaturePlot(SO.meso,c("Lhx1","Wnt4"),pt.size=1.2,sort=T) # Tubule differentiation
# From https://doi.org/10.1038/s41598-018-37485-8
FeaturePlot(SO.meso,c("Cdh16"),pt.size=1.2,sort=T) # Renal tubular cells. Extremely specific! Canonical marker used to check tubules generation from ES cells. 
FeaturePlot(SO.meso,c("Aqp1","Cdh16","Pax2","Lrp2"),pt.size=1.2,sort=T) # Proximal tubule markers
FeaturePlot(SO.meso,c("Vim","Snai2"),pt.size=1.2,sort=T) # EMT markers
FeaturePlot(SO.meso,c("Emx2","Ccnd1"),pt.size=1.2,sort=T) # Re-epithelialization markers
# From "The zebrafish pronephros: a model to study nephron segmentation"
FeaturePlot(SO.meso,c("Slc20a1","Slc6a14","Slc7a8","Slc4a4","Slc9a3","Slc5a2","Slc5a1","Slc13a3","Slc7a13")) # Mammalian metanephric -> not expressed
FeaturePlot(SO.meso,c("Apeh","Slc12a1","Slc12a3","Kcnj1","Ret","Gata3"),ncol=3,sort=T,pt.size=1.2) # Mammalian metanephric
# Pronephros
FeaturePlot(SO.meso,c("Lhx1","Hnf1b"),pt.size=1.2,sort=T)
# (Pro)Nephric duct (Expected structure at E8.5)
FeaturePlot(SO.meso,c("Lhx1","Ret","Pax2"),pt.size=1.2,sort=T)
# nice kidney markers summary:
FeaturePlot(SO.meso,c("Ret","Gfra1","Wnt11","Epcam","Wnt4","Pax8","Osr1","Pax2","Cdh16","Hnf1b","Gata3"),pt.size=1.2,sort=T) # kidney specification - good markers
# Markers for cluster 40 : it seems to be linked to the notochord, does it differ from sclerotome? -> seems it doesnt, only difference is a bit more COX signaling. 
FeaturePlot(SO.meso,c("Foxd1","Vtn","Pdlim4","Cox4i2","Arg1","Cxcl12"),pt.size=1.2,sort=T)
FeaturePlot(SO.meso,c("Ncam1","Sox9","Foxc2","Cdh11"),pt.size=1.2,sort=T)
FeaturePlot(SO.meso,c("Ndufa4l2","Igfbp2","Cox4i2","Anxa2","Pgm2","Col1a1"),pt.size=1.2,sort=T)
# Early/Nascent mesoderm definition
FeaturePlot(SO.meso,c("Mesp1","Fgf8","T","Fgf17","Wnt3a","Eomes","Tbx6","Msgn1","Isl1"),pt.size=1.2,sort=T) # Note the rostral (eomes) vs caudal (Wnt3a) specification. 
FeaturePlot(SO.meso,c("Hes7","Hoxb6","Gsc","Eomes"),pt.size=1.2,sort=T) # Hes7: caudal. Hoxb6: caudal lateral. Goosecoid: mid PS marker. Eomes: anterior PS.
# Anterior primitive streak, definitive endoderm, gut/foregut:
FeaturePlot(SO.meso,c("T","Eomes","Sox17","Epcam","Spink1","Mixl1","Foxg1","Nkx2-5","Nkx2-3"),pt.size=1.2,sort=T) # Hes7: caudal. Hoxb6: caudal lateral. Goosecoid: mid PS marker. Eomes: anterior PS.
# Lateral plate mesoderm
FeaturePlot(SO.meso,c("Fzd4","Pdgfra","Kdr"),ncol=3,sort=T,pt.size=1.2)
FeaturePlot(SO.meso,c("Foxf1","Pitx2","Isl1","Nkx2-5"),ncol=2,sort=T,pt.size=1.2) # Robust derivation of epicardium and its differentiated smooth muscle cell progeny from human pluripotent stem cells
# Lateral plate vs intermediate vs paraxial markers
FeaturePlot(SO.meso,c("Foxf1","Isl1","Nkx2-5","Osr1","Pax2","Pax8"),ncol=3,sort=T,pt.size=1.2) # Robust derivation of epicardium and its differentiated smooth muscle cell progeny from human pluripotent stem cells. Monitoring and robust induction of nephrogenic intermediate mesoderm from human pluripotent stem cells for Osr1. 
# Caudal lateral primitive mesoderm
FeaturePlot(SO.meso,c("Hoxb6"),sort=T,pt.size=1.2) # Region-specific Etv2 ablation revealed the critical origin of hemogenic capacity from Hox6-positive caudal-lateral primitive mesoderm
FeaturePlot(SO.meso,c("Gata4","Gata6"),sort=T,pt.size=1.2) # Region-specific Etv2 ablation revealed the critical origin of hemogenic capacity from Hox6-positive caudal-lateral primitive mesoderm
# Note that Tbx5, Pitx1 are limb bud markers, and there is a small part of the "pharyngeal mesoderm" on top which is Hox+ and starts to have a bit of those markers, which might correspond to very early hindlimb bud. 
# Limb buds # Cxcr7 not found, Alcam quite flat
FeaturePlot(SO.meso,c("Grem1","Sox9","Bmp4","Hand2","Tbx5","Wnt5a","Sall1","Spry4","Gria2","Irx3","Mab21l1","Pitx1"),sort=T,pt.size=1.3)
FeaturePlot(SO.meso,c("Irx3","Gli3","Asb4","Lhx2","Nr4a2","Hhex","Rarb","Meis2","Tbx1","Pkdcc","Cyp26b1"),sort=T,pt.size=1.3)
FeaturePlot(SO.meso,c("Fgf8","Shh","Osr1","Osr2","Nr2f2","Gbx2","Gsc","Eomes"),sort=T,pt.size=1.3)
FeaturePlot(SO.meso,c("Fgf10","Jag1","Tbx5"),sort=T,pt.size=1.3)
}

# Cluster naming in the mesendoderm #
tmp <- read.table(file = "InputTables/ClustersNamingMeso.tsv",sep="\t",header = T)
Idents(SO.meso) <- SO.meso$louvain_clusters2
levels(Idents(SO.meso)) <- tmp$celltype
SO.meso$celltype <- as.character(Idents(SO.meso))
DimPlot(SO.meso,pt.size=1,label=T,label.size = 4,repel = F,split.by = "Epithelialized") + NoLegend()

# Export plots of gene expression on UMAP of mesendoderm
# Custom list
gene_list <- readLines(con = "InputTables/Mesendoderm_Highlights.txt")
gene_list <- c("Noto")
# Make the plots
dir.plots <- "GenePlotsMesendoderm"
Idents(SO.meso) <- SO.meso$celltype
if(!dir.exists(dir.plots)){dir.create(dir.plots)}
for(gene in intersect(gene_list,rownames(SO.meso))){
  for(condition in c("Both")){ #,"Control","Epithelialized"
    if(condition=="Both"){
      cells.use <- colnames(SO.meso)
    }else{
      cells.use <- colnames(SO.meso)[SO.ecto$Epithelialized==condition]
    }
    png(filename = paste(dir.plots,"/",gene,"_",condition,".png",sep=""),width = 800,height = 800)
    print(FeaturePlot(SO.meso,gene,sort.cell = T,pt.size=1.5,label = T,repel=T,cells = cells.use))
    dev.off()
  }
}
dev.off()

##### ExE #####
SO.exe <- SO[,SO$GermLayer=="ExE"]
SO.exe <- FindVariableFeatures(SO.exe)
SO.exe <- ScaleData(SO.exe,features = VariableFeatures(SO.exe))
SO.exe <- RunPCA(SO.exe)
ElbowPlot(SO.exe)
dims.use.exe <- 1:10
set.seed(72)
SO.exe <- RunHarmony(object = SO.exe, group.by.vars = "Replicate", reduction = "pca", dims.use = dims.use.exe, nclust = 20)
res <- 4
SO.exe <- FindNeighbors(object = SO.exe, reduction = "harmony",dims = dims.use.exe)
SO.exe <- FindClusters(object = SO.exe, resolution = res,random.seed = 72)
SO.exe$louvain_clusters2 <- SO.exe$seurat_clusters
SO.exe@meta.data[,grep("RNA_snn_res",colnames(SO.exe@meta.data))] <- NULL
SO.exe$seurat_clusters <- NULL
SO.exe <- RunUMAP(SO.exe,dims = dims.use.exe,reduction="harmony",n.neighbors = 150,min.dist = 1,spread = 1)
DimPlot(SO.exe,label=T,pt.size=1.5)+NoLegend()

# Remove doublets contaminated with neural ectoderm
# Neuron signature is not particularly localized on any cell type - probably not too much mature neuron contamination
FeaturePlot(SO.exe, c("Nefl","Nefm","Myt1","Nav2","Kif5c","Klhl35","Tubb3","L1cam","Lhx1"),label=T,pt.size=2,sort=T,ncol=3)+NoLegend() 
# Some clusters are strongly positive for neural stem cells markers - exclude
FeaturePlot(SO.exe, c("Sox1","Bex4","Sox11","Crabp1","Crabp2","Bex2","Cd24a","Tspan13","Mex3a"),label=T,pt.size=2,sort=T,ncol=3)+NoLegend() 
SO.exe <- AddMetaData(SO.exe,metadata = colSums(SO.exe@assays$RNA@data[c("Sox1","Bex4","Sox11","Crabp1","Crabp2","Bex2","Cd24a"),]),"NeuralContamination")
FeaturePlot(SO.exe, "NeuralContamination",label=T,pt.size=2,sort=T)
#SO.exe <- SO.exe[,!(SO.exe$louvain_clusters2 %in% c(8,9,11)) & SO.exe@assays$RNA@counts["Sox1",]==0] # Earlier version of the filtering, update using the larger signature instead of just Sox1, which misses many cells. 
SO$celltype[colnames(SO.exe)[SO.exe$NeuralContamination>1.5]] <- "ExE - Neural Doublets"
SO.exe <- SO.exe[,!SO.exe$NeuralContamination>1.5]
# Re-run the analysis:
SO.exe <- FindVariableFeatures(SO.exe,500)
SO.exe <- ScaleData(SO.exe,features = VariableFeatures(SO.exe))
SO.exe <- RunPCA(SO.exe)
ElbowPlot(SO.exe)
dims.use.exe <- 1:10
set.seed(72)
SO.exe <- RunHarmony(object = SO.exe, group.by.vars = "Replicate", reduction = "pca", dims.use = dims.use.exe, nclust = 20)
res <- 5
SO.exe <- FindNeighbors(object = SO.exe, reduction = "harmony",dims = dims.use.exe)
SO.exe <- FindClusters(object = SO.exe, resolution = res,random.seed = 72)
SO.exe$louvain_clusters2 <- SO.exe$seurat_clusters
SO.exe@meta.data[,grep("RNA_snn_res",colnames(SO.exe@meta.data))] <- NULL
SO.exe$seurat_clusters <- NULL
SO.exe <- RunUMAP(SO.exe,dims = dims.use.exe,reduction="harmony",n.neighbors = 150,min.dist = 5,spread = 5)
DimPlot(SO.exe,label=T,pt.size=1.5)+NoLegend()

# Identify the cell types (interactive plots)
if(FALSE){
FeaturePlot(SO.exe,c("Cdx2","Elf5","Eomes","Bmp4"),pt.size=1.5,sort=T) # Extraembryonic ectoderm (Stem cells) -> Should be Esrrb negative, which they are not.
FeaturePlot(SO.exe,c("Cdx2","Elf5","Eomes","Bmp4","Esrrb"),ncol=3,pt.size=1.5,sort=T) # Chorion -> this is what we have
FeaturePlot(SO.exe,c("Ascl2","Hand1"),pt.size=1.5,sort=T) # Ectoplacental cone
FeaturePlot(SO.exe,c("Tpbpa","Ascl2","Tcf4","Cdx2","Flt1"),pt.size=1.5,sort=T) # Spongiotrophoblasts
FeaturePlot(SO.exe,c("Ncoa6","Hsf1","Ascl2","Krt8","Krt19"),pt.size=1.5,sort=T) # Spongiotrophoblasts - Trophoblast lineages 3 in histology folder
FeaturePlot(SO.exe,c("Prl3b1", "Prl2c", "Prl8a6", "Prl8a8", "Prl8a1", "Prl7a2", "Prl3a1", "Prl2b1", "Prl5a1"),pt.size=1.5,sort=T) # Spongiotrophoblasts signature from prolactin paper. 8a8 is the only one totally specific, but it is expressed later. Prl3b1 comes early, but is also expressed in TGCs. 
FeaturePlot(SO.exe,c("Gcm1","Syna","Synb","Ovol2","Dlx3","Cebpa","Tfeb","Esx1"),pt.size=1.5,sort=T) # Labyrinth progenitors (Gcm1 and Synb are type2 progenitors, Syna type1 progenitors). Note Esx1 is also expected in chorion, which we see. 
FeaturePlot(SO.exe,c("Prl7a1","Prl3d1","Prl3b1","Prl2c2","Cts7","Cts8","Ascl2","Ctsq"),pt.size=1.5,sort=T) # parietal trophoblast giant cells. Ascl2 and Ctsq should be negative, others pos. Ctsq marks only all the others. 
# Outhwaite 2016: A specialized subtype of trophoblast giant cells (TGCs) line the torturous sinusoids of the murine placental labyrinth, and can be distinguished from most other TGCs by the expression of Ctsq.  
# In many papers, Ctsq considered a sinusoidal TGC marker. 
FeaturePlot(SO.exe,c("Prl3b1","Ctsq"),pt.size=1.5,sort=T) # Double positive are Sinusoidal trophoblast giant cells
FeaturePlot(SO.exe,c("Prl2c2","Ctsq"),pt.size=1.5,sort=T) # Double positive are Channel trophoblast giant cells
FeaturePlot(SO.exe,c("Prl3d1","Prl3b1","Prl2c2","Prl2c5","Cts7","Cts8","Ctsq")) # All trophoblast giant cells (including diverse subpopulations). Should be Ascl2 negative.
# Suspected syncitiotrophoblasts (enriched genes rather than literature markers):
FeaturePlot(SO.exe,c("Prl8a9","Foxo4","Nupr1","Cited2","Klk6","Gadd45b","Hist1h2bc","Gcm1"),pt.size=1.5,sort=T)
# Syncitiotrophoblasts (Foxo1+, Cdh1-, Cdh1 marking villous cytotrophoblasts)
FeaturePlot(SO.exe,c("Foxo1","Cdh1"),pt.size=1.5,sort=T)
# Markers from trophoblast lineages 1 in histology folder: 
FeaturePlot(SO.exe,c("Cdx2","Esrrb","Eomes","Ascl2","Tpbpa","Flt1","Hand1","Cenpx","Mdfic"),ncol=3,pt.size=1.5,sort=T)
# Important secreted factors:
FeaturePlot(SO.exe,c("Bmp4","Wnt7b"),pt.size=1.5,sort=T) 
# Other QC and checks
FeaturePlot(SO.exe,c("TSCELLS","Tubb3","Sox1","Ctsq"),pt.size=1.5,sort=T) 
FeaturePlot(SO.exe,c("Meg3","Mest","Sox11","Cd24a","Bex2","Bex4","Crabp2"),pt.size=1.5,sort=T) 
FeaturePlot(SO.exe,c("Shroom1","Prl2c5","Prl3d3","Prl7a1"),pt.size=1.5,sort=T) 
# Genes enriched with the cluster identified as parietal trophoblast cells
FeaturePlot(SO.exe,c("Prl7a1","Tfpi","Prl2c5","Plac1","Procr","Serpine1","Prl2c3","Cd52","Ecm1","Shroom1","Plac8","Anxa2")) # Prl proteins are placental proteins, Simmons BMC genomics 2008
# Best short signature to define all cells: 
# Essrb: chorion, Ascl2 ectoplacental cone, Plac8 are trophoblast giant cells (Prl7a1 and Ctsq marck the two main subpopulations, parietal and others), Dlx3 are labyrinth, Tpbpa spongiotrophoblasts. 
FeaturePlot(SO.exe,c("Esrrb","Ascl2","Ctsq","Prl7a1","Dlx3","Tpbpa"),ncol=3,pt.size=1.5,sort=T) # ,"Plac8"
# For assembling gene overlay:
FeaturePlot(SO.exe,c("Esrrb","Ascl2","Ctsq","Prl7a1","Dlx3","Tpbpa"),ncol=3,pt.size=2.5,sort=T,cols = c("White","Red"),min.cutoff = 0,max.cutoff = "q98")
}

# Cluster naming in the extraembryonic cells
tmp <- read.table(file = "InputTables/ClustersNamingExE.tsv",sep="\t",header = T)
Idents(SO.exe) <- SO.exe$louvain_clusters2
levels(Idents(SO.exe)) <- tmp$celltype
SO.exe$celltype <- as.character(Idents(SO.exe))
DimPlot(SO.exe,pt.size=1,label=T,label.size = 4,repel = F) + NoLegend()

##### Manually pick up cells to find markers #####
if(FALSE){
plot1 <- FeaturePlot(SO.meso,"Ripply2",pt.size = 1,reduction = "umap",sort=T)
cells.tmp <- CellSelector(plot1)
markers <- FindMarkers(object = SO.meso, ident.1 = cells.tmp)
markers.sign <- markers[markers$p_val_adj<1e-3 & markers$avg_logFC > log(1.2),]
markers.sign <- markers.sign[order(-markers.sign$avg_logFC),]
head(markers.sign,50)
writeClipboard(rownames(markers.sign))
gene_list <- rownames(markers.sign)
write.table(x=markers.sign,file = "celltypemarkers.tsv",sep = "\t",row.names = F,col.names = T)
}

# Add into SO the info about the close-ups: 
SO$celltype[colnames(SO.exe)]  <- SO.exe$celltype 
SO$celltype[colnames(SO.ecto)] <- SO.ecto$celltype
SO$celltype[colnames(SO.meso)] <- SO.meso$celltype
SO$celltype[SO$GermLayer=="PGC"] <- "PGC"
SO$louvain_clusters_closeups <- "NA"
SO$louvain_clusters_closeups[colnames(SO.exe)]  <- paste0("ExE_",SO.exe$louvain_clusters2)
SO$louvain_clusters_closeups[colnames(SO.ecto)] <- paste0("Epiblast-Ectoderm_",SO.ecto$louvain_clusters2)
SO$louvain_clusters_closeups[colnames(SO.meso)] <- paste0("Mesendoderm_",SO.meso$louvain_clusters2)

##### Load custom colors for all cell types #####
colors.table <- read.table(file = "InputTables/ClusterColors.tsv", sep="\t", header = T, comment.char = "", as.is = T)
colors.use_celltypes   <- setNames(colors.table$color,colors.table$celltype)
colors.use_stages <- setNames(c('#5e5bf4','#f854ee','#f85554','#f8ea54'),c("Day5","Day6","Day7","Day8"))
colors.use_replicates <- c("#eed239","#7070fc")
colors.use_germlayers <- setNames(c("#ffd100","#474747","#9f9f9f","#6f8def","#f56969"),c("ExE","PGC","Epiblast","Ectoderm","Mesendoderm"))
  
# Make colored ident plots. Export 12x12 in for figure. 
Idents(SO) <- SO$celltype
DimPlot(SO,label=T,pt.size=1,cols = colors.use_celltypes)+NoLegend()

Idents(SO) <- SO$GermLayer
DimPlot(SO,label=T,pt.size=1,cols = colors.use_germlayers)+NoLegend()

Idents(SO) <- SO$Stage
DimPlot(SO,label=T,pt.size=1,cols = colors.use_stages)+NoLegend()

Idents(SO.ecto) <- SO.ecto$celltype
DimPlot(SO.ecto,label=T,pt.size=1.5,cols = colors.use_celltypes)+NoLegend()

Idents(SO.ecto) <- SO.ecto$Stage
DimPlot(SO.ecto,label=T,pt.size=1.5,cols = colors.use_stages)+NoLegend()

Idents(SO.meso) <- SO.meso$celltype
DimPlot(SO.meso,label=T,pt.size=1.3,cols = colors.use_celltypes)+NoLegend()

Idents(SO.meso) <- SO.meso$Stage
DimPlot(SO.meso,label=T,pt.size=1.3,cols = colors.use_stages)+NoLegend()

Idents(SO.exe) <- SO.exe$celltype
DimPlot(SO.exe,label=T,pt.size=3,cols = colors.use_celltypes)+NoLegend()

Idents(SO.exe) <- SO.exe$Stage
DimPlot(SO.exe,label=T,pt.size=3,cols = colors.use_stages)+NoLegend()

##### Proportions vs time #####
# By cell type
df.all <- DataFrame(matrix(nrow=0,ncol=7))
for(condition in c("Control","Epithelialized")){
  SO.tmp <- SO[,SO$Epithelialized==condition]
  # SO.tmp <- SO[,SO$Epithelialized==condition & SO$GermLayer%in%c("Epiblast","Ectoderm")] # To restrict to the ecto close-up
  # SO.tmp <- SO[,SO$Epithelialized==condition & SO$GermLayer=="Mesendoderm"] # To restrict to the meso close-up
  print(paste("Condition :",condition))
  df.rep <- list()
  for(rep in unique(SO.tmp$Replicate)){
    print(paste("Rep",rep))
    df <- matrix(data = 0,nrow = length(unique(SO.tmp$Day))*length(unique(SO.tmp$celltype)),ncol = 3)
    colnames(df) <- c("Day","Celltype","Percent")
    df <- as.data.frame(df)
    k=1
    for(day in unique(SO.tmp$Day)){
      print(paste("Day",day))
      tmp <- SO.tmp[,SO.tmp$Day==day & SO.tmp$Replicate==rep]$celltype
      total_cells <- length(tmp)
      for(celltype in unique(SO.tmp$celltype)){
        df[k,"Day"] <- day
        df[k,"Celltype"] <- celltype
        df[k,"Percent"] <- sum(tmp==celltype)/total_cells*100
        k <- k+1
      }
    }
    df.rep[[rep]] <- df
  }
  df$Percent <- (df.rep[[1]]$Percent+df.rep[[2]]$Percent)/2 # Recycle the temporary df for the two first identical columns, anyway identical in all dfs.
  df$SE <- abs(df.rep[[1]]$Percent-df.rep[[2]]$Percent)/2
  df$PercentReplicate1 <- df.rep[[1]]$Percent
  df$PercentReplicate2 <- df.rep[[2]]$Percent
  df$Celltype <- factor(df$Celltype,levels = intersect(names(colors.use_celltypes),unique(df$Celltype)))
  pdf(file = paste("ProportionsVsTime/Celltype_",condition,".pdf",sep=""),width = 10,height = 6)
  print(ggplot(as.data.frame(df), aes(x = Day, y = Percent, fill = Celltype)) + geom_area(color = "black") + labs(title = "Cell type proportions",subtitle = paste(condition,"gastruloids"),caption = "LSCB 2020",x = "Day",y = "Proportions of cells",fill = "Cell type") + theme_minimal() + scale_fill_manual(aesthetics = "fill",values = colors.use_celltypes))
  dev.off()
  
  df$condition <- condition
  df.all <- rbind(df.all,df)
}
write.table(df.all,"OutputTables/CellProportions.tsv",sep="\t",row.names = F)

# By Germ Layer
condition <- "Epithelialized"
condition <- "Control"
SO.tmp <- SO[,SO$Epithelialized==condition]
df <- matrix(data = 0,nrow = length(unique(SO.tmp$Day))*length(unique(SO.tmp$GermLayer)),ncol = 3)
colnames(df) <- c("Day","GermLayer","Number")
df <- as.data.frame(df)
k=1
for(i in unique(SO.tmp$Day)){
  print(paste("Day",i))
  tmp <- SO.tmp[,SO.tmp$Day==i]$GermLayer
  total_cells <- length(tmp)
  for(j in unique(SO.tmp[,]$GermLayer)){
    df[k,"Day"] <- i
    df[k,"GermLayer"] <- j
    df[k,"Number"] <- sum(tmp==j)/total_cells*100
    k <- k+1
  }
}
df$GermLayer <- factor(df$GermLayer,levels = intersect(names(colors.use_germlayers),unique(df$GermLayer)))
ggplot(as.data.frame(df), aes(x = Day, y = Number, fill = GermLayer)) + geom_area(color = "black") + labs(title = "Germ layer proportions",subtitle = paste(condition,"gastruloids"),caption = "LSCB 2020",x = "Day",y = "Proportions of cells",fill = "Cell type") + theme_minimal() + scale_fill_manual(aesthetics = "fill",values = colors.use_germlayers)

##### Cell proportions per celltype (collapse the days) in Epithelialized vs Ctrl shown on UMAP #####
# Compute the proportions
df.all2 <- matrix(nrow=length(unique(df.all$Celltype)),ncol=length(unique(df.all$condition)),dimnames = list(unique(df.all$Celltype),unique(df.all$condition)))
for(condition in unique(df.all$condition)){
  for(celltype in unique(df.all$Celltype)){
    df.all2[celltype,condition] <- sum(df.all[df.all$Celltype==celltype & df.all$condition==condition,"Percent"])
  }
}
df.all2 <- as.data.frame(df.all2/4)
df.all2$proportion <- df.all2$Epithelialized/(df.all2$Control+df.all2$Epithelialized)*100
# Assign them to a metadata column of the seurat objects
SO$proportion_epi <- 0
SO.exe$proportion_epi <- 0
SO.ecto$proportion_epi <- 0
SO.meso$proportion_epi <- 0
for(celltype in rownames(df.all2)){
  SO$proportion_epi[SO$celltype==celltype] <- df.all2[celltype,"proportion"]
  SO.exe$proportion_epi[SO.exe$celltype==celltype] <- df.all2[celltype,"proportion"]
  SO.meso$proportion_epi[SO.meso$celltype==celltype] <- df.all2[celltype,"proportion"]
  SO.ecto$proportion_epi[SO.ecto$celltype==celltype] <- df.all2[celltype,"proportion"]
}
# Make plots (interactive)
if(FALSE){
colors.tmp <- colorRampPalette(c("Red","#CB9898","Grey","#9797CB","Blue"),interpolate="spline")(100)
FeaturePlot(SO,"proportion_epi",cols = colors.tmp[round(min(SO$proportion_epi)):round(max(SO$proportion_epi))],cells = sample(colnames(SO)),pt.size=1)+NoLegend()
FeaturePlot(SO.exe,"proportion_epi",cols = colors.tmp[round(min(SO.exe$proportion_epi)):round(max(SO.exe$proportion_epi))],cells = sample(colnames(SO.exe)),pt.size=3)+NoLegend()
FeaturePlot(SO.ecto,"proportion_epi",cols = colors.tmp[round(min(SO.ecto$proportion_epi)):round(max(SO.ecto$proportion_epi))],cells = sample(colnames(SO.ecto)),pt.size=1.5)+NoLegend()
FeaturePlot(SO.meso,"proportion_epi",cols = colors.tmp[round(min(SO.meso$proportion_epi)):round(max(SO.meso$proportion_epi))],cells = sample(colnames(SO.meso)),pt.size=1.3)+NoLegend()
# Just for the purpose of getting the full scale bar, can set two cells to 0 and 100%: 
SO.exe$proportion_epi[1] <- 0
SO.exe$proportion_epi[2] <- 100
FeaturePlot(SO.exe,"proportion_epi",cols = colors.tmp[round(min(SO.exe$proportion_epi):round(max(SO.exe$proportion_epi)))],cells = sample(colnames(SO.exe)),pt.size=3)
}

# Add the close-ups UMAP coordinates, and the neuron close-up detailed annotations
SO$exe_closeup_UMAP_1 <- SO.exe[["umap"]]@cell.embeddings[,"UMAP_1"]
SO$exe_closeup_UMAP_2 <- SO.exe[["umap"]]@cell.embeddings[,"UMAP_2"]
SO$ecto_closeup_UMAP_1 <- SO.ecto[["umap"]]@cell.embeddings[,"UMAP_1"]
SO$ecto_closeup_UMAP_2 <- SO.ecto[["umap"]]@cell.embeddings[,"UMAP_2"]
SO$meso_closeup_UMAP_1 <- SO.meso[["umap"]]@cell.embeddings[,"UMAP_1"]
SO$meso_closeup_UMAP_2 <- SO.meso[["umap"]]@cell.embeddings[,"UMAP_2"]
SO$neuron_closeup_UMAP_1 <- SO.neurons[["umap"]]@cell.embeddings[,"UMAP_1"]
SO$neuron_closeup_UMAP_2 <- SO.neurons[["umap"]]@cell.embeddings[,"UMAP_2"]
SO$neuron_closeup_celltypes <- SO.neurons$celltype

# Re-export for sharing, including the close-ups annotations (go to the folder with bash and run "gzip *" to compress the files as is usual)
if(!dir.exists(data_export_folder)){dir.create(data_export_folder)}
tmp <- GetAssayData(object = SO,slot = "counts",assay = "RNA")
write(colnames(tmp), file = paste0(data_export_folder,"/barcodes.tsv"))
write(rownames(tmp), file = paste0(data_export_folder,"/features.tsv"))
writeMM(obj = tmp, file = paste0(data_export_folder,"/matrix.mtx"))
tmp <- cbind(SO@meta.data, Embeddings(SO[["harmony"]])[,dims.use])
colnames(tmp)[grep("harmony",colnames(tmp))] <- paste0("overview_harmony",dims.use)
tmp[rownames(Embeddings(SO.ecto[["harmony"]])[,dims.use.ecto]),paste0("ecto_harmony",dims.use.ecto)] <- Embeddings(SO.ecto[["harmony"]])[,dims.use.ecto]
tmp[rownames(Embeddings(SO.meso[["harmony"]])[,dims.use.meso]),paste0("meso_harmony",dims.use.meso)] <- Embeddings(SO.meso[["harmony"]])[,dims.use.meso]
tmp[rownames(Embeddings(SO.exe[["harmony"]])[,dims.use.exe]),paste0("exe_harmony",dims.use.exe)] <- Embeddings(SO.exe[["harmony"]])[,dims.use.exe]
tmp[rownames(Embeddings(SO.neurons[["harmony"]])[,dims.use.neurons]),paste0("neurons_harmony",dims.use.neurons)] <- Embeddings(SO.neurons[["harmony"]])[,dims.use.neurons]
tmp %>% head
write.table(x = tmp, file = paste0(data_export_folder,"/metadata.tsv"), sep = "\t", row.names = T, col.names = NA)

##### Prepare data for re-import into python/scanpy for scVelo RNA velocity analysis #####
# Name of the close-up / re-export
getwd()
if(!dir.exists("OutputTables")){dir.create("OutputTables")}
for(name.tmp in c("Overview","Mesendoderm","Ectoderm","Extraembryonic","Neurons")){
  if(name.tmp == "Overview"){SO.tmp <- SO}
  if(name.tmp == "Mesendoderm"){SO.tmp <- SO.meso}
  if(name.tmp == "Ectoderm"){SO.tmp <- SO.ecto}
  if(name.tmp == "Extraembryonic"){SO.tmp <- SO.exe}
  if(name.tmp == "Neurons"){SO.tmp <- SO.neurons}
  write.table(x = SO.tmp@meta.data,file = paste("./OutputTables/",name.tmp,"_metadata.tsv",sep=""),append = F,quote = F,sep = "\t",row.names = T,col.names = NA)
  write.table(x = cbind(Embeddings(SO.tmp[["umap"]]),Embeddings(SO.tmp[["harmony"]])), file = paste("./OutputTables/",name.tmp,"_DR.tsv",sep=""),append = F,quote = F,sep = "\t",row.names = T,col.names = NA)
  write.table(x = VariableFeatures(SO.tmp),file = paste("./OutputTables/",name.tmp,"_varfeats.tsv",sep=""),append = F,quote = F,sep = "\n",row.names = F,col.names = F)
}
rm(SO.tmp)

##### Ectoderm markers heatmap #####
max.cells.per.cluster <- 500
cluster.order <- intersect(names(colors.use_celltypes),unique(SO.ecto$celltype))
markers <- read.table(file = "./InputTables/Ectoderm_Highlights.txt",sep = "\n",header=F,as.is = T)$V1
Idents(SO.ecto) <- SO.ecto$celltype
Idents(SO.ecto) <- factor(Idents(SO.ecto),levels = cluster.order)
cells.plot <- c()
for(cluster in cluster.order){
  cells.tmp <- WhichCells(object = SO.ecto,idents = cluster)
  cells.plot <- c(cells.plot,sample(cells.tmp,size = min(max.cells.per.cluster,length(cells.tmp))))
}
cells.plot <- names(sort(-SO.ecto[["umap"]]@cell.embeddings[cells.plot,1]+SO.ecto[["umap"]]@cell.embeddings[cells.plot,2]))
if(!dir.exists("Heatmaps")){dir.create("Heatmaps")}
#png(filename = paste("./Heatmaps/Ectoderm_highlights1.png",sep=""),width = 1000,height = 1000)
pdf(file = paste("./Heatmaps/Ectoderm_highlights4.pdf",sep=""), width = 12, height = 10)
DoHeatmap(object = SO.ecto, cells = cells.plot, disp.max = 4, features = markers, slot = "data",raster = T,draw.lines = T, lines.width = 15,group.colors = colors.use_celltypes[levels(Idents(SO.ecto))])+ scale_fill_gradientn(colors = c("white","blue","black"))#+NoLegend()
dev.off()

##### Mesendoderm markers heatmap #####
max.cells.per.cluster <- 500
cluster.order <- intersect(names(colors.use_celltypes),unique(SO.meso$celltype))
markers <- read.table(file = "./InputTables/Mesendoderm_Highlights.txt",sep = "\n",header=F,as.is = T)$V1
Idents(SO.meso) <- SO.meso$celltype
Idents(SO.meso) <- factor(Idents(SO.meso),levels = cluster.order)
cells.plot <- c()
for(cluster in cluster.order){
  cells.tmp <- WhichCells(object = SO.meso,idents = cluster)
  cells.plot <- c(cells.plot,sample(cells.tmp,size = min(max.cells.per.cluster,length(cells.tmp))))
}
cells.plot <- names(sort(SO.meso[["umap"]]@cell.embeddings[cells.plot,2]))
if(!dir.exists("Heatmaps")){dir.create("Heatmaps")}
#png(filename = paste("./Heatmaps/Mesendoderm_highlights1.png",sep=""),width = 1000,height = 1000)
pdf(file = paste("./Heatmaps/Mesendoderm_highlights6.pdf",sep=""), width = 12, height = 10)
DoHeatmap(object = SO.meso, cells = cells.plot, disp.max = 4, features = markers, slot = "data",raster = T,draw.lines = T, lines.width = 15,group.colors = colors.use_celltypes[levels(Idents(SO.meso))])+ scale_fill_gradientn(colors = c("white","blue","black"))#+NoLegend()
dev.off()

FeaturePlot(SO.ecto,c("Vim","Ckb","Hoxb9","Hoxb9as","Gas1","Metrn","Zic1","Plagl1","Pantr1","Hoxb8","Igfbp2","Pax3"),sort=T,pt.size=1.2)

##### ExE markers heatmap #####
max.cells.per.cluster <- 100
cluster.order.exe <- c("Chorion","Ectoplacental Cone","Labyrinth progenitors","Parietal Trophoblast Giant Cell","Sinusoidal Trophoblast Giant Cell","Spongiotrophoblast")
Idents(SO.exe) <- SO.exe$celltype
Idents(SO.exe) <- factor(Idents(SO.exe),levels = cluster.order.exe)
# markers.all <- FindAllMarkers(SO.exe,min.pct = 0.2,min.diff.pct = 0.2,logfc.threshold = log(1.5),only.pos = T)
# markers.keep <- markers.all[markers.all$p_val_adj<0.01 & markers.all$avg_logFC>log(2) & abs(markers.all$pct.1-markers.all$pct.2)>0.2,]
# markers.keep$score <- 6*abs(markers.keep$pct.1-markers.keep$pct.2) # Make a custom measure for prioritizing, mix of logFC and dropout pct. 
# markers.keep %>% group_by(cluster) %>% top_n(12,score) %>% as.data.frame -> markers.keep
# markers <- markers.keep$gene
markers <- read.table(file = "./InputTables/Extraembryonic_Highlights.txt",sep = "\n",header=F,as.is = T)$V1
cells.plot <- c()
for(cluster in cluster.order.exe){
  cells.tmp <- WhichCells(object = SO.exe,idents = cluster)
  cells.plot <- c(cells.plot,sample(cells.tmp,size = min(max.cells.per.cluster,length(cells.tmp))))
}
cells.plot <- names(sort(-SO.exe[["umap"]]@cell.embeddings[cells.plot,1]-SO.exe[["umap"]]@cell.embeddings[cells.plot,2]))
if(!dir.exists("Heatmaps")){dir.create("Heatmaps")}
#png(filename = paste("./Heatmaps/Ectoderm_highlights1.png",sep=""),width = 1000,height = 1000)
pdf(file = paste("./Heatmaps/Extraembryonic_highlights.pdf",sep=""), width = 6, height = 10)
DoHeatmap(object = SO.exe, cells = cells.plot, disp.max = 4, features = markers, slot = "data",raster = T,draw.lines = T, lines.width = 3,group.colors = colors.use_celltypes[levels(Idents(SO.exe))])+ scale_fill_gradientn(colors = c("white","blue","black"))#+NoLegend()
dev.off()
