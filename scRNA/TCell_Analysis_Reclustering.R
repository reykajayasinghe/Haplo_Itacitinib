#########################################
##ANALYSIS AND CLUSTERING OF T CELLS#####
#########################################


library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
#library(SingleR)
#library(celldex)
library(dplyr)
library(tidyr)
library(ggrepel)
library(tidyverse)
library(dittoSeq)

obj<-readRDS("tcell_harmony_res0.5_dim10.rds")

# convert a v5 assay to a v3 assay
obj[["RNA3"]] <- as(object = obj[["RNA"]], Class = "Assay")
DefaultAssay(obj)<-"RNA3"

#https://satijalab.org/seurat/articles/integration_rpca
###Try anchoring to improve cluster quality

#Remove several samples with less than 50 T cells
Idents(obj)<-obj@meta.data$sample
subset_combined<-subset(obj,idents=c("Control_7_Day_28-5","Patient_29_Day_28-19","Patient-016-Day-28-4","Control-1-Day-28-13","Patient_29_Donor-21","Patient_33_Day_28-24"),invert=TRUE)

ifnb.list <- SplitObject(subset_combined, split.by = "sample")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	DefaultAssay(x)<-"RNA3"
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

print("step 1 done")

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    DefaultAssay(x)<-"RNA3"
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

print("step 2 done")

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list,anchor.features = features, k.filter = 200, reduction = "rpca")
print("step 3 done")

timmune.combined <- IntegrateData(anchorset = immune.anchors,k.weight = 40)
DefaultAssay(timmune.combined) <- "integrated"
print("step 4 done")

DefaultAssay(timmune.combined)<-"integrated"

# Run the standard workflow for visualization and clustering
timmune.combined <- ScaleData(timmune.combined, verbose = FALSE)
timmune.combined <- RunPCA(timmune.combined, npcs = 10, verbose = FALSE)
timmune.combined <- RunUMAP(timmune.combined, reduction = "pca", dims = 1:10)
timmune.combined <- FindNeighbors(timmune.combined, reduction = "pca", dims = 1:10)
timmune.combined <- FindClusters(timmune.combined, resolution = c(0.5))

DimPlot(timmune.combined,group.by=c("seurat_clusters"),raster=TRUE,label=TRUE)&coord_equal()

VlnPlot(timmune.combined,features=c("percent_mito"),group.by=c("seurat_clusters"))
VlnPlot(timmune.combined,features=c("nCount_RNA"),group.by=c("seurat_clusters"))
VlnPlot(timmune.combined,features=c("nFeature_RNA"),pt.size=0,group.by=c("seurat_clusters"))

FeaturePlot(timmune.combined,features="FOXP3",cols=c("lightblue","magenta4"),order=TRUE,pt.size=2)&coord_equal()


Tgenes<-c("TIGIT","CTLA4","FOXP3", #CD4Treg
         "PDCD1","LAG3", #exhausted
         "EOMES","HAVCR2",
         "GZMK","SELL","IL7R", #memory
          "GZMB","KLRD1","GZMH",
          "NCAM1", #CD56=NCAM1 NK
          "FCGR3A", #CD16=FCGR3A high NK
          "ITG2B",
          "CCR7","FOXO1","KLF2","LEF1","TCF7","ACTN1","FOXP1", #Naive
         "CD4","CD8A","CD8B")
DotPlot(timmune.combined,group.by=c("integrated_snn_res.0.5"),features=Tgenes)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(0,12))

  table(timmune.combined@meta.data$sample,timmune.combined@meta.data$seurat_clusters)

  ##Effector genes
goi<-c("FAS","FASLG","CD44","CD69","CD38","NKG7","KLRB1","KLRD1","KLRF1","KLRG1","KLRK1","FCGR3A","CX3CR1","CD300A","FGFBP2","ID2","ID3","PRDM1","RUNX3","TBX21","ZEB2","BATF","IRF4","NR4A1","NR4A2","NR4A3","PBX3","ZNF683","HOPX","FOS","FOSB","JUN","JUNB","JUND","STAT1","STAT2","STAT5A","STAT6","STAT4","EOMES")

DotPlot(timmune.combined,group.by=c("seurat_clusters"),features=goi)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(0,12))

  ###EXHAUSTION GENE SET
goi<-c("PDCD1","LAYN","HAVCR2","LAG3","CD244","CTLA4","LILRB1","TIGIT","TOX","VSIR","BTLA","ENTPD1","CD160","LAIR1")

DotPlot(tcell_nodonor,group.by=c("sample_summary"),features=goi)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(0,12))

Idents(timmune.combined)<-timmune.combined@meta.data$integrated_snn_res.0.5
clusterDEGs <- FindAllMarkers(timmune.combined, log2FC.threshold = 0.1, only.pos = TRUE, assay = "integrated",max.cells.per.ident=1000)
DEGs_sig_all_samples = subset(clusterDEGs, p_val_adj < 0.05 )

#Remove clusters with high mito content and low nFeatureRNA
Idents(timmune.combined)<-timmune.combined@meta.data$seurat_clusters
subset_combined<-subset(timmune.combined,idents=c("5","6"),invert=TRUE)

#Remove sample with less than 50 cells
Idents(subset_combined)<-subset_combined@meta.data$sample
subset_combined2<-subset(subset_combined,idents=c("Patient_29_Day_60-20"),invert=TRUE)

ifnb.list <- SplitObject(subset_combined2, split.by = "sample")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	DefaultAssay(x)<-"RNA3"
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

print("step 1 done")

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    DefaultAssay(x)<-"RNA3"
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

print("step 2 done")

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list,anchor.features = features, reduction = "rpca")
print("step 3 done")

all_genes <- lapply(ifnb.list, row.names) %>% Reduce(intersect, .) 

timmune.combined2 <- IntegrateData(anchorset = immune.anchors,k.weight = 40)
DefaultAssay(timmune.combined2) <- "integrated"
print("step 4 done")

# Run the standard workflow for visualization and clustering
timmune.combined2 <- ScaleData(timmune.combined2, verbose = FALSE)
timmune.combined2 <- RunPCA(timmune.combined2, npcs = 10, verbose = FALSE)
timmune.combined2 <- RunUMAP(timmune.combined2, reduction = "pca", dims = 1:10)
timmune.combined2 <- FindNeighbors(timmune.combined2, reduction = "pca", dims = 1:10)
timmune.combined2 <- FindClusters(timmune.combined2, resolution = 0.5)

saveRDS(timmune.combined,"tcell_with_highMT_hinFeature_022024.rds")
saveRDS(timmune.combined2,"tcell_subset_022024.rds")

DimPlot(timmune.combined2,group.by="seurat_clusters",label=TRUE)&coord_equal()

FeaturePlot(timmune.combined2,features=c("CD4","IL7R","CD8A","CD8B"),order=TRUE)&coord_equal()

Tgenes<-c("TIGIT","CTLA4","FOXP3", #CD4Treg
         "PDCD1","LAG3", #exhausted
         "EOMES","HAVCR2",
         "GZMK","SELL","IL7R", #memory
          "GZMB","KLRD1","GZMH",
          "NCAM1", #CD56=NCAM1 NK
          "FCGR3A", #CD16=FCGR3A high NK
          "ITG2B",
          "CCR7","FOXO1","KLF2","LEF1","TCF7","ACTN1","FOXP1", #Naive
         "CD4","CD8A","CD8B")
DotPlot(timmune.combined2,group.by=c("seurat_clusters"),features=Tgenes)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(0,12))

  goi<-c("GZMA","GZMB","GZMH","GZMK","GNLY","PRF1",
       "IFNG","TNF","SERPINB1","SERPINB6","SERPINB9","CTSA",
       "CTSB","CTSC","CTSD","CTSW","CST3","CST7","CSTB","LAMP1","LAMP3","CAPN2")

DotPlot(timmune.combined2,group.by=c("seurat_clusters"),features=goi)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(0,12))

  Idents(timmune.combined2)<-timmune.combined2@meta.data$integrated_snn_res.0.5
clusterDEGs <- FindAllMarkers(timmune.combined2, log2FC.threshold = 0.1, only.pos = TRUE, assay = "integrated",max.cells.per.ident=1000)
DEGs_sig_all_samples = subset(clusterDEGs, p_val_adj < 0.05 )

Idents(timmune.combined2)<-timmune.combined2@meta.data$seurat_clusters

timmune.combined2<-RenameIdents(timmune.combined2,
  '0'='NK',
  '1'='CD4',
  '2'='CD4',
  '3'='NK',
  '4'='NK',
  '5'='CD8',
  '6'='CD8',
  '7'='CD4',
  '8'='NK',
  '9'='MonocyteMacrophage',
  '10'='NK',
  '11'='Megakaryocytes')

timmune.combined2@meta.data$general_T_NK_annotation<-Idents(timmune.combined2)

###Subset T cell and redo anchoring
Tcellonly<-subset(timmune.combined2,idents=c("CD8","CD4"))
table(Tcellonly@meta.data$sample)

##Samples with less than 40 T cells removed
Idents(Tcellonly)<-Tcellonly@meta.data$sample
Tcellonly_2<-subset(Tcellonly,idents=c("Patient_27_Day_28-16","Patient-026-Day-28-10","Patient_11_Day_28-7","Patient_11_Day_60-8","Patient_38_Day_28-30","Patient_27_Day_60-17","Patient_31_Day_28-22","Patient_24_Day_28-13"),invert=TRUE)

ifnb.list <- SplitObject(Tcellonly_2, split.by = "sample")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	DefaultAssay(x)<-"RNA3"
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

print("step 1 done")

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    DefaultAssay(x)<-"RNA3"
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

print("step 2 done")

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list,anchor.features = features, k.filter = 200, reduction = "rpca")
print("step 3 done")

tcell_3 <- IntegrateData(anchorset = immune.anchors,k.weight = 40)

DefaultAssay(tcell_3) <- "integrated"
print("step 4 done")

# Run the standard workflow for visualization and clustering
tcell_3 <- ScaleData(tcell_3, verbose = FALSE)
tcell_3 <- RunPCA(tcell_3, npcs = 10, verbose = FALSE)
tcell_3 <- RunUMAP(tcell_3, reduction = "pca", dims = 1:10)
tcell_3 <- FindNeighbors(tcell_3, reduction = "pca", dims = 1:10)
tcell_3 <- FindClusters(tcell_3, resolution = c(0.5))

library(clustree)

###Add additional resolutions for resolving some clustering differences between naive and central memory T cells.

DefaultAssay(tcell_3)<-"integrated"
resolution.range <- seq(from = 0, to = 1, by = 0.1)
tcell_3 <- FindClusters(tcell_3, resolution = resolution.range)

clustree(tcell_3, prefix = "integrated_snn_res.")

DimPlot(tcell_3,group.by=c("integrated_snn_res.1"),label=TRUE)&coord_equal()
DimPlot(tcell_3,group.by=c("predicted.celltype.l2"),label=TRUE)&coord_equal()

table(tcell_3@meta.data$predicted.celltype.l2,tcell_3@meta.data$integrated_snn_res.1)

saveRDS(tcell_3,"tcell_subset_cd4cd8_032024.rds")

T_Naive<-list(c("IL7R","CCR7","SELL","FOXO1","KLF2","KLF3","LEF1","TCF7","ACTN1","FOXP1"))
Active_Effector<-list(c("FAS","FASLG","CD44","CD69","CD38","NKG7","KLRB1","KLRD1","KLRF1","KLRG1","KLRK1","FCGR3A","CX3CR1","CD300A","FGFBP2","ID2","ID3","PRDM1","RUNX3","TBX21","ZEB2","BATF","IRF4","NR4A1","NR4A2","NR4A3","PBX3","ZNF683","HOPX","FOS","FOSB","JUN","JUNB","JUND","STAT1","STAT2","STAT5A","STAT6","STAT4","EOMES"))
Exhaustion<-list(c("PDCD1","LAYN","HAVCR2","LAG3","CD244","CTLA4","LILRB1","TIGIT","TOX","VSIR","BTLA","ENTPD1","CD160","LAIR1"))
TCR_Signaling<-list(c("CALM1","CALM2","CALM3","CAST","CD247","CD3D","CD3E","CD3G","CSK","DOK1","DOK2","FYN","LCK","NFATC2","NFATC1","NFATC4","NFATC3","PLEK","PAG1","PTPN11","PTPN2","PTPN22","PTPN4","PTPN6","PTPN7","PTPRC","PTPRCAP","S100A10","S100A11","S100A13","S100A4","S100A6","ZAP70","DUSP1","DUSP2","DUSP4","DUSP5","DUSP10","LAT","PLCG1","PLCG2","PPP3CA","PPP3CC","FOS","FOSB","FOSL1","FOSL2","JUN","JUNB","JUND","NR4A1","NR4A2","NR4A3","BATF","IRF4","SH2D2A"))
cytotox<-list(c("GZMA","GZMB","GZMH","GZMK","GNLY","PRF1","IFNG","TNF","SERPINB1","SERPINB6","SERPINB9","CTSA","CTSB","CTSC","CTSD","CTSW","CST3","CST7","CSTB","LAMP1","LAMP3","CAPN2"))
cytokine<-list(c("CSF1","IL10RA","IL16","IL17RA","IL18RAP","IL21R","IL2RB","IL2RG","IL32","IL9R","ADAM10","ADAM8","METRNL","CD70"))
chemokine<-list(c("CCR4","CCR5","CCR7","CXCR3","CXCR4","CXCR5","CXCR6","CCL3","CCL4","CCL4L1","CCL4L2","CCL5","CXCL13","CXCL8","XCL1","XCL2"))
senesence<-list(c("KLRK1","IFNA1","IFNAR1","TAB1","SESN2"))
mapk<-list(c("MAP2K1","MAP2K2","MAP2K3","MAP3K4","MAP3K5","MAP3K8","MAP4K1","MAP4K5","MAPK1","MAPK3","MAPK11","MAPK13","MAPK14","MAPK1IP1L","MAPKAPK3","MAPK8","MAPK9","MAP2K7","MAPK10","MAP3K7"))
ifn_response<-list(c("IFIT1","IFIT2","IFIT3","IFIT5","STAT1","STAT2","MX1","IRF1","IRF4","IRF7","IRF8","IRF9","ISG15","ISG20","IFITM1","IFITM2","IFITM3","OAS1","OAS2","OAS3","JAK1","JAK2","SOCS1","SOCS3","TRIM14","TRIM21","TRIM22","APOL1","APOL2","APOL6","IFNGR1","GBP1","GBP2","GBP4","GBP5","GBP3","BST2","CMPK2","DDX58","DDX60","DDX60L","IFI30","IFI35","IFI44","IFI44L","IFI6","IFIH1","PARP10","PARP12","PARP14"))
antiapoptosis<-list(c("BCL11B","BCL2L1","BAG3","BAG4","BIRC3","GADD45B"))
pro_apoptosis<-list(c("BCL2L11","BID","BIN1","BIN2","BAX","BAK1","CASP8","CASP3"))
glycolosis<-list(c("SLC2A3","SLC2A8","PFKFB3","ALDOA","ENO1","GAPDH","GPI","PGAM1","PGK1","PKM","TPI1","GYG1","MDH1","MDH2","LDHA","LDHB","ADH5","IDH2"))

tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = T_Naive,name = 'T_Naive')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = Active_Effector,name = 'Active_Effector')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = Exhaustion,name = 'Exhaustion')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = TCR_Signaling,name = 'TCR_Signaling')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = cytotox,name = 'cytotox')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = cytokine,name = 'cytokine')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = chemokine,name = 'chemokine')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = senesence,name = 'senesence')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = mapk,name = 'mapk')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = ifn_response,name = 'ifn_response')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = antiapoptosis,name = 'antiapoptosis')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = pro_apoptosis,name = 'pro_apoptosis')
tcell_3 <- AddModuleScore(object = tcell_3,assay="SCT",features = glycolosis,name = 'glycolosis')

FeaturePlot(tcell_3,features="T_Naive1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="Active_Effector1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="Exhaustion1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="TCR_Signaling1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="cytotox1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="cytokine1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="senesence1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="mapk1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="ifn_response1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="antiapoptosis1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="pro_apoptosis1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()
FeaturePlot(tcell_3,features="glycolosis1",order=TRUE)&scale_colour_viridis(option="H")&coord_equal()

tcell_3<-PrepSCTFindMarkers(tcell_3)


DefaultAssay(tcell_3)<-"SCT"

Idents(tcell_3)<-tcell_3@meta.data$`integrated_snn_res.0.5`

tcell_3<-RenameIdents(tcell_3,
  '0'='CD4_Central_Memory',
  '1'='CD4_Central_Memory',
  '2'='CD4_Naive_Central_Memory_Mixed',
  '3'='CD8_Naive',
  '4'='CD8_Effector_Memory',
  '5'='CD4_Central_Memory',
  '6'='CD4_Treg',
  '7'='CD4_Central_Memory', #STAT1 expression
  '8'='CD8_Effector_Memory') #HAs GZMK expression

tcell_3@meta.data$cell_type_specific<-Idents(tcell_3)

tcell_3@meta.data$cell_type_specific <- factor(tcell_3@meta.data$cell_type_specific,levels = c("CD8_Effector_Memory","CD4_Treg","CD4_Central_Memory","CD4_Naive_Central_Memory_Mixed","CD8_Naive"))


Tgenes<-c("ZEB2","ZNF683", "CCL5", "GZMH", "KLRD1", "NKG7","GZMK", "GZMB",#cd8 effector memory"
                      "TIGIT","CTLA4","FOXP3", #CD4Treg
          "IL7R","TMSB10", "ITGB1", "LTB", "TRAC", "AQP3", #cd4 central mem
          "CCR7","FOXO1","LEF1","TCF7","ACTN1","FOXP1", #Naive
         "CD4","CD8A","CD8B")
DotPlot(tcell_3,assay="SCT",group.by=c("cell_type_specific"),features=Tgenes)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(0,12))+ theme(legend.position="bottom")

a<-DotPlot(tcell_3,assay="SCT",group.by=c("cell_type_specific"),features=Tgenes)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(0,12))
b<-DimPlot(tcell_3,group.by="cell_type_specific",label=TRUE)&coord_equal()& theme(legend.position="none")

a|b

