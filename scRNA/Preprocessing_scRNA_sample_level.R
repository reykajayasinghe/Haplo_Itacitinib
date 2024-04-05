#!/usr/bin/env Rscript --vanilla

# load required libraries
library(optparse)
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(DoubletFinder)
library(ggpubr)


# create user options
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL,
              help="path to data folder (e.g. cellranger output's raw matrices folder)",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./",
              help="output folder path",
              metavar="character"),
  make_option(c("-s","--sample_id"),
              type="character",
              default="single_cell_study",
              help="Name of your sample",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Path to data is required (--input).n", call.=FALSE)
}

###########################################################################################
## Run doublet finder
###########################################################################################

#Reference: https://github.com/chris-mcginnis-ucsf/DoubletFinder
RunDoubletFinder = function(seu){

  min.pc = 10
  ## Run Doublet finding algorithm
  sweep.list = paramSweep(seu, PCs = 1:min.pc, num.cores = 2, sct = TRUE)
  sweep.stats = summarizeSweep(sweep.list,GT=FALSE)
  bcmvn = find.pK(sweep.stats)
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  ## Homotypic doublet proportion estimate
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(seu@meta.data)) 
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  # run DoubletFinder
  seu <- doubletFinder(seu = seu,PCs = 1:min.pc,pK = optimal.pk,nExp = nExp.poi.adj,sct = T)
  colnames(seu@meta.data)[grep("DF.classifications",names(head(seu)))]="doublet_finder"
  #If you want to remove doublets then uncomment next line
  #seu = subset(seu, doublet_finder %in% "Singlet")
  seu@meta.data = seu@meta.data[,colnames(seu@meta.data)[-grep("pANN_",colnames(seu@meta.data))]]
  seu = RunPCA(seu,verbose=F)
  return(seu)
}

###################################################################
## Run SCTransform Seurat pipeline
###################################################################

RunSeurat = function(dat,project_name,out_path){
  seu = CreateSeuratObject(counts=dat,project=project_name)
  #seu = CreateSeuratObject(counts=expression_matrix,project="test")
  ## Compute %mito, %ribo and %RBC in the data
  seu = PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
  seu = PercentageFeatureSet(seu, "^RP[SL]", col.name = "percent_ribo")
  seu = PercentageFeatureSet(seu, "^HB[^(P)]", col.name = "percent_hb")

# plot pre-filter metadata
  pdf(paste(out_path,"/QC_in_sample_",project_name, ".pdf", sep=""), width=6, height=6)
    plot<-VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,raster=TRUE)
    #plot<-VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_hb", "percent_ribo"), ncol = 5,raster=TRUE)
	print(plot)
  dev.off()

  ## Data filtering 
  ##if using raw feature matrix use below filtering first
  seu = subset(x = seu, subset = nFeature_RNA > 200 & 
    nFeature_RNA < 10000 & 
    nCount_RNA > 1000 & 
    nCount_RNA < 80000 & percent_mito<20)
 #if use filtered feature matrix below filtering alone should suffice
 # seu = subset(seu, subset = percent_mito<20 & percent_ribo < 5 & percent_hb < 1 )

  ## Data normalization and cell cycle scoring
  seu = NormalizeData(seu,verbose=F)
  seu = CellCycleScoring(seu,g2m.features=cc.genes$g2m.genes,s.features=cc.genes$s.genes,verbose=F)

  #Plot pre-filter metadata
   pdf(paste(out_path,"/Post_QC_in_sample_",project_name, ".pdf", sep=""), width=12, height=7)
    plot1<-VlnPlot(object = seu, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_hb", "percent_ribo"), ncol = 5,raster=TRUE)
	print(plot1)
  dev.off()

  seu = SCTransform(seu, method = "glmGamPoi", vars.to.regress=c("nCount_RNA","S.Score","G2M.Score","percent_mito","percent_hb","percent_ribo"), verbose = FALSE)
  seu = RunPCA(seu,verbose=F)
  seu = RunUMAP(seu,dims=1:20,verbose=F)
  seu = FindNeighbors(seu, reduction = "pca", dims = 1:20,verbose=F)
  seu = FindClusters(seu,resolution=c(0.6),random.seed=123,verbose=F)
  ## Run Doublet finder
  seu = RunDoubletFinder(seu)
  ## SCTransform
  
  return(seu)
}

# read in initial arguments
sample_id <- opt$sample_id
out_path <- opt$output
matrix_dir = opt$input

# make output dir if it doesn't exist
dir.create(out_path)

barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "/features.tsv.gz")
matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")
expression_matrix <- ReadMtx(
  mtx = matrix.path, features = features.path,
cells = barcode.path
)

#expression_matrix <- Read10X(data.dir = matrix_dir)

panc=RunSeurat(expression_matrix,sample_id,out_path)


 pdf(paste(out_path,"/DimPlot_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
 	plot2<-DimPlot(objec = panc,group.by=c("seurat_clusters"),label=TRUE)&coord_equal()
	 print(plot2)
	plot3<-DimPlot(objec=panc,group.by=c("doublet_finder"),label=FALSE)&coord_equal()
	print(plot3)
 dev.off()

#Save sample name in metadata
panc@meta.data$sample_id<-sample_id

# save object so far
saveRDS(panc,file = paste(out_path,"/",sample_id, "_processed.rds", sep=""))
print("sample rds saved")

