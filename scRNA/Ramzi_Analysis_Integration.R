
library(Seurat)
library(ggplot2)
library(reshape2)
library(viridis)
library(dplyr)
library(tidyr)
library(ggrepel)
library(tidyverse)
library(reticulate)
library(Azimuth)

#
set.seed(1234)
options(Seurat.object.assay.version = 'v5')


###Read in RDS files
s1=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_11_Donor-9/MGI3382_TWJG-Patient_11_Donor-9_processed.rds"))
s2=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_27_Donor-18/MGI3382_TWJG-Patient_27_Donor-18_processed.rds"))
s3=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_36_Day_28-27/MGI3382_TWJG-Patient_36_Day_28-27_processed.rds"))
s4=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Control_5_Day_28-1/MGI3382_TWJG-Control_5_Day_28-1_processed.rds"))
s5=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_17_Day_28-10/MGI3382_TWJG-Patient_17_Day_28-10_processed.rds"))
s6=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_29_Day_28-19/MGI3382_TWJG-Patient_29_Day_28-19_processed.rds"))
s7=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_36_Day_60-28/MGI3382_TWJG-Patient_36_Day_60-28_processed.rds"))
s8=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Control_5_Day_60-2/MGI3382_TWJG-Control_5_Day_60-2_processed.rds"))
s9=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_17_Day_60-11/MGI3382_TWJG-Patient_17_Day_60-11_processed.rds"))
s10=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_29_Day_60-20/MGI3382_TWJG-Patient_29_Day_60-20_processed.rds"))
s11=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_36_Donor-29/MGI3382_TWJG-Patient_36_Donor-29_processed.rds"))
s12=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Control_6_Day_28-3/MGI3382_TWJG-Control_6_Day_28-3_processed.rds"))
s13=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_17_Donor-12/MGI3382_TWJG-Patient_17_Donor-12_processed.rds"))
s14=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_29_Donor-21/MGI3382_TWJG-Patient_29_Donor-21_processed.rds"))
s15=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_38_Day_28-30/MGI3382_TWJG-Patient_38_Day_28-30_processed.rds"))
s16=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Control_6_Day_60-4/MGI3382_TWJG-Control_6_Day_60-4_processed.rds"))
s17=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_24_Day_28-13/MGI3382_TWJG-Patient_24_Day_28-13_processed.rds"))
s18=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_31_Day_28-22/MGI3382_TWJG-Patient_31_Day_28-22_processed.rds"))
s19=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_38_Day_60-31/MGI3382_TWJG-Patient_38_Day_60-31_processed.rds"))
s20=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Control_7_Day_28-5/MGI3382_TWJG-Control_7_Day_28-5_processed.rds"))
s21=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_24_Day_60-14/MGI3382_TWJG-Patient_24_Day_60-14_processed.rds"))
s22=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_31_Day_60-23/MGI3382_TWJG-Patient_31_Day_60-23_processed.rds"))
s23=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_38_Donor-32/MGI3382_TWJG-Patient_38_Donor-32_processed.rds"))
s24=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Control_7_Day_60-6/MGI3382_TWJG-Control_7_Day_60-6_processed.rds"))


s25=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_24_Donor-15/MGI3382_TWJG-Patient_24_Donor-15_processed.rds"))
s26=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_33_Day_28-24/MGI3382_TWJG-Patient_33_Day_28-24_processed.rds"))
s27=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_11_Day_28-7/MGI3382_TWJG-Patient_11_Day_28-7_processed.rds"))
s28=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_27_Day_28-16/MGI3382_TWJG-Patient_27_Day_28-16_processed.rds"))
s29=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_33_Day_60-25/MGI3382_TWJG-Patient_33_Day_60-25_processed.rds"))
s30=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_11_Day_60-8/MGI3382_TWJG-Patient_11_Day_60-8_processed.rds"))
s31=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_27_Day_60-17/MGI3382_TWJG-Patient_27_Day_60-17_processed.rds"))
s32=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI3382_TWJG-Patient_33_Donor-26/MGI3382_TWJG-Patient_33_Donor-26_processed.rds"))
s33=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-1683278a_14-3/MGI2413_TWJG-1683278a_14-3_processed.rds"))
s34=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-1692941a_10-6/MGI2413_TWJG-1692941a_10-6_processed.rds"))
s35=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-1703980a_5-9/MGI2413_TWJG-1703980a_5-9_processed.rds"))
s36=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-1713398a_7-12/MGI2413_TWJG-1713398a_7-12_processed.rds"))
s37=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Control-1-Day-28-13/MGI2413_TWJG-Control-1-Day-28-13_processed.rds"))
s38=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Control-1-Day-60-14/MGI2413_TWJG-Control-1-Day-60-14_processed.rds"))
s39=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Control-2-Day-28-15/MGI2413_TWJG-Control-2-Day-28-15_processed.rds"))
s40=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Control-2-Day-60-16/MGI2413_TWJG-Control-2-Day-60-16_processed.rds"))
s41=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Patient-009-Day-28-1/MGI2413_TWJG-Patient-009-Day-28-1_processed.rds"))
s42=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Patient-009-Day-60-2/MGI2413_TWJG-Patient-009-Day-60-2_processed.rds"))
s43=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Patient-016-Day-28-4/MGI2413_TWJG-Patient-016-Day-28-4_processed.rds"))
s44=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Patient-016-Day-60-5/MGI2413_TWJG-Patient-016-Day-60-5_processed.rds"))
s45=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Patient-022-Day-60-8/MGI2413_TWJG-Patient-022-Day-60-8_processed.rds"))
s46=UpdateSeuratObject(readRDS("/diskmnt/Datasets/mmy_scratch/Dipersio/ramzi/Analysis_012024/MGI2413_TWJG-Patient-026-Day-28-10/MGI2413_TWJG-Patient-026-Day-28-10_processed.rds"))


alldata = merge(s1,y=c(x,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46),project = "Merged")
#Join layers before splitting
DefaultAssay(alldata)<-"RNA"
#alldata <- JoinLayers(alldata)


#Split layers by dataset
alldata[["RNA"]] <- split(alldata[["RNA"]], f = alldata$orig.ident)

#run sctransform
alldata <- SCTransform(alldata, assay = "RNA",vars.to.regress=c("nCount_RNA","S.Score","G2M.Score","percent_mito","percent_hb","percent_ribo"), vst.flavor = "v2",verbose=TRUE)
alldata <- RunPCA(alldata, npcs = 20, verbose = TRUE)

# one-liner to run Integration
alldata <- IntegrateLayers(object = alldata, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = 'harmony',
                       assay = "SCT", verbose = FALSE)
alldata <- FindNeighbors(alldata, reduction = "harmony", dims = 1:20)
alldata <- FindClusters(alldata, resolution = 0.6, cluster.name = "harmony_clusters")
alldata <- RunUMAP(alldata, reduction = "harmony", dims = 1:10, reduction.name = "umap10.harmony") ##reduce to 15/10
alldata <- RunUMAP(alldata, reduction = "harmony", dims = 1:15, reduction.name = "umap15.harmony") ##reduce to 15/10
alldata <- RunUMAP(alldata, reduction = "harmony", dims = 1:20, reduction.name = "umap20.harmony") ##reduce to 15/10

# one-liner to run Integration
alldata <- IntegrateLayers(object = alldata, method = RPCAIntegration,
                       orig.reduction = "pca", new.reduction = 'integrated.rpca',
                       normalization.method = "SCT", verbose = TRUE)
alldata <- FindNeighbors(alldata, reduction = "integrated.rpca", dims = 1:20)
alldata <- FindClusters(alldata, resolution = 0.6, cluster.name = "rpca_clusters")
alldata <- RunUMAP(alldata, reduction = "integrated.rpca", dims = 1:10, reduction.name = "umap10.rpca") ##reduce to 15/10
alldata <- RunUMAP(alldata, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap15.rpca") ##reduce to 15/10
alldata <- RunUMAP(alldata, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap20.rpca") ##reduce to 15/10

#Percent mito needs to be saved as percent.mt for RunAzimuth
alldata@meta.data$percent.mt=paste0(alldata@meta.data$percent_mito)
#DefaultAssay(s1)<-"SCT"

#Fixing object and rerunning SCT to help fix IntegrateLayers errors
#https://github.com/satijalab/seurat/issues/7542#issuecomment-1631534467
DefaultAssay(alldata)<-"RNA"
alldata <- JoinLayers(alldata)

#Run reference annotation
alldata <- RunAzimuth(alldata,query.modality = "SCT",reference = "pbmcref")

#https://satijalab.org/seurat/reference/prepsctfindmarkers
alldata<-PrepSCTFindMarkers(alldata, assay = "SCT", verbose = TRUE)

DefaultAssay(alldata)<-"SCT"

saveRDS(alldata,"combinedobject_v3_azimuth.rds")
