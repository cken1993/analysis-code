## 单细胞处理

library(Seurat)
library(SeuratObject)
library(future)
library(future.apply)
library(ggplot2)
library(tidyverse)
library(data.table)

plan(multisession, workers = 16)
options(future.globals.maxSize = 160000 * 1024^2)

### 读取数据和表型
{
  file_path = './data/3_GSE149614/GSE149614_count.txt.gz'
  scdata = fread(file_path,sep = '\t')
  scdata <- setDF(scdata) ;gc()
  
  rownames(scdata) <- scdata[,1]
  scdata <- scdata[,-1]
  
  scobj = CreateSeuratObject(counts = scdata, ### 创建seurat对象，分配给scobj
                             projtec = "scobj_3_GSE149614",
                             min.cells = 3,
                             min.features = 200)
  
  pheno = data.table::fread('./data/3_GSE149614/GSE149614_metadata.txt.gz')
  
  ## cellid加上数据集的名字
  identical(colnames(scobj),pheno$Cell) # T
  scobj = RenameCells(scobj,add.cell.id = "GSE149614")
  
  ### metadata中加入pheno
  scobj@meta.data = cbind(scobj@meta.data,pheno)
  
}

### 细胞分群
{
  scobj <- NormalizeData(scobj)
  
  scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 3000)
  scobj <- ScaleData(scobj, features = rownames(scobj))
  
  scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj),reduction.name = "pca")
  scobj <- RunUMAP(scobj,reduction = "pca", dims = 1:30, reduction.name = "umap")
  scobj <- RunTSNE(scobj,reduction = "pca", dims = 1:30, reduction.name = "tsne")
  scobj <- FindNeighbors(scobj, dims = 1:30)
  scobj <- FindClusters(scobj, resolution = 1)
}

{
  scobj$new_major_celltype = scobj$res.3
  
  scobj$new_major_celltype <- recode(scobj$new_major_celltype,
                                     '1'='T/NK cells',
                                     '2'='T/NK cells',
                                     '3'='Hepatocyte/bipotent cells',
                                     '4'='Hepatocyte/bipotent cells',
                                     '5'='Myeloid cells',
                                     '6'='Myeloid cells',
                                     '7'='Endothelial cells',
                                     '8'='T/NK cells',
                                     '9'='T/NK cells',
                                     '10'='Myeloid cells',
                                     '11'='T/NK cells',
                                     '12'='Hepatocyte/bipotent cells',
                                     '13'='T/NK cells',
                                     '14'='T/NK cells',
                                     '15'='Hepatocyte/bipotent cells',
                                     '16'='Myeloid cells',
                                     '17'='Hepatocyte/bipotent cells',
                                     '18'='T/NK cells',
                                     '19'='Hepatocyte/bipotent cells',
                                     '20'='T/NK cells',
                                     '21'='Myeloid cells',
                                     '22'='Hepatocyte/bipotent cells',
                                     '23'='Myeloid cells',
                                     '24'='Hepatocyte/bipotent cells',
                                     '25'='T/NK cells',
                                     '26'='Myeloid cells',
                                     '27'='Hepatocyte/bipotent cells',
                                     '28'='T/NK cells',
                                     '29'='Hepatocyte/bipotent cells',
                                     '30'='Endothelial cells',
                                     '31'='Fibroblasts',
                                     '32'='T/NK cells',
                                     '33'='Fibroblasts',
                                     '34'='B/Plasma cells',
                                     '35'='T/NK cells',
                                     '36'='B/Plasma cells',
                                     '37'='B/Plasma cells',
                                     '38'='Myeloid cells',
                                     '39'='Myeloid cells',
                                     '40'='B/Plasma cells',
                                     '41'='Myeloid cells',
                                     '42'='Hepatocyte/bipotent cells',
                                     '43'='Hepatocyte/bipotent cells',
                                     '44'='Myeloid cells',
                                     '45'='Hepatocyte/bipotent cells',
                                     '46'='Myeloid cells',
                                     '47'='Hepatocyte/bipotent cells',
                                     '48'='Endothelial cells',
                                     '49'='Hepatocyte/bipotent cells',
                                     '50'='B/Plasma cells',
                                     '51'='T/NK cells',
                                     '52'='Myeloid cells',
                                     '53'='Myeloid cells')
  
  scobj$new_major_celltype <- factor(scobj$new_major_celltype, levels = c("T/NK cells", "Myeloid cells", "B/Plasma cells", "Fibroblasts",
                                                                          "Endothelial cells", "Hepatocyte/bipotent cells")) #因子顺序
  
  ### 构造color ###
  color_vec = c('#7d9ac0','#ecab44','#59975d','#d47464','#9ea2a1','#133e63')
  
  
  ### 绘制tsne图
  Idents(scobj) = 'new_major_celltype'
  DimPlot(scobj, reduction = "umap", label = F,cols = color_vec)
}

{
  ### 构造marker -- 原文注释
  epithelial_marker_genes <- c('ALB','SERPINA1','HNF4A','EPCAM')
  B_marker_genes <- c('IGHG1','JCHAIN','CD79A')
  TNK_marker_genes <- c('CD3D','CD3E','NKG7')
  Myeloid_marker_genes <- c('CD68','CD14','CD163','CD1C','CLEC4C','KIT')
  Endothelial_marker_genes <- c('VWF','PECAM1','FCGR2B')
  Fibroblast_marker_genes <- c('ACTA2','COL1A1','COL1A2')
  
  marker_all_v1 = c(epithelial_marker_genes,
                    B_marker_genes,
                    TNK_marker_genes,
                    Myeloid_marker_genes,
                    Endothelial_marker_genes,
                    Fibroblast_marker_genes)
  
  DotPlot(scobj, features = marker_all_v1,group.by = "new_major_celltype",) + coord_flip()
}

saveRDS(scobj,'./resource/annoRDS/3_GSE149614_scobj.rds')