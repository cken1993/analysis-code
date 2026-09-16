options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

rm(list = ls())
suppressPackageStartupMessages(library(dplyr));library(data.table);library(magrittr);library(tibble);library(stringr)

library(GSVA)
library(ComplexHeatmap)
library(clusterProfiler)


load('project_TLS/1.for_bulk/TLS_ornot_bulk_data2.0.rdata')
Tpm_TLS_ornot <- as.data.frame(apply(Tpm_TLS_ornot[,setdiff(colnames(Tpm_TLS_ornot), "gene_id")], 2, function(x) tapply(x, INDEX=factor(Tpm_TLS_ornot$gene), FUN=median, na.rm=TRUE))) 

Tpm_TLS_ornot = log2(Tpm_TLS_ornot+1)




library(org.Hs.eg.db)
gmt = list()
cluster_index = top10_gmt$cluster %>% unique()

top10_gmt = read.table('project_TLS/1.for_bulk/3.GSVA/markers_for_22_fine_celltypes.txt',
                header = T)
for(i in seq_len(length(cluster_index))){
    gmt[[cluster_index[i]]] = top10_gmt %>% filter(cluster == cluster_index[i]) %>% .[,'gene']
}

# 备选gmt
if(T){
    load('project_TLS/1.for_bulk/3.GSVA/easy_input_immunity.rdata')
    str(immunity)
    for(i in seq_len(length(immunity))){
        gmt[[names(immunity[i])]] = bitr(immunity[[i]],fromType = 'ENTREZID',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db') %>% .[,2] %>% unlist()
    }
}

gmt2 = list()
gmt2[['9-gene signature']] = c('CD79B', 'CD1D', 'CCR6', 'LAT', 'SKAP1', 'CETP', 'EIF1AY', 'RBP5', 'PTGDS')
gmt2[['12-gene signature']] = c('CCL2', 'CCL3', 'CCL4', 'CCL5', 'CCL8', 'CCL18', 'CCL19', 'CCL21', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL13')
gmt2[['M1_sig']] = c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','IDO1','KYNU','CCR7')
gmt2[['M2_sig']] = c('IL4R','CCL4','CCL13','CCL20','CCL17','CCL18','CCL22','CCL24','LYVE1','VEGFA','VEGFB','VEGFC','VEGFD','EGF','CTSA','CTSB','CTSC','CTSD','TGFB1','TGFB2','TGFB3','MMP14','MMP19','MMP9','CLEC7A','WNT7B','FASL','TNFSF12','TNFSF8','CD276','VTCN1','MSR1','FN1','IRF4')

tcga_gsva1 <- as.data.frame(t(gsva(as.matrix(Tpm_TLS_ornot), gmt, method = "ssgsea")))
tcga_gsva2 <- as.data.frame(t(gsva(as.matrix(Tpm_TLS_ornot), gmt2, method = "ssgsea")))

tcga_gsva = cbind.data.frame(tcga_gsva1,tcga_gsva2)

#compare
outTab0 <- NULL
tcga_gsva$group = Clinical_TLS_ornot[match(rownames(tcga_gsva),Clinical_TLS_ornot$sample),'group.y']


for (i in setdiff(colnames(tcga_gsva), "group")) {
  kt <- kruskal.test(as.numeric(tcga_gsva[,i]) ~ tcga_gsva$group)
  outTab0 <- rbind.data.frame(outTab0,
                             data.frame(feature = i,
                                        p = kt$p.value,
                                        mean_TLS = median(tcga_gsva[tcga_gsva$group=='TLS',i]),
                                        mean_No = median(tcga_gsva[tcga_gsva$group=='No',i]),
                                        dif = median(tcga_gsva[tcga_gsva$group=='TLS',i]) - median(tcga_gsva[tcga_gsva$group=='No',i]),
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}
outTab0$fdr <- p.adjust(outTab0$p, method = "fdr") # 矫正p值
outTab0$txt <- formatC(outTab0$fdr, digits = 2, format = "e") # 将FDR改为科学计数法
outTab0$txt <- gsub("e"," × 10", outTab0$txt) # 修改文本


#write.csv(outTab0,'project_TLS/1.for_bulk/3.GSVA/outTab.csv')
#write.csv(tcga_gsva,'project_TLS/1.for_bulk/3.GSVA/tcga_gsva.csv')



outTab1 <- NULL
id = "12-gene signature"
for (i in setdiff(colnames(tcga_gsva), "group")) {
  r <- cor.test(as.numeric(tcga_gsva[,id]),as.numeric(tcga_gsva[,i]))
  outTab1 <- rbind.data.frame(outTab1,
                             data.frame(feature = i,
                                        r = as.numeric(r$estimate),
                                        p = as.numeric(r$p.value),
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}
outTab1$fdr <- p.adjust(outTab1$p, method = "fdr") # 矫正p值
outTab1$txt <- formatC(outTab1$fdr, digits = 2, format = "e") # 将FDR改为科学计数法
outTab1$txt <- gsub("e"," × 10", outTab1$txt) # 修改文本
outTab1$main = id


outTab2 = NULL
id = "9-gene signature"
for (i in setdiff(colnames(tcga_gsva), "group")) {
  r <- cor.test(as.numeric(tcga_gsva[,id]),as.numeric(tcga_gsva[,i]))
  outTab2 <- rbind.data.frame(outTab2,
                             data.frame(feature = i,
                                        r = as.numeric(r$estimate),
                                        p = as.numeric(r$p.value),
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}
outTab2$fdr <- p.adjust(outTab2$p, method = "fdr") # 矫正p值
outTab2$txt <- formatC(outTab2$fdr, digits = 2, format = "e") # 将FDR改为科学计数法
outTab2$txt <- gsub("e"," × 10", outTab2$txt) # 修改文本
outTab2$main = id



# -------------------------
for_cox = Clinical_TLS_ornot
tcga_gsva = tcga_gsva[match(Clinical_TLS_ornot$sample,
                        rownames(tcga_gsva)),]
for_cox$"9-gene signature" = tcga_gsva$"9-gene signature"
for_cox$"12-gene signature" = tcga_gsva$"12-gene signature"
write.csv(for_cox,'for_cox.csv')
save.image(file = 'gsva.rdata')
