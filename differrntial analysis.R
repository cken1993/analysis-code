options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

rm(list = ls())
suppressPackageStartupMessages(library(dplyr));library(data.table);library(magrittr);library(tibble);library(stringr)

load('project_TLS/1.for_bulk/TLS_ornot_bulk_data2.0.rdata')


Count_TLS_ornot <- as.data.frame(apply(Count_TLS_ornot[,setdiff(colnames(Count_TLS_ornot), "gene_id")], 2, function(x) tapply(x, INDEX=factor(Count_TLS_ornot$gene_id), FUN=median, na.rm=TRUE))) 
Count_TLS_ornot = round(Count_TLS_ornot)

library(DESeq2)

index = intersect(colnames(Count_TLS_ornot),Clinical_TLS_ornot$sample)
Count_TLS_ornot = Count_TLS_ornot[,index]
Clinical_TLS_ornot = Clinical_TLS_ornot[match(index,Clinical_TLS_ornot$sample),]
identical(colnames(Count_TLS_ornot),Clinical_TLS_ornot$sample) #TRUE

condition = Clinical_TLS_ornot$group.y %>% as.factor()
coldata <- data.frame(row.names = colnames(Count_TLS_ornot), condition)
dds <- DESeqDataSetFromMatrix(countData = Count_TLS_ornot,
                                  colData = coldata,
                                  design = ~condition)

dds$condition = relevel(dds$condition,ref = "No")
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","TLS","No"))



resData <- as.data.frame(res[order(res$padj),])
resData$id <- rownames(resData)
resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
colnames(resData) <- c("id","baseMean","log2FC","lfcSE","stat","PValue","FDR")

write.csv(resData, 'project_TLS/1.for_bulk/2.DEG/deg.csv',row.names = F)

geneup <- resData %>%
    filter(log2FC > 1 & FDR < 0.05)
genedown <- resData %>%
    filter(log2FC < (-1) & FDR < 0.05)
write.csv(geneup, 'project_TLS/1.for_bulk/2.DEG/geneup.csv',row.names = F)
write.csv(genedown, 'project_TLS/1.for_bulk/2.DEG/genedown.csv',row.names = F)

resData %>% arrange(desc(log2FC)) %>% select(id,log2FC) ->gsea_df
rownames(gsea_df) = seq_len(nrow(gsea_df))
write.csv(gsea_df, 'project_TLS/1.for_bulk/2.DEG/gsea_df.csv',row.names = F)
