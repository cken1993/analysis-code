library(Seurat)
library(rhdf5)
library(Matrix)
library(Seurat)
library(tidyverse)


{
  ### 构造color ###
  color_vec = c('#7d9ac0','#ecab44','#59975d','#d47464','#9ea2a1','#133e63')
  
  ### 构造marker 
  Hepatocytes = c('HP', 'ALDOP', 'KRT8', 'KRT18', 'EPCAM'); names(Hepatocytes) = rep('Hepatocytes',length(Hepatocytes))
  B = c('CD79A', 'CD19', 'MS4A1', 'MZB1', 'IGHA1', 'IGHG3'); names(B) = rep('B/Plasma',length(B))
  TNK = c('NKG7', 'KLRF1', 'CD3D', 'IL7R', 'CD8B', 'ICOS'); names(TNK) = rep('TNK',length(TNK))
  Myeloid = c('S100A8', 'S100A9', 'CD68', 'VCAN', 'C1QC', 'C1QA', 'TPSB2', 'LILRA4'); names(Myeloid) = rep('Myeloid',length(Myeloid))
  Endothelial = c('FCN3', 'VWF', 'GNG11', 'CRHBP'); names(Endothelial) = rep('Endothelial',length(Endothelial))
  Fibroblast = c('ACTA2', 'RGS5', 'COL1A1', 'TAGLN'); names(Fibroblast) = rep('Fibroblast',length(Fibroblast))
  
  marker_all_v2 = c(Hepatocytes,
                    Endothelial,
                    Fibroblast,
                    B,
                    Myeloid,
                    TNK
  )
}

seu = readRDS('GSE149614_scobj.rds')


DimPlot(seu, reduction = "umap", label = F,cols = color_vec)

celltype_column <- "new_major_celltype"

seu <- subset(
  seu,
  cells = colnames(seu)[
    !is.na(seu@meta.data[[celltype_column]])
  ]
)

if ("Assay5" %in% class(seu[["RNA"]])) {
  
  data.input <- GetAssayData(
    seu,
    assay = "RNA",
    layer = "data"
  )
  
} else {
  
  ## Seurat v4
  
  data.input <- GetAssayData(
    seu,
    assay = "RNA",
    slot = "data"
  )
}

# ============================================================
# 7. 创建 CellChat metadata
# ============================================================

meta <- data.frame(
  celltype = seu@meta.data[[celltype_column]],
  row.names = colnames(seu)
)

# ============================================================
# 8. 创建 CellChat 对象
# ============================================================

## 注意：
##
## 这里的 group.by = "celltype"
##
## 不是 responder / non-responder 分组。
##
## 它代表：
## CellChat 按“细胞类型”计算细胞之间的通讯。

cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "celltype"
)

{
  # ============================================================
  # 9. 设置物种数据库
  # ============================================================
  
  ## 如果你的数据是人 HCC：
  
  CellChatDB <- CellChatDB.human
  
  ## 如果是小鼠，请改为：
  ##
  # CellChatDB <- CellChatDB.mouse
  
  
  # 查看数据库
  
  showDatabaseCategory(CellChatDB)
  
  
  # 查看数据库中的 interaction
  
  head(CellChatDB$interaction)
  
  colnames(CellChatDB$interaction)
}

{
  # ============================================================
  # 10. 设置 CellChat 数据库
  # ============================================================
  
  cellchat@DB <- CellChatDB
}

{
  # ============================================================
  # 11. 首先分析 Secreted Signaling
  # ============================================================
  
  ## 趋化因子 CCL / CXCL 主要属于分泌型信号。
  
  CellChatDB.use <- subsetDB(
    CellChatDB,
    search = "Secreted Signaling"
  )
  
  cellchat@DB <- CellChatDB.use
}

{
  # ============================================================
  # 12. 查看数据库中的趋化因子通路
  # ============================================================
  
  ## 查看 pathway 名称
  
  unique(
    cellchat@DB$interaction$pathway_name
  )
  
  
  # 查看 CCL / CXCL
  
  chemokine_DB <- cellchat@DB$interaction %>%
    filter(
      pathway_name %in% c("CCL", "CXCL")
    )
  
  head(chemokine_DB)
  
  dim(chemokine_DB)
  
}

{
  # ============================================================
  # 13. 检查 CCL19 / CCR7 是否存在
  # ============================================================
  
  chemokine_DB %>%
    filter(
      grepl("CCL19|CCR7", interaction_name, ignore.case = TRUE)
    )
  
}

{
  # ============================================================
  # 14. CellChat 标准分析流程
  # ============================================================
  
  ## Step 1
  ## 提取数据库中存在的基因
  
  cellchat <- subsetData(cellchat)
  
  
  ## Step 2
  ## 寻找过表达基因
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  
  ## Step 3
  ## 寻找过表达 ligand-receptor pairs
  
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  
  ## Step 4
  ## 计算通讯概率
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean"
  )
  
  
  ## Step 5
  ## 过滤细胞数量过少的细胞群
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = 1
  )
  
  
  ## Step 6
  ## 计算 signaling pathway 层面的通讯
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  
  ## Step 7
  ## 汇总网络
  
  cellchat <- aggregateNet(cellchat)
}

# -----------------------------------
{
  # ============================================================
  # 15. 查看总体细胞通讯网络
  # ============================================================
  
  groupSize <- as.numeric(
    table(cellchat@idents)
  )
  
  
  # ------------------------------------------------------------
  # 15.1 通讯数量网络
  # ------------------------------------------------------------
  
  pdf(
    "01_CellChat_total_number.pdf",
    width = 10,
    height = 10
  )
  
  netVisual_circle(
    cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of cell-cell interactions"
  )
  
  dev.off()
  
  
  # ------------------------------------------------------------
  # 15.2 通讯强度网络
  # ------------------------------------------------------------
  
  pdf(
    "02_CellChat_total_weight.pdf",
    width = 10,
    height = 10
  )
  
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction strength"
  )
  
  dev.off()
  
}


# ============================================================
# 15. 查看总体细胞通讯网络
# ============================================================

groupSize <- as.numeric(
  table(cellchat@idents)
)


# ------------------------------------------------------------
# 15.1 通讯数量网络
# ------------------------------------------------------------

pdf(
  "01_CellChat_total_number.pdf",
  width = 10,
  height = 10
)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of cell-cell interactions"
)

dev.off()


# ------------------------------------------------------------
# 15.2 通讯强度网络
# ------------------------------------------------------------

pdf(
  "02_CellChat_total_weight.pdf",
  width = 10,
  height = 10
)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction strength"
)

dev.off()


# ============================================================
# 16. 提取全部 ligand-receptor interaction
# ============================================================

all.communication <- subsetCommunication(
  cellchat
)

head(all.communication)

colnames(all.communication)

write.xlsx(
  all.communication,
  file = "03_All_CellChat_communication.xlsx",
  rowNames = FALSE
)


# ============================================================
# 17. 提取 CCL / CXCL 趋化因子通讯
# ============================================================

chemokine.communication <- subsetCommunication(
  cellchat,
  signaling = c("CCL", "CXCL")
)


# 查看结果

head(chemokine.communication)

dim(chemokine.communication)


# ============================================================
# 18. 导出 CCL / CXCL 通讯结果
# ============================================================

write.xlsx(
  chemokine.communication,
  file = "04_CCL_CXCL_communication.xlsx",
  rowNames = FALSE
)


# ============================================================
# 19. 查看 CCL 通路
# ============================================================

CCL.communication <- subsetCommunication(
  cellchat,
  signaling = "CCL"
)

head(CCL.communication)

write.xlsx(
  CCL.communication,
  file = "05_CCL_communication.xlsx",
  rowNames = FALSE
)


# ============================================================
# 20. 查看 CXCL 通路
# ============================================================

CXCL.communication <- subsetCommunication(
  cellchat,
  signaling = "CXCL"
)

head(CXCL.communication)

write.xlsx(
  CXCL.communication,
  file = "06_CXCL_communication.xlsx",
  rowNames = FALSE
)


# ============================================================
# 21. 筛选 CCL19–CCR7
# ============================================================

CCL19_CCR7 <- all.communication %>%
  filter(
    grepl("CCL19", interaction_name, ignore.case = TRUE) |
      grepl("CCR7", interaction_name, ignore.case = TRUE)
  )

CCL19_CCR7


# 导出

write.xlsx(
  CCL19_CCR7,
  file = "07_CCL19_CCR7.xlsx",
  rowNames = FALSE
)


# ============================================================
# 22. 查看所有与 CCL19 有关的 interaction
# ============================================================

CCL19_interactions <- all.communication %>%
  filter(
    grepl("CCL19", interaction_name, ignore.case = TRUE) |
      grepl("CCL19", ligand, ignore.case = TRUE)
  )

CCL19_interactions


# ============================================================
# 23. 导出 CCL19 相关通讯
# ============================================================

write.xlsx(
  CCL19_interactions,
  file = "08_CCL19_related_interactions.xlsx",
  rowNames = FALSE
)


# ============================================================
# 24. CCL pathway 网络图
# ============================================================

pdf(
  "09_CCL_network.pdf",
  width = 10,
  height = 10
)

netVisual_aggregate(
  cellchat,
  signaling = "CCL",
  layout = "circle"
)

dev.off()


# ============================================================
# 25. CXCL pathway 网络图
# ============================================================

pdf(
  "10_CXCL_network.pdf",
  width = 10,
  height = 10
)

netVisual_aggregate(
  cellchat,
  signaling = "CXCL",
  layout = "circle"
)

dev.off()


# ============================================================
# 26. CCL 通路气泡图
# ============================================================

pdf(
  "11_CCL_bubble.pdf",
  width = 12,
  height = 8
)

netVisual_bubble(
  cellchat,
  signaling = "CCL",
  remove.isolate = FALSE
)

dev.off()


# ============================================================
# 27. CXCL 通路气泡图
# ============================================================

pdf(
  "12_CXCL_bubble.pdf",
  width = 12,
  height = 8
)

netVisual_bubble(
  cellchat,
  signaling = "CXCL",
  remove.isolate = FALSE
)

dev.off()


# ============================================================
# 28. CCL / CXCL sender-receiver 网络
# ============================================================

pdf(
  "13_CCL_CXCL_bubble.pdf",
  width = 14,
  height = 10
)

netVisual_bubble(
  cellchat,
  signaling = c("CCL", "CXCL"),
  remove.isolate = FALSE
)

dev.off()


# ============================================================
# 29. 分析 CCL pathway 的信号贡献
# ============================================================

pdf(
  "14_CCL_contribution.pdf",
  width = 10,
  height = 8
)

netAnalysis_contribution(
  cellchat,
  signaling = "CCL"
)

dev.off()


# ============================================================
# 30. 分析 CCL / CXCL signaling role
# ============================================================

pdf(
  "15_CCL_CXCL_signaling_role.pdf",
  width = 10,
  height = 8
)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netAnalysis_signalingRole_network(
  cellchat,
  signaling = c("CCL", "CXCL")
)

dev.off()


# ============================================================
# 31. 分析主要发送者和接收者
# ============================================================

pdf(
  "16_CCL_CXCL_signaling_role_scatter.pdf",
  width = 10,
  height = 8
)

netAnalysis_signalingRole_scatter(
  cellchat,
  signaling = c("CCL", "CXCL")
)

dev.off()


# ============================================================
# 32. 进一步：提取发送者和接收者
# ============================================================

## 对每一个 chemokine interaction：
##
## source = ligand 来源细胞
## target = receptor 接收细胞

chemokine_sender_receiver <- chemokine.communication %>%
  select(
    source,
    target,
    ligand,
    receptor,
    interaction_name,
    pathway_name,
    prob,
    pval
  ) %>%
  arrange(
    desc(prob)
  )


head(
  chemokine_sender_receiver,
  30
)


# 导出

write.xlsx(
  chemokine_sender_receiver,
  file = "17_Chemokine_sender_receiver.xlsx",
  rowNames = FALSE
)


# ============================================================
# 33. 找出最强的 chemokine interaction
# ============================================================

top_chemokine <- chemokine.communication %>%
  arrange(
    desc(prob)
  ) %>%
  select(
    source,
    target,
    ligand,
    receptor,
    interaction_name,
    pathway_name,
    prob,
    pval
  ) %>%
  head(30)


top_chemokine


write.xlsx(
  top_chemokine,
  file = "18_Top30_Chemokine_interactions.xlsx",
  rowNames = FALSE
)


# ============================================================
# 34. 专门分析 CCL19
# ============================================================

CCL19_network <- all.communication %>%
  filter(
    grepl(
      "CCL19",
      interaction_name,
      ignore.case = TRUE
    )
  ) %>%
  arrange(
    desc(prob)
  )


CCL19_network


# ============================================================
# 35. CCL19 来源细胞
# ============================================================

CCL19_sender <- CCL19_network %>%
  select(
    source,
    ligand,
    receptor,
    target,
    interaction_name,
    prob,
    pval
  ) %>%
  arrange(
    desc(prob)
  )

CCL19_sender


# ============================================================
# 36. 导出 CCL19 来源和靶细胞
# ============================================================

write.xlsx(
  CCL19_sender,
  file = "19_CCL19_sender_receiver.xlsx",
  rowNames = FALSE
)


# ============================================================
# 37. 如果数据库中存在 CCL19–CCR7
# ============================================================

CCL19_CCR7_network <- all.communication %>%
  filter(
    grepl(
      "CCL19",
      interaction_name,
      ignore.case = TRUE
    ) &
      grepl(
        "CCR7",
        interaction_name,
        ignore.case = TRUE
      )
  )


CCL19_CCR7_network


write.xlsx(
  CCL19_CCR7_network,
  file = "20_CCL19_CCR7_network.xlsx",
  rowNames = FALSE
)


# ============================================================
# 38. 保存 CellChat 对象
# ============================================================

saveRDS(
  cellchat,
  file = "CellChat_HCC_chemokine.rds"
)


# ============================================================
# 39. 最终分析结果汇总
# ============================================================

cat("\n")
cat("============================================\n")
cat(" CellChat analysis completed\n")
cat("============================================\n")
cat("\n")

cat("Number of cells:",
    ncol(seu),
    "\n")

cat("Number of cell types:",
    length(unique(meta$celltype)),
    "\n")

cat("\nCell types:\n")
print(unique(meta$celltype))

cat("\n")

cat("Total communication pairs:",
    nrow(all.communication),
    "\n")

cat("Chemokine communication pairs:",
    nrow(chemokine.communication),
    "\n")

cat("CCL19-related interactions:",
    nrow(CCL19_network),
    "\n")

cat("\nResults saved in current working directory.\n")

---
  
  tmp_chemokine_DB = as.data.frame(chemokine_DB)




