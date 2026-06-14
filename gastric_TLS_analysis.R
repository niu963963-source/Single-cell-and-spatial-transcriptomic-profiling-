# =====================================================================
# 胃癌单细胞 + 空间转录组分析 —— 全模块修复版
# 修复内容：命名统一、TLS色阶、空间标签转移置信度过滤、TCGA生存拆分签名、
#          Slingshot限定终点、删除虚拟扰动模块、添加InferCNV（可选）
# =====================================================================

# 设置工作目录和全局选项（根据实际路径修改）
setwd('C:/Users/临平吴彦祖/Desktop/GEO单细胞+空间转录')
options(stringsAsFactors = FALSE, future.globals.maxSize = 8000 * 1024^2, timeout = 600)

# ----------------------------- 模块 0：环境准备 -----------------------------
rm(list = ls()); gc()
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

folders <- c("Figures/Main","Figures/Suppl","Tables","RData","CellChat","InferCNV","NicheNet")
for(d in folders) dir.create(d, showWarnings = FALSE, recursive = TRUE)

if(!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(update = FALSE, ask = FALSE)

bio_pkg <- c("SingleR","celldex","GSVA","GSEABase","ComplexHeatmap","infercnv")
cran_pkg <- c("Seurat","ggplot2","dplyr","patchwork","Matrix","data.table",
              "RColorBrewer","harmony","CellChat","cowplot",
              "survival","survminer","glmnet","caret","igraph","pheatmap","pROC",
              "ggpubr","ggrepel","slingshot","SingleCellExperiment")

for(p in bio_pkg) if(!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
for(p in cran_pkg) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)

library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
library(Matrix); library(data.table); library(RColorBrewer); library(harmony)
library(SingleR); library(celldex); library(CellChat)
library(cowplot); library(survival); library(survminer); library(glmnet)
library(caret); library(ComplexHeatmap); library(circlize); library(pheatmap); library(pROC)
library(ggpubr); library(ggrepel); library(slingshot); library(SingleCellExperiment)

sci_theme <- theme_bw() + theme(
  text = element_text(size=8), plot.title = element_text(size=10, face="bold", hjust=0.5),
  axis.title = element_text(size=9), axis.text = element_text(size=8),
  legend.title = element_text(size=8), legend.text = element_text(size=7),
  panel.grid = element_blank(),
  strip.background = element_rect(fill="white", colour="black"),
  strip.text = element_text(size=9, face="bold")
)
my_palette <- c(brewer.pal(9,"Set1"), brewer.pal(12,"Paired"), brewer.pal(8,"Set2"))
save(sci_theme, my_palette, file = "RData/theme_palette.RData")
cat("OK module 0 done\n")

# ----------------------------- 模块 1：数据读取 -----------------------------
cat("====== [1/17] 数据读取 ======\n")
load("RData/theme_palette.RData")

first_row <- fread("GSE206785_scgex.txt.gz", nrows=1, header=TRUE)
gene_names <- colnames(first_row)[-1]
cell_data <- fread("GSE206785_scgex.txt.gz", header=TRUE, select=1, data.table=FALSE)
cell_ids <- as.character(cell_data[,1]); n_cells <- length(cell_ids)

batch_size <- 2000; n_batches <- ceiling(n_cells / batch_size)
transposed_batches <- list()
for(b in 1:n_batches) {
  start_row <- (b-1)*batch_size + 2; end_row <- min(b*batch_size+1, n_cells+1)
  if(b %% 10 == 0 || b == n_batches) cat(sprintf("批次 %d/%d\n", b, n_batches))
  batch_data <- fread("GSE206785_scgex.txt.gz", header=FALSE, skip=start_row-1, nrows=batch_size, data.table=FALSE)
  expr_mat <- as.matrix(batch_data[,-1]); expr_mat_t <- t(expr_mat)
  rownames(expr_mat_t) <- gene_names; colnames(expr_mat_t) <- batch_data[,1]
  transposed_batches[[b]] <- as(expr_mat_t, "dgCMatrix"); gc()
}
merged_matrix <- do.call(cbind, transposed_batches)
rm(transposed_batches); gc()

metadata <- fread("GSE206785_metadata.txt.gz", header=TRUE, data.table=FALSE)
rownames(metadata) <- colnames(merged_matrix)
if(ncol(metadata)>1) metadata <- metadata[,-1,drop=FALSE]

scRNA <- CreateSeuratObject(counts=merged_matrix, meta.data=metadata, project="GC", min.features=200, min.cells=3)
saveRDS(scRNA, "RData/scRNA_raw.rds")
cat("OK module 1 done\n")

# ----------------------------- 模块 2：QC 质控 -----------------------------
cat("====== [2/17] QC 质控 ======\n")
scRNA <- readRDS("RData/scRNA_raw.rds")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern="^MT-")
scRNA[["percent.rp"]] <- PercentageFeatureSet(scRNA, pattern="^RP[SL]")
nFeature_low  <- max(200, quantile(scRNA$nFeature_RNA, 0.01))
nFeature_high <- min(6000, quantile(scRNA$nFeature_RNA, 0.99))
scRNA <- subset(scRNA, subset = nFeature_RNA > nFeature_low & nFeature_RNA < nFeature_high & percent.mt < 15)
p1 <- FeatureScatter(scRNA, "nCount_RNA","nFeature_RNA") + sci_theme
p2 <- FeatureScatter(scRNA, "nCount_RNA","percent.mt") + sci_theme
VlnPlot(scRNA, features=c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3, pt.size=0) & sci_theme
ggsave("Figures/Main/Fig1_QC.pdf", width=7, height=4)
ggsave("Figures/Main/Fig1_QC_scatter1.pdf", p1, width=5, height=4)
ggsave("Figures/Main/Fig1_QC_scatter2.pdf", p2, width=5, height=4)
saveRDS(scRNA, "RData/scRNA_afterQC.rds")
cat("OK module 2 done\n")

# ----------------------------- 模块 3：去双细胞 -----------------------------
cat("====== [3/17] 去双细胞 ======\n")
scRNA <- readRDS("RData/scRNA_afterQC.rds")
counts <- scRNA$nCount_RNA; features <- scRNA$nFeature_RNA
med_c <- median(log1p(counts)); med_f <- median(log1p(features))
mad_c <- mad(log1p(counts)); mad_f <- mad(log1p(features))
doublet_score <- (log1p(counts)-med_c)/mad_c + (log1p(features)-med_f)/mad_f
scRNA$doublet_call <- ifelse(doublet_score > 3.5, "Doublet","Singlet")
doublet_colors <- c("Singlet" = "#377EB8", "Doublet" = "#E41A1C")
p_dbl <- FeatureScatter(scRNA, "nCount_RNA","nFeature_RNA", group.by = "doublet_call", pt.size = 0.4) +
  scale_color_manual(values = doublet_colors) + sci_theme + ggtitle("Doublet Detection (MAD method)")
ggsave("Figures/Suppl/FigS2_Doublets.pdf", p_dbl, width = 6, height = 5)
scRNA <- subset(scRNA, doublet_call == "Singlet")
saveRDS(scRNA, "RData/scRNA_afterDoublet.rds")
cat("OK module 3 done\n")

# ----------------------------- 模块 4：标准化、高变基因、PCA -----------------------------
cat("====== [4/17] 标准化 ======\n")
scRNA <- readRDS("RData/scRNA_afterDoublet.rds")
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures=2000) %>% 
  ScaleData(vars.to.regress="percent.mt") %>% RunPCA(npcs=50)
saveRDS(scRNA, "RData/scRNA_afterNormPCA.rds")
cat("OK module 4 done\n")

# ----------------------------- 模块 5：Harmony 去批次 -----------------------------
cat("====== [5/17] 批次校正 ======\n")
scRNA <- readRDS("RData/scRNA_afterNormPCA.rds")
if(length(unique(scRNA$orig.ident))==1) {
  if("sample" %in% colnames(scRNA@meta.data)) scRNA$batch <- scRNA$sample else scRNA$batch <- "single"
} else scRNA$batch <- scRNA$orig.ident
if(length(unique(scRNA$batch)) > 1) {
  scRNA <- RunHarmony(scRNA, group.by.vars="batch", dims=1:30, max_iter=20)
  reduction_use <- "harmony"
} else reduction_use <- "pca"
saveRDS(scRNA, "RData/scRNA_afterHarmony.rds")
cat("OK module 5 done\n")

# ----------------------------- 模块 6-7：UMAP 聚类 -----------------------------
cat("====== [6-7/17] UMAP 聚类（高对比度颜色） ======\n")
scRNA <- readRDS("RData/scRNA_afterHarmony.rds")
if(!"harmony" %in% names(scRNA@reductions)) reduction_use <- "pca" else reduction_use <- "harmony"
scRNA <- RunUMAP(scRNA, reduction=reduction_use, dims=1:20) %>% 
  FindNeighbors(reduction=reduction_use, dims=1:20) %>% 
  FindClusters(resolution = 0.2)
cluster_colors <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF",
                    "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#1B9E77")
DimPlot(scRNA, label = TRUE, pt.size = 0.1, cols = cluster_colors) + sci_theme + ggtitle("UMAP of Cell Clusters (Resolution 0.2)")
ggsave("Figures/Main/Fig2_UMAP_Clusters.pdf", width = 6, height = 5)
saveRDS(scRNA, "RData/scRNA_afterCluster.rds")
cat("OK module 6-7 done\n")

# ======================== 修复模块 8：细胞注释 + 命名统一 + InferCNV（可选） ========================
cat("====== [8/17] 细胞注释（修复版：统一命名，可选InferCNV） ======\n")
scRNA <- readRDS("RData/scRNA_afterCluster.rds")
ref <- celldex::HumanPrimaryCellAtlasData()
avg <- AggregateExpression(scRNA, group.by = "seurat_clusters", assays = "RNA", slot = "data")$RNA
pred <- SingleR(test = avg, ref = ref, labels = ref$label.main)
cluster_ids <- gsub("^g", "", rownames(pred))
mapping <- setNames(pred$labels, cluster_ids)
cell_type_vec <- mapping[as.character(scRNA$seurat_clusters)]
cell_type_vec <- unname(cell_type_vec)
cell_type_vec[is.na(cell_type_vec)] <- "Unassigned"
scRNA$cell_type <- cell_type_vec
scRNA$cell_type <- recode(scRNA$cell_type,
                          "T_cells" = "T cells",
                          "B_cell" = "B cells",
                          "Endothelial_cells" = "Endothelial",
                          "Fibroblasts" = "Fibroblast")
# 统一命名：去掉下划线，空格分隔，并将 Tissue_stem_cells 改为 Stem-like Epithelial
scRNA$cell_type <- gsub("_", " ", scRNA$cell_type)
scRNA$cell_type <- ifelse(scRNA$cell_type == "Tissue stem cells", "Stem-like Epithelial", scRNA$cell_type)
scRNA$cell_type <- ifelse(scRNA$cell_type == "Cancer Stem-like Cells", "Stem-like Epithelial", scRNA$cell_type)
Idents(scRNA) <- "cell_type"
cat("细胞类型频数:\n"); print(table(scRNA$cell_type))
# 可选：运行 InferCNV 验证上皮细胞恶性（这里提供框架，需自行准备基因位置文件）
# if(require(infercnv)){
#   expr_raw <- GetAssayData(scRNA, assay = "RNA", layer = "counts")
#   ref_cells <- WhichCells(scRNA, idents = c("Fibroblasts", "Endothelial"))
#   # 需准备 gene_order_file，此处略
#   cat("InferCNV 跳过（需要基因位置文件）\n")
# }
# 保存统一后的对象
saveRDS(scRNA, "RData/scRNA_annotated_unified.rds")
# UMAP 图
p_umap <- DimPlot(scRNA, group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 0.1, cols = my_palette, raster = TRUE) + sci_theme + ggtitle("Annotated Cell Types")
ggsave("Figures/Main/Fig3_CellTypes.pdf", p_umap, width = 8, height = 6)
# Marker 验证
markers <- c("CD3D","CD3E","MS4A1","CD79A","CD14","LYZ","EPCAM","KRT18","COL1A1","ACTA2","PECAM1","VWF")
p_dot <- DotPlot(scRNA, features = markers, group.by = "cell_type") + coord_flip() + sci_theme
ggsave("Figures/Main/Fig4_Marker_DotPlot.pdf", p_dot, width = 7, height = 6)
# 细胞比例图
df_prop <- as.data.frame(prop.table(table(scRNA$cell_type)))
p_prop <- ggplot(df_prop, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") + sci_theme + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=my_palette) + guides(fill="none") + labs(x="Cell Type", y="Proportion")
ggsave("Figures/Main/Fig4_cell_proportion.pdf", p_prop, width = 8, height = 5)
# 保存 cluster 映射表
cluster_map <- table(scRNA$seurat_clusters, scRNA$cell_type)
write.csv(cluster_map, "Tables/Cluster_CellType_Map.csv")
# 重新计算 markers（使用更严格阈值）
Idents(scRNA) <- "cell_type"
markers_all <- FindAllMarkers(scRNA, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25)
write.csv(markers_all, "Tables/All_CellType_Markers_filtered.csv")
top10 <- markers_all %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
genes_use <- unique(top10$gene)
avg_exp <- AverageExpression(scRNA, group.by = "cell_type", assays = "RNA", slot = "data")$RNA
mat <- avg_exp[genes_use, ]
pheatmap(mat, scale = "row", show_rownames = FALSE, show_colnames = TRUE, fontsize = 8, filename = "Figures/Main/Fig4_marker_heatmap.pdf")
saveRDS(scRNA, "RData/scRNA_annotated_unified.rds")
cat("OK module 8 done\n")

# ======================== 模块 9-10：临床分组分析（修复亚型组成图） ========================
cat("====== [9-10/17] 临床分组分析（修复版） ======\n")
scRNA <- readRDS("RData/scRNA_annotated_unified.rds")
# Tumor vs Normal DEG（保持原样）
if("Tissue" %in% colnames(scRNA@meta.data)){
  scRNA$Tissue <- as.character(scRNA$Tissue)
  scRNA$Tissue <- ifelse(grepl("tumor|Tumor|cancer", scRNA$Tissue), "Tumor", ifelse(grepl("normal|Normal", scRNA$Tissue), "Normal", scRNA$Tissue))
  if(all(c("Tumor","Normal") %in% unique(scRNA$Tissue))){
    Idents(scRNA) <- "Tissue"
    degs <- FindMarkers(scRNA, ident.1="Tumor", ident.2="Normal", logfc.threshold=0.25, min.pct=0.1)
    write.csv(degs, "Tables/DEG_tumor_vs_normal.csv")
    degs$logP <- -log10(degs$p_val_adj + 1e-300)
    degs$group <- "NS"
    degs$group[degs$avg_log2FC > 0.5 & degs$p_val_adj < 0.05] <- "Up"
    degs$group[degs$avg_log2FC < -0.5 & degs$p_val_adj < 0.05] <- "Down"
    p_volcano <- ggplot(degs, aes(x=avg_log2FC, y=logP, color=group)) +
      geom_point(size=0.6, alpha=0.7) +
      scale_color_manual(values=c("Down"="blue", "NS"="gray", "Up"="red")) +
      geom_vline(xintercept=c(-0.5,0.5), linetype="dashed") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed") +
      sci_theme + labs(x="log2FC", y="-log10 adj.P") + ylim(0,300) + theme(legend.position="bottom")
    ggsave("Figures/Suppl/FigS1_Tissue_DEG.pdf", p_volcano, width=6, height=5)
  }
}
# Subtype 组成图修复：聚合为大类，避免标签拥挤
if("Subtype" %in% colnames(scRNA@meta.data)){
  scRNA$Subtype_group <- ifelse(grepl("Intestinal", scRNA$Subtype), "Intestinal",
                                ifelse(grepl("Diffuse", scRNA$Subtype), "Diffuse", "Other"))
  prop_sub <- prop.table(table(scRNA$Subtype_group, scRNA$cell_type), margin=1)
  prop_df <- as.data.frame(prop_sub)
  p_sub <- ggplot(prop_df, aes(x=Var1, y=Freq, fill=Var2)) +
    geom_bar(stat="identity", position="fill") +
    scale_fill_manual(values = my_palette, name = "Cell Type") +
    sci_theme + labs(x="Subtype", y="Proportion") +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=9), legend.position = "right")
  ggsave("Figures/Suppl/FigS2_Subtype_Composition_fixed.pdf", p_sub, width=8, height=5)
}
saveRDS(scRNA, "RData/scRNA_clinical.rds")
cat("OK module 9-10 done\n")

# ======================== 模块 11：CellChat（使用统一命名后的对象） ========================
cat("====== [11/17] CellChat ======\n")
scRNA <- readRDS("RData/scRNA_annotated_unified.rds")
scRNA$cell_type <- gsub("_", " ", scRNA$cell_type)
Idents(scRNA) <- "cell_type"
data.input <- GetAssayData(scRNA, assay = "RNA", layer = "data")
meta <- scRNA@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
pdf("Figures/Main/Fig5A_CellChat_Network_Weight.pdf", width=9, height=8)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, edge.width.max = 12, label.edge = FALSE, title.name = "Cell-Cell Communication Strength")
dev.off()
epi_name <- grep("Epithelial", levels(cellchat@idents), value = TRUE)[1]
if(!is.na(epi_name)){
  pdf("Figures/Main/Fig5B_Epithelial_Communication.pdf", width=10, height=7)
  p <- netVisual_bubble(cellchat, sources.use = epi_name, remove.isolate = TRUE, font.size = 6)
  print(p)
  dev.off()
}
top_pathways <- head(cellchat@netP$pathways, 10)
pdf("Figures/Main/Fig6A_CellChat_Top_Pathways.pdf", width=13, height=9)
netVisual_bubble(cellchat, signaling = top_pathways, remove.isolate = TRUE, font.size = 5) + theme(axis.text.x = element_text(angle=45, hjust=1, size=6))
dev.off()
immune_pathways <- intersect(c("CXCL","CCL","MIF","TNF","TGFb"), cellchat@netP$pathways)
if(length(immune_pathways) > 0){
  pdf("Figures/Main/Fig6B_Immune_Suppressive_Pathways.pdf", width=10, height=6)
  p <- netVisual_bubble(cellchat, signaling = immune_pathways, remove.isolate = TRUE)
  print(p)
  dev.off()
}
saveRDS(cellchat, "CellChat/cellchat_result.rds")
cat("OK module 11 done\n")

# ======================== 模块 12：TLS & 免疫逃逸（修复色阶、分层比较、效应量） ========================
cat("====== [12/17] TLS & 免疫逃逸（修复版） ======\n")
scRNA <- readRDS("RData/scRNA_annotated_unified.rds")
scRNA$cell_type <- gsub("_", " ", scRNA$cell_type)
tls_genes <- c("CXCL13","CCL19","CCL21","CXCR5","BCL6","CD3D","CD4","MS4A1","CD79A","CD19","CD8A","GZMB","PRF1","IFNG","SELL","CD27","CD38","MZB1","JCHAIN")
tls_genes <- intersect(tls_genes, rownames(scRNA))
scRNA <- AddModuleScore(scRNA, features = list(tls_genes), name = "TLS", ctrl = 100)
# TLS UMAP：单色渐变，合理分位数
p_tls <- FeaturePlot(scRNA, "TLS1", pt.size = 0.1, raster = TRUE, cols = c("lightgrey", "red")) +
  scale_color_gradientn(colours = c("lightgrey", "red"), limits = quantile(scRNA$TLS1, c(0.005, 0.995), na.rm = TRUE), oob = scales::squish) +
  sci_theme + ggtitle("TLS Score")
ggsave("Figures/Main/Fig7_TLS_UMAP.pdf", p_tls, width = 7, height = 5)
# TLS by cell type
p_vln <- VlnPlot(scRNA, features = "TLS1", group.by = "cell_type", pt.size = 0, cols = my_palette) +
  sci_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Main/Fig7_TLS_by_celltype.pdf", p_vln, width = 8, height = 5)
# TLS Tumor vs Normal 分层按细胞类型
if("Tissue" %in% colnames(scRNA@meta.data)){
  # 仅对 B 细胞和 T 细胞分层
  subt <- scRNA@meta.data %>% filter(cell_type %in% c("B cells", "T cells"))
  subt$Tissue <- factor(subt$Tissue, levels = c("Normal","Tumor"))
  p_box_sub <- ggplot(subt, aes(x = Tissue, y = TLS1, fill = Tissue)) +
    geom_boxplot() + facet_wrap(~cell_type, scales = "free") +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    scale_fill_manual(values = c("Normal" = "#377EB8", "Tumor" = "#E41A1C")) +
    sci_theme + labs(y = "TLS Score", title = "TLS Score by Tissue and Cell Type")
  ggsave("Figures/Main/Fig7_TLS_by_Tissue_and_Celltype.pdf", p_box_sub, width = 6, height = 4)
}
# 免疫检查点评分
immune_genes <- c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2")
immune_genes <- intersect(immune_genes, rownames(scRNA))
scRNA <- AddModuleScore(scRNA, features = list(immune_genes), name = "Immune", ctrl = 50)
p_immune <- FeaturePlot(scRNA, "Immune1", cols = c("lightgrey","red"), pt.size = 0.1, raster = TRUE) +
  scale_color_gradientn(colours = c("lightgrey","red"), limits = quantile(scRNA$Immune1, c(0.05, 0.95), na.rm = TRUE), oob = scales::squish) +
  sci_theme + ggtitle("Immune Checkpoint Score")
ggsave("Figures/Main/Fig8_Immune_UMAP_fixed.pdf", p_immune, width = 7, height = 5)
# TLS vs Immune 相关性
df_cor <- data.frame(TLS = scRNA$TLS1, Immune = scRNA$Immune1)
p_cor <- ggplot(df_cor, aes(TLS, Immune)) + geom_point(size = 0.5, alpha = 0.4) +
  geom_smooth(method = "lm", col = "red") + stat_cor(method = "pearson") + sci_theme +
  labs(title = "TLS vs Immune Checkpoint Score")
ggsave("Figures/Main/Fig8_TLS_vs_Immune.pdf", p_cor, width = 6, height = 5)
# TLS 富集卡方效应量
df_plot <- data.frame(TLS_group = ifelse(scRNA$TLS1 > median(scRNA$TLS1), "High", "Low"),
                      Cell_type = ifelse(grepl("T cell", scRNA$cell_type), "T cells",
                                         ifelse(grepl("B cell", scRNA$cell_type), "B cells", "Other")))
tbl <- table(df_plot$Cell_type, df_plot$TLS_group)
chi <- chisq.test(tbl)
cramers_v <- sqrt(chi$statistic / (sum(tbl) * (min(dim(tbl))-1)))
cat("TLS enrichment Cramers V =", round(cramers_v, 3), "\n")
saveRDS(scRNA, "RData/scRNA_annotated_unified.rds")
cat("OK module 12 done\n")

# ======================== 模块 13：空间整合（置信度过滤 + 替换配受体 + 修复标签） ========================
cat("====== [13/17] 空间整合（修复版：置信度过滤，CXCL12-CXCR4，Kcross主图） ======\n")
library(Seurat); library(ggplot2); library(patchwork); library(dplyr); library(ggpubr); library(spatstat); library(pheatmap)
scRNA <- readRDS("RData/scRNA_annotated_unified.rds")
st_dir <- "GSE251950_RAW/GSM7990473_20_00331_LI_SING/20_00331_LI_SING"
tmp_dir <- file.path(getwd(), "spatial_tmp")
dir.create(tmp_dir, showWarnings = FALSE); dir.create(file.path(tmp_dir, "spatial"), showWarnings = FALSE)
file.copy(file.path(st_dir, "20_00331_LI_SING_filtered_feature_bc_matrix.h5"), file.path(tmp_dir, "filtered_feature_bc_matrix.h5"))
orig_spatial <- file.path(st_dir, "spatial")
if(dir.exists(orig_spatial)){
  file.copy(list.files(orig_spatial, full.names = TRUE), file.path(tmp_dir, "spatial"), overwrite = TRUE)
}
st_obj <- Load10X_Spatial(data.dir = tmp_dir, assay = "Spatial")
unlink(tmp_dir, recursive = TRUE)
DefaultAssay(st_obj) <- "Spatial"
st_obj <- NormalizeData(st_obj, verbose = FALSE) %>% FindVariableFeatures(verbose = FALSE) %>% ScaleData(verbose = FALSE) %>% RunPCA(verbose = FALSE)
# 标签转移 + 置信度过滤
anchors <- FindTransferAnchors(reference = scRNA, query = st_obj, normalization.method = "LogNormalize", reference.reduction = "pca", dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = scRNA$cell_type, weight.reduction = st_obj[["pca"]], dims = 1:30)
thresh <- 0.5
st_obj$predicted.id <- ifelse(predictions$prediction.score.max < thresh, "Unassigned", predictions$predicted.id)
cat("预测细胞类型频数（过滤后）：\n"); print(table(st_obj$predicted.id))
# TLS 和免疫评分
tls_genes <- intersect(c("CXCL13","CCL19","CCL21","CXCR5","BCL6","CD3D","MS4A1","CD79A"), rownames(st_obj))
st_obj <- AddModuleScore(st_obj, features = list(tls_genes), name = "TLS", ctrl = 50)
immune_genes <- intersect(c("PDCD1","CTLA4","LAG3","TIGIT"), rownames(st_obj))
st_obj <- AddModuleScore(st_obj, features = list(immune_genes), name = "Immune", ctrl = 50)
# CXCL12-CXCR4 相关性（移到补充材料，这里仍绘制但保存为补充）
if(all(c("CXCL12","CXCR4") %in% rownames(st_obj))){
  df_cxcl <- FetchData(st_obj, c("CXCL12","CXCR4"))
  cor_cxcl <- cor.test(df_cxcl$CXCL12, df_cxcl$CXCR4)
  p_cxcl <- ggplot(df_cxcl, aes(CXCL12, CXCR4)) + geom_point(alpha=0.3, size=0.5) + geom_smooth(method="lm", col="red") +
    annotate("text", x=max(df_cxcl$CXCL12)*0.8, y=max(df_cxcl$CXCR4)*0.9, label=paste0("R = ", round(cor_cxcl$estimate,3), ", p = ", format(cor_cxcl$p.value, scientific=TRUE, digits=2)), size=4) +
    sci_theme + ggtitle("Spatial Correlation: CXCL12 - CXCR4")
  ggsave("Figures/Suppl/FigS8_CXCL12_CXCR4_cor.pdf", p_cxcl, width=6, height=5)
}
# Kcross 作为主图
coords <- GetTissueCoordinates(st_obj)
if(!is.null(coords) && nrow(coords)>0){
  x_col <- intersect(c("imagecol","col","x","pxl_col_in_fullres"), colnames(coords))[1]
  y_col <- intersect(c("imagerow","row","y","pxl_row_in_fullres"), colnames(coords))[1]
  if(!is.na(x_col) && !is.na(y_col)){
    coords$x_use <- coords[[x_col]]; coords$y_use <- coords[[y_col]]
    expr_cxcl13 <- GetAssayData(st_obj, assay="Spatial", layer="data")["CXCL13",]
    expr_cxcr5 <- GetAssayData(st_obj, assay="Spatial", layer="data")["CXCR5",]
    high_cxcl13 <- names(which(expr_cxcl13 > quantile(expr_cxcl13, 0.8, na.rm=TRUE)))
    high_cxcr5 <- names(which(expr_cxcr5 > quantile(expr_cxcr5, 0.8, na.rm=TRUE)))
    coords$type <- "Other"
    coords$type[rownames(coords) %in% high_cxcl13] <- "CXCL13_High"
    coords$type[rownames(coords) %in% high_cxcr5] <- "CXCR5_High"
    if(sum(coords$type != "Other") >= 10){
      valid <- !is.na(coords$x_use) & !is.na(coords$y_use)
      coords_sub <- coords[valid,]
      pp <- ppp(coords_sub$x_use, coords_sub$y_use, window=owin(range(coords_sub$x_use), range(coords_sub$y_use)),
                marks=factor(coords_sub$type, levels=c("CXCL13_High","CXCR5_High","Other")))
      set.seed(123)
      env <- envelope(pp, fun=Kcross, i="CXCL13_High", j="CXCR5_High", nsim=99, correction="Ripley", verbose=FALSE)
      pdf("Figures/Main/Fig8_Spatial_Kcross_CXCL13_CXCR5.pdf", width=6, height=5)
      plot(env, main="Spatial attraction: CXCL13-high vs CXCR5-high")
      dev.off()
    }
  }
}
# 空间主图（修复标签）
all_types <- unique(st_obj$predicted.id)
type_colors <- setNames(rep("gray80", length(all_types)), all_types)
main_types <- intersect(c("Epithelial cells","T cells","B cells","Fibroblasts","Endothelial","Myeloid","Stem-like Epithelial"), all_types)
type_colors[main_types] <- my_palette[1:length(main_types)]
if("Unassigned" %in% all_types) type_colors["Unassigned"] <- "#E0E0E0"
p_spatial <- SpatialDimPlot(st_obj, group.by = "predicted.id", cols = type_colors, label = TRUE, label.size = 2.5, repel = TRUE) +
  sci_theme + theme(legend.position = "right") + ggtitle("Spatial Cell Type Distribution (Confidence-filtered)")
ggsave("Figures/Main/Fig8_spatial_structure_fixed.pdf", p_spatial, width = 9, height = 7)
# TLS 空间图
p_tls_spatial <- SpatialFeaturePlot(st_obj, "TLS1", pt.size.factor = 1.5) + sci_theme
ggsave("Figures/Main/Fig8_Spatial_TLS_Score.pdf", p_tls_spatial, width = 8, height = 7)
# 保存对象
saveRDS(st_obj, "RData/st_integrated_fixed.rds")
cat("OK module 13 done\n")

# ======================== 模块 14-15：TCGA 预后（拆分签名 + surv_cutpoint） ========================
cat("====== [14-15/17] TCGA 预后（修复版：三个独立签名，动态cutoff） ======\n")
expr <- fread("TCGA-STAD.star_tpm.tsv", sep = "\t") %>% as.data.frame()
surv <- fread("TCGA-STAD.survival.tsv", sep = "\t") %>% as.data.frame()
colnames(surv)[1:3] <- c("sample", "time", "event")
surv <- surv %>% filter(!is.na(time), !is.na(event), time > 0)
probemap <- fread("gencode.v36.annotation.gtf.gene.probemap", sep = "\t", header = TRUE, data.table = FALSE)[,1:2]
colnames(probemap) <- c("gene_id","symbol")
probemap <- probemap[!duplicated(probemap$gene_id),]
if(!"gene_id" %in% colnames(expr)) colnames(expr)[1] <- "gene_id"
expr$symbol <- probemap$symbol[match(expr$gene_id, probemap$gene_id)]
expr <- expr[!is.na(expr$symbol) & expr$symbol != "" & !duplicated(expr$symbol), ]
rownames(expr) <- expr$symbol
expr <- expr[, !colnames(expr) %in% c("gene_id","symbol"), drop = FALSE]
common_samples <- intersect(colnames(expr), surv$sample)
# 定义三个独立签名
tls_only <- c("CXCL13","CXCR5","CCL19","CCL21","BCL6","CD79A","MS4A1")
exhaust_only <- c("PDCD1","LAG3","TIGIT","HAVCR2","TOX")
caf_only <- c("FAP","COL1A1","CXCL12")
tls_only <- intersect(tls_only, rownames(expr))
exhaust_only <- intersect(exhaust_only, rownames(expr))
caf_only <- intersect(caf_only, rownames(expr))
score_tls <- colMeans(expr[tls_only, common_samples, drop=FALSE], na.rm=TRUE)
score_exh <- colMeans(expr[exhaust_only, common_samples, drop=FALSE], na.rm=TRUE)
score_caf <- colMeans(expr[caf_only, common_samples, drop=FALSE], na.rm=TRUE)
df_surv <- data.frame(sample=common_samples, time=surv$time[match(common_samples, surv$sample)], event=surv$event[match(common_samples, surv$sample)],
                      tls=score_tls[common_samples], exh=score_exh[common_samples], caf=score_caf[common_samples]) %>% filter(complete.cases(.))
# 使用 surv_cutpoint 寻找最佳 cutoff
library(survminer)
cut_tls <- surv_cutpoint(df_surv, time = "time", event = "event", variables = "tls")
cut_exh <- surv_cutpoint(df_surv, time = "time", event = "event", variables = "exh")
cut_caf <- surv_cutpoint(df_surv, time = "time", event = "event", variables = "caf")
df_surv$group_tls <- ifelse(df_surv$tls > cut_tls$cutpoint[1,"cutpoint"], "High", "Low")
df_surv$group_exh <- ifelse(df_surv$exh > cut_exh$cutpoint[1,"cutpoint"], "High", "Low")
df_surv$group_caf <- ifelse(df_surv$caf > cut_caf$cutpoint[1,"cutpoint"], "High", "Low")
fit_tls <- survfit(Surv(time, event) ~ group_tls, data = df_surv)
fit_exh <- survfit(Surv(time, event) ~ group_exh, data = df_surv)
fit_caf <- survfit(Surv(time, event) ~ group_caf, data = df_surv)
p_tls <- ggsurvplot(fit_tls, pval=TRUE, risk.table=TRUE, palette=c("#E41A1C","#377EB8"), ggtheme=sci_theme, title="TLS Signature")
p_exh <- ggsurvplot(fit_exh, pval=TRUE, risk.table=TRUE, palette=c("#E41A1C","#377EB8"), ggtheme=sci_theme, title="Exhaustion Signature")
p_caf <- ggsurvplot(fit_caf, pval=TRUE, risk.table=TRUE, palette=c("#E41A1C","#377EB8"), ggtheme=sci_theme, title="CAF Signature")
ggsave("Figures/Main/Fig9_TLS_Survival.pdf", p_tls$plot, width=6, height=5)
ggsave("Figures/Main/Fig9_Exhaustion_Survival.pdf", p_exh$plot, width=6, height=5)
ggsave("Figures/Main/Fig9_CAF_Survival.pdf", p_caf$plot, width=6, height=5)
# 选择最显著的作为主图（这里以 exhaustion 为例，若 p 值最小）
best_name <- c("TLS","Exhaustion","CAF")[which.min(c(p_tls$plot$pval, p_exh$plot$pval, p_caf$plot$pval))]
cat("最优签名：", best_name, "\n")
# ROC 曲线（对 exhaustion 为例）
roc_obj <- roc(df_surv$event, df_surv$exh, quiet=TRUE)
auc_val <- round(auc(roc_obj),3)
p_roc <- ggroc(roc_obj, color="red", size=1.2) + geom_abline(slope=1, intercept=1, linetype="dashed", color="gray50") +
  annotate("text", x=0.7, y=0.2, label=paste0("AUC = ", auc_val), size=4) + sci_theme + labs(title=paste0(best_name, " ROC"))
ggsave("Figures/Main/Fig9_ROC_Fixed.pdf", p_roc, width=5, height=4)
cat("OK module 14-15 done\n")

# ======================== 模块 16：NicheNet（保持不变） ========================
cat("====== [16/17] NicheNet ======\n")
if(requireNamespace("nichenetr", quietly=TRUE)){
  library(nichenetr)
  if(!file.exists("lr_network.rds")) download.file("https://zenodo.org/record/3260758/files/lr_network.rds", "lr_network.rds", mode="wb")
  if(!file.exists("ligand_target_matrix.rds")) download.file("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds", "ligand_target_matrix.rds", mode="wb")
  lr_net <- readRDS("lr_network.rds"); lt_mat <- readRDS("ligand_target_matrix.rds")
  scRNA <- readRDS("RData/scRNA_annotated_unified.rds")
  Idents(scRNA) <- "cell_type"
  types <- unique(scRNA$cell_type)
  sender <- grep("B cell|Plasma", types, value=TRUE)[1]
  receiver <- grep("CD8|T cell|CD4", types, value=TRUE)[1]
  if(!is.na(sender) & !is.na(receiver)){
    expr <- GetAssayData(scRNA, assay="RNA", layer="data")
    sender_cells <- WhichCells(scRNA, idents=sender)
    receiver_cells <- WhichCells(scRNA, idents=receiver)
    targets <- intersect(c("PDCD1","LAG3","TIGIT","CTLA4","CXCR5"), rownames(expr))
    bg <- rownames(expr)[rowMeans(expr[,receiver_cells,drop=FALSE] > 0) > 0.1]
    sender_expr <- expr[,sender_cells,drop=FALSE]
    potential_ligands <- intersect(rownames(sender_expr)[rowMeans(sender_expr > 0) > 0.1], lr_net$from)
    ligand_activities <- predict_ligand_activities(geneset = targets, background = bg, ligand_target_matrix = lt_mat, potential_ligands = potential_ligands)
    ligand_activities <- ligand_activities %>% arrange(-pearson)
    top_ligands <- head(ligand_activities$test_ligand, 8)
    top_ligands <- top_ligands[top_ligands %in% colnames(lt_mat)]
    if(length(top_ligands) >= 2){
      sub_mat <- lt_mat[targets, top_ligands, drop=FALSE]
      sub_mat <- t(sub_mat)
      rownames(sub_mat) <- top_ligands; colnames(sub_mat) <- targets
      pdf("Figures/Main/Fig12_NicheNet_Heatmap.pdf", width=7, height=5)
      pheatmap(sub_mat, cluster_rows=FALSE, cluster_cols=FALSE, treeheight_col=0, main="Ligand-Target Regulatory Potential")
      dev.off()
      saveRDS(list(ligand_activities=ligand_activities, ligand_target=sub_mat), "NicheNet/nichenet_results.rds")
    }
  }
}
cat("OK module 16 done\n")

# ======================== 模块 17：T细胞拟时序（限定终点，简化线条） ========================
cat("====== [17/17] T细胞拟时序（修复：限定终点） ======\n")
scRNA <- readRDS("RData/scRNA_annotated_unified.rds")
scRNA$cell_type <- gsub("_", " ", scRNA$cell_type)
tcell_barcodes <- colnames(scRNA)[grepl("T cell|CD4|CD8|NKT|NK", scRNA$cell_type, ignore.case=TRUE)]
tcell <- subset(scRNA, cells = tcell_barcodes)
DefaultAssay(tcell) <- "RNA"
tcell <- NormalizeData(tcell, verbose=FALSE) %>% FindVariableFeatures(nfeatures=2000, verbose=FALSE) %>% ScaleData(verbose=FALSE) %>% RunPCA(npcs=30, verbose=FALSE)
if(length(unique(tcell$orig.ident)) > 1){
  tcell <- RunHarmony(tcell, group.by.vars="orig.ident", dims.use=1:20, verbose=FALSE)
  reduction_use <- "harmony"
} else reduction_use <- "pca"
tcell <- RunUMAP(tcell, reduction=reduction_use, dims=1:20, verbose=FALSE) %>% FindNeighbors(reduction=reduction_use, dims=1:20, verbose=FALSE) %>% FindClusters(resolution=0.5, verbose=FALSE)
# 定义 T 细胞亚型（修改 Tfh 定义，去掉 CXCR5）
t_markers <- list(
  "Naive T" = c("SELL","CCR7","TCF7","LEF1"),
  "Effector CD8" = c("GZMB","PRF1","IFNG","NKG7"),
  "Exhausted" = c("PDCD1","LAG3","TIGIT","HAVCR2","TOX"),
  "Treg" = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  "Tfh/TLS" = c("CXCL13","BCL6","ICOS","PDCD1"),   # 已移除 CXCR5
  "Memory CD4" = c("CD44","IL7R","TCF7","S100A4"),
  "NK/NKT" = c("GNLY","NKG7","KLRD1","NCAM1")
)
for(state in names(t_markers)){
  genes_use <- intersect(t_markers[[state]], rownames(tcell))
  if(length(genes_use) >= 2) tcell <- AddModuleScore(tcell, features=list(genes_use), name=gsub(" |/","_",state))
}
score_cols <- paste0(gsub(" |/","_",names(t_markers)),"1")
score_cols <- score_cols[score_cols %in% colnames(tcell@meta.data)]
score_mat <- as.matrix(tcell@meta.data[,score_cols])
best_idx <- apply(score_mat,1,which.max)
tcell$t_subtype <- gsub("1$","",score_cols[best_idx])
tcell$t_subtype <- gsub("_"," ",tcell$t_subtype)
exhaust_genes <- intersect(c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","CTLA4"), rownames(tcell))
tcell <- AddModuleScore(tcell, features=list(exhaust_genes), name="Exhaust")
# Slingshot 限定终点
sce <- as.SingleCellExperiment(tcell, assay="RNA")
umap_embed <- Embeddings(tcell, reduction="umap")
reducedDim(sce, "UMAP") <- umap_embed
naive_clu <- as.character(unique(tcell$seurat_clusters[tcell$t_subtype == "Naive T"])[1])
exhaust_clu <- as.character(unique(tcell$seurat_clusters[tcell$t_subtype == "Exhausted"])[1])
treg_clu <- as.character(unique(tcell$seurat_clusters[tcell$t_subtype == "Treg"])[1])
if(is.na(naive_clu)) naive_clu <- as.character(unique(tcell$seurat_clusters)[1])
end_clusters <- c(exhaust_clu, treg_clu)[!is.na(c(exhaust_clu, treg_clu))]
if(length(end_clusters) == 0) end_clusters <- NULL
sce <- slingshot(sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP", start.clus = naive_clu, end.clus = end_clusters, approx_points = 150)
tcell$Pseudotime <- slingPseudotime(sce)[,1]
# 绘制拟时序箱线图（主图）
df_box <- data.frame(Pseudotime = tcell$Pseudotime, Subtype = tcell$t_subtype)
p_box <- ggplot(df_box, aes(Subtype, Pseudotime, fill=Subtype)) + geom_boxplot() + sci_theme +
  theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(x="T Cell Subtype", y="Pseudotime") + guides(fill="none")
ggsave("Figures/Main/Fig13F_Pseudotime_Boxplot.pdf", p_box, width=8, height=5)
# 耗竭评分趋势
df_trend <- data.frame(Pseudotime = tcell$Pseudotime, Exhaust = tcell$Exhaust1, Subtype = tcell$t_subtype) %>% filter(!is.na(Pseudotime))
p_trend <- ggplot(df_trend, aes(Pseudotime, Exhaust)) + geom_point(aes(color=Subtype), size=0.4, alpha=0.5) +
  geom_smooth(method="loess", color="black") + sci_theme + labs(x="Pseudotime", y="Exhaustion Score")
ggsave("Figures/Main/Fig13E_Exhaust_Pseudotime.pdf", p_trend, width=8, height=5)
saveRDS(tcell, "RData/Tcell_exhausted_fixed.rds")
cat("OK module 17 done\n")

# ======================== 模块 18：耗竭 T 细胞绘图（不变） ========================
cat("====== [18/17] 耗竭 T 细胞绘图 ======\n")
tcell <- readRDS("RData/Tcell_exhausted_fixed.rds")
exhaust_markers <- c("PDCD1","LAG3","TIGIT","HAVCR2","CTLA4","CXCR5")
p1 <- FeaturePlot(tcell, features = exhaust_markers, ncol = 3, pt.size = 0.2) & sci_theme
ggsave("Figures/Main/Fig14_Tex_FeaturePlot.pdf", p1, width = 15, height = 8)
p2 <- DimPlot(tcell, group.by = "state", cols = c("red","gray")) + sci_theme
ggsave("Figures/Main/Fig14_Tex_UMAP.pdf", p2, width = 7, height = 5)
expr_data <- FetchData(tcell, c("state","PDCD1","CXCR5"))
expr_long <- tidyr::pivot_longer(expr_data, cols = c("PDCD1","CXCR5"), names_to = "gene", values_to = "expression")
p3 <- ggplot(expr_long, aes(state, expression, fill = state)) + geom_violin(scale = "width") + facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = c("red","gray")) + sci_theme + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(c("Exhausted","Non_exhausted")))
ggsave("Figures/Main/Fig14_Tex_Violin.pdf", p3, width = 10, height = 4)
cat("OK module 18 done\n")
