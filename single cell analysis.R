## =========================
## Libraries
## =========================
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(pheatmap)
library(ggplot2)
# If you want gene symbols:
ibrary(AnnotationDbi); library(org.Hs.eg.db)

set.seed(123)

## =========================
## 0) Column names in colData(sce)
##    
## =========================
ct_col <- grep("(cell.?type|cluster|annotation|label|cell_class|celltype)",
               colnames(colData(sce)), ignore.case = TRUE, value = TRUE)[1]
reg_col <- grep("^region$", colnames(colData(sce)), ignore.case = TRUE, value = TRUE)[1]
dx_col  <- grep("(diagnos|dx|status)", colnames(colData(sce)), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(ct_col), !is.na(reg_col), !is.na(dx_col))

## =========================
## 1) Subset to our 8 strata (Micro/Oligo × ACC/PFC × ASD/CTRL)
## =========================
ctype  <- tolower(as.character(colData(sce)[[ct_col]]))
region <- toupper(as.character(colData(sce)[[reg_col]]))
dx     <- toupper(as.character(colData(sce)[[dx_col]]))

is_micro <- grepl("microg", ctype)
is_oligo <- grepl("oligo",  ctype)
is_accpfc <- region %in% c("ACC","PFC")
is_asdctrl<- dx %in% c("ASD","CONTROL")

keep8 <- (is_micro | is_oligo) & is_accpfc & is_asdctrl
sce8  <- sce[, keep8]

CellType8 <- ifelse(is_micro[keep8], "Microglia", "Oligodendrocyte")
Region8   <- region[keep8]
Dx8       <- ifelse(dx[keep8] == "CONTROL", "Control", "ASD")

grp8 <- factor(paste(CellType8, Region8, Dx8, sep = "_"),
               levels = c("Microglia_ACC_ASD","Microglia_ACC_Control",
                          "Microglia_PFC_ASD","Microglia_PFC_Control",
                          "Oligodendrocyte_ACC_ASD","Oligodendrocyte_ACC_Control",
                          "Oligodendrocyte_PFC_ASD","Oligodendrocyte_PFC_Control"))

## =========================
## 2) Helper: run DE per stratum (ASD vs Control *within* CT × Region)
## =========================
run_stratum <- function(target_ct = c("Micro","Oligo"), target_region = c("ACC","PFC"),
                        fdr_cut = 1) {
  target_ct <- match.arg(target_ct)
  ct_vec <- tolower(as.character(colData(sce)[[ct_col]]))
  reg_vec<- toupper(as.character(colData(sce)[[reg_col]]))
  dx_vec <- toupper(as.character(colData(sce)[[dx_col]]))
  
  is_ct  <- if (target_ct=="Micro") grepl("microg", ct_vec) else grepl("oligo", ct_vec)
  keep   <- is_ct & reg_vec == toupper(target_region) & dx_vec %in% c("ASD","CONTROL")
  sce_sub <- sce[, keep]
  grp_dx  <- factor(ifelse(dx_vec[keep]=="ASD","ASD","Control"), levels=c("ASD","Control"))
  
  mk <- scran::findMarkers(sce_sub, groups = grp_dx, direction = "up", lfc = 0)
  # Tidy the ASD table (ASD vs Control)
  df <- as.data.frame(mk$ASD)
  df$gene <- rownames(df)
  # find an LFC column (pairwise or summary)
  lfc_col <- if ("summary.logFC" %in% names(df)) "summary.logFC" else grep("logFC", names(df), value = TRUE)[1]
  out <- df[, c("gene","FDR", lfc_col)]
  names(out) <- c("gene","FDR","logFC")
  out$CellType <- ifelse(target_ct=="Micro","Microglia","Oligodendrocyte")
  out$Region   <- toupper(target_region)
  dplyr::filter(out, FDR < fdr_cut)
}

## =========================
## 3) Gather DE per stratum & pick genes by DE (NOT variance)
## =========================
fdr_cut <- 0.05
de_ACC_Micro <- run_stratum("Micro","ACC", fdr_cut)
de_PFC_Micro <- run_stratum("Micro","PFC", fdr_cut)
de_ACC_Oligo <- run_stratum("Oligo","ACC", fdr_cut)
de_PFC_Oligo <- run_stratum("Oligo","PFC", fdr_cut)

de_all <- dplyr::bind_rows(de_ACC_Micro, de_PFC_Micro, de_ACC_Oligo, de_PFC_Oligo)
stopifnot(nrow(de_all) > 0)

topN_per <- 25  # tweak per panel if you want fewer/more
pick_top <- function(df) df %>%
  dplyr::arrange(FDR, dplyr::desc(abs(logFC))) %>%
  dplyr::slice_head(n = topN_per) %>%
  dplyr::pull(gene)

genes_to_plot <- unique(c(
  pick_top(de_ACC_Micro), pick_top(de_PFC_Micro),
  pick_top(de_ACC_Oligo), pick_top(de_PFC_Oligo)
))

## =========================
## 4) Build mean-expression matrix (genes × 8 groups)
## =========================
logm <- assay(sce8, "logcounts")
genes_plot <- intersect(genes_to_plot, rownames(logm))
stopifnot(length(genes_plot) > 0)

mean_by_group <- function(m, g)
  sapply(levels(g), function(L) Matrix::rowMeans(m[, g==L, drop=FALSE]))

MEAN <- mean_by_group(logm[genes_plot, , drop=FALSE], grp8)  # genes × 8
M_z  <- t(scale(t(as.matrix(MEAN)))); M_z[is.na(M_z)] <- 0

(Optional) Map Ensembl IDs → gene symbols for nicer row labels
ens2sym <- function(ens_ids){
core <- sub("\\.\\d+$", "", ens_ids)
sym  <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = core, keytype="ENSEMBL",
                                 column="SYMBOL", multiVals="first")
   sym  <- sym[core]
   out  <- ifelse(is.na(sym) | sym=="", ens_ids, sym)
   make.unique(out)
 }
 rownames(MEAN) <- rownames(M_z) <- ens2sym(rownames(MEAN))

## =========================
## 5) Column annotations & order (ASD LEFT of Control)
## =========================
parts <- do.call(rbind, strsplit(colnames(M_z), "_"))
ann_col <- data.frame(
  CellType  = parts[,1],
  Region    = parts[,2],
  Diagnosis = parts[,3],
  row.names = colnames(M_z),
  check.names = FALSE
)

# Order: CellType → Region; within each block, ASD then Control
ord <- with(ann_col,
            order(factor(CellType,  c("Microglia","Oligodendrocyte")),
                  factor(Region,    c("ACC","PFC")),
                  factor(Diagnosis, c("ASD","Control"))))
M_z     <- M_z[, ord, drop = FALSE]
ann_col <- ann_col[ord, , drop = FALSE]

# Vertical gaps to "split the middle" between blocks (great on iPad)
block <- paste(ann_col$CellType, ann_col$Region, sep = "_")
gaps_col <- cumsum(rle(block)$lengths); gaps_col <- gaps_col[-length(gaps_col)]

## =========================
## 6) Palettes (nicer diverging colors, centered at 0)
## =========================
hm_cols <- colorRampPalette(c("#2b8cbe","#f7fbff","#ef8a62","#b2182b"))(255)

ann_colors <- list(
  CellType  = c(Microglia = "#2ca25f", Oligodendrocyte = "#f39c12"),
  Region    = c(ACC = "#e41a1c", PFC = "#377eb8"),
  Diagnosis = c(ASD = "#7b3294", Control = "#9e9e9e")
)

## =========================
## 7) Heatmap (poster ready)
## =========================
hp <- pheatmap::pheatmap(
  M_z,
  color = hm_cols,
  show_rownames = TRUE,
  fontsize_row  = 7,
  main = "Significant DEGs (ASD left of Control within each CT × Region)",
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,     # keep our ASD-first order
  gaps_col = gaps_col,
  border_color = NA
)

## =========================
## 8) Dot plot in the SAME order as the heatmap
## =========================
det_mat <- (logm[genes_plot, ] > 0) * 1
PCT <- sapply(levels(grp8), function(L) Matrix::rowMeans(det_mat[, grp8==L, drop=FALSE]))

# Keep the same row/column order as the heatmap
ord_genes  <- rownames(M_z)[hp$tree_row$order]
ord_groups <- colnames(M_z)  # cols not clustered, already ordered

df_dot <- data.frame(
  Gene  = rep(rownames(MEAN), times = ncol(MEAN)),
  Group = rep(colnames(MEAN), each  = nrow(MEAN)),
  mean  = as.vector(MEAN),
  pct   = as.vector(PCT)
)
df_dot$Gene  <- factor(df_dot$Gene,  levels = ord_genes)
df_dot$Group <- factor(df_dot$Group, levels = ord_groups)

ggplot(df_dot, aes(x = Group, y = Gene)) +
  geom_point(aes(size = pct, color = mean)) +
  scale_size(range = c(1, 6), breaks = c(0.25, 0.5, 0.75, 1.0), name = "% detected") +
  scale_color_gradientn(colors = hm_cols, name = "Mean log-expression") +
  labs(title = "Dotplot: Significant DEGs (ASD vs Control)",
       x = "Group (CellType_Region_Dx)", y = "Gene") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
