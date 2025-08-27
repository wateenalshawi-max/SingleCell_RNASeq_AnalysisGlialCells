############################################################
# Single-cell RNA-seq QC, Normalization, and Embeddings
# Dataset: ASD (Velmeshev et al. 2019, Science)
# Merged script (yours + fixes)
############################################################

## --- 0) Packages ------------------------------------------------------------
 install.packages("data.table")
 install.packages("BiocManager")

suppressPackageStartupMessages({
  library(data.table)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(scuttle)
  library(AnnotationHub)
  library(AnnotationDbi)
  library(ensembldb)
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
})

bpp <- MulticoreParam(max(1, parallel::detectCores() - 1))

## --- 1) Load metadata & counts ---------------------------------------------
m <- fread("/home/participant/Downloads/metaasd/meta.tsv")

# Clean metadata column names (spaces/punct -> underscores)
setnames(m, gsub("[^A-Za-z0-9_]+","_", names(m)))

cat("Samples:\n"); print(table(m$sample))
cat("Unique samples:", length(unique(m$sample)), "\n")

# Load 10x counts (expects rawMatrix/ with matrix.mtx, features.tsv/genes.tsv, barcodes.tsv)
sce <- read10xCounts("rawMatrix", col.names = TRUE)
print(list.files("rawMatrix"))

## --- 2) Quick expression landscape -----------------------------------------
obj <- sce
genesPerCell <- colSums(counts(obj) > 0)

total_genes_detected <- sum(rowSums(counts(obj) > 0) > 0)
ref_genes <- nrow(obj)
cat("Total genes detected across cells:", total_genes_detected, "out of", ref_genes, "\n")
print(summary(genesPerCell))
rv <- range(genesPerCell)
cat("Range (min–max):", rv[1], "to", rv[2], "\n")
cat("Median:", median(genesPerCell), "\n")

plot(density(genesPerCell), main = "", xlab = "Genes per cell")

## --- 3) Merge metadata onto SCE --------------------------------------------
# Ensure we have a barcode column; fall back to colnames if needed
if (!"Barcode" %in% colnames(colData(sce))) {
  colData(sce)$Barcode <- colnames(sce)
}
# Expect m$cell to match the barcodes/colnames used
stopifnot(all(colnames(sce) %in% m$cell))
m_ord <- m[match(colnames(sce), m$cell)]
stopifnot(all(m_ord$cell == colnames(sce)))

# Bind all metadata except the duplicate "cell" column
colData(sce) <- cbind(colData(sce), DataFrame(m_ord[, !c("cell")]))

# Sanity checks
print(colSums(is.na(as.data.frame(colData(sce)))[, c("cluster","diagnosis","sample","individual")]))
print(with(as.data.frame(colData(sce)), table(cluster, diagnosis))[1:10, ])
print(with(as.data.frame(colData(sce)), table(sample, diagnosis))[1:10, ])

cat("dim(meta): ", paste(dim(m), collapse=" x "), "\n")
cat("dim(sce): ", paste(dim(sce), collapse=" x "), "\n")
cat("dim(colData): ", paste(dim(colData(sce)), collapse=" x "), "\n")
print(rownames(counts(sce))[1:6])
print(dim(counts(sce)))

## --- 4) Gene-level quick views ---------------------------------------------
# Total UMI for a gene vs. detection proportion
plot(rowSums(counts(sce)) / pmax(rowSums(counts(sce) > 0), 1),
     rowMeans(counts(sce) > 0),
     log = "x",
     xlab="Mean UMIs per expressing cell",
     ylab="Proportion of cells expressing the gene")

# Distribution of counts for most expressed genes
rel_expression <- t( t(counts(sce)) / colSums(counts(sce)) ) * 100
# try to use Symbol if present
if (!"Symbol" %in% colnames(rowData(sce))) rowData(sce)$Symbol <- rownames(sce)
rownames(rel_expression) <- rowData(sce)$Symbol
most_expressed <- sort(rowSums(rel_expression), decreasing = TRUE)[20:1]
plot_data <- as.matrix(t(rel_expression[names(most_expressed),]))
boxplot(plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)

# Keep only detected genes
detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)
sce <- sce[detected_genes,]

## --- 5) Gene annotation (EnsDb v107) ---------------------------------------
ah <- AnnotationHub()
ensdb <- query(ah, c("Homo sapiens", "EnsDb", "107"))[[1]]

ids <- rowData(sce)$ID
if (is.null(ids)) ids <- rownames(sce)
ids_nov <- sub("\\.\\d+$", "", ids)  # strip version suffixes

ann <- AnnotationDbi::select(
  ensdb,
  keys    = unique(ids_nov),
  keytype = "GENEID",
  columns = c("GENEID","SEQNAME","SYMBOL","GENEBIOTYPE")
)

idx <- match(ids_nov, ann$GENEID)
rowData(sce)$Chromosome   <- ann$SEQNAME[idx]
rowData(sce)$Gene_biotype <- ann$GENEBIOTYPE[idx]
need_sym <- is.na(rowData(sce)$Symbol) | rowData(sce)$Symbol == ""
rowData(sce)$Symbol[need_sym] <- ann$SYMBOL[idx][need_sym]

# Flag mitochondrials & compute per-cell QC
chrM_alias <- c("MT","chrM","chrMT","M")
rowData(sce)$is_mito <- rowData(sce)$Chromosome %in% chrM_alias

qc <- perCellQCMetrics(sce, subsets = list(Mt = rowData(sce)$is_mito))
colData(sce) <- cbind(colData(sce), qc)

## --- 6) QC visualization & thresholds --------------------------------------
# Distributions
plotColData(sce, x="sample", y="sum") + scale_y_log10() + ggtitle("Total count")
plotColData(sce, x="sample", y="detected") + scale_y_log10() + ggtitle("Detected features")
plotColData(sce, x="sample", y="subsets_Mt_percent") + ggtitle("Mito percent")

# UMI vs genes colored by mito
as.data.frame(colData(sce)) %>%
  arrange(subsets_Mt_percent) %>%
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = subsets_Mt_percent > 10)) +
  facet_wrap(vars(diagnosis)) +
  labs(x = "Total UMIs per cell", y = "Number of detected genes", colour = "Mito > 10%") +
  theme_minimal()

summary(colData(sce)$subsets_Mt_percent)
sum(colData(sce)$subsets_Mt_percent > 10, na.rm = TRUE)
table(rowData(sce)$is_mito, useNA = "ifany")

# Data-driven mito threshold (95th percentile)
thr <- quantile(colData(sce)$subsets_Mt_percent, 0.95, na.rm = TRUE)
colData(sce)$mito_high <- colData(sce)$subsets_Mt_percent > thr

as.data.frame(colData(sce)) %>%
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = mito_high), alpha = 0.5) +
  facet_wrap(vars(diagnosis)) +
  labs(x = "Total UMIs per cell", y = "Number of detected genes",
       colour = paste0("Mito >", round(thr, 2), "%")) +
  theme_minimal()

# Outlier flags (batch-aware)
low_lib_size <- isOutlier(colData(sce)$sum,      log = TRUE, type = "lower", nmads = 3)
low_n_genes  <- isOutlier(colData(sce)$detected, log = TRUE, type = "lower", nmads = 3)
high_mito    <- isOutlier(colData(sce)$subsets_Mt_percent,    type = "higher", nmads = 3)

table(low_lib_size); table(low_n_genes); table(high_mito)
attr(high_mito, "thresholds")
table(sce$Seqbatch, high_mito)

# Summary per sample
tibble(
  sample = sce$sample,
  low_lib_size = low_lib_size,
  low_n_features = low_n_genes,
  high_Mito_percent = high_mito
) %>%
  mutate(discard = low_lib_size | low_n_features | high_Mito_percent) %>%
  group_by(sample) %>%
  summarise(
    n_low_lib_size    = sum(low_lib_size,    na.rm = TRUE),
    n_low_n_features  = sum(low_n_features,  na.rm = TRUE),
    n_high_mito       = sum(high_Mito_percent, na.rm = TRUE),
    n_discard         = sum(discard,         na.rm = TRUE),
    n_cells           = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_discard)) %>%
  print(n = 50)

## quickPerCellQC: all cells vs per-sample thresholds
cell_qc_results <- quickPerCellQC(sce, sub.fields = TRUE)
batch.cell_qc_results <- quickPerCellQC(colData(sce), sub.fields = TRUE, batch = sce$sample)

all.thresholds <- tibble(
  Sample = "All",
  `Library Size` = attr(cell_qc_results$low_lib_size, "thresholds")[1],
  `Genes detected` = attr(cell_qc_results$low_n_features, "thresholds")[1],
  `Mitochondrial UMIs` = attr(cell_qc_results$high_subsets_Mt_percent, "thresholds")[2]
)

batch.thresholds <- tibble(
  Sample = names(attr(batch.cell_qc_results$low_lib_size, "thresholds")[1, ]),
  `Library Size` = attr(batch.cell_qc_results$low_lib_size, "thresholds")[1, ],
  `Genes detected` = attr(batch.cell_qc_results$low_n_features, "thresholds")[1, ],
  `Mitochondrial UMIs` = attr(batch.cell_qc_results$high_subsets_Mt_percent, "thresholds")[2, ]
)

print(bind_rows(batch.thresholds, all.thresholds))

# Store batch-aware flags on sce
sce$low_lib_size      <- batch.cell_qc_results$low_lib_size
sce$low_n_features    <- batch.cell_qc_results$low_n_features
sce$high_Mito_percent <- batch.cell_qc_results$high_subsets_Mt_percent
sce$discard           <- batch.cell_qc_results$discard

plotColData(sce, x = "sample", y = "sum", colour_by = "discard") +
  scale_y_log10() +
  labs(y = "Total count", title = "Total count by sample")

## --- 7) Size factors & normalization ---------------------------------------
# Quick clustering (for deconvolution pools)
set.seed(100)
clust <- quickCluster(sce, BPPARAM = bpp)
table(clust)

# Pooled size factors
sce <- computePooledFactors(sce, clusters = clust, min.mean = 0.1, BPPARAM = bpp)
deconv.sf <- sizeFactors(sce)
lib.sf    <- librarySizeFactors(sce)

data.frame(
  LibrarySizeFactors = lib.sf,
  DeconvolutionSizeFactors = deconv.sf,
  Sample = sce$sample
) %>%
  ggplot(aes(x = LibrarySizeFactors, y = DeconvolutionSizeFactors)) +
  geom_point(aes(col = Sample)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Library size factors", y = "Deconvolution size factors",
       title = "Comparison of size factors") +
  theme_classic() + theme(legend.position = "bottom")

# Log-normalize & store on object (adds assay 'logcounts')
sce <- logNormCounts(sce)

# If you want non-log normalized counts per cell (for barplot), do it explicitly:
norm_mat <- normalizeCounts(counts(sce), size_factors = sizeFactors(sce))  # non-log
norm_per_cell <- colSums(norm_mat)
norm_df <- tibble(Barcode = colnames(sce), normCounts = log2(norm_per_cell + 1))


## --- 7) HVGs, PCA, t-SNE, UMAP ---------------------------------------------
set.seed(1)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 2000)

sce <- runPCA(sce, subset_row = hvg, ncomponents = 50, BSPARAM = BiocSingular::ExactParam())
sce <- runTSNE(sce, dimred = "PCA", perplexity = 30)
sce <- runUMAP(sce, dimred = "PCA", n_neighbors = 30, min_dist = 0.3)

# Quick embedding plots
plotPCA(sce, colour_by = "diagnosis")
plotTSNE(sce, colour_by = "cluster")
plotUMAP(sce, colour_by = "sample")

## --- 8) Save if desired -----------------------------------------------------
saveRDS(sce, file = "sce_qc_norm_embeddings.rds")
