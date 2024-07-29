
library(SpotSweeper)

library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir = here('plots',"new_qc_stuff")
processed_dir = here('raw-data')


# large DLPFC dataset
load(here("raw-data","dlPFC_raw.RData"))
spe.pfc <- spe_raw
spe.pfc
# class: SpatialExperiment
# dim: 36601 149757
# metadata(0):
#   assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817 ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(149757): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
# colData names(27): sample_id in_tissue ... sample_id_complete count
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor



# subset to hangnail artifact



# ========= Testing MALAT1 expression as a new QC feature ==========

spe.pfc <- spe.pfc[,colData(spe.pfc)$in_tissue == TRUE]
rownames(spe.pfc) <- rowData(spe.pfc)$gene_name
logcounts(spe.pfc) <- log2(counts(spe.pfc) + 1)


# get random other sample
spe.hangnail <- spe.pfc[,colData(spe.pfc)$sample_id == "Br8325_ant"]
spe.other <- spe.pfc[,colData(spe.pfc)$sample_id == "Br8325_post"]

spe.hangnail$counts_MALAT1 <- log2(counts(spe.hangnail)[which(rowData(spe.hangnail)$gene_name=="MALAT1"),]+1)
spe.other$counts_MALAT1 <- log2(counts(spe.other)[which(rowData(spe.other)$gene_name=="MALAT1"),]+1)

make_escheR(spe.other) |>
  add_fill("counts_MALAT1")

make_escheR(spe.hangnail) |>
  add_fill("counts_MALAT1")

spe.pfc$counts_MALAT1 <- log2(counts(spe.pfc)[which(rowData(spe.pfc)$gene_name=="MALAT1"),]+1)

plotQCpdf(spe.pfc, metric="counts_MALAT1", outliers=NULL, fname=here(plot_dir, "MALAT1_expression.pdf"), point_size=1.5)





# plot mito ratio
plotQC(spe.hangnail, metric="expr_chrM_ratio")

# ====== New way to normalize =====

# house keeping genes
hk.genes <- c("GAPDH", "ACTB", "B2M", "RPL13A", "TUBA1A", "PGK1", "PPIA", "YWHAG")
# List of heat shock protein genes
hk.genes <- c(
  "HSPA1A", "HSPA1B", "HSPA2", "HSPA5", "HSPA6", "HSPA8",
  "HSPB1", "HSPB2", "HSPB3", "HSPB6", "HSPB7", "HSPB8",
  "HSP90AA1", "HSP90AB1", "HSP90B1",
  "HSPH1", "HSPD1", "HSPE1",
  "DNAJB1", "DNAJB2", "DNAJB4", "DNAJB6",
  "DNAJC1", "DNAJC5", "DNAJC6", "DNAJC7", "DNAJC10",
  "HSPH1", "HSPH2"
)

# Print the list
print(heatshock_genes)


# get proportion hk.genes
hk.genes <- hk.genes[hk.genes %in% rownames(spe.hangnail)]

hk.df <- perCellQCMetrics(spe.hangnail, subsets=list(House=hk.genes))


spe.hangnail$subsets_House_detected <- hk.df$subsets_House_detected
spe.hangnail$subsets_House_sum <- hk.df$subsets_House_sum
spe.hangnail$subsets_House_percent <- hk.df$subsets_House_percent


plotQC(spe.hangnail, metric="subsets_House_sum")


# ====== Quantile normalization =====
library(dplyr)

total_counts <- rowSums(counts(spe),1)
names(total_counts) <- rownames(spe)
total_counts_df <- data.frame(gene = names(total_counts), count = total_counts, stringsAsFactors = FALSE)


bottom_25_percent_genes <- total_counts_df %>%
  arrange(count) %>%
  slice(1:round(0.95 * n())) %>%
  pull(gene)


quantile.df <- perCellQCMetrics(spe, subset=list(Q25=c("MOBP","MBP", "SLC17A7","SNAP25")))
quantile.df

spe$subsets_Q25_percent <- quantile.df$subsets_Q25_detected


make_escheR(spe) |>
  add_fill(var="subsets_Q25_percent")




# ======== Geometric mean ==========

# Load necessary libraries
library(SpatialExperiment)
library(dplyr)
library(matrixStats)

# Function to calculate the geometric mean
geometric_mean <- function(x) {
  exp(mean(log(x[x > 0]))) # use only positive values to avoid -Inf in log
}

# Apply the geometric mean function to each column (spot)
geometric_means <- apply(assay(spe), 2, geometric_mean)

# Calculate the total counts per spot
total_counts <- spe$sum_umi

# Normalize the total counts by the geometric means
normalized_counts <- total_counts / geometric_means

# Add the normalized counts to the SPE object (optional)
colData(spe)$normalized_counts <- normalized_counts

# View the normalized counts
head(colData(spe)$normalized_counts)

make_escheR(spe) |>
  add_fill(var="normalized_counts")





# ======== using crtlGene =========
library("TREG")

spe <- logNormCounts(spe)


sum(assays(spe)$counts == 0) / (nrow(spe) * ncol(spe))
# 1] 0.9651388

row_means <- rowMeans(assays(spe)$logcounts)
(median_row_means <- median(row_means))

spe<- spe[row_means > median_row_means, ]
dim(spe)
# [1] 18300  3529

sum(assays(spe)$counts == 0) / (nrow(spe) * ncol(spe))
# [1] 0.9303263

spe$group4treg <- TRUE
rownames(spe) <- make.unique(rownames(spe))

prop_zeros <- get_prop_zero(spe, group_col="group4treg")
head(prop_zeros)

ggplot( data = prop_zeros, aes(x=prop_zeros)) +
  geom_histogram()

hist(prop_zeros$TRUE)

filtered_genes <- filter_prop_zero(prop_zeros, cutoff = 0.75)
length(filtered_genes)

spe <- spe[filtered_genes, ]

sum(assays(spe)$counts == 0) / (nrow(spe) * ncol(spe))
# [1] 0.5354088

# == Rank invariance
## Get the rank of the gene in each group
group_rank <- rank_group(spe, group_col = "group4treg")

## Get the rank of the gene for each cell
spot_rank <- rank_cells(spe, group_col = "group4treg")

## Use both rankings to calculate rank_invariance()
rank_invar <- rank_invariance(group_rank, spot_rank)
head(sort(rank_invar,decreasing=TRUE))



top_20 <- sort(rank_invar,decreasing=TRUE)[1:1000]
# MT-CO3    MT-ND2    MT-CO2    MT-ND4   MT-ATP6    MT-CO1    MT-ND3    MT-ND1    MT-CYB MTRNR2L12    MT-ND5    SNAP25       MT3       CLU     RPL41      FTH1      GNAS
# 1167      1166      1165      1164      1163      1162      1161      1160      1159      1158      1157      1156      1155      1154      1153      1152      1151
# MTRNR2L1    NDUFA4    RPL37A
# 1150      1149      1148

# NOTE: MT genes are among the top rank invariant genes, which is a good gut check!

rank_invar_no_mito <- rank_invar[!grepl("^MT-", names(rank_invar))]

top_20 <- sort(rank_invar_no_mito,decreasing=TRUE)[5:10]
top_20
# MTRNR2L12    SNAP25       MT3       CLU     RPL41      FTH1      GNAS  MTRNR2L1    NDUFA4    RPL37A     GAPDH    TUBA1B    TMSB10     CALM2      NRGN     CALM1       CKB
#      1158      1156      1155      1154      1153      1152      1151      1150      1149      1148      1147      1146      1145      1144      1143      1142      1140
#      CST3     CALM3      ACTB
#      1139      1138      1137

# === Now using our new rank invariant genes ===
invar_genes <- names((top_20))

quantile.df <- perCellQCMetrics(spe, subset=list(RI=invar_genes))
quantile.df

spe$subsets_RI_percent <- quantile.df$subsets_RI_percent


plotQC(spe, metric="subsets_RI_percent")


plotColData(spe, x="subsets_RI_sum", y="sum_umi")
plotColData(spe, x="expr_chrM", y="sum_umi")




# =========== Find genes highly correlated with mito expression ==========


# Extract expression data from SPE object
counts <- counts(spe)

# Identify mitochondrial genes
mito_genes <- grep("^MT-", rownames(counts), value = TRUE)

# Calculate mean expression of mitochondrial genes
mito_expr <- colMeans(counts[mito_genes, , drop = FALSE])

# Calculate correlations
correlation_vector <- apply(counts, 1, function(gene_expr) cor(gene_expr, mito_expr, use = "pairwise.complete.obs"))
correlation_vector <- correlation_vector[!grepl("^MT-", names(correlation_vector))]
top_20 <- names(sort(correlation_vector,decreasing=TRUE))[1:20]

# === Now using our new rank invariant genes ===
invar_genes <- top_20

quantile.df <- perCellQCMetrics(spe, subset=list(MtCorr=invar_genes))
quantile.df

spe$subsets_MtCorr_percent <- quantile.df$subsets_MtCorr_percent
spe$subsets_MtCorr_sum <- quantile.df$subsets_MtCorr_sum


plotQC(spe, metric="subsets_MtCorr_percent")

plotColData(spe, x="subsets_MtCorr_sum", y="sum_umi")
