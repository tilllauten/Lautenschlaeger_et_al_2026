
  library(tidyverse)
  library(Seurat)
  library(openxlsx)
  library(anndataR)
  library(SummarizedExperiment)
  library(DOtools)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplotify)
  library(ggrastr)


koenig_robj <- "GSE183852_DCM_Nuclei.Robj"
koenig_metadata <- "Metadata_koenig2022.xlsx"
chaffin_rds <- "human_dcm_hcm_scportal_03.17.2022.rds"
kuppe_h5ad <- "54d24dbe-d39a-4844-bb21-07b5f4e173ad.h5ad"
kuppe_h5ad_raw <- "RAWCOUNTS_54d24dbe-d39a-4844-bb21-07b5f4e173ad.h5ad"
dimmeler_rds <- "251112_Human_HF_complete.rds"
scvi_model_dir <- "scvi_model"
scanvi_model_dir <- "scanvi_model"
out_dir <- "output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Koenig et al. 2022
load(koenig_robj)
DefaultAssay(nuclei) <- "RNA"
meta_koenig <- read.xlsx(koenig_metadata)
nuclei$sex <- plyr::revalue(meta_koenig$Sex[match(nuclei$orig.ident, meta_koenig$Samples)],
                            c("male" = "M", "female" = "F"))
nuclei$age <- meta_koenig$Age[match(nuclei$orig.ident, meta_koenig$Samples)]
nuclei$chemistry <- "5p_v1_1"
nuclei <- DO.Subset(nuclei, ident = "orig.ident",
                    ident_name = c("TWCM-13-84", "TWCM-13-285", "TWCM-LVAD2", "TWCM-10-5", "TWCM-13-17",
                                   "TWCM-13-208", "TWCM-LVAD3", "TWCM-13-102", "H_ZC-LVAD"))
nuclei$dataset <- "koenig"
nuclei <- NormalizeData(nuclei)

# Chaffin et al. 2022
Seu_Chaffin <- readRDS(chaffin_rds)
Seu_Chaffin$nCount_RNA <- Seu_Chaffin$cellbender_ncount
Seu_Chaffin$nFeature_RNA <- Seu_Chaffin$cellbender_ngenes
Seu_Chaffin$sex <- plyr::revalue(Seu_Chaffin$sex, c("female" = "F", "male" = "M"))
Seu_Chaffin <- DO.Subset(Seu_Chaffin, ident = "donor_id",
                         ident_name = c("P1290", "P1300", "P1304", "P1358", "P1371", "P1430",
                                        "P1437", "P1472", "P1504", "P1606", "P1617"))
Seu_Chaffin$chemistry <- "3p_v3"
Seu_Chaffin$dataset <- "chaffin"

# Kuppe et al. 2022
reduced_dims <- list("PCA" = c(sampleFactors = "X_pca"),
                     "UMAP" = c(sampleFactors = "X_umap"),
                     "harmony" = c(sampleFactors = "X_harmony"))
Kuppe_sce <- read_h5ad(kuppe_h5ad)$as_SingleCellExperiment(x_mapping = "data", reducedDims_mapping = reduced_dims)
Kuppe_sce_raw <- read_h5ad(kuppe_h5ad_raw)$as_SingleCellExperiment(x_mapping = "counts", reducedDims_mapping = reduced_dims)
counts(Kuppe_sce) <- counts(Kuppe_sce_raw)
Kuppe_sce <- Kuppe_sce[!grepl("^ENSG", rowData(Kuppe_sce)$feature_name), ]
rownames(Kuppe_sce) <- make.names(rowData(Kuppe_sce)$feature_name, unique = TRUE)

kuppe_Seu <- as.Seurat(Kuppe_sce, data = "data", counts = "counts")
kuppe_Seu <- RenameAssays(kuppe_Seu, assay.name = "originalexp", new.assay.name = "RNA")
kuppe_Seu[["RNA"]] <- as(kuppe_Seu[["RNA"]], Class = "Assay5")
kuppe_Seu$sex <- plyr::revalue(kuppe_Seu$sex, c("female" = "F", "male" = "M"))
kuppe_Seu <- DO.Subset(kuppe_Seu, ident = "donor_id", ident_name = c("P19", "P14", "P4", "P5", "P18", "P20"))
kuppe_Seu$chemistry <- "3p_v3"
kuppe_Seu$dataset <- "kuppe"

# Ruz Jurado et al. 2025 (Dimmeler)
Human_HF_complete <- readRDS(dimmeler_rds)
Human_HF_complete <- DO.Subset(Human_HF_complete, ident = "disease", ident_name = c("CTRL.LV", "HFrEF"))
Human_HF_complete$disease <- factor(Human_HF_complete$disease, levels = c("CTRL.LV", "HFrEF"))
v2_chems <- c("Control_LV_CBTM362C", "Control_LV_CBTM364B", "Control_LV_CBTM390C",
              "Control_LV_CBTM417C", "Control_LV_CBTM423C")
Human_HF_complete$chemistry <- ifelse(Human_HF_complete$orig.ident %in% v2_chems, "3p_v2", "3p_v3")
Human_HF_complete$dataset <- "dimmeler"

# Merge and QC
Seu_merged <- merge(
  x = Human_HF_complete,
  y = list(kuppe_Seu, Seu_Chaffin, nuclei),
  add.cell.ids = c("dimmeler", "kuppe", "chaffin", "koenig"),
  project = "Till_HFrEF_Atlas"
)
Seu_merged <- JoinLayers(Seu_merged)
Seu_merged <- subset(Seu_merged, subset = nFeature_RNA >= 300 & nFeature_RNA <= 6000 &
                       nCount_RNA >= 500 & nCount_RNA <= 15000)
Seu_merged@assays$SCT <- NULL

Seu_merged$orig.ident <- ifelse(Seu_merged$orig.ident == "SeuratProject", NA, Seu_merged$orig.ident)
Seu_merged$orig.ident <- coalesce(Seu_merged$orig.ident, Seu_merged$donor_id)
Seu_merged$percent.mito <- coalesce(Seu_merged$percent.mito, Seu_merged$percent_mito,
                                    Seu_merged$cellranger_percent_mito)
Seu_merged$annotation <- ifelse(is.na(Seu_merged$annotation), "Unknown", Seu_merged$annotation)

Seu_merged <- readRDS()

# scVI integration
Seu_merged <- DO.scVI(
  sce_object = Seu_merged,
  batch_key = "orig.ident",
  layer_counts = "counts",
  layer_logcounts = "logcounts",
  categorical_covariates = list("chemistry", "dataset"),
  continuos_covariates = list("nCount_RNA", "nFeature_RNA", "percent.mito"),
  get_model = TRUE,
  model_save_path = scvi_model_dir
)

# scANVI integration and label transfer from the Dimmeler annotation
DO.scANVI <- function(sce_object, model_path, labels_key, unlabeled_category,
                      batch_size = 128, save_model = NULL) {
  args <- list(model_path = model_path, labels_key = labels_key,
               unlabeled_category = unlabeled_category, batch_size = as.integer(batch_size),
               save_model = save_model)
  res <- basilisk::basiliskRun(
    env = DOtools:::DOtoolsEnv,
    fun = function(args) {
      scvi <- reticulate::import("scvi")
      model <- scvi$model$SCVI$load(args$model_path)
      model_scanvi <- scvi$model$SCANVI$from_scvi_model(scvi_model = model,
                                                        unlabeled_category = args$unlabeled_category,
                                                        labels_key = args$labels_key)
      model_scanvi$train(batch_size = args$batch_size)
      if (!is.null(args$save_model)) model_scanvi$save(args$save_model, overwrite = TRUE, save_anndata = TRUE)
      list(embedding = model_scanvi$get_latent_representation(), labels_pred = model_scanvi$predict())
    },
    args = args
  )
  emb <- res$embedding
  rownames(emb) <- colnames(sce_object)
  sce_object[["scANVI"]] <- CreateDimReducObject(embeddings = emb, key = "scANVI_",
                                                 assay = DefaultAssay(sce_object))
  sce_object$scANVI_labels <- res$labels_pred
  sce_object
}

Seu_merged <- DO.scANVI(Seu_merged, model_path = scvi_model_dir, labels_key = "annotation",
                        unlabeled_category = "Unknown", save_model = scanvi_model_dir)
Seu_merged <- FindNeighbors(Seu_merged, reduction = "scANVI", dims = 1:30)
Seu_merged <- RunUMAP(Seu_merged, reduction = "scANVI", reduction.name = "umap", dims = 1:30)

# Endothelial / immune reclustering: cluster 20 = immune-like EC
immune <- c("T_cells", "MP_recruited", "MP_resident", "Lym_EC", "Mast")
ec_broad <- c("Art_EC", "Cap_EC", "Ven_EC", "PC_EC")
Seu_merged$broad <- ifelse(Seu_merged$scANVI_labels %in% immune, "Immune",
                           ifelse(Seu_merged$scANVI_labels %in% ec_broad, "EC", "rest"))

Seu_merged_EC_IM <- DO.Subset(Seu_merged, ident = "broad", ident_name = c("EC", "Immune"))
Seu_merged_EC_IM <- FindNeighbors(Seu_merged_EC_IM, reduction = "scANVI", dims = 1:30)
Seu_merged_EC_IM <- FindClusters(Seu_merged_EC_IM, resolution = 1, algorithm = 4, random.seed = 1234)
Seu_merged_EC_IM$annotation_new <- as.character(Seu_merged_EC_IM$scANVI_labels)
Seu_merged_EC_IM$annotation_new[Seu_merged_EC_IM$seurat_clusters == "20"] <- "Immune like EC"

Seu_merged$annotation <- Seu_merged$scANVI_labels
Seu_merged <- DO.TransferLabel(sce_object = Seu_merged, Subset_obj = Seu_merged_EC_IM,
                               annotation_column = "annotation", subset_annotation = "annotation_new")

# Endothelial reclustering: cluster 4 = endocardial EC
Seu_merged_EC <- DO.Subset(Seu_merged, ident = "annotation", ident_name = c(ec_broad, "Immune like EC"))
Seu_merged_EC <- FindNeighbors(Seu_merged_EC, reduction = "scANVI", dims = 1:30)
Seu_merged_EC <- FindClusters(Seu_merged_EC, resolution = 0.3, algorithm = 4, random.seed = 1234)
Seu_merged_EC$annotation <- as.character(Seu_merged_EC$annotation)
Seu_merged_EC$annotation[as.character(Seu_merged_EC$seurat_clusters) == "4"] <- "End_EC"

Seu_merged <- DO.TransferLabel(sce_object = Seu_merged, Subset_obj = Seu_merged_EC,
                               annotation_column = "annotation", subset_annotation = "annotation")
saveRDS(Seu_merged, file.path(out_dir, "Seu_firstversion_annotated.rds"))

# Harmonise condition labels across datasets
m <- Seu_merged[[]]
blank <- function(z) is.na(z) | trimws(as.character(z)) %in% c("", "NA", "N/A")
source_label <- as.character(m$disease)
source_label[blank(source_label)] <- as.character(m$Condition[blank(source_label)])
condition <- unname(c("CTRL.LV" = "ctrl", "DCM" = "HFrEF", "HFrEF" = "HFrEF",
                      "myocardial infarction" = "HFrEF")[source_label])
stopifnot(!anyNA(condition))

clean <- m[, names(m)[vapply(m, function(z) !any(blank(z)), logical(1))], drop = FALSE]
clean$condition <- factor(condition, levels = c("ctrl", "HFrEF"))
clean$sample_id <- paste(clean$dataset, clean$orig.ident, sep = "::")
Seu_merged@meta.data <- clean
saveRDS(Seu_merged, file.path(out_dir, "Seu_clean_harmonized.rds"), compress = "gzip")

rm(Seu_merged, Seu_merged_EC, Seu_merged_EC_IM, nuclei, Seu_Chaffin, kuppe_Seu, Kuppe_sce, Kuppe_sce_raw, Human_HF_complete)
invisible(gc())

seu <- readRDS('./Seu_clean_harmonized.rds')

font <- "Arial"
set.seed(1234)

ctrl_values <- c("CTRL.LV", "CTRL", "ctrl", "Ctrl", "Control", "control", "NF", "Donor", "donor", "non-failing")
hfref_values <- c("HFrEF", "DCM", "myocardial infarction", "HCM", "ICM")
dataset_names <- c(chaffin = "chaffin", dimmeler = "ruz jurado", koenig = "koenig", kuppe = "kuppe")

celltype_cols <- c(Act_FB = "#C7C7C7", Adip = "#B0C4DE", Art_EC = "#EE8262", Cap_EC = "#AEC7E8",
                   CM = "#CD661D", End_EC = "#DBDB8D", FB = "#C49C94", Lym_EC = "#FF7F0E",
                   Mast = "#C5B0D5", MP_recruited = "#228B22", MP_resident = "#9467BD",
                   NC_glial = "#BCBD22", PC = "#F7B6D2", PC_EC = "#1F77B4", SMC = "#FF9896",
                   T_cells = "#DA7970", Ven_EC = "#FFE4B5", Monocytes = "#17BECF")

markers <- c("COL1A1", "COL1A2", "ADIPOQ", "PLIN1", "GJA5", "EFNB2", "TNNT2", "MYH6", "NPR3", "NRG1",
             "CD74", "HLA-DRA", "LYVE1", "PROX1", "PDPN", "KIT", "CPA3", "PLP1", "SOX10", "RGS5",
             "PDGFRB", "KDR", "VWF", "MYH11", "TAGLN", "PTPRC","CD3E", "CD3D",  "ACKR1", "SELE")
celltype_order <- c("Act-FB", "Adip", "Art-EC", "CM", "End-EC", "MP-recruited", "MP-resident", "Lym-EC",
                    "Mast", "NC-glial", "PC", "PC-EC", "Cap-EC", "SMC", "T-cells", "Monocytes", "Ven-EC")

ec_subtypes <- c("Art_EC", "Cap_EC", "Ven_EC", "End_EC", "Lym_EC", "PC_EC")
ec_cols <- c(Art_EC = "#EE8262", Cap_EC = "#AEC7E8", Ven_EC = "#FFE4B5",
             End_EC = "#DBDB8D", Lym_EC = "#FF7F0E", PC_EC = "#1F77B4")
ec_markers <- c("PECAM1", "CDH5", "VWF", "GJA5", "EFNB2", "SEMA3G", "HEY1", "CA4", "RGCC", "BTNL9",
                "ACKR1", "SELE", "NPR3", "NRG1", "CDH11", "PROX1", "LYVE1", "PDPN", "CCL21", "RGS5", "PDGFRB")

target_gene <- "CIITA"
exclude_fine_clusters <- c("12")
min_cells <- 5
min_ctrl_patients <- 2
paper_cols <- c(Ctrl = "#D4D4D4", HFrEF = "#F94040")

save_all <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
}

rasterise_points <- function(p) {
  p$layers[[1]]$aes_params$stroke <- 0
  p$layers[[1]] <- ggrastr::rasterise(p$layers[[1]], dpi = 600)
  p
}

make_heatmap <- function(mat) {
  col_fun <- circlize::colorRamp2(seq(-2, 2, length.out = 7),
                        c("#083269", "#8DC1DA", "#D4E6EE", "white", "#F9DCCC", "#EC9977", "#BF0027"))
  ht <- ComplexHeatmap::Heatmap(mat, col = col_fun, name = "Z-score",
                cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_side = "left", column_names_side = "bottom", column_names_rot = 90,
                row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 10),
                rect_gp = grid::gpar(col = "grey80", lwd = 0.5))
  ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(ht)))
}

scaled_avg <- function(obj, features, group_by, col_order) {
  avg <- AverageExpression(obj, features = features, group.by = group_by, assays = "RNA")[[1]]
  colnames(avg) <- gsub("_", "-", colnames(avg))
  z <- t(scale(t(as.matrix(avg))))
  z <- z[features[features %in% rownames(z)], col_order[col_order %in% colnames(z)], drop = FALSE]
  z <- z[apply(z, 1, function(x) all(is.finite(x))), , drop = FALSE]
  z[z > 2] <- 2
  z[z < -2] <- -2
  z
}

umap_label <- function(emb, n) {
  x0 <- min(emb[, 1]); y0 <- min(emb[, 2])
  lx <- diff(range(emb[, 1])) * 0.15; ly <- diff(range(emb[, 2])) * 0.15
  list(annotate("text", x = x0 + lx * 0.1, y = y0 + ly * 0.35,
                label = paste0("n = ", format(n, big.mark = ",")), hjust = 0, size = 5),
       theme(axis.line = element_blank(),
             axis.title.x = element_text(hjust = 0, size = 14),
             axis.title.y = element_text(hjust = 0, size = 14)))
}

DefaultAssay(seu) <- "RNA"

ct <- as.character(seu$annotation)
ct[ct == "Immune like EC"] <- "Monocytes"
seu$annotation <- ct
seu$annotation_plot <- droplevels(factor(ct, levels = union(names(celltype_cols), sort(unique(ct)))))

# EC subset
ec <- subset(seu, cells = colnames(seu)[seu$annotation %in% ec_subtypes])
ec <- NormalizeData(ec, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

dis <- as.character(ec$condition)
condition <- rep(NA_character_, length(dis))
condition[dis %in% ctrl_values] <- "Ctrl"
condition[dis %in% hfref_values] <- "HFrEF"
stopifnot(!anyNA(condition))
ec$condition <- factor(condition, levels = c("Ctrl", "HFrEF"))
ec$sample_id <- paste(ec$dataset, ec$orig.ident, sep = "::")

ec <- RunUMAP(ec, reduction = "scANVI", dims = 1:30, reduction.name = "umap_ec", seed.use = 1234, verbose = FALSE)
ec <- FindNeighbors(ec, reduction = "scANVI", dims = 1:30, verbose = FALSE)
ec <- FindClusters(ec, resolution = 0.5, algorithm = 4, random.seed = 1234, verbose = FALSE)
ec <- subset(ec, cells = colnames(ec)[!(as.character(ec$seurat_clusters) %in% exclude_fine_clusters)])
ec <- RunUMAP(ec, reduction = "scANVI", dims = 1:30, reduction.name = "umap_ec", seed.use = 1234, verbose = FALSE)

ec_subtypes_present <- intersect(ec_subtypes, unique(ec$annotation))
ec$subcluster <- factor(ec$annotation, levels = ec_subtypes_present)

# Fig. 3f
pf <- rasterise_points(
  DimPlot(seu, group.by = "annotation_plot", cols = unname(celltype_cols[levels(seu$annotation_plot)]),
          pt.size = 0.4, raster = FALSE, label = FALSE) +
    theme(legend.position = "right") + labs(color = NULL, title = NULL)
) + umap_label(Embeddings(seu, "umap"), ncol(seu))
pf
save_all(pf, "Fig3f_UMAP_Human_celltypes", 9, 7)

# Fig. 3g
pg <- make_heatmap(scaled_avg(seu, markers, "annotation_plot", celltype_order))
pg
save_all(pg, "Fig3g_Heatmap_Human_celltypes", 6, 7)

# Fig. 3h
ph <- rasterise_points(
  DimPlot(ec, group.by = "subcluster", reduction = "umap_ec", cols = unname(ec_cols[ec_subtypes_present]),
          pt.size = 0.8, raster = FALSE, label = FALSE) +
    theme(legend.position = "right") + labs(color = NULL, title = NULL)
) + umap_label(Embeddings(ec, "umap_ec"), ncol(ec))
ph
save_all(ph, "Fig3h_UMAP_Human_EC_subclusters", 8, 6.5)

# Fig. 3i
pi <- make_heatmap(scaled_avg(ec, ec_markers, "subcluster", gsub("_", "-", ec_subtypes_present)))
pi
save_all(pi, "Fig3i_Heatmap_Human_EC_subclusters", 3.5, 6.5)

# Fig. 3j
pj <- DO.Dotplot(sce_object = ec, Feature = target_gene, group.by.x = "subcluster", group.by.y = "condition",
                 sort_x = ec_subtypes_present, stats_y = FALSE, point_stroke = 0.1) +
  theme(text = element_text(family = font),
        axis.text.x = element_text(family = font), axis.text.y = element_text(family = font),
        strip.text.x = element_text(family = font, face = "bold"))
pj
save_all(pj, "Fig3j_Dotplot_CIITA_EC_subclusters", 5.2, 3.6)

# Fig. 3k
cell_df <- data.frame(sample_id = ec$sample_id,
                      dataset = unname(dataset_names[as.character(ec$dataset)]),
                      condition = ec$condition,
                      subcluster = ec$subcluster,
                      expr = LayerData(ec, assay = "RNA", layer = "data")[target_gene, ])

keep_sub <- cell_df %>%
  dplyr::count(sample_id, condition, subcluster, name = "n_cells") %>%
  filter(n_cells >= min_cells, condition == "Ctrl") %>%
  dplyr::count(subcluster) %>%
  filter(n >= min_ctrl_patients) %>%
  pull(subcluster) %>%
  as.character()
keep_sub <- intersect(ec_subtypes_present, keep_sub)

pat <- cell_df %>%
  group_by(sample_id, dataset, condition, subcluster) %>%
  dplyr::summarise(n_cells = n(), mean_log_expression = mean(expr), .groups = "drop") %>%
  filter(n_cells >= min_cells, subcluster %in% keep_sub) %>%
  mutate(subcluster = factor(subcluster, levels = keep_sub),
         dataset = factor(dataset, levels = dataset_names))

stats <- pat %>%
  group_by(subcluster) %>%
  dplyr::summarise(n_ctrl = sum(condition == "Ctrl"), n_HFrEF = sum(condition == "HFrEF"),
            p_wilcox = wilcox.test(mean_log_expression ~ condition, exact = FALSE)$p.value,
            .groups = "drop") %>%
  mutate(q_wilcox_BH = p.adjust(p_wilcox, "BH"))

write.csv(pat, file.path(out_dir, "Fig3k_CIITA_patient_values.csv"), row.names = FALSE)
write.csv(stats, file.path(out_dir, "Fig3k_CIITA_patient_level_stats.csv"), row.names = FALSE)

fmt_q <- function(x) ifelse(x < 0.001, "q < 0.001", paste("q =", formatC(signif(x, 3), format = "fg")))

dat <- pat %>%
  group_by(subcluster, condition) %>%
  mutate(x_label = paste0(condition, "\n(n = ", n(), ")")) %>%
  ungroup() %>%
  arrange(condition) %>%
  mutate(x_label = factor(x_label, levels = unique(x_label)))

br <- dat %>%
  group_by(subcluster) %>%
  dplyr::summarise(top = max(mean_log_expression), .groups = "drop") %>%
  mutate(top = ifelse(top > 0, top, 0.01), y = top * 1.08, tick = top * 0.03) %>%
  left_join(stats, by = "subcluster") %>%
  mutate(label = fmt_q(q_wilcox_BH))

pk <- ggplot(dat, aes(x_label, mean_log_expression)) +
  geom_boxplot(aes(fill = condition), colour = "black", width = .6, outlier.shape = NA, linewidth = .5) +
  geom_point(aes(shape = dataset), size = 1.8, colour = "black", fill = "white", stroke = .6,
             position = position_jitter(width = .13, height = 0, seed = 4)) +
  geom_segment(data = br, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = .5) +
  geom_segment(data = br, aes(x = 1, xend = 1, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = .5) +
  geom_segment(data = br, aes(x = 2, xend = 2, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = .5) +
  geom_text(data = br, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE,
            vjust = -0.5, size = 3.6, family = font) +
  facet_wrap(~subcluster, nrow = 1, scales = "free") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .15))) +
  scale_fill_manual(values = paper_cols, guide = "none") +
  scale_shape_manual(values = c(21, 24, 22, 23), drop = FALSE) +
  labs(x = NULL, y = paste0(target_gene, " expression\n(mean log-normalized counts)"), shape = "Dataset") +
  theme_classic(base_size = 12, base_family = font) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 14, family = font),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(face = "bold"),
        axis.ticks = element_line(colour = "black"),
        panel.spacing = unit(1.2, "lines"))
pk
save_all(pk, "Fig3k_CIITA_per_patient_Ctrl_vs_HFrEF", 2.3 * length(keep_sub) + 1.5, 4.2)

# Suppl. Fig. 2f
dataset_labels <- c(chaffin = "Chaffin et al., 2022", dimmeler = "Ruz Jurado et al., 2025",
                    koenig = "Koenig et al., 2022", kuppe = "Kuppe et al., 2022")
dataset_cols <- c(chaffin = "#FF9896", dimmeler = "#98DF8A", koenig = "#1F77B4", kuppe = "#CD661D")
seu$dataset_plot <- factor(unname(dataset_labels[as.character(seu$dataset)]), levels = dataset_labels)

psf <- rasterise_points(
  DimPlot(seu, group.by = "dataset_plot", cols = unname(dataset_cols[names(dataset_labels)]),
          pt.size = 0.4, raster = FALSE, label = FALSE) +
    theme(legend.position = "right", legend.text = element_text(face = "plain")) + labs(color = NULL, title = NULL)
) + umap_label(Embeddings(seu, "umap"), ncol(seu))
save_all(psf, "SupplFig2f_UMAP_Human_dataset", 9, 7)

# Suppl. Fig. 2g
seu <- FindNeighbors(seu, reduction = "scANVI", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, algorithm = 4, random.seed = 1234, verbose = FALSE)
seu$unbiased_cluster <- seu$seurat_clusters
print(table(seu$annotation_plot, seu$unbiased_cluster))

psg <- rasterise_points(
  DimPlot(seu, group.by = "unbiased_cluster", label = TRUE, pt.size = 0.4, raster = FALSE) +
    theme(legend.position = "right") + labs(color = NULL, title = NULL)
) + umap_label(Embeddings(seu, "umap"), ncol(seu))
save_all(psg, "SupplFig2g_UMAP_Human_unbiased_clusters", 9, 7)

# Suppl. Fig. 2h
wholedata_markers <- c(
  "UPK3B", "BNC1", "MSLN", "FCN1", "HLA-DRA", "PTPRC", "KLRD1", "NPPA", "NPR3",
  "LYVE1", "FOLR2", "TREM2", "KIT", "TPSAB1", "TPSB2", "CD1C",
  "ADIPOQ", "PLIN1", "PROX1", "PDPN", "PLP1", "SOX10", "NRG1", "HBB", "HBA1", "HBA2", "CD74",
  "ACKR1", "SELE", "CXCR2", "NKG7", "GNLY", "CD3D", "CD3E", "CD79A", "FCGR3A", "CD68",
  "MYH11", "GJA5", "EFNB2", "PECAM1", "VWF", "KDR", "CDH5", "CD14", "LYZ", "ITGAM",
  "RGS5", "PDGFRB", "TAGLN", "TOP2A", "PDGFRA", "COL1A1", "COL1A2", "TNNT2", "MYL2", "MYH6"
)
wholedata_markers <- intersect(wholedata_markers, rownames(seu))

psh <- DotPlot(seu, features = wholedata_markers, group.by = "unbiased_cluster", cols = c("lightgrey", "red")) +
  RotatedAxis() + labs(title = "Marker panel across all unsupervised clusters")
save_all(psh, "SupplFig2h_Dotplot_Human_unbiased_clusters", 14, 8)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_Fig3_human.txt"))
