library(Seurat)
library(tidyverse)
library(DOtools)
library(scDblFinder)
library(decontX)
library(SingleCellExperiment)
library(pheatmap)
library(ggplotify)

options(future.globals.maxSize = 16 * 1024^3)
set.seed(42)

input_rds <- "SeuratObject_Mice_HFrEF_Ctrl.rds"
out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- readRDS(input_rds) %>% subset(condition %in% c("Mice-CTRL", "Mice-HFrEF"))
data$orig.ident <- droplevels(data$orig.ident)
DefaultAssay(data) <- "RNA"

sce <- as.SingleCellExperiment(data, assay = "RNA")
colnames(sce) <- colnames(data)
sce <- decontX(sce)
sce <- scDblFinder(sce, samples = sce$orig.ident)

data$decontX_contamination <- colData(sce)$decontX_contamination
data$scDblFinder_class <- colData(sce)$scDblFinder.class
data[["decontX"]] <- CreateAssayObject(counts = decontXcounts(sce))
DefaultAssay(data) <- "decontX"
data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)

data_filtered <- subset(data, subset = scDblFinder_class == "singlet" & decontX_contamination < 0.2)

# Suppl. Fig. 2a
svg(file.path(out_dir, "SupplFig2a_DecontX_hist.svg"), width = 6, height = 5)
hist(data$decontX_contamination, breaks = 50,
     main = "DecontX Contamination Scores", xlab = "Contamination Fraction")
abline(v = 0.2, col = "red", lty = 2, lwd = 2)
dev.off()

# Suppl. Fig. 2b
svg(file.path(out_dir, "SupplFig2b_DecontX_boxplot.svg"), width = 8, height = 8)
boxplot(decontX_contamination ~ orig.ident, data = data@meta.data, las = 2,
        main = "Contamination per Sample", xlab = "", ylab = "decontX_contamination")
abline(h = 0.2, col = "red", lty = 2, lwd = 2)
dev.off()

# Suppl. Fig. 2c
doublets <- as.data.frame(table(sample = data$orig.ident, class = data$scDblFinder_class))
write.csv(doublets, file.path(out_dir, "SupplFig2c_Doublet_counts.csv"), row.names = FALSE)

ggplot(doublets, aes(sample, Freq, fill = factor(class, levels = c("doublet", "singlet")))) +
  geom_col(position = "fill", colour = "black", width = 0.8) +
  scale_fill_manual(values = c(doublet = "red", singlet = "grey80")) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(title = "Doublet Frequency", x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("SupplFig2c_Doublet_frequency.pdf", path = out_dir, width = 8, height = 5)

# Fig. 3a
DimPlot(data_filtered, group.by = "annotation_final") + labs(title = NULL)
ggsave("Fig3a_UMAP_Mice_Allcells.pdf", path = out_dir, width = 8, height = 8)

# Fig. 3b
markers <- c(
  "Col1a1", "Col1a2", "Adipoq", "Plin1", "Tnnt2", "Myh6", "Npr3", "Nrg1", "Pecam1", "Vwf",
  "Cd74", "H2-Aa", "Runx1", "Ptprc", "Kit", "Kdr", "Gja5", "Efnb2", "Lyve1", "Prox1",
  "Wt1", "Msln", "Plp1", "Sox10", "Rgs5", "Pdgfrb", "Myh11", "Tagln", "Cd3e", "Cd3d"
)
cell_order <- c(
  "FB", "Adip", "CM", "VeinEC", "MP-recruited", "MP-resident", "IMEC", "CapEC",
  "ArtEC", "LymphEC", "Meso", "NC", "PC", "SMC", "B-cells"
)

avg.exp <- AverageExpression(data_filtered, features = markers, group.by = "annotation_final",
                             assays = DefaultAssay(data_filtered))[[1]]
avg.exp_scaled <- t(scale(t(avg.exp)))
avg.exp_scaled <- avg.exp_scaled[markers[markers %in% rownames(avg.exp)],
                                 cell_order[cell_order %in% colnames(avg.exp)], drop = TRUE]
avg.exp_scaled[avg.exp_scaled > 2] <- 2
avg.exp_scaled[avg.exp_scaled < -2] <- -2

as.ggplot(pheatmap(
  avg.exp_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#083269", "#8DC1DA", "#D4E6EE", "white",
                             "#F9DCCC", "#EC9977", "#BF0027"))(100),
  breaks = seq(-2, 2, length.out = 101),
  border_color = "grey80",
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 90
))
ggsave("Fig3b_Heatmap_Mouse_Allcells.pdf", path = out_dir, width = 4)

# Fig. 3c
endo_immune <- subset(
  data_filtered,
  subset = annotation_final %in% c("ArtEC", "B_cells", "CapEC", "LymphEC",
                                   "MP_recruited", "MP_resident", "VeinEC", "IMEC")
)

imec_markers <- c(
  "Cdh5", "Pecam1", "Kdr", "Esam", "Tie1", "Eng", "Cd36", "Fabp4", "Flt1", "Emcn", "Cldn5", "Lyve1",
  "Ptprc", "Runx1", "Lcp1", "Cd44", "Irf8", "Cd68", "Spi1", "Cxcl16", "Adgre5", "Cd83",
  "Ciita", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1", "H2-DMb2", "H2-Ob"
)
celltype_order <- c("CapEC", "ArtEC", "VeinEC", "LymphEC", "IMEC", "MP_resident", "MP_recruited", "B_cells")

DO.Dotplot(endo_immune, imec_markers, group.by.x = "annotation_final",
           sort_x = imec_markers, sort_y = rev(celltype_order), scale_gene = TRUE)
ggsave("Fig3c_DotPlot_IMECs_EndoImmune.pdf", path = out_dir, width = 11, height = 5)

# Fig. 3d
cell_numbers <- as.data.frame.matrix(table(data_filtered$annotation_final, data_filtered$orig.ident))
ec_types <- c("ArtEC", "CapEC", "VeinEC", "LymphEC", "IMEC")
imec_fraction <- tibble(
  sample = colnames(cell_numbers),
  condition = if_else(grepl("HFrEF", colnames(cell_numbers)), "MI 28", "Control"),
  IMEC = unlist(cell_numbers["IMEC", ]),
  EC_total = colSums(cell_numbers[ec_types, ]),
  cells_total = colSums(cell_numbers)
) %>%
  mutate(IMEC_per_EC_pct = 100 * IMEC / EC_total,
         IMEC_per_cells_pct = 100 * IMEC / cells_total)

write.csv(cell_numbers, file.path(out_dir, "Fig3d_Cell_numbers_per_sample.csv"))
write.csv(imec_fraction, file.path(out_dir, "Fig3d_IMEC_fractions.csv"), row.names = FALSE)
