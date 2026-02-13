#Import data and libaries
setwd("/Users/tilllautenschlager/OneDrive - Johann Wolfgang Goethe Universität/Uni/PhD/Bioinformatic/23_LongTermAMI/For Publication/")
```

## Human

```
{library(tidyverse)}
library(Seurat)
library(openxlsx2)
library(SCpubr)
library(DOtools)

#Import Human dataset
Human <- readRDS("/Volumes/Storage+/Objects/snRNA_Human.Mice_HFpEF,HFrEF,AS_Mariano/SeuratObject_Human_AS_HFpEF_HFrEF_Ctrl_Annotated_20_10_22.rds")

Human <- UpdateSeuratObject(Human)

Human_ECs <- subset(Human, subset = cell_type == "Endothelial") %>%
    subset(subset = condition %in% c("Human-CTRLlv", "Human-HFrEF"))
DefaultAssay(Human_ECs) <- "RNA"

genes_vector_human <- c(
  "SMURF1","SMAD3","CIITA", "PLEKHG5",  "NDST1",  "SAP30BP", "LRRTM3"
)

htmap_human <- do_ExpressionHeatmap(Human_ECs, genes_vector_human, group.by = "condition", diverging.palette = "Spectral", sequential.palette = "Reds", flip = T, features.order = rev(genes_vector_human), legend.width = 1,legend.length = 8)

ggsave( filename = "Heatmap_expression_Human.pdf", plot = htmap_human, path = "./Figure 1/", units = "cm", limitsize = FALSE)
```

## Mouse

```
{library(tidyverse)}
library(Seurat)
library(openxlsx2)
library(SCpubr)
library(DOtools)

#Import snRNA Seq. Mouse dataset
Mice <- readRDS("/Volumes/Storage+/Objects/snRNA_Human.Mice_HFpEF,HFrEF,AS_Mariano/240712_SeuratObject_Mouse_ECs_Ctrl_HFpEF_HFrEF.rds")

Mice_ECs <- subset(Mice, subset = cell_type == "Endothelial") %>%
    subset(subset = disease %in% c("CTRL", "HFrEF"))

genes_vector_mouse <- c("Smurf1", "Smad3", "Ciita","Plekhg5","Ndst1", "Sap30bp", "Lrrtm3")

DefaultAssay(Mice_ECs) <- "RNA"

htmap_mouse <- do_ExpressionHeatmap(
  Mice_ECs,
  genes_vector_mouse,
  group.by = "condition",
  sequential.palette = "Reds",
  flip = T,
  features.order = rev(genes_vector_mouse),
  legend.width = 1,
  legend.length = 8
)

htmap_mouse <- do_ExpressionHeatmap(
  Mice_ECs,
  overlap_genes,
  group.by = "condition",
  min.cutoff = 0.1,
  sequential.palette = "Reds",
  flip = F,
  legend.width = 1,
  legend.length = 8
)

htmap_mouse

ggsave( filename = "Heatmap_expression_Mouse.pdf", plot = htmap_mouse, path = "./Figure 1/", units = "cm", limitsize = FALSE)
```

# 
