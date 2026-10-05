import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read in ATAC Dataset object
adata = sc.read_h5ad('/Volumes/Storage+/Objects/snATAC_Mice_AMI_d28_Till/ext442_ML_v6.h5ad')

adata_AMI = adata[adata.obs["condition"] == "AMI"].copy()

#Define marker genes
genes =["Dcn", "Pdgfra", # Fibroblasts
        'Tie1','Jag2', 'Kdr', # Endothelial
        "Cd14", "Cd68", "C1qa", # Macrophages
        "Ryr2", "Myh7", "Ttn", # Cardiomyoytes
        "Rgs5",'Pdgfrb',  # Pericyte
        "Cd8a",'Cd27', "Il7r", # T cells
        "Cd19", "Cd20", 'Ms4a1', 'Aff3',  # B cells
        ]


genes = [g for g in genes if g in adata.var_names]¯
print("Using genes:", genes)

#Average GeneActivity signal per cell type
df = pd.DataFrame(
    adata[:, genes].X.toarray(),
    index=adata.obs_names,
    columns=genes
)

avg = df.groupby(adata.obs["cell_types_leiden_0.1"]).mean().T

#Z-score per gene (row-wise)
avg_z = (avg - avg.mean(axis=1).values[:, None]) / avg.std(axis=1).values[:, None]

#Plot (seaborn heatmap)
plt.figure(figsize=(4, 6))

sns.heatmap(
    avg_z,
    linewidths=1,
    cmap= "Reds",
    linecolor="black",
    cbar_kws={"label": "Z-score"},
    square=True
)

plt.title("snATAC-seq")

plt.tight_layout()
plt.savefig("Marker_Heatmap.pdf", format="pdf")
plt.show()
