#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File:    Analysis
Author:  davidr
Created: 24.09.26 21:05

Quality control and processing of the dataset
"""

from prelude_py import do, plt, sns

main_path = do.utility.set_path(
    path=#"/Volumes/Sequencing_Storage/20260924_DRM006_snRNA_26056_mouse_heart_shCiitaVsMock_Till"
    "/mnt/davidr/Sequencing_Storage/20260924_DRM006_snRNA_26056_mouse_heart_shCiitaVsMock_Till"
)

results_path = do.utility.set_path(path=main_path / "results")
object_path = do.utility.set_path(path=results_path / "Objects")
figure_path = do.utility.set_path(path=results_path / "Figures")
raw_data_path = do.utility.set_path(path=main_path / "raw_data" / "cellranger")


# ! Quality Control
# region
samples_path = [
    str(raw_data_path / "Mock_1/outs/filtered_feature_bc_matrix.h5"),
    str(raw_data_path / "Mock_2/outs/filtered_feature_bc_matrix.h5"),

    str(raw_data_path / "shCiita_1/outs/filtered_feature_bc_matrix.h5"),
    str(raw_data_path / "shCiita_2/outs/filtered_feature_bc_matrix.h5"),
]
metadata = {
    "condition": ["Mock", "Mock", "shCiita", "shCiita"],
    "sex": ["female", "female", "female", "female"],
    "strain": ["C57BL6", "C57BL6", "C57BL6", "C57BL6"],
    "tissue": ["heart", "heart", "heart", "heart"],
    "kit_10x":["3v4", "3v4", "3v4", "3v4"],
    "modality":["snrna", "snrna", "snrna", "snrna"],
    "intervention": ["MI", "MI", "MI", "MI"]
}

adata = do.pp.importer_py(
    paths=samples_path,
    ids=["Mock1", "Mock2", "shCiita1", "shCiita2"],
    metadata=metadata,
    min_genes_in_cell=300,
    min_cells_with_genes=5,
    cut_mt=20,
    min_counts=250,
    n_reads=10_000,
    high_quantile=98,
    remove_doublets=True,
    doublet_tool="scDblFinder",
    normalisation_method="LogNormalisation"
)
adata.write(object_path / "snRNA_heart_mouse_shCiitaVsMock.h5ad")
adata = do.io.read_h5ad(object_path / "snRNA_heart_mouse_shCiitaVsMock.h5ad")
# endregion

# ! Integration
# region
do.tl.integrate_data(
    adata=adata,
    batch_key="batch",
    integration_method="harmony",
    resolution=0.5,
)

# * Check Integration
do.pl.umap(adata, color="leiden", labels="leiden", figsize=(6, 5), path=figure_path, filename="UMAP_leiden.pdf")
do.pl.split_embedding(adata, split_by="batch", figsize=(6, 5), ncols=2, path=figure_path, filename="UMAPSplit_batch.pdf")

adata.write(object_path / "snRNA_heart_mouse_shCiitaVsMock.h5ad")
# endregion

# ! Annotation
# region
do.tl.auto_annot(adata, cluster_key="leiden")
do.pl.umap(adata, color="autoAnnot", labels="leiden", figsize=(6, 5), path=figure_path, filename="UMAP_autoAnnot.pdf")

markers = {
    "CM": ["Tnni3", "Ttn"],
    "ECs":["Cdh5", "Sema3g", "Npr3", "Rgcc", "Vcam1", "Lyve1", "Prox1"],
    "ICs": ["Ptprc", "Cd3d", "Cd4", "Cd8a", "Kit", "C1qb", "Ccr2", "Bank1", "Igkc", "Themis", "Cd163", "Lyz2"],
    "Mural":["Rgs5", "Pdgfrb", "Acta2"],
    "Fibroblasts": ["Pdgfra", "Postn", "Negr1"],
    "Misc": ["Pax3", "Cdh19", "Tubb3", "Adipoq", "C3", "Slc39a8"]
}

do.pl.dotplot(adata, x_axis="leiden", features=markers, path=figure_path, filename="Dotplot_markers_on_leiden.pdf")

adata.obs["annotation"] = adata.obs["leiden"].map({
    "0":"EndoEC",
    "1":"ArtEC",
    "2":"CapEC",
    "3":"Fibro_activ",
    "4":"Fibroblasts",
    "5":"CM",
    "6":"MP",
    "7":"Mast",
    "8":"LymphEC",
    "9":"MP",
    "10":"Pericytes",
    "11":"CM",
    "12":"MP_recruited",
    "13":"T_cells",
    "14":"Epi_cells",
    "15":"SMC",
    "16":"MP",
    "17":"B_cells",
    "18":"Adip",
})
do.pl.umap(adata, color="annotation", figsize=(6, 5), path=figure_path, filename="UMAP_annotation.pdf")


markers = {
    "CM": ["Tnni3"],
    "ECs":["Cdh5", "Sema3g", "Npr3", "Prox1", "Rgcc", ],
    "ICs": ["Ptprc", "Cd3d", "Themis","Cd4", "Bank1", "Igkc",  "Lyz2", "Cd163", "Ccr2", "Kit"],
    "Mural":["Rgs5", "Pdgfrb", "Acta2"],
    "Fibro": ["Pdgfra", "Negr1", "Postn",],
    "Misc": ["Adipoq", "Slc39a8"]
}

do.pl.dotplot(
    adata, x_axis="annotation",
    features=markers, path=figure_path,
    x_categories_order=["CM", "ArtEC", "EndoEC", "LymphEC", "CapEC",
                        "T_cells", "B_cells", "MP", "MP_recruited", "Mast",
                        "Pericytes", "SMC",
                        "Fibroblasts", "Fibro_activ",
                        "Adip", "Epi_cells"],
    filename="Dotplot_markers_on_annotation.pdf"
)



# endregion

# ! Save Objects
# region
do.io.save_rds(
    path_rds=object_path / "snRNA_heart_mouse_shCiitaVsMock.rds",
    adata=adata
)
adata.write(object_path / "snRNA_heart_mouse_shCiitaVsMock.h5ad")
# endregion


# ! Project Requests
# ? Ciita expression in ECs
# region
ecs = do.get.subset(adata, obs_key="annotation", obs_groups=["CapEC", "ArtEC", "EndoEC"])
do.pl.barplot(
    ecs, x_axis="annotation", feature="Ciita", figsize=(4, 5), hue="condition",
    edgecolor="k", linewidth=1.5, path=figure_path, filename="Barplot_ECs_Ciita.pdf"
)

do.pl.boxplot(
    ecs, x_axis="annotation", feature="Ciita", figsize=(4, 5), hue="condition",
     path=figure_path, filename="Boxplot_ECs_Ciita.pdf",scatter=True,
    batch_key="batch", pseudobulk=True
)

# endregion

# ? Do we see changes in immune cells?
# region
do.pl.cell_composition(
    adata, annot_key="annotation", condition_key="condition", batch_key="batch",
    figsize=(6, 5), xticks_rotation=45, path=figure_path, filename="CellComposition.pdf"
)
# endregion























