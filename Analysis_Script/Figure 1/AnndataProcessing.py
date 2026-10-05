from matplotlib import pyplot as plt
import pandas as pd
import scanpy as sc
import copy
import HelperScripts.Various as Various
import HelperScripts.Heatmaps as Heatmaps
import HelperScripts.ReadAnndata_Helpers as ReadAnndata_Helpers

tag = 'subEC'
h5_file = "AMI-scATAC/ext442_ML_v6.h5ad"
fragments_file = 'AMI-scATAC/all.fragments_snapatac2.tsv.gz'
peaks_out = f"AMI-scATAC/ATACPeaks_{tag}_CPM.bed"
plots_out = "plots/0626_RevisionReproduce_"
sample_sums_out = f"AMI-scATAC/ATAC_{tag}_SampleSums.txt"
annotation = "base_data/mouse/gencode.vM27.annotation.gtf.gz"
tmp_dir = '/tmpdir/'

# AnnData object with n_obs × n_vars = 71632 × 168364
scdata = sc.read(h5_file)
X = copy.deepcopy(scdata.X)  # CARE if raw or the already normalized one (whatever that means).

clustering = "subEC_leiden0.2"

# ---------------------------------------------------------------------------------------------------------------------
# Add a column with EC subclusters
# ---------------------------------------------------------------------------------------------------------------------
# Clustering of EC subtypes, after normalization.
sc.pp.normalize_total(scdata)
sc.pp.log1p(scdata)
ec_data = scdata[scdata.obs['cell_types_leiden_0.1'] == 'Endothelial cells'].copy()
sc.pp.neighbors(ec_data)

for res in [0.2]:
    # Get the cluster labels for ECs from the sub-clustering, but plot the UMAP
    # with the original cell type labels for all cells.
    res_key = 'subEC_leiden' + str(res)
    if res_key not in scdata.obs.columns:
        sc.tl.leiden(ec_data, n_iterations=12, resolution=res, key_added=res_key)
        ec_data.obs[res_key] = ['EC'+i for i in ec_data.obs[res_key].values]
        scdata.obs = scdata.obs.join(ec_data.obs[[res_key]], how='left')
    # Fill the NAs with the labels from cell_types_leiden_0.1
    scdata.obs[res_key] = scdata.obs[res_key].astype('object').fillna(scdata.obs['cell_types_leiden_0.1'])

    # Reviewer Response b - Make a UMAP with the ECs only (excluding EC3), recomputing the UMAP coordinates.
    non3_ec_data = ec_data[ec_data.obs[res_key] != 'EC3'].copy()
    sc.tl.umap(non3_ec_data)
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.umap(non3_ec_data, color=res_key, show=False, ax=ax)
    fig.savefig(plots_out+"EC_subclustering_leiden"+str(res)+"_ECs_only.pdf", bbox_inches='tight')


# ---------------------------------------------------------------------------------------------------------------------
# Pseudo-bulked activity in peaks
# ---------------------------------------------------------------------------------------------------------------------
# Peaks should be var, they already have a gene mapped to them.
peak_df_out = ReadAnndata_Helpers.peak_pseudobulk_cpm(scobj=scdata, x_matrix=X, cell_col=clustering, condition_col='condition')

# Separately also add the bulk ECs again.
all_pseudo_peak_cpm = ReadAnndata_Helpers.peak_pseudobulk_cpm(scobj=scdata, x_matrix=X, cell_col='cell_types_leiden_0.1', condition_col='condition')
peak_df_out = pd.concat([peak_df_out, all_pseudo_peak_cpm[[c for c in all_pseudo_peak_cpm.columns if 'Endothelial_cells' in c]]], axis=1)
peak_df_out.to_csv(peaks_out, sep="\t", header=True, index=False)

# ---------------------------------------------------------------------------------------------------------------------
# Sum the counts per sample for DESeq2 to shrink the log2FC
# ---------------------------------------------------------------------------------------------------------------------
sample_sums = ReadAnndata_Helpers.sample_pseudobulk_counts(scobj=scdata, cell_col=clustering, condition_col='condition', sample_col='Sample')
# Add bulk ECs again.
all_pseudo_sums = ReadAnndata_Helpers.sample_pseudobulk_counts(scobj=scdata, cell_col='cell_types_leiden_0.1', condition_col='condition', sample_col='Sample')
sample_sums = pd.concat([sample_sums, all_pseudo_sums[[c for c in all_pseudo_sums.columns if 'Endothelial_cells' in c]]], axis=1)
sample_sums.to_csv(sample_sums_out, sep="\t", header=True, index=False)


# ---------------------------------------------------------------------------------------------------------------------
# Reviewer Response c - Marker genes at the EC subcluster
# ---------------------------------------------------------------------------------------------------------------------
pseudo_bulk_bws = Various.fn_patternmatch('AMI-scATAC/split_bigwigs_pseudobulk//output_*.bw')
ec_bws = Various.fn_patternmatch('AMI-scATAC/split_bigwigs_ecs/output_*.bw')
ct_marker_genes = {"Fibroblasts": ["Dcn", 'Pdgfra'],
                "Endothelial cells": ["Tie1", 'Jag2', 'Kdr'],
                'Macrophages': ['Cd14', 'Cd68', 'C1qa'],
                'Cardiomyocytes': ['Ryr2', 'Myh7', 'Ttn'],
                'Pericytes': ['Rgs5', 'Pdgfrb'],
                'T cells': ['Cd8a', 'Cd27', 'Il7r'],
                'B cells': ['Cd19', 'Ms4a1', 'Aff3']
               }


def ec_marker_heatmap(gene_list, list_tag, include_all=False):
    ec_marker_counts = ReadAnndata_Helpers.get_gene_activities(gene_list, {val: k for k, val in ec_bws.items()}, annotation)
    ec_marker_counts = ec_marker_counts.loc[:, ['EC0', 'EC1', 'EC2', 'EC3']]
    all_marker_counts = ReadAnndata_Helpers.get_gene_activities(gene_list, {val: k for k, val in pseudo_bulk_bws.items()}, annotation)
    ec_bw_counts = pd.concat([all_marker_counts[[c for c in all_marker_counts if c!='Endothelial cells']], ec_marker_counts], axis=1)
    if not include_all:
        plot_cols = ['EC0', 'EC1', 'EC2', 'EC3']
    else:
        plot_cols = [c for c in ct_marker_genes.keys() if c != 'Endothelial cells'] + ['EC0', 'EC1', 'EC2', 'EC3']
    Heatmaps.clustermap(ec_bw_counts, columns=plot_cols, row_column='index', cbar_label='z-scored gene activity',
                        title=f"EC marker {list_tag}", plot_out=f"{plots_out}ECMarker{list_tag}_AllCTs{include_all}", centre=0, cmap=['#133E77', '#7CB8DB', '#ffffff', '#F96849', '#6D040F'],
                        x_size=8, y_size=max(int(len(gene_list)*0.9), 5), y_dendro=False, x_dendro=False, row_cluster=False, col_cluster=False,
                        z_score=1, x_rotation=90, y_rotation=0, main_space=0.75)


ec_marker2g = ['Cdh5', 'Pecam1', 'Sema3g', 'Hey1', 'Nr2f2', 'Vcam1', 'Rgcc', 'Kdr', 'Lyve1', 'Mmrn1', 'Vwf', 'Nrg1']
ec_marker_heatmap(ec_marker2g, list_tag='Rev', include_all=True)


