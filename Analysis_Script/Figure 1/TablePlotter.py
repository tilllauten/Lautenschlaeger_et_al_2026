import numpy as np
import pandas as pd
import json
import copy
import HelperScripts.Heatmaps as Heatmaps
import HelperScripts.GOEnrichment as GOEnrichment

"""Exploratory plotting with trying to visualise the changes across conditions."""

tag = '1026_RevisionReproduce'
info_dict = json.load(open(f'/projects/instituteassociated/work/ami_mouse/InfoFiles/{tag}_Meta_JSON.txt'))
plot_path = '/projects/instituteassociated/work/ami_mouse/plots/' + tag
gene_info_file = info_dict['gene_info_out'].replace(".txt", '_ECOnly.txt.gz')
interaction_info_file = info_dict['interaction_out'].replace(".txt", '_ECOnly.txt.gz')
interaction_folder = info_dict['abc_folder']
diff_atac_folder = info_dict['diff_atac_folder']
annotation = info_dict['annotation']

gene_df = pd.read_table(gene_info_file, header=0, sep="\t")
gene_df.columns = list(gene_df.columns)
gene_name_map = {val[0]: val[1] for val in gene_df[['Ensembl ID', 'Gene name']].values}

interaction_df = pd.read_table(interaction_info_file, header=0, sep="\t")
conditions = sorted(
    [c.split(':')[-1] for c in interaction_df.columns if 'ATAC:CPM:control-' in c or "ATAC:CPM:AMI-" in c],
    key=lambda x: (x.split('-')[1], x.split('-')[0]), reverse=True)
comparisons = sorted([c.split(':')[-1] for c in interaction_df.columns if '_vs_' in c], reverse=True)

ec_tags = ['Endothelial_cells', 'EC0', 'EC1', 'EC2']
ec_conditions = [c for c in conditions if any([e in c for e in ec_tags])]
ec_comparisons = [c for c in comparisons if "Endothelial" in c or 'EC' in c]

hfref_df = gene_df[gene_df['RNA:ECs_HFrEF_avg_log2FC'].notna()]

# -------------------------------------------------------------------------------------------------------------------
# Fig.1 c - GO of the genes with most enhancer changes.
# -------------------------------------------------------------------------------------------------------------------
open_genes_df = gene_df[gene_df['ATAC:PromoterPeakLogFC:Endothelial_cells_AMI_vs_control'].notna()]
# Make a version for EC only where we pool the databases.
ec_go_dict = {'MergedDBs': open_genes_df.set_index("Gene name")[[f'ATAC:sumLogFCAll:Endothelial_cells_AMI_vs_control']]}
ec_gsea_df = GOEnrichment.gsea_prerank(ec_go_dict, weight=0, godf_only=True, max_terms='all', nes_sign='both',
                                       permutation_num=1000,
                                       gmt_path_pattern="/projects/abcpp/work/base_data/GSEA_gmt/mouse//*.v2026.1.Mm.symbols.gmt",
                                       cores=30,
                                       annotation=annotation)

# Change the format so we can have a plot merged across databases.
merged_ec_gsea_df = pd.concat(ec_gsea_df['MergedDBs'].values(), axis=0)
GOEnrichment.plot_go(mode='prerank', wanted_sources=['MergedDBs'],
                     term_fetcher={"ECs": {"MergedDBs": merged_ec_gsea_df}},
                     keywords={}, cmap='plasma', fig_width=None, fig_height=None, go_genes=ec_go_dict,
                     dict_keys=['ECs'],
                     numerate=True, font_s=16, rotation=45, title_tag='Endothelial cells', legend_out=False,
                     max_terms='all',
                     out_tag=plot_path + 'sumLogFCAll_EndothelialCells_MergedDBs', formats=['pdf'])


# -------------------------------------------------------------------------------------------------------------------
# Fig. 1d, Supplementary Fig. 1b&e Plot the genes with the highest summed changes in enhancers.
# -------------------------------------------------------------------------------------------------------------------
def plot_atac_change(top_df, rank_col, rank_tag):
    # Append column for promoter changes.
    top_df = copy.deepcopy(top_df).reset_index()
    for ec in ec_tags:
        max_proms = int(np.nanmax(top_df[f'ATAC:PromoterPeakLogFC:{ec}_AMI_vs_control'].str.count(',')) + 1)
        expanded_proms = pd.DataFrame(top_df[f'ATAC:PromoterPeakLogFC:{ec}_AMI_vs_control'].apply(
            lambda x: map(float, str(x).split(','))).tolist()).rename(
            {i: f'{ec} Promoter {i + 1}' for i in range(max_proms)}, axis=1)
        top_df = pd.concat([top_df, expanded_proms], axis=1)
    cmap_cols = {0: {'cols': [f'ATAC:sumLogFCAll:{ec}_AMI_vs_control' for ec in ec_tags],
                     'centre': 0,
                     'cmap': ['#133E77', '#7CB8DB', '#ffffff', '#F96849', '#6D040F'],
                     'cbar_label': 'summed ATAC logFC'},
                 1: {'cols': ['ATAC:PromoterMaxPeakabsLogFC:Endothelial_cells_AMI_vs_control'],
                     'centre': 0,
                     'cmap': ['#133E77', '#7CB8DB', '#ffffff', '#F96849', '#6D040F'],
                     'cbar_label': 'maximum promoter ATAC logFC'}, }
    for ec in ec_tags:
        cmap_cols[ec] = {'cols': [c for c in top_df if c.startswith(f'{ec} Promoter ')],
                         'centre': 0,
                         'cmap': ['#133E77', '#7CB8DB', '#ffffff', '#F96849', '#6D040F'],
                         'cbar_label': 'ATAC promoter log2FC'}

    Heatmaps.heatmap_cols(top_df, cmap_cols, x_rotation=90,
                          plot_out=plot_path + "_" + rank_col + "_" + rank_tag + str(25),
                          row_label_col="Gene name", class_col=None, x_size=int(10 + len(ec_tags) * 3.5), y_size=10,
                          wspace=1.8)


# Get the top for each comparison with the highest change.
for ec in ec_tags:
    rank_col = f'ATAC:sumLogFCAll:{ec}_AMI_vs_control'
    col_sorted = copy.deepcopy(hfref_df).sort_values(by=rank_col, ascending=False).reset_index()
    for rank_side, top_df in [['neg', col_sorted.iloc[-25:][::-1]], ['pos', col_sorted.iloc[:25]]]:
        plot_atac_change(top_df, rank_col, f"{rank_side}")


