import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from matplotlib import cm
from itertools import chain
import ColoursAndShapes
import matplotlib.colors
import matplotlib.gridspec
from matplotlib import cm, colors, colorbar
import pandas as pd
import gzip
from pybedtools import BedTool
import os
import itertools
from Various import sanitize_path


def heatmap_cols(plot_df, cmap_cols, plot_out, row_label_col=None, column_labels=None, class_col=None,
                 x_size=20, y_size=40, title="", annot_cols=None, width_ratios=None, wspace=0.4, rasterized=True,
                 annot_s=10, ticksize=14, heat_ticksize=14, square=False, x_rotation=70, y_rotation=0,
                 ax_fontweight='normal', row_label_first=False, x_label_pos='top', formats=['pdf']):
    """
    Multiple heatmaps side-by-side but the same rows. Allows to show several metrics for the same rows with different
    colourmaps etc. E.g. for a list of top differential genes first a heatmap of baseline expression coloured by TPM,
    followed by a separate heatmap-block with the log2FC for the same genes.

    Args:
        cmap_cols: Dictionary with one entry for each block. The keys don't matter as long as they are unique. E.g.
            {0: {'cols': ['Mean_Control_FM', 'Mean_FM_Mock_Ctrl', 'Mean_Tcf15_FM', 'Mean_FM_Tcf15_OE'],
                     'centre': 0, (optional)
                     'cmap': 'mako',  # Allows a list of colours to create a custom colormap.
                     'cbar_label': 'TPM',
                     'vmax': 200, (optional)
                     'vmin': 0, (optional)
                     }
        row_label_col: Column where to fetch the row-strings from. Set to None to use the index.
        column_labels: Alternative to using the column names as indicated in cmap_cols.
        class_col: Column that should be added as separate first heatmap-block, should be categorical.
        annot_cols: Dictionary of {"column": "column with annotation string"} to write the strings in the value into the cells of columns.
        width_ratios: Ratios of the widths of each heatmap-block.
        wspace: Additional horizontal space between blocks.
        rasterized: Whether to draw thin white lines around cells.
        square: Whether cells should be squares.
        row_label_first: Only write the row names for the first entry and skip for the others.
    """
    all_cols = list(chain(*[c['cols'] for c in cmap_cols.values()]))
    f, axes = plt.subplots(nrows=1, ncols=len(cmap_cols)+bool(class_col), figsize=(x_size, y_size),
                           gridspec_kw={'width_ratios': [0.1] * bool(class_col) + [0.9*len(c['cols'])/len(all_cols) for c in cmap_cols.values()] if not width_ratios else width_ratios})
    # Since we might have different colourmaps we define two matrices for each, one for the values
    # and one for the annotation.
    for n, c_attrs in enumerate(cmap_cols.values()):
        if class_col:
            n += 1
        if len(cmap_cols) == 1 and not class_col:  # If we only had one heatmap we can't index axes.
            this_ax = axes
        else:
            this_ax = axes[n]
        # this_cmap = cm.get_cmap(c_attrs['cmap'])

        if isinstance(c_attrs['cmap'], list):
            c_attrs['cmap'] = matplotlib.colors.LinearSegmentedColormap.from_list("custom_diverging", c_attrs['cmap'])

        value_mat = np.zeros([len(plot_df), len(c_attrs['cols'])])
        annot_mat = np.full([len(plot_df), len(c_attrs['cols'])], '', dtype=object)  # Numpy complains otherwise.
        for c, col in enumerate(c_attrs['cols']):
            value_mat[:, c] = plot_df[col].values
            if annot_cols:
                if col in annot_cols:
                    annot_mat[:, c] = plot_df[annot_cols[col]].values
        heat = sns.heatmap(value_mat, ax=this_ax,  rasterized=rasterized,
                           yticklabels=plot_df.index if not row_label_col else plot_df[row_label_col].values,
                           xticklabels=c_attrs['cols'] if not column_labels else column_labels, cbar=True,
                           cmap=c_attrs['cmap'], fmt='', annot=annot_mat,
                           annot_kws={'size': annot_s}, cbar_kws={'label': c_attrs['cbar_label'], 'shrink': 0.7},
                           center=None if 'centre' not in c_attrs or c_attrs['centre'] is False else c_attrs['centre'],
                           vmin=None if 'vmin' not in c_attrs else c_attrs['vmin'],
                           vmax=None if 'vmax' not in c_attrs else c_attrs['vmax'], square=square)
        if row_label_first and n > 0:
            heat.tick_params(left=False, labelleft=False)
        heat.set_xticklabels(heat.get_xmajorticklabels(), fontsize=ticksize, rotation=x_rotation, fontweight=ax_fontweight)
        heat.set_yticklabels(heat.get_ymajorticklabels(), fontsize=ticksize, rotation=y_rotation, fontweight=ax_fontweight)
        if 'row_labels' in c_attrs and not c_attrs['row_labels']:
                heat.axes.get_yaxis().set_visible(False)
        if x_label_pos == 'top':
            this_ax.tick_params(axis='x', labeltop=True, top=True, labelbottom=False, bottom=False)
        else:
            this_ax.tick_params(axis='x', labeltop=False, top=False, labelbottom=True, bottom=True)
        heat_cbar = heat.collections[0].colorbar
        heat_cbar.ax.tick_params(labelsize=heat_ticksize)
        heat_cbar.ax.yaxis.label.set_fontsize(heat_ticksize)

        if class_col and n == 1:
            # Add the class bar as separate one-column heatmap to the left.
            class_to_int = {c: i for i, c in enumerate(set(plot_df[class_col]))}
            if len(class_to_int) == 2:
                class_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("two_contrast", ColoursAndShapes.two_contrasts[0], N=2)
            else:
                class_cmap = cm.get_cmap("tab20", len(class_to_int))
            sns.heatmap([[class_to_int[c]] for c in plot_df[class_col].values], cmap=class_cmap, ax=axes[0],
                        xticklabels=False, yticklabels=False, square=square,
                        cbar_kws={'label': class_col, 'location': 'left', 'shrink': 1})  # Shrink is somehow ignored here.
            colorbar = axes[0].collections[0].colorbar
            r = colorbar.vmax - colorbar.vmin
            colorbar.set_ticks([colorbar.vmin + r / len(class_to_int) * (0.5 + i) for i in range(len(class_to_int))])
            colorbar.set_ticklabels(list(class_to_int.keys()))
            colorbar.ax.yaxis.set_label_position('left')
            colorbar.ax.yaxis.set_ticks_position('left')
            colorbar.ax.yaxis.label.set_fontsize(14)
            colorbar.ax.tick_params(labelsize=14)
    if title:
        plt.title(title, size=18, fontweight='bold')
    plt.subplots_adjust(wspace=0.4 if not wspace else wspace)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig(sanitize_path(plot_out + "_MultiColHeatmap."+form), bbox_inches='tight', format=form)
    plt.close()


def clustermap(plot_df, columns, row_column, cbar_label, class_col='', class_row='', title="", plot_out="", vmin=None,
               vmax=None, annot_cols=None, cmap='viridis', x_size=12, y_size=10, y_dendro=False, x_dendro=True,
               column_labels=None, row_cluster=True, col_cluster=True, centre=None, tick_size=12, mask=None,
               metric='euclidean', z_score=None, class_col_colour=None, formats=['pdf'], hlines=[], vlines=[],
               main_space=0.82, col_colours=None, x_rotation=0, y_rotation=0):#, return_linkage=False):
    # TODO allow col_colors and also allow to give indices instead of columns.
    """
    Create a heatmap that can be additionally clustered with seaborn. CARE: the class_col_order and class_row parameters
    are not properly tested.

    Args:
        columns: List of the columns, should be present in the plot_df.
        row_column: Column with which rows should be taken, set to 'index' to take the df.index.
        class_col: Add a column into the plot with a categorical value that will be coloured and gets a separate colourbar.
        class_col_colour: Allows a list that will be taken iteratively, or a dict with {label: colour}.
        class_row: Same as class_col but add a row instead.
        annot_cols: Dictionary {col: other-col} to add text into the entries from col taken from other_col.
        cmap: Colourmap for the heatmap, or list of colour from which a custom colormap will be created.
        y_dendro: Whether to plot the dendrogram on y.
        x_dendro: Whether to plot the dendrogram on x.
        column_labels: List that will replace the names from columns if given.
        row_cluster: Whether to cluster the rows.
        col_cluster: Whether to cluster the columns.
        centre: Centre for the colormap, e.g. 0 for bwr.
        mask: Must match the dimensions of the plot_df. If given will not show data where entries are True.
        metric: Metric for clustering for the scipy function, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist.
        z_score: Whether to do z-score normalization before clustering. If None won't do z-scoring, otherwise takes 0 or 1 for the axis along which the normalization should be done.
        hlines: List of indices where white horizontal lines will be added, e.g. for grouping subsets.
        vlines: Same but vertical.
        main_space: How much space the main heatmap will be given. Reduce if the colourbar overlaps the heatmap.
    """
    if z_score is not None:  # Not using seaborns z_score flag as it messes with the separate colourbar.
        plot_part = zscore(plot_df[columns].astype(float), axis=z_score)
    else:
        plot_part = plot_df[columns].astype(float)

    if isinstance(cmap, list):
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_diverging", cmap)

    if class_col and class_row:
        print("ERROR: having colours at both rows and columns is not implemented, because colourbars are fun")
        return
    if class_col:
        class_list = list(reversed(list(dict.fromkeys(plot_df[class_col]))))  # Maintain the list order.
    if class_row:
        class_list = list(reversed(list(dict.fromkeys(plot_df.loc[class_col]))))  # Maintain the list order.
    if class_col or class_row:
        if type(class_col_colour) == dict:
            class_colours = class_col_colour
        else:
            if class_col_colour and type(class_col_colour) == list:
                class_cmap = class_col_colour
            elif len(set(plot_df[class_col if class_col else class_row])) == 2:
                class_cmap = ColoursAndShapes.two_contrasts[0]
            elif len(set(plot_df[class_col if class_col else class_row])) <= 7:
                class_cmap = ColoursAndShapes.tol_vibrant[:len(set(plot_df[class_col if class_col else class_row]))]
            else:
                class_cmap = ColoursAndShapes.glasbey_palettes['glasbey']
            class_colours = {c: class_cmap[i] for i, c in enumerate(class_list)}

    row_colors = None if not class_col else pd.Series([class_colours[cl] for cl in plot_df[class_col].values])
    annot_mat = np.full([len(plot_df), len(columns)], '', dtype=object)  # Numpy complains otherwise.
    if annot_cols:
        for c, col in enumerate(columns):
            if col in annot_cols:
                annot_mat[:, c] = plot_df[annot_cols[col]].values
    if row_column == 'index':
        yticklabels = plot_df.index
    elif row_column:
        yticklabels = plot_df[row_column].values
    else:  # If row_columns is set to None, meaning no labels should be shown.
        yticklabels = False

    clustermap = sns.clustermap(plot_part, vmin=vmin, vmax=vmax, rasterized=True, cmap=cmap, center=centre,
                                figsize=(x_size, y_size), annot=annot_mat, annot_kws={'size': 10}, fmt='',
                                yticklabels=yticklabels,
                                xticklabels=columns if not column_labels else column_labels, row_cluster=row_cluster,
                                col_cluster=col_cluster, row_colors=None if not class_col else pd.Series(
                                    [class_colours[cl] for cl in plot_df[class_col].values], index=plot_df.index),
                                col_colors=None if not col_colours else col_colours,# else pd.Series(
                                #     [class_colors[cl] for cl in plot_df.loc[class_col].values], index=plot_df.columns),
                                dendrogram_ratio=0.1, mask=mask, metric=metric)
    clustermap.ax_col_dendrogram.set_visible(y_dendro)
    clustermap.ax_row_dendrogram.set_visible(x_dendro)
    clustermap.ax_heatmap.tick_params(labelsize=tick_size)
    plt.setp(clustermap.ax_heatmap.get_yticklabels(), rotation=y_rotation)
    plt.setp(clustermap.ax_heatmap.get_xticklabels(), rotation=x_rotation)
    clustermap.cax.set_visible(False)
    clustermap.gs.update(right=main_space)
    # Add a manual colormap to be able to edit it.
    gs2 = matplotlib.gridspec.GridSpec(1, 1, right=0.95, left=0.9, bottom=0.2, top=0.8)
    cbar_ax = clustermap.fig.add_subplot(gs2[0])
    if not vmin:
        vmin = plot_part.min().min()
    if not vmax:
        vmax = plot_part.max().max()
    if centre is not None and centre is not False:
        norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=centre, vmax=vmax)
    else:
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
    if z_score == 0:
        cbar_label = 'column zscore ' + cbar_label
    elif z_score == 1:
        cbar_label = 'row zscore ' + cbar_label
    sep_cbar = colorbar.ColorbarBase(cbar_ax, cmap=cm.get_cmap(cmap), norm=norm, label=cbar_label)
    sep_cbar.ax.yaxis.label.set_fontsize(12)

    if hlines:
        clustermap.ax_heatmap.hlines(hlines, *clustermap.ax_heatmap.get_xlim(), colors='white', linewidth=5)
    if vlines:
        clustermap.ax_heatmap.vlines(vlines, *clustermap.ax_heatmap.get_xlim(), colors='white', linewidth=5)

    if class_col:
        # Add the class bar.
        clustermap.gs.update(left=0.15)
        gs3 = matplotlib.gridspec.GridSpec(1, 1, right=0.04, bottom=0.3, top=0.7)
        class_bar_ax = clustermap.fig.add_subplot(gs3[0])
        cbar_map = colors.ListedColormap([class_colours[c] for c in class_list])
        class_bar_norm = colors.BoundaryNorm([x for x in range(len(class_list) + 1)], cbar_map.N)
        sep_classbar = colorbar.ColorbarBase(class_bar_ax, cmap=cbar_map, norm=class_bar_norm,
                                             ticks=[x + 0.5 for x in range(len(class_list))])
        sep_classbar.ax.set_yticklabels(class_list, size=12)

    clustermap.fig.suptitle(title, y=0.95, size=20, fontweight='bold')
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        clustermap.savefig(sanitize_path(plot_out + "_Clustermap."+form), bbox_inches='tight', format=form)
    plt.close()
    # if return_linkage:
    #     return clustermap.linkage.dendrogram_row or clustermap.linkage.dendrogram_col

# TODO currently not doing what intended?
def interaction_intersection_diff_enhancer(abc_folder, tag_order=None, plot_path='', x_size=16, y_size=10, annot_s=15):
    """Plot an asymmetric heatmap of the ABC-interactions, for cases where we started with different enhancer sets.
    Takes care of cases with multiple overlaps.
    @param abc_folder: Folder from which all files will be taken for the heatmap.
    @param tag_order: Order of the ABCpp tags from the files TAG_ABCpp_scoredInteractions..., if None take them as
    they come from the folder."""

    abc_files = {x.split('_ABCpp_scoredInteractions')[0]: abc_folder+'/'+x for x in os.listdir(abc_folder) if 'ABCpp_scoredInteractions' in x and x.endswith('.gz')}
    if not tag_order:
        tag_order = abc_files.keys()

    enh_gene_interactions = {}
    # Fetch an enhancer_gene dict for each.
    for tag, file in abc_files.items():
        print("reading", tag)
        abc_head = {x: i for i, x in enumerate(gzip.open(file, 'rt').readline().strip().split('\t'))}
        enh_gene = {}
        with gzip.open(file, 'rt') as abc_in:
            abc_in.readline()
            for entry in abc_in:
                enh = '\t'.join([entry.strip().split('\t')[abc_head[c]] for c in ['#chr', 'start', 'end']])
                if enh not in enh_gene:
                    enh_gene[enh] = set()
                enh_gene[enh].add(entry.strip().split('\t')[abc_head['Ensembl ID']])
        enh_gene_interactions[tag] = enh_gene
    number_interactions = {t: sum([len(val) for val in enh_gene_interactions[t].values()]) for t in abc_files.keys()}

    # Now go through each pair of files and count the number of matches.
    # We need to do a pairwise comparison for each, to find which enhancers intersect and if they map to the same gene.
    mat_idx = {t: i for i, t in enumerate(tag_order)}
    shared_frac = np.ones([len(abc_files), len(abc_files)])
    shared_inter = np.zeros([len(abc_files), len(abc_files)], dtype=object)
    for tag in tag_order:
        shared_inter[mat_idx[tag]][mat_idx[tag]] = str(number_interactions[tag])
    tag_pairs = list(itertools.combinations(abc_files.keys(), 2))
    for pair in tag_pairs:
        print(pair)
        pair0_bed = BedTool('\n'.join(enh_gene_interactions[pair[0]]), from_string=True)
        pair1_bed = BedTool('\n'.join(enh_gene_interactions[pair[1]]), from_string=True)
        # First we need to map enhancers to each other, having the enhancers from each unique, to not count multiple
        # overlaps wrongly.
        enh0_inter = {}
        enh1_inter = {}
        for inter in pair0_bed.intersect(pair1_bed, wo=True):
            enh0 = '\t'.join(inter.fields[:3])
            enh1 = '\t'.join(inter.fields[3:6])
            if enh0 not in enh0_inter:
                enh0_inter[enh0] = set()
            enh0_inter[enh0].add(enh1)
            if enh1 not in enh1_inter:
                enh1_inter[enh1] = set()
            enh1_inter[enh1].add(enh0)

        shared0 = 0
        shared1 = 0
        for enh0, other_enh in enh0_inter.items():
            other_enh_genes = set.union(*[enh_gene_interactions[pair[1]][e] for e in other_enh])
            shared0 += len(enh_gene_interactions[pair[0]][enh0] & other_enh_genes)
        for enh1, other_enh in enh1_inter.items():
            other_enh_genes = set.union(*[enh_gene_interactions[pair[0]][e] for e in other_enh])
            shared1 += len(enh_gene_interactions[pair[1]][enh1] & other_enh_genes)
        shared_frac[mat_idx[pair[0]]][mat_idx[pair[1]]] = shared0 / number_interactions[pair[0]]
        shared_frac[mat_idx[pair[1]]][mat_idx[pair[0]]] = shared1 / number_interactions[pair[1]]
        shared_inter[mat_idx[pair[0]]][mat_idx[pair[1]]] = str(shared0)
        shared_inter[mat_idx[pair[1]]][mat_idx[pair[0]]] = str(shared1)

    total_inter = pd.DataFrame([number_interactions[t] for t in tag_order], index=tag_order)

    f, axes = plt.subplots(nrows=1, ncols=3, figsize=(x_size, y_size), gridspec_kw={'width_ratios': [0.2, 0.4, 8]})
    sns.heatmap(total_inter, cmap="Blues", ax=axes[1], xticklabels=False,
                rasterized=True, yticklabels=False, cbar=False, fmt='',
                annot=total_inter, annot_kws={'size': 12})
    cobar = f.colorbar(axes[1].get_children()[0], cax=axes[0], orientation="vertical", label='total number of interactions')
    cobar.ax.yaxis.label.set_fontsize(14)
    cobar.ax.tick_params(labelsize=14)
    cobar.ax.yaxis.set_label_position('left')
    cobar.ax.yaxis.set_ticks_position('left')
    shared_heat = sns.heatmap(shared_frac, cmap="Blues", ax=axes[2], square=True, vmin=0, vmax=1, rasterized=True,
                              yticklabels=tag_order, xticklabels=['shared w/ '+c for c in tag_order],
                              fmt='', annot=shared_inter, annot_kws={'size': annot_s},
                              cbar_kws={'label': "Fraction of interactions shared interactions", 'pad': 0.01})
    shared_heat.set_xticklabels(shared_heat.get_xmajorticklabels(), fontsize=14, rotation=90)
    shared_heat.set_yticklabels(shared_heat.get_ymajorticklabels(), fontsize=14, rotation=0)
    shared_heat_cbar = shared_heat.collections[0].colorbar
    shared_heat_cbar.ax.tick_params(labelsize=14)
    shared_heat_cbar.ax.yaxis.label.set_fontsize(14)
    plt.suptitle("Shared interactions", y=0.93, size=20, fontweight='bold')
    plt.subplots_adjust(wspace=0.02)
    plt.savefig(sanitize_path(plot_path + "_ABCMultiIntersectHeat.pdf"), bbox_inches='tight')
    plt.close()







