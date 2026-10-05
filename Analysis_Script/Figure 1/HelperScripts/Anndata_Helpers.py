import numpy as np
import GTF_Processing
import BigWig_Counter


def peak_pseudobulk_cpm(scobj, x_matrix, cell_col, condition_col, peak_cols=['peak_chr', 'peak_start', 'peak_end', 'peak_id']):
    """Pseudo-bulked activity in peaks."""

    peak_df = scobj.var
    peak_df_out = peak_df[peak_cols]
    peak_df_out.columns = ["#chr", "start", "end", "peak_id"]

    # For each cell type and condition sum up the vectors of the individual cells and take as CPM.
    for cell_tag in set(scobj.obs[cell_col]):
        for cond in set(scobj.obs[condition_col]):
            sc_indices = np.where((scobj.obs[cell_col] == cell_tag) & (scobj.obs[condition_col] == cond))[0]
            print(cell_tag, cond, len(sc_indices))
            summed_counts = np.asarray(x_matrix[sc_indices, :].sum(axis=0))[0]
            cpms = summed_counts * 10**6 / summed_counts.sum()
            peak_df_out[cond + '-' + cell_tag.replace(" ", "_")] = cpms
        
    return peak_df_out


def sample_pseudobulk_counts(scobj, cell_col, condition_col, sample_col, peak_cols=['peak_chr', 'peak_start', 'peak_end', 'peak_id']):
    """Write a matrix where we sum the raw counts in a peak per sample per cell type for DESeq2."""

    sample_sums = scobj.var
    sample_sums = sample_sums[peak_cols]
    sample_sums.columns = ["#chr", "start", "end", "peak_id"]
    
    if scobj.raw is None:
        print("WARNING: No raw count slot for X, using scobj.X")
        x_mat = scobj.X
    else: 
        x_mat = scobj.raw.X

    for cell_tag in set(scobj.obs[cell_col]):
        for cond in set(scobj.obs[condition_col]):
            for sample in set(scobj.obs[sample_col]):
                sc_indices = np.where((scobj.obs[cell_col] == cell_tag) & (scobj.obs[condition_col] == cond) & (scobj.obs[sample_col] == sample))[0]
                if len(sc_indices) > 0:
                    print(cell_tag, cond, sample, len(sc_indices))
                    # Get the sum across cells per peak.
                    summed_counts = np.asarray(x_mat[sc_indices, :].sum(axis=0))[0]
                    sample_sums[cond + '-' + cell_tag.replace(" ", "_") + '-' + sample] = summed_counts
    return sample_sums


def get_gene_activities(gene_list, bigwig_files, annotation):
    """Separate function, since we also pull them for the comparison to the EC subcluster below."""
    # Similar to Seurat, get 2kb upstream of the TSS plus the gene body.
    marker_gbs_up = GTF_Processing.gene_body_bed(gtf_file=annotation, gene_set=gene_list, add_upstream=2000)
    tss_dict = GTF_Processing.gene_window_bed(gtf_file=annotation, gene_set=gene_list, dict_only=True)  # To get the gene names back.
    region_ecgene_map = {'\t'.join(x.fields[:3]): x.fields[3] for x in marker_gbs_up}

    bw_counts, ctrl_errors = BigWig_Counter.bigwig_counts(bed_file=marker_gbs_up, bigwigs=bigwig_files, n_cores=4)
    bw_counts.index = [tss_dict[region_ecgene_map['\t'.join(x)]]['name'] for x in bw_counts[['#chr', 'start', 'end']].values]
    bw_counts = bw_counts.loc[gene_list, [c for c in bw_counts.columns if c not in ['#chr', 'start', 'end']]]
    bw_counts.columns = [c.replace('_', ' ') for c in bw_counts.columns]

    return bw_counts