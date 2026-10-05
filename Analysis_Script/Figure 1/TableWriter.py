from pybedtools import BedTool
import gzip
import pandas as pd
import json
import os
import random
import re
import fnmatch
from collections import Counter
import subprocess
from itertools import chain
from timeit import default_timer as clock
import HelperScripts.GTF_Processing as GTF_Processing
random.seed(1234)

"""Create joint information files, one for the interactions and one on gene level."""


tag = '1026_RevisionReproduce'
out_path = "InfoFiles/"
plot_out = 'plots/' + tag
info_dict = {"out_tag": tag,
            'abc_folder': "ABC_subEC/",
            'enhancer_file': "AMI-scATAC/ATACPeaks_subEC_CPM.bed",
            'diff_atac_folder': "DESeq2_subEC_ATAC/",
            'diff_atac_pattern': "Diff_*_AMI_vs_control.bed",
            'interaction_out': out_path + "MouseAMI_InteractionInfo_"+tag+".txt",
            "gene_info_out": out_path + "MouseAMI_GeneInfo_"+tag+".txt",
            'rna_files': { # CARE the log2FC will be flipped
            'HFrEF': '/RNA/HFrEF/DEGs_Mouse_snRNA_ECs_CTRL.vs.HFrEF.csv'},
            'annotation': "base_data/mouse/gencode.vM27.annotation.gtf.gz"}
open(out_path + info_dict['out_tag'] + '_Meta_JSON.txt', 'w').write(json.dumps(info_dict))

promoter_bed = GTF_Processing.gene_window_bed(gtf_file=info_dict['annotation'], extend=200, tss_type='all', merge=True)
promoter_bed = BedTool('\n'.join(['\t'.join(x.fields[:4]) for x in promoter_bed]), from_string=True)
gene_name_map = {g.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].split('.')[0]:
                 g.strip().split('\t')[8].split('gene_name "')[-1].split('";')[0] for g in gzip.open(info_dict['annotation'], 'rt').readlines() if not g.startswith('#') and g.split('\t')[2] == 'gene'}

gene_dict_genes = list(gene_name_map.keys())

# Already get the Ensembl IDs for the RNA files.
# And grab the genes from the HFrEF RNA table, the gene_df will be limited to it.
rna_dfs = {}
for rna_tag, rna_file in info_dict['rna_files'].items():
    rna_df = pd.read_table(rna_file, sep=',')
    name_hits, missed = GTF_Processing.match_gene_identifiers(rna_df['gene_name'].values, gtf_file=info_dict['annotation'],
                                                              species='mouse', scopes="symbol,alias,uniprot",
                                                              fields="ensembl,symbol", ensemblonly=True)
    rna_df['Ensembl ID'] = [None if s not in name_hits else name_hits[s]['ensembl'] for s in rna_df['gene_name'].values]
    rna_df = rna_df.drop_duplicates(subset=['Ensembl ID'], keep='first')
    rna_dfs[rna_tag] = rna_df


# Get it from the header of the enhancer file we already generated.
conditions = open(info_dict['enhancer_file']).readline().strip().split('\t')[4:]
main_condition_labels = ['AMI', 'control']

re_pattern = re.compile(info_dict['diff_atac_pattern'].replace('*', '(.*)'))
diff_atac_files = {info_dict['diff_atac_folder'] + '/' + x: re_pattern.search(x).group(1) + '_AMI_vs_control'
                   for x in os.listdir(info_dict['diff_atac_folder']) if fnmatch.fnmatch(x, info_dict['diff_atac_pattern'])}

enhancer_bed = BedTool(info_dict['enhancer_file'])
enhancer_dict = {'\t'.join([str(y) for y in k]): vals for k, vals in
                 pd.read_table(info_dict['enhancer_file'], header=0, dtype=str, index_col=['#chr', 'start', 'end']).to_dict(orient="index").items()}
for enh in enhancer_dict:
    enhancer_dict[enh].update({'ovPromoter': set()})

# ---------------------------------------------------------------------------------------------------
# Interaction-Info
# ---------------------------------------------------------------------------------------------------
# We write one file pooling all interactions from all cell types and conditions with the score if the interaction
# was present and 0 otherwise. Similarly, we add columns with the logFC if the enhancer intersected a differential
# ATAC peak, otherwise 0. Also add overlap with any 5' promoter.
# We collect all enhancer-centric features for each enhancer, for that we gather all from the interaction files.
interaction_files = {info_dict['abc_folder'] + x: x.split('ABCpp_scoredInteractions_')[-1].split('.')[0]
                     for x in os.listdir(info_dict['abc_folder']) if "ABCpp_scoredInteractions" in x and x.endswith(".gz")}

interaction_set = {}
start = clock()
for file in interaction_files:
    inter_head = {x: i for i, x in enumerate(gzip.open(file, 'rt').readline().strip().split('\t'))}
    for entry in [x.strip().split('\t') for x in gzip.open(file, 'rt').readlines()[1:]]:
        this_enhancer = '\t'.join([entry[inter_head[c]] for c in ['#chr', 'start', 'end']])
        this_inter = entry[inter_head['Ensembl ID']] + '#' + this_enhancer
        if this_inter not in interaction_set:
            interaction_set[this_inter] = {}
        interaction_set[this_inter][interaction_files[file]] = entry[inter_head['ABC-Score']]
print(clock() - start)


# Store whether an enhancer was differential in any condition.
for diff_file, comp in diff_atac_files.items():
    diff_head = {x: i for i, x in enumerate(open(diff_file).readline().strip().split('\t'))}
    for entry in open(diff_file).readlines()[1:]:
        entry = entry.strip().split('\t')
        if len(entry) == len(diff_head):  # Some lines are not complete.
            enhancer_dict['\t'.join([entry[diff_head[c]] for c in ['#chr', 'start', 'end']])][comp] = float(entry[diff_head['log2FC']])
        else:
            enhancer_dict['\t'.join([entry[diff_head[c]] for c in ['#chr', 'start', 'end']])][comp] = 0

# Get the enhancer intersection with 5' promoters.
promoter_inter = enhancer_bed.intersect(promoter_bed, wo=True)
for inter in promoter_inter:
    enhancer_dict['\t'.join(inter.fields[:3])]['ovPromoter'].add(inter.fields[-2])

# Define which columns we want to have in the output file.
abc_header = ['Peak_Chr', 'Peak_Start', 'Peak_End', 'Upeak id', 'Ensembl ID', 'Gene Name', 'ovPromoter'] + \
             ['ATAC:CPM:'+c for c in conditions] + list(chain(*[['ATAC:log2FC:'+c] for c in list(diff_atac_files.values())])) + \
             ["ABC-Score:" + c for c in conditions]

with open(info_dict['interaction_out'], 'w') as output:
    output.write('\t'.join(abc_header) + '\n')
    for inter, vals in interaction_set.items():
        enh = 'chr' + inter.split('#')[1]
        if enh in enhancer_dict:
            enh_vals = enhancer_dict[enh]
            output.write(enh + '\t' + enh_vals['peak_id'] + '\t' + inter.split('#')[0].split('.')[0] + '\t' +
                         gene_name_map[inter.split('#')[0].split('.')[0]])
            output.write('\t' + (','.join(enh_vals['ovPromoter']) if enh_vals['ovPromoter'] else '-'))
            output.write('\t' + '\t'.join([enh_vals[c] for c in conditions]))
            # Add the log2FC and if the CRE is differential.
            for comp in diff_atac_files.values():
                output.write('\t' + str(enh_vals[comp]))

            # Add the ABC-scores.
            output.write('\t' + '\t'.join(['0' if c not in vals else vals[c] for c in conditions]))
            output.write('\n')

subprocess.call("gzip -f " + info_dict['interaction_out'], shell=True)

# ---------------------------------------------------------------------------------------------------
# Gene-Info
# ---------------------------------------------------------------------------------------------------
# Base it on the genes that have an open promoter.
gene_df = pd.DataFrame()
gene_df['Ensembl ID'] = list(gene_dict_genes)
gene_df['Gene name'] = [gene_name_map[g] for g in gene_df['Ensembl ID']]
promoter_bed = BedTool(''.join([str(x) for x in promoter_bed if x.fields[3] in gene_dict_genes]), from_string=True)

# Get the max openness across TSS per condition.
max_openness = {g: {c: 0 for c in conditions} for g in gene_dict_genes}
enhancer_inter = promoter_bed.intersect(info_dict['enhancer_file'], wo=True)
for inter in enhancer_inter:
    for c in conditions:
        if float(enhancer_dict['\t'.join(inter.fields[4:7])][c]) > max_openness[inter.fields[3]][c]:
            max_openness[inter.fields[3]][c] = float(enhancer_dict['\t'.join(inter.fields[4:7])][c])

# Get the log2FC of all peaks overlapping promoters of a gene.
promoter_changes = {g: {comp: [] for comp in diff_atac_files.values()} for g in gene_dict_genes}
# Restructure the intersect object to keep only the intersection for a promoter with the maximum overlap.
inter_df = pd.DataFrame([x.fields for x in enhancer_inter])
inter_df.index = ['\t'.join(x) for x in inter_df[[0, 1, 2, 3]].values]
inter_df = inter_df.sort_values(inter_df.columns[-1], ascending=False)
inter_df = inter_df[~inter_df.index.duplicated(keep='first')]
for i, inter in inter_df.iterrows():
    for comp in diff_atac_files.values():
        promoter_changes[inter[3]][comp].append(enhancer_dict['\t'.join(inter[4:7])][comp])

for c in conditions:
    gene_df['ATAC:MaxPromCPM:' + c] = [max_openness[g][c] for g in gene_dict_genes]
for comp in diff_atac_files.values():
    gene_df['ATAC:PromoterPeakLogFC:'+comp] = [','.join([str(x) for x in promoter_changes[g][comp]]) for g in gene_dict_genes]
    gene_df['ATAC:PromoterMaxPeakabsLogFC:'+comp] = [None if not promoter_changes[g][comp] else max(promoter_changes[g][comp], key=abs) for g in gene_dict_genes]

# Collect the interaction information on gene level.
interaction_df = pd.read_table(info_dict['interaction_out']+".gz", header=0, sep='\t')
# CRE counts per condition and cell type.
for abc_cond in [c for c in interaction_df.columns if c.startswith('ABC-Score')]:
    enhancer_counter = Counter(interaction_df[interaction_df[abc_cond] > 0]['Ensembl ID'])
    gene_df['ATAC:Enhancer:'+abc_cond.replace("ABC-Score:", '')] = [0 if g not in enhancer_counter else enhancer_counter[g] for g in gene_df['Ensembl ID']]

for comp in diff_atac_files.values():
    print(comp, 'collecting interaction info per gene')
    # Keep only the interactions matching to the comparison. CARE: assumes strict naming scheme.
    match_cols = main_condition_labels
    match_ct = comp.split('_'+main_condition_labels[0])[0]
    match_interaction_df = interaction_df[(interaction_df[f"ABC-Score:{match_cols[0]}-{match_ct}"] > 0) |
                                          (interaction_df[f"ABC-Score:{match_cols[1]}-{match_ct}"] > 0)]

    for direction in ['Down', 'Up']:
        if direction == 'Down':
            dir_tag = 'Neg'
            all_fc = match_interaction_df[match_interaction_df[f'ATAC:log2FC:{comp}'] < 0].groupby("Ensembl ID")[f'ATAC:log2FC:{comp}'].sum().to_dict()
        else:
            dir_tag = 'Pos'
            all_fc = match_interaction_df[match_interaction_df[f'ATAC:log2FC:{comp}'] > 0].groupby("Ensembl ID")[f'ATAC:log2FC:{comp}'].sum().to_dict()

        gene_df[f'ATAC:sum{dir_tag}LogFCAll:{comp}'] = [0 if g not in all_fc else all_fc[g] for g in gene_df['Ensembl ID']]

    gene_df['ATAC:sumAbsLogFCAll:'+comp] = gene_df['ATAC:sumNegLogFCAll:'+comp].abs() + gene_df['ATAC:sumPosLogFCAll:'+comp].abs()
    gene_df['ATAC:sumLogFCAll:'+comp] = gene_df['ATAC:sumNegLogFCAll:'+comp] + gene_df['ATAC:sumPosLogFCAll:'+comp].abs()


# Add the RNA values, and also if they are differentially expressed.
for rna_tag, rna_df in rna_dfs.items():
    rna_df = rna_df[(rna_df['Ensembl ID'].isin(gene_dict_genes)) & (rna_df['p_val_adj'] != 'NA')].set_index("Ensembl ID")
    rna_df["DE"] = (rna_df['p_val_adj'] <= 0.05) & (rna_df['avg_log2FC'].abs() >= 0.25)
    rna_df['avg_log2FC'] = rna_df['avg_log2FC'] * -1
    rna_df = rna_df[['pct.1', 'pct.2', 'avg_log2FC', 'p_val_adj', 'DE']]
    rna_df.columns = ['RNA:ECs_'+rna_tag+'_'+c for c in rna_df.columns]
    gene_df = gene_df.join(rna_df, on='Ensembl ID', how='left')

gene_df.to_csv(info_dict['gene_info_out'], sep='\t', header=True, index=False)
subprocess.call("gzip -f " + info_dict['gene_info_out'], shell=True)

# ---------------------------------------------------------------------------------------------------
# Subset tables for EC information
# ---------------------------------------------------------------------------------------------------
full_gene_df = pd.read_table(info_dict['gene_info_out'] + ".gz", header=0, sep="\t")
gene_cols = [c for c in full_gene_df.columns if "EC" in c or 'Endothelial' in c or c in ['Ensembl ID', 'Gene name']]
full_gene_df[gene_cols].to_csv(info_dict['gene_info_out'].replace('.txt', '_ECOnly.txt.gz'), header=True, sep='\t', index=False)

full_interaction_df = pd.read_table(info_dict['interaction_out'] + '.gz', header=0, sep='\t')
inter_cols = [c for c in full_interaction_df.columns if ('ATAC' not in c and 'ABC' not in c) or 'Endothelial' in c or 'EC' in c]
ec_abc_cols = [c for c in full_interaction_df.columns if 'ABC' in c and ('Endothelial' in c or 'EC' in c)]
ec_interaction_df = full_interaction_df[(full_interaction_df[ec_abc_cols] > 0).any(axis=1)][inter_cols]
ec_interaction_df.to_csv(info_dict['interaction_out'].replace('.txt', '_ECOnly.txt.gz'), header=True, sep='\t', index=False)

