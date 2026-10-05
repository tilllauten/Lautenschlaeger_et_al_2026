import pandas as pd
import pickle
from timeit import default_timer as clock
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

"""DESeq2 for calling differential peaks to get the shrunken log2FC values."""

plot_out = "/plots/DESeq2_subEC_"
pickle_out = "DESeq2_subEC_Tests.pkl"

counts_file = "AMI-scATAC/ATAC_subEC_SampleSums.txt"
counts_df = pd.read_table(counts_file, sep='\t', header=0).set_index(['#chr', 'start', 'end'])
counts_df.index = ['-'.join([str(x) for x in i]) for i in counts_df.index]
counts_df = counts_df.drop('peak_id', axis=1)

cell_results = {}
for cell in set([c.split('-')[1] for c in counts_df.columns]):
    cell_df = counts_df[[c for c in counts_df.columns if cell in c]].T  # Has to be samples * peaks
    # The metadata is samples * metafields, to match the peak counts.
    meta_df = pd.DataFrame([i.split('-') for i in cell_df.index], columns=['condition', 'cell', 'sample']).set_index(cell_df.index)

    start_d = clock()
    dds = DeseqDataSet(
        counts=cell_df,
        metadata=meta_df,
        design_factors="condition",
        refit_cooks=True,
    )
    print(clock() - start_d, 'dds built')

    start_de = clock()
    dds.deseq2()
    print(clock() - start_de, 'DESeq2 run')

    start_stat = clock()
    stat_res = DeseqStats(dds, contrast=['condition', 'AMI', 'control'])
    stat_res.summary()
    print(clock() - start_stat, 'Tests run')

    results = stat_res.results_df
    cell_results[cell] = results
    print(cell, sum(results['padj'] <= 0.05))

pickle.dump(cell_results, open(pickle_out, 'wb'))

# Write the DESeq2 results to output and match the format to the older ones from BN.
cell_results = pickle.load(open(pickle_out, 'rb'))
out_folder = "DESeq2_subEC_ATAC/"
for cell in cell_results:
    cell_df = cell_results[cell]
    cell_df[['#chr', 'start', 'end']] = [x.split('-') for x in cell_df.index]
    cell_df = cell_df[['#chr', 'start', 'end', 'baseMean', 'log2FoldChange', 'pvalue', 'padj']]
    open(out_folder + "Diff_"+cell+"_AMI_vs_control.bed", 'w').write('\t'.join(['#chr', 'start', 'end', 'baseMean', 'log2FC', 'pvalue', 'padj'])+'\n')
    cell_df.to_csv(out_folder + "Diff_"+cell+"_AMI_vs_control.bed", sep='\t', header=False, index=False, mode='a')
