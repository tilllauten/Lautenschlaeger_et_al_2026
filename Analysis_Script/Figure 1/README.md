### Analysis and plotting scripts

Collection of analysis and plotting scripts. Please note that some scripts produce multiple panels that are not exclusive to Figure 1.

- bashcommands.txt: Command how STARE was run to link ATAC peaks to genes.
- TableWriter.py: Script that collects information on interaction and on gene level.
- TablePlotter.py: Uses the tables written above to generate plots.
- AnndataProcessing.py: Processes the h5ad file, produces intermediary files and generates plots.
- DifferentialATAC.py: Getting the log2FC ATAC values in the peaks, uses DESeq2 for log2FC shrinkage.
- Fig. 1b/g - Heatmap.py: Scripts specific for those panels.
- HelperScripts: Collection of scripts with functions that are imported by the scripts described above.
