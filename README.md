# Mito_scRNA

The scripts were used to generage the variant calls in the single cell data. 

varReadDepth.py is used to produce the variant call and base counts of target sites.
singleCellVar.py depends on the output of varReadDepth.py to generate summary of variants identified.
singleCellVar_collapse2pos.R, singleCellVar_plot.R and seuratDataExc.2.R were used for generating visualization of the result.
