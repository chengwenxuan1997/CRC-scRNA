# CRC-scRNA
Codes for a single cell RNA sequencing data analysis

**CellChat.R**: Infer the cell-cell interaction with [CellChat](https://github.com/sqjin/CellChat), and draw a dotplot to show strength of interaction.  
**cellcycle.R**: Calculate cellcyle score using [Seurat](https://satijalab.org/seurat/articles/cell_cycle_vignette.html), and draw a barplot to describe the cell-cycle compostition for each celltype.  
**cellphoneDB.R**: Infer the cell-cell interaction using [cellphoneDB](https://github.com/Teichlab/cellphonedb) for a compare with CellChat.  
**CIBERSORTx_survival.R**: Prepare the Input for [CIBERSORTx](https://cibersortx.stanford.edu/), and evaluate the prognostic value of cell fraction.  
**inferCNV.R**: Infer the copy number variation using [InferCNV](https://github.com/broadinstitute/infercnv), and draw a heatmap and violin plot to show malignancy of each group.  
**monocle.R**: Do trajectory analysis on CD8T subtypes, to check whether it exists a biological process along the trajectory.  
**SCENIC**: Infer transciption factor activity using [SCENIC](https://pyscenic.readthedocs.io/en/latest/), and show the distribution of TFs AUCell scores.  
**scMetabolism.R**: Infer the metabolism activity using [scMetabolism](https://github.com/wu-yc/scMetabolism), and draw a heatmap.
