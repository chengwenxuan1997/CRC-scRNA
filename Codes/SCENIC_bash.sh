arboreto_with_multiprocessing.py \
  loom.loom \
  cisTarget_databases/hs_hgnc_tfs.txt \
  --method grnboost2 \
  --output output/adj.tsv \
  --num_workers 20 \
  --seed 777
  
pyscenic ctx \
  output/adj.tsv \
  cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
  --annotations_fname cisTarget_databases/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname loom.loom \
  --mode "dask_multiprocessing" \
  --output output/reg.csv \
  --num_workers 20 \
  --mask_dropouts
  
pyscenic aucell \
  loom.loom \
  output/reg.csv \
  --output output/auc_mtx.csv \
  --num_workers 20 \
  --seed 777