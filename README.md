# sc_seurat
Scripts for single-cell 'omics analysis


Full pipeline in seurat 3.0 including
* filter cellranger raw BC matrices
* removes ambient RNA (SoupX)
* removes doublets (doubletFinder)
* merges or integrates data (Seurat 3.0)
* tSNE feature plots and diagnostic QC plots

## Usage
e.g. 

`time Rscript ./seurat_pipeline_seurat3.0.R  --dirs_project_10x "c('/nfsdata/data/sc-10x/data-runs/181119-serup-pancreas/1-5000_cells/','/nfsdata/data/sc-10x/data-runs/181119-serup-pancreas/2-5000_cells/','/nfsdata/data/sc-10x/data-runs/181119-serup-pancreas/3-5000_cells/','/nfsdata/data/sc-10x/data-runs/181119-serup-pancreas/4-5000_cells/','/nfsdata/data/sc-10x/data-runs/181119-serup-pancreas/5-5000_cells/')" --dir_out /projects/jonatan/pub-perslab/181119-serup-pancreas/ --flag_datatype  sc --flag_organism hsapiens --prefix_data serup_panc_2 --prefix_run  seurat_2 --n_cells_loaded 9000 --n_cells_recovered NULL --use_filtered_feature_bc_matrix T --nCount_RNA_min 2000 --nCount_RNA_max  50000 --run_SoupX F   --nFeature_RNA_min 0 --nFeature_RNA_max 25000 --percent.mito_max 1 --percent.ribo_max 1 --rm_sc_multiplets T --vars.to.regress 'c("percent.mito","nCount_RNA")' --nAnchorFeatures  2000 --n_comp 50 --use_jackstraw T --res_primary 1.2 --res_to_calc  "c(0.8,1.2,1.5,2.0,2.5)" --feats_to_plot "c('nCount_RNA','nFeature_RNA','percent.mito','percent.ribo')" --feats_plot_separate T --RAM_Gb_max 250`

## Help

Rscript ./seurat_pipeline_seurat3.0.R --help



