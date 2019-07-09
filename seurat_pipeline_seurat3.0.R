# Seurat pipeline
# Usage: 
# export R_MAX_NUM_DLLS=999 # optional
# time Rscript /projects/jonatan/tools/seurat-src/seurat_pipeline.R --dirs_project_10x 'c("/nfsdata/data/sc-10x/data-runs/180511-perslab-immunometab/247L-5000_cells/","/nfsdata/data/sc-10x/data-runs/180511-perslab-immunometab/249L-5000_cells/","/nfsdata/data/sc-10x/data-runs/180511-perslab-immunometab/255L-5000_cells/","/nfsdata/data/sc-10x/data-runs/180511-perslab-immunometab/260L-5000_cells/")' --dir_out /projects/jonatan/tmp-liver/ --flag_datatype sc --flag_organism hsapiens --prefix_data liver  --prefix_run seurat_3_align  --use_filtered_feature_bc_matrix F --nCount_RNA_min 500  --nCount_RNA_max 15000 --nFeature_RNA_min 100 --nFeature_RNA_max 10000 --percent.mito_max Inf --percent.ribo_max Inf --rm_sc_multiplets T --n_comp 75 --align_group_IDs "c('all')" --feats_to_plot "c('percent.mito','percent.ribo','nCount_RNA','nFeature_RNA','MALAT1')"
# time Rscript /projects/jonatan/tools/seurat-src/seurat_pipeline.R --dirs_project_10x 'c("/nfsdata/data/sc-10x/data-runs/180511-perslab-immunometab/")' --dir_out /projects/jonatan/tmp-liver/ --flag_organism hsapiens --prefix_data liver  --prefix_run seurat_3_merge  --use_filtered_feature_bc_matrix F --nCount_RNA_min 500  --nCount_RNA_max 15000 --nFeature_RNA_min 100 --nFeature_RNA_max 10000 --percent.mito_max Inf --percent.ribo_max Inf --rm_sc_multiplets T --n_comp 75 --merge_specify "list(sc=c('245L',247L','249L'),sn=c('255L',259L','260L'))" --feats_to_plot "c('percent.mito','percent.ribo','nCount_RNA','nFeature_RNA','MALAT1')" 

# Sources
## Alignment tutorial: PBMCs - 10x and seqwell: https://satijalab.org/seurat/Seurat_AlignmentTutorial.html
## Alignment tutorial: stimulated vs. control PBMCs: https://satijalab.org/seurat/immune_alignment.html
## Aligning multiple dataset workflow: https://www.dropbox.com/s/aji4ielg8gc70vj/multiple_pancreas_workflow.R?dl=1

# SoupX 
#...

## Left Truncated Mixture Gaussian
## Biorxiv https://www.biorxiv.org/content/early/2018/09/29/430009
## Github https://github.com/zy26/LTMGSCA/blob/master/vignettes/gcr_vignette.pdf

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--dirs_project_10x", type="character", default = NULL,
              help = "Corresponds to cellranger's $DATARUNS_DIR/$PROJECT_ID or e.g. $SAMPLE_NAME subfolders with /outs/filtered_feature_bc_matrix/mm10 subdirs. Takes a quoted vector, [default %default]"),  
  make_option("--paths_data", type="character", default = NULL,
              help = "Path(s) to single counts data matrix stored in standard (gz compressed) delimited text format, Seurat 3.0 or loom. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g. ''c('<path1>','<path2>')'',  [default %default]"),  
  make_option("--dir_out", type="character",
              help = "Project directory to which to write files. Should include subdirectories /tables, /RObjects, /plots, /log"),  
  make_option("--flag_datatype", type="character", default = "sc",
              help = "Accepts arguments sc or bulk, [default %default]"), 
  make_option("--flag_organism", type="character", default = "mmusculus",
              help = "One of mmusculus or hsapiens, [default %default]"), 
  make_option("--prefix_data", type="character", default = "Seurat_out",
              help = "Dataset prefix for output files, [default %default]"), 
  make_option("--prefix_run", type="character", default=substr(gsub("-","",as.character(Sys.Date())),3,1000),
              help = "Run prefix for output files, [default %default]"), 
  make_option("--paths_metadata", type="character", default = NULL,
              help = "Path(s) to metadata file(s) in one of the standard (compressed) character separated formats. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g. ''c('<path1>','<path2>')''. Sample (cell) IDs should be in the first column. If paths_metadata is not given, takes any metadata stored within the path_data object(s). [default %default]"),  
  make_option("--paths_cellCycleGenes", type="character", default = NULL,#'c("s.genes"="/projects/jonatan/genesets/cell_cycle_vignette_files/regev_lab_cell_cycle_genes_s.txt","cc.genes"="/projects/jonatan/genesets/cell_cycle_vignette_files/regev_lab_cell_cycle_genes_g2m.txt")',
              help = "Path(s) to file(s) with cell cycle genes in one of the standard (compressed) character separated formats, for scoring cells on and, if included in vars.to.regress, regressing out, cell cycle related expression variation. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g. ''c('<path1>','<path2>')''. [default %default]"),  
  make_option("--n_cells_loaded", type="character", default = 'c(9000)',
              help = "Approximately how many cells loaded in each sample? Takes a vector, in quotes, with one value per group of samples or a single common value. Used by doubletFinder to estimate doublet rate, [default %default]"),
  make_option("--n_cells_recovered", type="character", default = NULL,
              help = "Approximately how many cells recovered per sample? Takes a vector, in quotes, with one value per group of samples or a single common value. Used to cut off cells in QC. If left as NULL, a sensible value is computed on basis on n_cells_loaded [default %default]"),
  make_option("--use_filtered_feature_bc_matrix", type="logical", default = F,
              help = "Use the filtered_feature_bc_matrix generated by 10x Genomics cell ranger? If F, uses the raw, [default %default]"),
  make_option("--nCount_RNA_min", type="integer", default = 1000L,
              help = "Minimum number of Unique Molecular Identifiers needed to keep a cell [default %default]"),
  make_option("--nCount_RNA_max", type="integer", default = 25000L,
              help = "Max number of Unique Molecular Identifiers tolerated in a cell [default %default]"),
  make_option("--run_SoupX", type="logical", default = TRUE,
              help = "Run the data through SoupX?, [default %default]"),
  make_option("--SoupX_genes", type="character", default = NULL,
              help = "If using SoupX to filter out ambient RNA, optionally provide a custom set of genes for each sample. Takes a quoted list of vectors of gene names, each vector named by sample. Alternatively, a single vector to use for all samples. [default %default]"),
  make_option("--nFeature_RNA_min", type="integer", default = 500L,
              help = "Minimum number of genes needed to keep a cell [default %default]"),
  make_option("--nFeature_RNA_max", type="integer", default = 10000L,
              help = "Max number of genes tolerated in a cell [default %default]"),
  make_option("--percent.mito_max", type="numeric", default = 0.15,
              help = "Maximum proportion of mitochrondrial genes tolerated in a cell, [default %default]"),
  make_option("--percent.ribo_max", type="numeric", default = 1,
              help = "Maximum proportion of ribosomal genes tolerated in a cell, [default %default]"),
  make_option("--rm_sc_multiplets", type="logical", default = T,
              help = "Use DoubletFinder to remove suspected multiplets? [default %default]"),
  make_option("--vars.to.regress", type="character", default='c("nCount_RNA", "percent.mito")',
              help="Provide arguments to Seurat's ScaleData function in the form of a vector in quotes, e.g.''c('nCount_RNA', 'percent.mito', 'percent.ribo')''. For regressing out cell cycle, use ''c('S.Score', 'G2M.Score')'' or ''c('CC.Difference')'', and path_cellcycleGenes must be provided. See https://satijalab.org/seurat/cell_cycle_vignette.html [default %default]"),
  make_option("--merge_group_IDs", type="character", default = NULL,
              help = "Merge samples by some group identifiers? Precedes any sample alignment. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g. ''c('GF','CR')'' to merge all expression data containing the string 'GF' or 'CR', respectively. The group identifier must be part of the names of the sample sub-directories. Use 'c('all')' to merge all samples, [default %default]"), 
  make_option("--merge_specify", type="character", default = NULL,
              help = "Merge samples by individual sample identifiers? Precedes any sample alignment. Similar to merge_group_IDs, but specifying samples in each group by matching strings. Takes a list of vectors of sample identifiers, named by desired group_IDs. The whole list should be in quotes (single if using double inside, and vice versa). E.g. ''list(colon=c('1_GF', '2_CR'), ileum=c('3_GF', '4_CR'))'', [default %default]"), 
  make_option("--align_group_IDs", type="character", default = NULL,
              help = "Align samples by some identifiers? Succeeds any sample merge. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g.''c('tissue.1','tissue.2')''. Use ''c('all')'' to align all samples, [default %default]"), 
  make_option("--align_specify", type="character", default = NULL,
              help = "align samples by individual sample identifiers? Succeeds any sample merge. Similar to align_group_IDs, but specifying samples in each group by matching strings. Takes a list of vectors of sample identifiers, named by desired group identifier. The whole list should be in quotes (single if using double inside, and vice versa). E.g. ''list(colon=c('1_GF', '2_CR'), ileum=c('3_GF', '4_CR'))'', [default %default]"), 
  make_option("--nAnchorFeatures", type="integer", default = 2000,
              help = "How many features (e.g. genes) to use when integrating/aligning datasets [default %default]"), 
  make_option("--n_comp", type="integer", default = 40L,
              help = "How many principal and canonical correlation components?  [default %default]"),
  make_option("--use_jackstraw", type="logical", default = TRUE,
              help = "Use Jackstraw resampling to select subset of PCs with significant gene loadings? May improve results but slow [default %default]"),
  make_option("--res_primary", type="double", default = NULL,
              help = "resolution parameter to pass to Seurat::FindClusters. The usual default is 0.8. Leave NULL to skip clustering [default %default]"),
  make_option("--res_to_calc", type="character", default = NULL,
              help = "multiple resolution parameters to pass to Seurat::FindClusters. Must include res_primary. Format as a quoted character with a vector of values, without whitespace. [default %default]"),
  make_option("--path_transferLabelsRef", type="character", default = NULL,
              help = "If provided, Seurat transfers labels from path_transferLabelsRef, which must be a seurat object with labels[default %default]"),
  make_option("--colLabels", type="character", default = NULL,
              help = "If path_transferLabelsRef is given, give the column name for the labels to transfer [default %default]"),
  # make_option("--minPredictionScore", type="numeric", default = 0.6,
  #             help = "If path_transferLabelsRef is given, minimum prediction score (0-1), [default %default]"),
  make_option("--feats_to_plot", type="character", default = 'c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo","Malat1")',
              help = "Features to plot. Format as a quoted character with a vector of values, without whitespace. [default %default]"),
  make_option("--feats_plot_separate", type="logical", default = T,
              help = "Plot features separately or in one plot? If FALSE, only the first two features will be plotted [default %default]"),
  make_option("--RAM_Gb_max", type="integer", default=200,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]")
  # make_option("--path_runLog", type="character", default=NULL,
  #             help = "Path to file to log the run and the git commit. If left as NULL, write to a file called runLog.text in the dir_log [default %default]")
)

######################################################################
######################### GET CURRENT DIRECTORY ######################
######################################################################

# TODO: replace with here() package

LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  #' @usage returns the current location of the script
  #' @value directory of the script, character
  
  if (interactive()) {
    stop("LocationOfThisScript does not work in interactive sessions")
  }
  
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}

dir_current = LocationOfThisScript()
#dir_current = getwd()

######################################################################
######################### UTILITY FUNCTIONS ##########################
######################################################################

source(file=paste0(dir_current, "/perslab-sc-library/utility_functions.R"))
source(file=paste0(dir_current, "/perslab-sc-library/functions_sc.R"))

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("devtools", "optparse", "Matrix", "Matrix.utils", "Seurat", "ggplot2", "scales", "dplyr", "parallel", "reshape", "reshape2", "cowplot"))#, "pSI", "loomR", "doubletFinder")

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

dirs_project_10x <- opt$dirs_project_10x ; if (!is.null(dirs_project_10x)) dirs_project_10x <- eval(parse(text=dirs_project_10x))
paths_data <- opt$paths_data ; if (!is.null(paths_data)) paths_data <- eval(parse(text=paths_data))
dir_out <- opt$dir_out
flag_datatype <- opt$flag_datatype
flag_organism <- opt$flag_organism
prefix_data <- opt$prefix_data
prefix_run <- opt$prefix_run
paths_metadata <- opt$paths_metadata ; if (!is.null(paths_metadata)) paths_metadata <- eval(parse(text=paths_metadata))
paths_cellCycleGenes <- opt$paths_cellCycleGenes ; if (!is.null(paths_cellCycleGenes)) paths_cellCycleGenes <- eval(parse(text=paths_cellCycleGenes))
n_cells_loaded <- opt$n_cells_loaded ; if (!is.null(n_cells_loaded)) n_cells_loaded <- eval(parse(text=n_cells_loaded))
n_cells_recovered <- opt$n_cells_recovered ; if (!is.null(n_cells_recovered)) n_cells_recovered <- eval(parse(text=n_cells_recovered))
use_filtered_feature_bc_matrix = opt$use_filtered_feature_bc_matrix 
nFeature_RNA_min = opt$nFeature_RNA_min
nFeature_RNA_max = opt$nFeature_RNA_max
nCount_RNA_min = opt$nCount_RNA_min
nCount_RNA_max = opt$nCount_RNA_max
run_SoupX <- opt$run_SoupX
SoupX_genes <- opt$SoupX_genes
if (!is.null(SoupX_genes)) SoupX_genes <- eval(parse(text=SoupX_genes))
percent.mito_max = opt$percent.mito_max
percent.ribo_max = opt$percent.ribo_max
rm_sc_multiplets = opt$rm_sc_multiplets
vars.to.regress <- opt$vars.to.regress ; if (!is.null(vars.to.regress)) vars.to.regress <- eval(parse(text=vars.to.regress))
merge_group_IDs <- opt$merge_group_IDs ; if (!is.null(merge_group_IDs)) merge_group_IDs <- eval(parse(text=merge_group_IDs))
merge_specify <- opt$merge_specify ; if (!is.null(merge_specify)) merge_specify <- eval(parse(text=merge_specify))
align_group_IDs <- opt$align_group_IDs ; if (!is.null(align_group_IDs)) align_group_IDs <- eval(parse(text=align_group_IDs))
align_specify <- opt$align_specify ; if (!is.null(align_specify)) align_specify <- eval(parse(text=align_specify))
nAnchorFeatures <- opt$nAnchorFeatures
n_comp <- opt$n_comp
use_jackstraw <- opt$use_jackstraw
#n_CC <- opt$n_CC
res_primary <- opt$res_primary
res_to_calc <- opt$res_to_calc ; if (!is.null(res_to_calc)) res_to_calc <- eval(parse(text = res_to_calc))
path_transferLabelsRef <- opt$path_transferLabelsRef
colLabels <- opt$colLabels
#minPredictionScore <- opt$minPredictionScore
feats_to_plot <- opt$feats_to_plot ; if (!is.null(feats_to_plot)) feats_to_plot <- eval(parse(text = feats_to_plot))
feats_plot_separate <- opt$feats_plot_separate
#n_cores <- opt$n_cores
RAM_Gb_max <- opt$RAM_Gb_max
#path_runLog <- opt$path_runLog

######################################################################
######################## CONDITIONED PACKAGES  #######################
######################################################################

#if (!is.null(res_primary)) ipak(c("xlsx"))
if (!is.null(res_primary)) ipak(c("openxlsx"))

if (run_SoupX) {
  devtools::install_github(repo="constantAmateur/SoupX", dependencies=NA, upgrade = "never")
  library(SoupX)
  if (is.null(SoupX_genes)){
    devtools::install_github("zy26/LTMGSCA", dependencies=NA, upgrade = "never")
    library(LTMGSCA)
  }
}

######################################################################
############################ OPTIONS #################################
######################################################################

options(stringsAsFactors = F, use = "pairwise.complete.obs", warn=1)

######################################################################
############################ CONSTANTS ###############################
######################################################################

as.character(Sys.time()) %>% gsub("\\ ", "_",.) %>% gsub("\\:", ".", .) ->tStart

# if specified output directory doesn't exist, create it 
if (!file.exists(dir_out)) {
  dir.create(dir_out) 
  message("Project directory not found, new one created")
}

dir_plots = paste0(dir_out,"plots/")
if (!file.exists(dir_plots)) dir.create(dir_plots) 

dir_tables = paste0(dir_out,"tables/")
if (!file.exists(dir_tables)) dir.create(dir_tables)

dir_RObjects = paste0(dir_out,"RObjects/")
if (!file.exists(dir_RObjects)) dir.create(dir_RObjects)

dir_log = paste0(dir_out,"log/")
if (!file.exists(dir_log)) dir.create(dir_log)

dir_scratch = "/scratch/tmp-seurat/"

#dir_current <- "/projects/jonatan/tools/seurat-src/"

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

pvalThreshold=0.05

randomSeed <- 12345
set.seed(randomSeed)

# DoubletFinder
proportion.artificial = 0.25
proportion.NN = 0.02

# Seurat FindIntegrationAnchors and FindTransferAnchors 
# https://satijalab.org/seurat/pancreas_integration_label_transfer.html
k.anchor=5
k.filter=200
k.score=30
max.features = 200
k.weight = 50
sd.weight = 1

######################################################################
########################### VERIFY INPUT #############################
######################################################################

if (!xor(is.null(dirs_project_10x), is.null(paths_data))) stop("Provide dirs_project_10x or path_data arg")
if (!is.null(merge_group_IDs) & !is.null(merge_specify)) stop("To merge datasets, pass merge_group_IDs or merge_specify but not both")
if (!is.null(align_group_IDs) & !is.null(align_specify)) stop("To align datasets, pass merge_group_IDs or merge_specify but not both")
if (!flag_organism %in% c("mmusculus", "hsapiens")) stop("flag_organism must be set to mmusculus or hsapiens")

######################################################################
############# LOAD 10x DATA AND CREATE SEURAT OBJECT #################
######################################################################

message("Loading raw data")
# Get dirs_sample and sample_IDs by searching dirs_project_10x
if (!is.null(dirs_project_10x)) {
  
  ref_transcript_cell <- if (flag_organism == "mmusculus") "mm10" else if (flag_organism=="hsapiens") "hg19"
  ref_transcript_nuclei <- if (flag_organism == "mmusculus") "mm10-1\\.2\\.0_premrna" else if (flag_organism=="hsapiens") "GRCh38"

  # Do initial rough search
  fun <- function(dir_project_10x) {
    dir(path = dir_project_10x, full.names = T, recursive = T, include.dirs = T)
  }
  list_iterable <- list("X" = dirs_project_10x)
  dirs_all <- safeParallel(fun=fun, list_iterable=list_iterable)
  dirs_all <- unlist(dirs_all, use.names=F)
  
  # For some reason, giving the pattern argument to dir() returns nothing
  pattern0 = if (!use_filtered_feature_bc_matrix | run_SoupX) { 
    paste0(".*outs/raw_feature_bc_matrix$|.*/raw_gene_bc_matrices/", ref_transcript_cell, "$|.*outs/raw_genes_bc_matrices/", ref_transcript_nuclei, "$")
    #paste0(".*outs/filtered_feature_bc_matrix/", ref_transcript_cell, "$|.*outs/filtered_feature_bc_matrix, "/", ref_transcript_nuclei, "$")
  } else { 
    paste0(".*outs/filtered_feature_bc_matrix$|.*/filtered_gene_bc_matrices/", ref_transcript_cell, "$|.*outs/filtered_genes_bc_matrices/", ref_transcript_nuclei, "$")
    #paste0(".*outs/raw_feature_bc_matrix/", ref_transcript_cell, "$|.*outs/raw_feature_bc_matrix", "/", ref_transcript_nuclei, "$")
  }
  
  # successively filter down vector of directories
  dirs_sample <- grep(pattern = pattern0, x = dirs_all, value=T)
  
  pattern_merge_group_IDs <- pattern_align_group_IDs <- pattern_align_specify <- pattern_merge_specify <- NULL
  ## Added : take into account merge_group_IDs, merge_specify, align_group_IDs, align_specify
  if (!is.null(merge_group_IDs)) {
    if (merge_group_IDs[1]!="all") merge_group_IDs[!sapply(merge_group_IDs, is.null)] %>% paste0(., collapse = "|") -> pattern_merge_group_IDs 
  }
  if (!is.null(align_group_IDs)) {
    if (align_group_IDs[1]!="all") align_group_IDs[!sapply(align_group_IDs, is.null)] %>% paste0(., collapse = "|") -> pattern_align_group_IDs 
  }
  
  if (!is.null(merge_specify)) {
  merge_specify[!sapply(merge_specify, is.null)] %>% 
    unlist(., use.names = F) %>% paste0(., collapse = "|") -> pattern_merge_specify 
  }
  
  if (!is.null(align_specify)) {
    align_specify[!sapply(align_specify, is.null)] %>% 
      unlist(., use.names = F) %>% paste0(., collapse = "|") -> pattern_align_specify 
  } 
  
  # continue filtering vector of directories
  
  if (!is.null(pattern_merge_group_IDs)) {
    dirs_sample <- grep(pattern = pattern_merge_group_IDs, x = dirs_sample, value=T)
  }
  if (!is.null(pattern_align_group_IDs)) {
    dirs_sample <- grep(pattern = pattern_align_group_IDs, x = dirs_sample, value=T)
  }
  if(!is.null(pattern_merge_specify)) {
    dirs_sample <- grep(pattern = pattern_merge_specify, x = dirs_sample, value=T)
  }
  if(!is.null(pattern_align_specify)) {
    dirs_sample <- grep(pattern = pattern_align_specify, x = dirs_sample, value=T)
  }
  
  dirs_sample <- paste0(dirs_sample, "/")
  
  # strip down the directory strings to get sample_IDs
  #sample_IDs <- gsub(paste0("/scratch|/nfsdata|/data|/sc-10x|/data-runs|-\\d{4}_cells/|outs|/raw_feature_bc_matrix|/filtered_feature_bc_matrix|/", ref_transcript_cell ,"/|", ref_transcript_nuclei ,"/"), "", dirs_sample)
  sample_IDs <- gsub(paste0("/scratch|/nfsdata|/data|/sc-10x|/data-runs|-\\d{4}_cells/|outs|/raw_feature_bc_matrix/|/raw_gene_bc_matrices|/filtered_feature_bc_matrix|/filtered_gene_bc_matrices", "|",ref_transcript_cell , "|", ref_transcript_nuclei), "", dirs_sample)
  sample_IDs <- gsub("^/[^/]*", "", sample_IDs)
  sample_IDs <- gsub("/", "", sample_IDs)
  
  rm(dirs_all)
  
} else {
  dirs_sample <- paths_data 
  sample_IDs <- gsub(".*/|\\.RData|\\.rda|\\.RDS|\\.loom|\\.csv|\\.tab|\\.txt", "", dirs_sample, ignore.case = T)
}

names(dirs_sample) <- sample_IDs

if (is.null(n_cells_recovered)) {
  train_cells_loaded <- c(870, 1700, 3500, 5300, 7000, 8700, 10500, 12200, 14000, 15700, 17400)
  train_cells_recovered <- c(500,seq(from=1000,to=10000, by=1000))
  n_cells_recovered_lm <- lm(train_cells_recovered~train_cells_loaded)
  n_cells_recovered <- round(n_cells_recovered_lm$coefficients[1]+n_cells_recovered_lm$coefficients[2]*n_cells_loaded, 0) %>% na.omit %>% as.integer
} 

if (length(n_cells_loaded)==1) n_cells_loaded <- rep(x = n_cells_loaded, times=length(sample_IDs))
if (length(n_cells_recovered)==1) n_cells_recovered <- rep(x = n_cells_recovered, times=length(sample_IDs))

# Check input param n_cells_loaded
if (length(n_cells_loaded) != length(dirs_sample)) stop("n_cells_loaded must be either a single integer or an integer per sample")
if (length(n_cells_recovered)!=length(dirs_sample)) stop("n_cells_recovered must be either a single integer or an integer per sample")
if (any(n_cells_loaded<n_cells_recovered)) stop("n_cells_recovered cannot exceed n_cells_loaded")

if (!is.null(dirs_project_10x)) {
  # call Read10X
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_Read10X.txt")
  list_data_tmp <- safeParallel(list_iterable=list("dir_sample"=dirs_sample), fun=Read10X, outfile=outfile)

} else {
  # if expression data files are in a different format than that output by 10x Genomics cellranger, e.g. as matrices 
  fun = function(dir_sample, sample_ID) {
    data_tmp <- load_obj(f=dir_sample)
    if (all(c(1,2,3) %in% rownames(data_tmp))) {
      if (!is.null(data_tmp[["X"]]))  {
        rownames(data_tmp) <- data_tmp[["X"]] 
        data_tmp[["X"]] <- NULL
      } else {
        stop("cannot identify expression data gene names. Try to add as a column 'X'")
      }
    }
    return(data_tmp)
  }

  list_iterable = list("dir_sample"=dirs_sample, "sample_ID"= sample_IDs)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_loadData.txt")
  list_data_tmp <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
}

# Get filtered cell idx - unless we loaded cellranger's filtered matrix and are not using SoupX

if (run_SoupX & use_filtered_feature_bc_matrix)  {
  fun = function(dir_sample, data_tmp) {
    #Get the barcodes that cell ranger considers to contain cells
    dir_up <- gsub(paste0("raw_feature_bc_matrix/$|raw_gene_bc_matrices/", ref_transcript_cell, "/$|raw_gene_bc_matrices/", ref_transcript_nuclei, "/$"), "", dir_sample)
    path_barcodes <- dir(path=dir_up, pattern ="barcodes.tsv", full.names = T, recursive = T)
    path_barcodes <- grep(pattern="filtered", x= path_barcodes,  value=T)
    #path_barcodes <- gsub(paste0("raw_feature_bc_matrix/$|raw_gene_bc_matrices/", ref_transcript_cell, "/$|raw_gene_bc_matrices/", ref_transcript_nuclei, "/$"), "filtered_feature_bc_matrix/barcodes.tsv.gz", dir_sample)
    cells = read.delim(path_barcodes,sep='\t',header=FALSE)
    cells = gsub('-1','',cells[,1])
    
    #Get the index in the big table
    cellIdx <- match(cells,colnames(data_tmp))
  }
  
  list_iterable = list("dir_sample"=dirs_sample, "data_tmp"=list_data_tmp)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_get_filtered_cell_idx.txt")
  
  list_cellIdx <- safeParallel(fun=fun, 
                               list_iterable=list_iterable, 
                               outfile=outfile)
  
} else if (!use_filtered_feature_bc_matrix) {
  # whether or not we use SoupX
  fun = function(n_cells, data_tmp, sample_ID) {
    
    data_tmp %>% colSums -> nCount_RNA_sums
    nCount_RNA_sums %>% rank -> rank_tmp
    which(rank_tmp >= (length(rank_tmp)-n_cells)) -> cellIdx
    #which(rank_tmp <= (length(rank_tmp)-5000) & rank_tmp > (length(rank_tmp)-10000)) -> cellIdx
    nCount_RNA_sums[cellIdx] %>% sort(decreasing=T) -> nCount_RNA_sums_sort
    # Plot nCount_RNA histogram 
    ggplot(data.frame(nCount_RNA=nCount_RNA_sums_sort, barcode=1:length(nCount_RNA_sums_sort)), aes(barcode,nCount_RNA)) + geom_line() + 
      scale_x_continuous(trans='log10', limits=c(1,n_cells+n_cells%/%5), breaks = c(sapply(1:4, function(x) 10^x))) + 
      scale_y_continuous(trans='log10', limits=c(1,nCount_RNA_sums_sort[1]), breaks = c(sapply(1:4, function(x) 10^x))) + 
      geom_hline(yintercept=nCount_RNA_min) + geom_hline(yintercept=nCount_RNA_max) + 
      ggtitle(paste0(sample_ID, ": UMI counts vs. barcodes, top ", as.character(n_cells), " barcodes"))
    saveMeta(savefnc=ggsave, filename=paste0(dir_plots, prefix_data,"_", prefix_run,"_", sample_ID,"_nCount_RNA_vs_barcode.pdf"), w=12, h=8)
    return(cellIdx)
  } 
  
  list_iterable = list("n_cells"=n_cells_recovered, "data_tmp"=list_data_tmp, "sample_ID"=sample_IDs)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_get_top_cell_idx.txt")
  
  list_cellIdx <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

} else {
  list_cellIdx <- NULL
}

# SoupX ambient RNA filtering
if (run_SoupX) {
  
  fun= function(data_tmp, cellIdx, channelName) {
    SoupChannel(tod=data_tmp,
                toc=data_tmp[,cellIdx,drop=FALSE],
                channelName=channelName,
                #ref=ref,
                #path=dataDir,
                dataType='10X')
    }
  list_iterable = list("data_tmp" = list_data_tmp, "cellIdx"=list_cellIdx, "channelName" = sample_IDs)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_SoupChannel.txt")
  list_channels <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
  
  SoupChannelList <- function (channels, ...)
  {
    if (any(duplicated(sapply(channels, function(e) e$channelName))))
      stop("Duplicate channel names found.  Please give each channel a unique name before continuing.")
    names(channels) = sapply(channels, function(channel) channel$channelName)
    scl = list(channels = channels)
    fun1= function(mat1,mat2) {
      Matrix.utils::merge.Matrix(x = mat1,
                                 y=mat2,
                                 by.x =rownames(mat1),
                                 by.y=rownames(mat2),
                                 all.x=T,
                                 all.y=T)
    }
    scl$toc = Reduce(f = fun1, x=lapply(channels, function(e) e$toc))
    scl$nUMIs = do.call(c, lapply(channels, function(e) e$nUMIs))
    scl = c(scl, list(...))
    class(scl) = c("list", "SoupChannelList")
    
    list_df <- mapply(function(channel, colname) {
      df <- data.frame("est"=channel$soupProfile[,"est"],
                       "gene"=rownames(channel$soupProfile))
      colnames(df) <- c(colname, "gene")
      return(df)}, channel=scl$channels, colname=names(list_channels), SIMPLIFY=F)
    
    fun2 <- function(df1, df2) {
      full_join(df1,df2, by="gene")
    }
    df_join <- Reduce(f=fun2, 
                      x=list_df)
    rownames(df_join) <- df_join[["gene"]]; df_join[["gene"]] <- NULL
    scl$soupMatrix <- df_join
    #rownames(scl$soupMatrix) = rownames(scl$toc)
    
    scl
  }
  scl = SoupChannelList(list_channels)

  if (!is.null(SoupX_genes)) {
    list_topgenes <- lapply(1:length(sample_IDs), function(i) SoupX_genes)
  }  else {
    # since no suitably specifically expressed genes provided by user
    # find bimodally distributed genes 
    scl <- inferNonExpressedGenes(scl)
    
    # Find genes with a bimodal distribution
    fun <- function(sample_ID) {
      if (any(scl$channels[[sample_ID]]$nonExpressedGenes$isUseful)) {
        nonExpressedGenes_useful <- scl$channels[[sample_ID]]$nonExpressedGenes[scl$channels[[sample_ID]]$nonExpressedGenes$isUseful,]
        nonExpressedGenes_useful <- nonExpressedGenes_useful[order(nonExpressedGenes_useful$extremity, decreasing = T),]
        
        if (sum(nonExpressedGenes_useful$extremity>0.7)<=1) return(NA) 
        
        topgenes <- rownames(nonExpressedGenes_useful)[nonExpressedGenes_useful$extremity>0.7]
        topgenes <- topgenes[1:min(10, length(topgenes))]
        
        # tag log of sparse matrix, avoiding error due to large size
        fun <- function(bigSpMat){
          submat1 <- bigSpMat[,1:(ncol(bigSpMat)%/%2)]
          submat1 <- submat1 %>% as.matrix %>% log 
          submat2 <- bigSpMat[,(ncol(bigSpMat)%/%2+1):ncol(bigSpMat)]
          submat2 <- submat2 %>% as.matrix %>% log 
          mat_out <- cbind(submat1, submat2)
          mat_out[is.infinite(mat_out)] <- 0
          return(mat_out)
          }
        toc_log <- fun(scl$toc)#log(as.matrix(scl$toc))
        
        list_fits <- lapply(topgenes, function(genes) LTMGSCA::SeparateKRpkmNew(x=toc_log[genes,, drop=F], n=100,q=0,k=2,err = 1e-10))
        idx_topgenes_ok <- sapply(list_fits, function(fit) fit[2,1]>0.4)
        topgenes <- topgenes[idx_topgenes_ok]
        
        if (length(topgenes)<=1) return(NA)
        
        return(topgenes) 
  
      } else {
        return(NA)
      }
    }
    list_iterable= list("sample_ID"=sample_IDs)
    outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_SoupX_findTopGenes.txt")
    list_topgenes <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
  }
  
  # If we have some suitable genes to estimate ambient RNA contamination, filter the expression matrix
  if (!any(sapply(list_topgenes, function(vec_topgenes) all(is.na(vec_topgenes)), simplify = T))) {
    fun = function(data_tmp, sample_ID, topgenes) {
      if (all(is.na(topgenes))) return(data_tmp)
      tryCatch({
        #         # This plot shows:
        #         # 1. for each selected gene sample 100 cells
        #         # 2. for each cell, compute the expected normalized expression based on soup
        #         # 3. for each cell, compute the observed normalized expression
        #         # 4. for each cell, compute the log ratio of observed to expected
        #         # 5. where the ratio is high, call the gene as actually expressed (not just by soup)
        #         #       Ideally the 'actually expressed' should look like a mode
        plotMarkerDistribution(scl, sample_ID, topgenes, maxCells=250)
        saveMeta(savefnc=ggsave, filename = paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_SoupX_plotMarkerDistribution.pdf"), w=15, h=8)
        scl = calculateContaminationFraction(scl, sample_ID, list(genes = topgenes), excludeMethod = 'pCut', exCut=0.5)
        # Cell level contamination fraction
        scl = interpolateCellContamination(scl, sample_ID, interpolationMethod = c("lowess"))
        scl = adjustCounts(scl) # adds adjusted (raw) counts matrix in scl$atoc. 
        # Alternative is strainCells, which modifies expression fraction. Columns (cells all sum to 1) 
        # What genes were set to zero min most cells? (output to log file)
        cntSoggy = rowSums(scl$toc > 0)
        cntAdjusted = rowSums(scl$atoc > 0)
        mostZeroed <- head(sort(x=(cntSoggy - cntAdjusted)/cntSoggy, decreasing=T), n= 20)
        # Continue with adjusted counts matrix
        data_tmp <- scl$atoc
        # Save topgenes temporarily for plotting later
        saveMeta(savefnc=saveRDS, object = topgenes, 
                file = paste0(dir_RObjects, prefix_data, "_", prefix_run, "_", sample_ID, "_SoupX_topgenes.RDS"))
        return(data_tmp)
      }, error= function(err) {
        data_tmp
      })
    }
      
    outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_SoupXAdjustcounts.txt")
    list_iterable=list("data_tmp"=list_data_tmp, "sample_ID" = sample_IDs, topgenes = list_topgenes)
    list_data_tmp <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile, scl=scl) 
  }
} else if (!use_filtered_feature_bc_matrix){
  fun = function(data_tmp, cellIdx) {
    data_tmp[,cellIdx, drop=F]
  }
  list_iterable=list(data_tmp=list_data_tmp, cellIdx = list_cellIdx)
  outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_filterDataTmp.txt")
  list_data_tmp <- safeParallel(fun=fun,list_iterable=list_iterable, outfile = outfile)
}

# Now we have filtered on expected number of cells and number of UMI or genes, we compute percent.mito, percent ribo
# and filter on nCount_RNA and nFeature_RNA min and max
fun = function(data_tmp, sample_ID, n_cells) {
  # Compute and add percent mito and percent ribo as metadata
  mito.genes <- grepl(pattern = "^mt-", x = rownames(data_tmp), ignore.case=T)
  ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(data_tmp), ignore.case=T)
  colSums_tmp <- colSums(x = data_tmp)
  
  metadata <- data.frame(percent.mito=colSums(x = data_tmp[mito.genes,])/colSums_tmp, 
                         percent.ribo = colSums(x = data_tmp[ribo.genes,])/colSums_tmp, 
                         row.names=colnames(data_tmp))
  # Add nCount_RNA
  metadata[["nCount_RNA"]] <- colSums(data_tmp)
  
  # add idents
  #if (!is.null(idents)) metadata[["idents"]] <- idents
  
  # Add nFeature_RNA
  metadata[["nFeature_RNA"]] <-  colSums(data_tmp>0)
  metadata[["sample_ID"]] <-  rep(sample_ID, nrow(metadata))
  
  # Add sample IDs
  rownames(metadata) <- colnames(data_tmp)
  # Plot QC metrics
  ## nCount_RNA_nFeature_RNA
  ggplot(metadata, aes(nCount_RNA, nFeature_RNA)) + 
    geom_point(shape=1) + 
    scale_x_continuous(trans='log10', limits=c(1,max(metadata[["nCount_RNA"]])), breaks = c(sapply(1:4, function(x) 10^x))) +
    scale_y_continuous(trans='log2', limits=c(1, max(metadata[["nFeature_RNA"]])), breaks = c(sapply(1:4, function(x) 10^x))) + 
    geom_vline(xintercept=nCount_RNA_min) + 
    geom_vline(xintercept=nCount_RNA_max) +
    geom_hline(yintercept=nFeature_RNA_min) + 
    geom_hline(yintercept=nFeature_RNA_max) +
    ggtitle(paste0(sample_ID, ": nFeature_RNA vs nCount_RNA counts vs. barcodes, top ", as.character(n_cells), " barcodes"))
  
  saveMeta(savefnc=ggsave, filename =  paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_nCount_RNA_nFeature_RNA.pdf"), w=12, h=8)
  ## percent.mito, percent.ribo
  ggplot(metadata, aes(percent.mito, percent.ribo)) + 
    geom_point(shape=1) + 
    geom_vline(xintercept=percent.mito_max) +
    geom_hline(yintercept = percent.ribo_max) + 
    #scale_y_continuous(breaks = seq(from=0.1, to=0.9, by=0.1)) +
    #scale_x_continuous(breaks = seq(from=0.1, to=0.9, by=0.1)) + 
    ggtitle(paste0(sample_ID, ": prop. ribo vs. prop. mito, top ", as.character(n_cells), " barcodes"))
  saveMeta(savefnc=ggsave, filename =  paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_percent.mito_percent.ribo.pdf"), w=12, h=8)
  
  # Filter data matrix on QC metrics
  idx_QC_ok <- metadata[["nCount_RNA"]] >= nCount_RNA_min & 
    metadata[["nCount_RNA"]] <= nCount_RNA_max & 
    metadata[["nFeature_RNA"]] >= nFeature_RNA_min & 
    metadata[["nFeature_RNA"]] <= nFeature_RNA_max & 
    metadata[["percent.mito"]] <= percent.mito_max & 
    metadata[["percent.ribo"]] <= percent.ribo_max
  
  if (sum(idx_QC_ok)>=50) {
    data_tmp <- data_tmp[,idx_QC_ok]
    metadata <- metadata[idx_QC_ok,]
    rownames(metadata) <- colnames(data_tmp)
    
    seurat_obj <- CreateSeuratObject(counts = data_tmp, 
                                     project = paste0(prefix_data, "_", prefix_run), 
                                     #min.cells = -Inf, 
                                     min.features = -Inf, 
                                     #normalization.method = "LogNormalize", 
                                     #scale.factor = 1e4, 
                                     #do.scale=F, 
                                     #do.center=F, 
                                     meta.data = metadata)
    
    # We scale and center here as we do not wish to regress out technical confounders before doubletFinder - see preprint p. 9: https://www.biorxiv.org/content/biorxiv/early/2018/06/20/352484.full.pdf 
    return(seurat_obj)
  } else {
    warning (paste0(sample_ID, " had fewer than 50 cells after QC and was discarded"))
    return(NULL)
  }
}

list_iterable=list("sample_ID"=sample_IDs, "data_tmp" = list_data_tmp, "n_cells"=n_cells_recovered)
outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_QCandMakeSeuratObj.txt")

list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

names(list_seurat_obj) <- sample_IDs

idx_sample_ok <- !sapply(list_seurat_obj, is.null)

samples_dropped <- paste0(sample_IDs[!idx_sample_ok], collapse=" ")

if (nchar(samples_dropped)>0) {
  samples_dropped_df <- data.frame(samples = samples_dropped, data = prefix_data, run = prefix_run, percent.mito_max=percent.mito_max, percent.ribo_max=percent.ribo_max)
  saveMeta(savefnc=write.csv, x=samples_dropped_df, file=paste0(dir_log, prefix_data, "_", prefix_run, "_samples_dropped.csv"), quote = F, row.names = F)
}

# Filter files 
list_seurat_obj <- list_seurat_obj[idx_sample_ok]
sample_IDs <- sample_IDs[idx_sample_ok]
n_cells_loaded <- n_cells_loaded[idx_sample_ok]
n_cells_recovered <- n_cells_recovered[idx_sample_ok]

rm(list_data_tmp)

######################################################################
####################### ADD META DATA  ###############################
######################################################################

#Add any additional metadata specified by the user
if (!is.null(paths_metadata)) {
  outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_load_metadata.txt")
  fun = function(path_metadata) {
    if (!is.na(path_metadata)) {
      df_metadata <- load_obj(path_metadata) 
      # Set the first column as rownames
      rownames(df_metadata) <- df_metadata[,1] 
      df_metadata[[1]] <- NULL
      return(df_metadata)
      } else { return(NA)}
  }
  list_iterable = list("X"=paths_metadata)
  list_metadata <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

  # add metadata to seurat objects
  if (length(list_metadata)==length(list_seurat_obj)) {
    list_seurat_obj <- mapply(function(metadata, seurat_obj) {
      if (!is.na(df_metadata)) AddMetaData(object=seurat_obj, metadata=df_metadata) else seurat_obj
    }, seurat_obj =  list_seurat_obj, df_metadata=list_metadata)
  } else {
    warning("list_metadata length did not match number of datasets, could not add metadata")
  }
}

######################################################################
###################### LOAD CELL CYCLE DATA  ###############################
######################################################################

#Add any additional metadata specified by the user
if (!is.null(paths_cellCycleGenes)) {
  outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_load_cellCycleGenes.txt")
  fun = load_obj
  list_iterable = list("X"=paths_cellCycleGenes)
  list_cellCycleGenes <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
  # convert to list
  list_cellCycleGenes <- lapply(list_cellCycleGenes, FUN = '[[', (1))
  names(list_cellCycleGenes) <- names(paths_cellCycleGenes)
  
}

if (rm_sc_multiplets) {

  ######################################################################
  ######################### FIND DOUBLETS ##############################
  ######################################################################
  
  message("Detecting single cell multiplets")
  
  ######################################################################
  ############################### DOUBLETFINDER ########################
  ######################################################################
  
  # changes to original:
  # * instead of using expected.doublets as cut-off, uses outlier criterion pANN > Q3+1.5*(Q3-Q1)
  # * rename artificial doublets "thisisfake" rather than "X" to avoid issues e.g. with 10X cells!
  # * fixed a bug with PCA
  
  list_real.cells <- lapply(list_seurat_obj, colnames)
  
  fun = function(seurat_obj, real.cells) {
    data <- GetAssayData(seurat_obj, slot="counts")[, colnames(seurat_obj)]
    #real.cells <- colnames(seurat_obj)
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - proportion.artificial) -
                          n_real.cells) # these will be artificial doublets
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE) # draw samples of real cells of size of n artificial doublets to make
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2 # make artificial doublets
    colnames(doublets) <- paste("thisisfake", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    seu_wdoublets <- Seurat::CreateSeuratObject(counts = data_wdoublets)
  }
  list_iterable = list("seurat_obj"=list_seurat_obj, "real.cells"=list_real.cells)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_addFakeDoublets.txt")
  list_seurat_obj_wfakeDoub <- safeParallel(fun=fun, 
                                            list_iterable=list_iterable, 
                                            outfile=outfile)
  
  # NormalizeData
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_NormalizeData_wFakeDoub.txt")    
  list_seurat_obj_wfakeDoub <- safeParallel(fun=Seurat::NormalizeData, list_iterable=list("X"=list_seurat_obj_wfakeDoub), outfile=outfile)
  
  # FindVariableFeatures
  fun = function(seurat_obj) tryCatch({
    FindVariableFeatures(object = seurat_obj,
                        nfeatures = nAnchorFeatures,
                        selection.method="vst", verbose=T)}, 
  error= function(err) {
    message("FindVariableFeatures failed with selection.method='vst', trying with 'mean.var.plot'")
    FindVariableFeatures(object = seurat_obj,
                         nfeatures = nAnchorFeatures,
                         selection.method="mean.var.plot",verbose=T)
  })
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_FindVariableFeatures_wFakeDoub.txt")    
  list_seurat_obj_wfakeDoub <- safeParallel(fun=fun, list_iterable=list("seurat_obj"=list_seurat_obj_wfakeDoub), outfile=outfile)#,

  # ScaleData
  list_seurat_obj_wfakeDoub <- lapply(X=list_seurat_obj_wfakeDoub, 
                            FUN = function(seurat_obj) {
                              ScaleData(object = seurat_obj,
                                        block.size=15000, 
                                        # Don't regress out, this is fake data, not to be used later
                                        min.cells.to.block=5000,
                                        verbose=T) 

                            })
  
  # PCA 
  fun <- function(seu_wdoublets) {tryCatch({
      pcs.compute1 <- min(40, min(ncol(GetAssayData(seu_wdoublets, slot="scale.data")), length(VariableFeatures(seu_wdoublets)))%/%2)
      Seurat::RunPCA(seu_wdoublets,
                     #pc.genes = seu_wdoublets@var.genes,
                     npcs = pcs.compute1,
                     verbose=F)
    }, error= function(err) {
      pcs.compute2 <- min(20, min(ncol(GetAssayData(seu_wdoublets, slot="scale.data")), length(VariableFeatures(seu_wdoublets))%/%3))
      warning(paste0("RunPCA failed with ", pcs.compute1," components, trying again with ", pcs.compute2, " components"))
      Seurat::RunPCA(seu_wdoublets,
                     #pc.genes = seu_wdoublets@var.genes,
                     do.print = F,
                     npcs = pcs.compute2,
                     verbose=F)})}
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_PCA_wFakeDoub.txt")     
  list_seurat_obj_wfakeDoub <- safeParallel(fun=fun, list_iterable=list("X"=list_seurat_obj_wfakeDoub), outfile=outfile)
  
  # Find make cell-cell distance matrix, find nearest neighbours, predict whether cells are doublets or singlets
  fun <- function(seu, seu_wdoublets, real.cells) {
    n_real.cells = length(real.cells)
    cell.names <- colnames(seu_wdoublets)
    nCells <- length(cell.names)
    PCs <- 1:ncol(Embeddings(object=seu_wdoublets, reduction="pca"))
    pca.coord <- Embeddings(object=seu_wdoublets, reduction="pca")#seu_wdoublets@dr$pca@cell.embeddings[, PCs]
    rm(seu_wdoublets)
    gc()
    dist.mat <- as.matrix(dist(pca.coord)) # distances in PCA space
    dist.mat <- dist.mat[, -grep("thisisfake", colnames(dist.mat))] # so keep fake cell rows but not columns
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * proportion.NN) # k is how many neighbours to check
    for (i in 1:n_real.cells) { # i goes up to the last real cell in the distance matrix
      neighbors <- order(dist.mat[, i]) # all cell numbers ranked by distance (i.e. closest first)
      neighbors <- neighbors[2:(k + 1)] # get k nearest neighbours
      neighbor.names <- rownames(dist.mat)[neighbors] # these will be cell names. Artificial cells have "X" in them
      pANN[i, 1] <- length(grep("thisisfake", neighbor.names))/k # for real cell i, how big a proportion of k nearest neighbours are artificial?
    }
    seu[["pANN"]] <- pANN # proportion of fake neighbours to original (real) seurat obj
    predictions <- as.data.frame(rep("Singlet", n_real.cells),
                                 ncol = 1, row.names = real.cells , stringsAsFactors = FALSE) # initialise prediction metadata column
    ### modified @author Jonatan Thompson, jjt3f2188@gmail.com  @date 181108 ###
    Q1 <- quantile(x = pANN[,1], probs=0.25)
    Q3 <- quantile(x = pANN[,1], probs=0.75)
    pANN_inter_Q_range <- as.numeric(Q3-Q1)
    doublet.predictions <- colnames(seu)[pANN>Q3+1.5*(Q3-Q1)]
    # doublet.predictions <- rownames(seu@meta.data)[order(seu@meta.data$pANN,
    #                                                      decreasing = TRUE)] # order real cells by prop fake neighbours
    # doublet.predictions <- doublet.predictions[1:expected.doublets] # take the top of the list, based on absolute number
    #
    predictions[doublet.predictions, ] <- "Doublet"
    colnames(predictions) <- "pANNPredictions"
    seu[["pANNPredictions"]] <-  predictions
    return(seu)
  }
  list_iterable=list("seu" = list_seurat_obj, "seu_wdoublets"=list_seurat_obj_wfakeDoub, "real.cells"=list_real.cells)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_doubletNNandPredict.txt")   
  list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable, timeout = 12000, outfile=outfile)
  
  # doublet rate data source: 10x https://pdfs.semanticscholar.org/presentation/0c5c/ad2edbb8f0b78cdba1c78b4324213ac20ab3.pdf
  # Save a log 
  log_prop_singlets <- sapply(X= list_seurat_obj, FUN = function(seurat_obj){
    tryCatch({
      round(sum(seurat_obj$pANNPredictions=="Singlet")/ncol(seurat_obj),2)}, error = function(err) NA)
  }, simplify = T)
  
  log_prop_singlets_df <- matrix(log_prop_singlets, nrow=1) %>% data.frame(.,stringsAsFactors = F, row.names = NULL)
  colnames(log_prop_singlets_df) <- names(list_seurat_obj)
  saveMeta(savefnc=write.table, x=log_prop_singlets, file=paste0(dir_log, prefix_data, "_", prefix_run, "_prop_singlets.tab"))
  
  message("Removing suspected doublets")
  # Filter out doublets
  outfile = paste0(dir_log, prefix_data, "_", prefix_run, "_doubletFinder_subset.txt")
  fun = function(seurat_obj, obj_name) {
    tryCatch({
      subset(x = seurat_obj, subset= pANNPredictions=="Singlet")
      #SubsetData(object=seurat_obj, subset.name = "pANNPredictions", accept.value="Singlet", subset.raw=T)
    }, error = function(err) {
      warning(paste0(obj_name, ": filtering out doublets failed with error: ", err, ". Maybe doubletFinder failed, possibly due to an upstream PCA or t-SNE failure."))
      seurat_obj
    })
  }
  list_iterable=list(seurat_obj=list_seurat_obj, obj_name = names(list_seurat_obj))  
  list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
}

######################################################################
############################# MERGE ##################################
######################################################################

if ((!is.null(merge_group_IDs) | !is.null(merge_specify)) & length(list_seurat_obj)>1) {
  
  list_seurat_merged = NULL 
  
  message("Merging samples")
  
  if (!is.null(merge_group_IDs)) {
    if (merge_group_IDs[1] == 'all') {
      list_vec_merge_idx <- list('all'=rep(TRUE, length=length(list_seurat_obj)))
      names(list_vec_merge_idx[[1]])=names(list_seurat_obj)
    } else {
      list_vec_merge_idx <- lapply(merge_group_IDs, function(id) {
        vec_idx <- grep(pattern=paste0("^",id), 
                              x = names(list_seurat_obj),
                              ignore.case=T)
        if (length(vec_idx)>0) names(vec_idx) <- names(list_seurat_obj)[vec_idx]
        vec_idx
        })
      names(list_vec_merge_idx) = merge_group_IDs
      list_vec_merge_idx <- Filter(f = length, x=list_vec_merge_idx)
      if (length(list_vec_merge_idx)==0) stop("merge_group_IDs: no matches")
    }
  } else if (!is.null(merge_specify)) {
    list_vec_merge_idx <- lapply(merge_specify, function(vec_sample_ids) {
      sapply(vec_sample_ids, function(id) {
        vec_idx <- grep(pattern=paste0("^",id), x = names(list_seurat_obj), ignore.case = T)
        if (length(vec_idx)>0) names(vec_idx) <- names(list_seurat_obj)[vec_idx]
        vec_idx
      }, simplify=T)
    }) 
    merge_group_IDs <- names(list_vec_merge_idx)  <- names(merge_specify)
    list_vec_merge_idx <- Filter(f=length, x = list_vec_merge_idx)
    if (length(list_vec_merge_idx)==0) stop("merge_specify: no matches")
  }

  list_list_seurat_to_merge <- list()
  for (name in names(list_vec_merge_idx)) {
    list_list_seurat_to_merge[[name]] <- list_seurat_obj[unlist(list_vec_merge_idx[[name]])] 
  }
  
  rm(list_seurat_obj)
  
  list_seurat_merged <- lapply(list_list_seurat_to_merge, function(list_seurat_to_merge) {
    
    # for (i in 1:length(list_seurat_to_merge)) {
    #   colnames(x=list_seurat_to_merge[[i]]@assays$RNA@scale.data) <- colnames(x=list_seurat_to_merge[[i]]@assays$RNA@data) <- 
    #     colnames(x=list_seurat_to_merge[[i]]@assays$RNA@counts) <- paste0(i, "_", colnames(list_seurat_to_merge[[i]]))
    #   if (i == 1) {
    #     seurat_merged <- list_seurat_to_merge[[i]]
    #   } else {
        seurat_merged <- merge(x=list_seurat_to_merge[[1]],#x = seurat_merged, 
                               y = list_seurat_to_merge[2:length(list_seurat_to_merge)],
                               add.cell.ids=paste0(names(list_seurat_to_merge)),
                               project = paste0(prefix_data, "_", prefix_run),
                               merge.data=F,
                               min.cells = -Inf,
                               min.genes = -Inf,
                               do.normalize = F,
                               do.scale=F,
                               do.center=F)
                                     #add.cell.id2 = names(list_seurat_to_merge)[i])
    #   }
    # }
    return(seurat_merged) 
  })

  rm(list_list_seurat_to_merge)
  
  names(list_seurat_merged) <- merge_group_IDs
  list_seurat_obj <- list_seurat_merged
  rm(list_seurat_merged)
  
  list_seurat_obj <- mapply(function(seurat_obj, merge_ID)
  { seurat_obj$merge_ID <- rep(merge_ID,  times=ncol(GetAssayData(object=seurat_obj)))
  return(seurat_obj)
  }, 
  seurat_obj = list_seurat_obj,
  merge_ID = names(list_seurat_obj), SIMPLIFY=F)
}


######################################################################
########################## NORMALISE #################################
######################################################################

fun = function(seurat_obj) NormalizeData(object=seurat_obj)
list_iterable = list("X"=list_seurat_obj)
list_seurat_obj <- lapply(FUN=fun, "X"=list_iterable[[1]])

######################################################################
#################### FIND HIGHLY VAR GENES ###########################
######################################################################

outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_log_FindVariableFeatures.txt")

fun = function(seurat_obj) tryCatch({FindVariableFeatures(object = seurat_obj,
                                                          nfeatures = nAnchorFeatures,
                                                          selection.method="vst", verbose=T)}, 
                                    error= function(err) {
                                      message("FindVariableFeatures failed with selection.method='vst', trying with 'mean.var.plot'")
                                      FindVariableFeatures(object = seurat_obj,
                                                           nfeatures = nAnchorFeatures,
                                                           selection.method="mean.var.plot",verbose=T)
                                    })
list_iterable = list("X"=list_seurat_obj)
list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

######################################################################
######################### DO INTEGRATION #############################
######################################################################

if ((!is.null(align_group_IDs) | !is.null(align_specify)) & length(list_seurat_obj)>1) {
  
  message("Integrating datasets")
  
  list_seurat_align = NULL
  
  if (!is.null(align_group_IDs)) {
    if ('all' %in% align_group_IDs) {
      list_vec_align_idx <- list('all'=1:length(list_seurat_obj))
    } else {
      # a list of vectors giving the idx within the list_seurat_obj of seurat_obj to merge
      list_vec_align_idx <- lapply(align_group_IDs, function(id) grep(pattern=paste0("^",id), 
                                                                      x = names(list_seurat_obj),
                                                                      ignore.case=T))
      names(list_vec_align_idx) <- align_group_IDs
      
      list_vec_align_idx <- Filter(f=length,x=list_vec_align_idx)
      if (length(list_vec_align_idx)==0) stop("align_group_IDs: no matches")
    }
  } else if (!is.null(align_specify)) {
    list_vec_align_idx <- lapply(align_specify, function(vec_sample_ids) {
      sapply(vec_sample_ids, function(id) {
        grep(pattern=paste0("^",id), x = names(list_seurat_obj), ignore.case = T)
      })
    })
    align_group_IDs <- names(list_vec_align_idx) <- names(align_specify)
    list_vec_align_idx <- Filter(f=length, x=list_vec_align_idx)
    if (length(list_vec_align_idx)==0) stop("align_specify: no matches")
  }
  
  # a list of lists of seurat objects to align
  group_list_seurat_obj_tmp <- lapply(list_vec_align_idx, function(vec_align_idx) {list_seurat_obj[vec_align_idx]}) 
  names(group_list_seurat_obj_tmp) <- names(list_vec_align_idx)
    
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_log_FindAchors.txt")
  fun = function(list_seurat_obj_tmp, name) {
    tryCatch({
      FindIntegrationAnchors(object.list = list_seurat_obj_tmp, 
                             anchor.features = nAnchorFeatures,
                             dims = 1:n_comp, 
                             k.anchor=k.anchor,
                             k.filter=k.filter,
                             k.score=k.score,
                             verbose = T)
                             
      
    }, error=function(err) {
      message(paste0(name, ": FindIntegrationAnchors failed with error ", err))
      idx_ok <- sapply(list_seurat_obj_tmp, function(seu){ncol(seu)>=max(k.anchor, k.filter, k.score)})
      names_rm <- names(list_seurat_obj_tmp)[!idx_ok]
      message(paste0("possible cause: cell count in ", paste0(names_rm, collapse=", "), " is less than k, the number of neighbours"))
      message("trying to remove samples with insufficient number of cells")
      list_seurat_obj_tmp <- list_seurat_obj_tmp[idx_ok]
      tryCatch({
        FindIntegrationAnchors(object.list = list_seurat_obj_tmp,
                               anchor.features = nAnchorFeatures,
                               dims = 1:n_comp, 
                               k.anchor=k.anchor,
                               k.filter=k.filter,
                               k.score=k.score,
                               verbose=T)}, 
               error = function(err1) {
                  message(paste0("FindIntegrationAnchors failed again with error ", err1, ", returning NULL"))
                  return(NULL)
        })
    })
  }
  
  #list_iterable = list(list_seurat_obj_tmp = group_list_seurat_obj_tmp, name = names(group_list_seurat_obj_tmp))
  message("Finding achors")  
  # parallelising seems to make it freeze at "Scaling features for provided objects"
  list_anchors <- mapply(FUN = fun, list_seurat_obj_tmp = group_list_seurat_obj_tmp, name = names(group_list_seurat_obj_tmp), SIMPLIFY=F)#safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
  
  names(list_anchors) <- names(group_list_seurat_obj_tmp)
 
  # filter out those groups that failed
  list_anchors <- list_anchors[!all(sapply(list_anchors, is.null))]
  
  message("Integrating data")
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_log_IntegrateData.txt")
  fun = function(anchors, name) {
    tryCatch({
      IntegrateData(anchorset = anchors, 
                    dims = 1:n_comp, 
                    k.weight=k.weight,
                    verbose = T)
    }, error=function(err) {
      message(paste0(name, ": dropped because IntegrateData failed with error ", err))
      return(NULL)
    })
  }
  #list_iterable = list(anchors = list_anchors, name = names(list_anchors))
  
  list_seurat_align <- mapply(FUN = fun, anchors = list_anchors, name = names(list_anchors), SIMPLIFY=F)
  #list_seurat_align <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile, n_comp=n_comp)
  names(list_seurat_align) <- align_group_IDs
  
  # Filter again
  list_seurat_obj <- list_seurat_align[!sapply(list_seurat_align, is.null)]
  rm(list_seurat_align, group_list_seurat_obj_tmp)
  
}

if (!is.null(align_group_IDs) | !is.null(align_specify)) {
  list_seurat_obj <- lapply(list_seurat_obj, function(seu){DefaultAssay(object = seu) <- "integrated"; return(seu)})
}

######################################################################
######################## COMPUTE CELL CYCLE SCORE ####################
######################################################################
# See https://satijalab.org/seurat/cell_cycle_vignette.html

if (!is.null(paths_cellCycleGenes)) {
  if (all(!is.na(list_cellCycleGenes)) & 
      any(c('S.Score', 'G2M.Score','CC.Difference') %in% vars.to.regress)) {
  
  # Compute cell cycle scores using provided genesets
  fun = function(seurat_obj, name) {
    tryCatch({
      seurat_obj <- CellCycleScoring(object = seurat_obj, 
                                     s.genes = if (!is.null(list_cellCycleGenes[["s.genes"]])) list_cellCycleGenes[["s.genes"]] else NULL, 
                                     g2m.genes = if (!is.null(list_cellCycleGenes[["g2m.genes"]])) list_cellCycleGenes[["g2m.genes"]] else NULL, 
                                     set.ident = F)
      
      if (!is.null(list_cellCycleGenes[["s.genes"]]) & is.null(list_cellCycleGenes[["g2m.genes"]])) {
        seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
        }
        
      return(seurat_obj)
      }, error=function(err) {
          warning(paste0("CellCycleScoring failed for", name))
          seurat_obj$CC.Difference <- seurat_obj$S.Score <- seurat_obj$G2M.Score <- rep(NA, nrow(seurat_obj@meta.data))
          return(seurat_obj)
        })
  }
  list_iterable=list("seurat_obj"=list_seurat_obj, name=names(list_seurat_obj))
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_log_CellCycleScoring.txt")
  list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable, timeout=1200, outfile=outfile)
  
  # Visualize the distribution of cell cycle markers 
  invisible(sapply(function(seurat_obj,name) {
    try({
      
      features.plot <- c(if (!is.null(list_cellCycleGenes[["s.genes"]])) list_cellCycleGenes[["s.genes"]][1:min(3, length(list_cellCycleGenes[["s.genes"]]))] else NULL, 
                         if (!is.null(list_cellCycleGenes[["g2m.genes"]])) list_cellCycleGenes[["g2m.genes"]][1:min(3, length(list_cellCycleGenes[["g2m.genes"]]))] else NULL) 
      
      p<-RidgePlot(object = marrow, 
              features.plot = features.plot,
              nCol = min(1,length(features.plot)%/%2))
      
      saveMeta(savefnc=ggsave, plot=p,filename = paste0(dir_plots, prefix_data, "_", prefix_run, "_", name, "_cellcycle_RidgePlot.pdf"), width=15, height=10)
    
    })
  }, seurat_obj=list_seurat_obj, name=names(list_seurat_obj)))
  
  }
}
######################################################################
########################## SCALE DATA ################################
######################################################################

message("scaling data")
fun = function(seurat_obj, name) {
  
  # Check that the vars.to.regress are in the seurat object and that they are not all NA
  vec_logicVarsOK <- sapply(vars.to.regress, function(eachVar){
    if (!is.null(seurat_obj[[eachVar]])) !all(sapply(seurat_obj[[eachVar]], is.na)) else F
  })
  if (!all(vec_logicVarsOK)) warning(paste0("could not regress out ", vars.to.regress[!vec_logicVarsOK]))
  seurat_obj <- ScaleData(object = seurat_obj, 
                          vars.to.regress = vars.to.regress[vec_logicVarsOK],
                          block.size=15000,
                          min.cells.to.block = 10000,
                          verbose = T)
}
list_iterable=list("seurat_obj" = list_seurat_obj, 
          "name" = names(list_seurat_obj))
outfile = paste0(dir_log, prefix_data, "_", prefix_run, "_log_ScaleData.txt")
#list_seurat_obj<-safeParallel(fun=fun, list_iterable=list_iterable,outfile=outfile)
list_seurat_obj <- mapply(FUN=fun, "seurat_obj"=list_iterable[[1]], "name" = list_iterable[[2]],SIMPLIFY=F)

######################################################################
#################################### PCA #############################
######################################################################

#if (!is.null(align_group_IDs) | !is.null(align_specify) | !is.null(merge_group_IDs) | !is.null(merge_specify)) {
message("Computing PCA")
fun =  function(seurat_obj, obj_name) {
    tryCatch({
      # precaution for numerical stability in case of very few variable features or cells
      pcs.compute1 = min(n_comp, min(ncol(seurat_obj), length(VariableFeatures(object=seurat_obj)))%/%2)
      seurat_obj<-RunPCA(object = seurat_obj, 
             weight.by.var = F,
             npcs = pcs.compute1, 
             seed.use=randomSeed,
             verbose=F)
      # project to get loadings for all genes, not just variable
      seurat_obj <- ProjectDim(object = seurat_obj)
    },   error = function(err1) {
      pcs.compute2=min(n_comp, min(ncol(seurat_obj), length(VariableFeatures(object=seurat_obj)))%/%3)
      warning(paste0("RunPCA failed using ", pcs.compute1, " components with error ", err1, ", maybe due to cell count: ", ncol(seurat_obj), ". Trying with ", pcs.compute2, " components"))
      tryCatch({
        seurat_obj<-RunPCA(object = seurat_obj, 
                       weight.by.var = F,
                       npcs = pcs.compute2, 
                       seed.use=randomSeed,
                       verbose=F)
        # project so we get loadings for all genes, not just variable 
        ProjectDim(object = seurat_obj)
      }, error = function(err2) {
                         warning(paste0(obj_name, ": PCA failed again with ", pcs.compute2, " components with error: ", err2, ". Maybe due to cell count: ", ncol(seurat_obj), ". Returning seurat object without PCA"))
                         seurat_obj})
    })}
list_iterable = list("seurat_obj"=list_seurat_obj, "obj_name"=names(list_seurat_obj))
outfile=paste0(dir_log,prefix_data,"_",prefix_run,"_PCA2.txt")
list_seurat_obj<- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)
#}  

######################################################################
############################ JACKSTRAW ###############################
#####################################################################

# TODO: what random seed does JackStraw() use? the system seed?

if (use_jackstraw) {
  message("Using Jackstraw resampling to determine significant principal components")
  pvalThreshold = 0.05
  list_PC_signif_idx <- mapply(function(seurat_obj, name) {
    tryCatch({
    PC_signif_idx <- rep(TRUE, ncol(Loadings(object=seurat_obj, reduction="pca"))) # default
      
    prop.freq <- max(0.016, round(4/length(VariableFeatures(object=seurat_obj)),3)) # to ensure we have at least 3 samples so the algorithm works well
    # see https://github.com/satijalab/seurat/issues/5
    pcs.compute = ncol(Loadings(object=seurat_obj, reduction="pca"))
    seurat_obj <- JackStraw(object = seurat_obj, 
                            reduction="pca",
                            dims = pcs.compute,
                            num.replicate = 300, 
                            #display.progress = T,
                            #do.par = T,
                            verbose=T,
                            #num.cores = n_cores,
                            prop.freq = prop.freq,
                            maxit=2000)

    seurat_obj <- ScoreJackStraw(object=seurat_obj, reduction="pca", dims=1:pcs.compute, do.plot=F) # do we need to run JackStrawPlot to get the ggplot object?
    p <- JackStrawPlot(object=seurat_obj, dims=1:pcs.compute, reduction="pca")
    saveMeta(savefnc=ggsave, plot=p, filename = paste0(dir_plots, prefix_data, "_", prefix_run, "_", name, "_JackStrawPlot.pdf"), width=10, height=10)
    PC_signif_idx <- seurat_obj@reductions$pca@jackstraw$overall.p.values[,2]< pvalThreshold
    return(PC_signif_idx)
    }, error= function(err) {
      message("JackStraw failed on ", name)
      message("Returning TRUE for all PCs")
      rep(TRUE, n_comp)
    })
  }, seurat_obj=list_seurat_obj, name=names(list_seurat_obj), SIMPLIFY = F)
} else {
  list_PC_signif_idx <- lapply(1:length(list_seurat_obj), function(seurat_obj){
    rep(TRUE, n_comp)
  })
} 
######################################################################
############################ T-SNE ###################################
######################################################################

message("Computing t-SNE")

fun =  function(seurat_obj, obj_name, PC_signif_idx) {
  
  perplexity1 = min(30, min(n_comp,ncol(GetAssayData(object=seurat_obj)))%/%1.5)
  tryCatch({
    RunTSNE(object = seurat_obj, 
                    tsne.method="Rtsne",
                    reduction = "pca",#if (is.null(align_group_IDs) & is.null(align_specify)) "pca" else "cca.aligned", 
                    dims= (1:n_comp)[PC_signif_idx][1:min(30,sum(PC_signif_idx))], # no need to use all PCs for t-SNE
                    seed.use = randomSeed,
                    #do.fast=T,
                    perplexity=perplexity1,
                    check_duplicates=F)
  }, error = function(err1) {
    perplexity2 = min(30, min(n_comp,ncol(GetAssayData(object=seurat_obj)))%/%2.5)
    warning(paste0(obj_name, ": RunTSNE failed with error: ", err1, ", using perplexity ", perplexity1, ". Maybe due to cell count: ", dim(GetAssayData(object=seurat_obj))[2], ". Trying with perplexity ", perplexity2))
    tryCatch({RunTSNE(object = seurat_obj, 
                      reduction = "pca",#if (is.null(align_group_IDs) & is.null(align_specify)) "pca" else "cca.aligned", 
                      dims= (1:n_comp)[PC_signif_idx][1:min(30,sum(PC_signif_idx))], # no need to use all PCs for t-SNE
                      #do.fast=T,
                      seed.use = randomSeed,
                      perplexity = perplexity2,
                      check_duplicates=F)
    }, error= function(err2) {
      warning(paste0("RunTSNE failed again with error ", err2, ", returning original seurat object"))
      seurat_obj})
    })
  }
outfile=paste0(dir_log,prefix_data,"_",prefix_run,"_TSNE2.txt")
list_iterable=list("seurat_obj"=list_seurat_obj,
          "obj_name" = names(list_seurat_obj), 
          "PC_signif_idx" = list_PC_signif_idx)
list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile)

######################################################################
############################ FIND CLUSTERS ###########################
######################################################################

if (!is.null(res_primary)) {
  
  message("Finding neighbors")
  fun = function(seurat_obj, PC_signif_idx) {
    FindNeighbors(object=seurat_obj, 
                  dims = (1:n_comp)[PC_signif_idx][1:min(30,sum(PC_signif_idx))],
                  verbose=T)
  }
  list_iterable = list("seurat_obj"=list_seurat_obj, "PC_signif_idx"=list_PC_signif_idx)
  list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable)
  
  message("Computing clusters")
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_FindClusters.txt")
  reduction <- "pca"#if (is.null(align_group_IDs) & is.null(align_specify)) "pca" else "cca.aligned"

  if (!is.null(res_to_calc)) {
    
    list_seurat_obj <- mapply(function(seurat_obj,
                                       PC_signif_idx) {
      # Run FindClusters in parallel using different resolution parameter values
      fun = function(res) { 
        FindClusters(object = seurat_obj,
                     #reduction.type = reduction.type,
                     #dims.use = if (reduction.type=="cca.aligned") 1:n_CC else (1:ncol(Embeddings(object = seurat_obj, reduction = reduction.type)))[PC_signif_idx],
                     #save.SNN = T,
                     random.seed = randomSeed,
                     verbose=T,
                     resolution = res)}
      
      #list_iterable= list("X"=res_to_calc)
      
      # list_seurat_obj_res <- safeParallel(fun=fun, 
      #                                      list_iterable=list_iterable, 
      #                                      outfile=outfile, 
      #                                      seurat_obj=seurat_obj,  
      #                                      randomSeed=randomSeed,
      #                                      timeout= 600)
      
      list_seurat_obj_res <- lapply(X = res_to_calc,FUN = fun)
                                 
      # Save all the cluster assignments as meta data in the original seurat object
      for (i in seq(1:length(res_to_calc))) {
        seurat_obj[[paste0("RNA_snn_res.", res_to_calc[i])]] <- Idents(object=list_seurat_obj_res[[i]])
      }
      # No need to have this massive list in session!
      rm(list_seurat_obj_res)
      
      Idents(object = seurat_obj) <- paste0("RNA_snn_res.", res_primary)
      
      return(seurat_obj)
    }, seurat_obj = list_seurat_obj, PC_signif_idx=list_PC_signif_idx, SIMPLIFY=F)

  } else {
    fun = function(seurat_obj, PC_signif_idx) {
     
      FindClusters(object = seurat_obj,
                  #reduction.type = reduction.type,
                  #dims = (1:ncol(Embeddings(object = seurat_obj, reduction = reduction.type)))[PC_signif_idx],
                  #dims.use = if (reduction.type=="cca.aligned") 1:n_CC else (1:ncol(Embeddings(object = seurat_obj, reduction = reduction.type)))[PC_signif_idx],
                  #print.output = 0,
                  verbose=T,
                  random.seed= randomSeed,
                  #save.SNN = T,
                  resolution = res_primary)

      }
    list_iterable = list("seurat_obj"=list_seurat_obj, "PC_signif_idx"=list_PC_signif_idx)
    list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable, outfile=outfile, res_primary=res_primary, randomSeed=randomSeed)
  }
}

######################################################################
######################### TRANSFER LABELS ############################
######################################################################

if (!is.null(path_transferLabelsRef)) {
  seuratObjRef <- load_obj(path_transferLabelsRef)
  
  if (length(seuratObjRef@assays$RNA@var.features)==0) {
    seuratObjRef <- tryCatch({
      FindVariableFeatures(object = seuratObjRef,
                           nfeatures = nAnchorFeatures,
                           selection.method="vst", 
                           verbose=T)}, 
                          error= function(err) {
                            message("FindVariableFeatures failed with selection.method='vst', trying with 'mean.var.plot'")
                            FindVariableFeatures(object = seuratObjRef,
                                                 nfeatures = nAnchorFeatures,
                                                 selection.method="mean.var.plot", 
                                                 verbose=T)
                          })
  }
  
  if (is.null(seuratObjRef@meta.data$percent.mito) & is.null(seuratObjRef@meta.data$percent.ribo)) {
    # Compute and add percent mito and percent ribo as metadata
    mito.genes <- grepl(pattern = "^mt-", x = rownames(seuratObjRef), ignore.case=T)
    ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(seuratObjRef), ignore.case=T)
    colSums_tmp <- colSums(x = seuratObjRef@assays$RNA@counts)
    
    metadata <- data.frame(percent.mito=colSums(x = seuratObjRef@assays$RNA@counts[mito.genes,])/colSums_tmp, 
                           percent.ribo = colSums(x = seuratObjRef@assays$RNA@counts[ribo.genes,])/colSums_tmp, 
                           rownames=colnames(seuratObjRef))
    
    seuratObjRef <- AddMetaData(object=seuratObjRef, metadata=metadata)
  }
  
  if (all(dim(seuratObjRef@assays$RNA@data)==0)) {
    seuratObjRef <- NormalizeData(seuratObjRef)
  }
  
  # to avoid errors, scale, find var features and compute PCA and JackStraw afresh
  seuratObjRef <- FindVariableFeatures(seuratObjRef,
                                       nfeatures = nAnchorFeatures,
                                       selection.method="vst", 
                                       verbose=T)
  
  seuratObjRef <- ScaleData(seuratObjRef,
                            #vars.to.regress=vars.to.regress,
                            block.size=15000,
                            min.cells.to.block=5000,
                            verbose=T)
  
  seuratObjRef <- RunPCA(object= seuratObjRef,
                         npcs= n_comp,
                         weight.by.var=F,
                         verbose=F,
                         seed.use = randomSeed)
  
  if (use_jackstraw) {
    tryCatch({
      
      #TODO: debug
      # |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 51s
      # Error in fake.vals.raw[[y]][, x] : incorrect number of dimensions
      # In addition: There were 23 warnings (use warnings() to see them)
      # warnings()
      # 3: In r1.use:r2.use : numerical expression has 30 elements: only the first used
      # [...]
      # 13: In 1:dims : numerical expression has 30 elements: only the first used
      
      seuratObjRef <- JackStraw(seuratObjRef,
                              reduction="pca", 
                              dims = 1:n_comp, 
                              num.replicate = 300, 
                              prop.freq = 0.01, 
                              verbose=T, 
                              maxit=1000)
      
      seuratObjRef <- ScoreJackStraw(seuratObjRef, dims=1:n_comp, do.plot=F)
      
      PC_signif_idx_ref <- seuratObjRef@reductions$pca@jackstraw$overall.p.values[,2]< pvalThreshold
      
    }, error = function(err1) {
      PC_signif_idx_ref <- rep(TRUE, n_comp)
    })
  } else {
    PC_signif_idx_ref <- rep(TRUE, n_comp)
  }
   
  fun = function(seuratObjQuery, seuratObjName, PC_signif_idx) {
    anchors <- tryCatch({
      # https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html
      # In data transfer, Seurat has an option (set by default) to project the 
      # PCA structure of a reference onto the query, instead of learning a joint structure with CCA. 
      # We generally suggest using this option when projecting data between scRNA-seq datasets.
      # We use this option 
      FindTransferAnchors(reference = seuratObjRef, 
                         query = seuratObjQuery,
                         reduction = "pcaproject",
                         project.query = F, # i.e. project reference
                         npcs = NULL,
                         l2.norm=TRUE,
                         dims= (1:ncol(Embeddings(object=seuratObjRef, reduction="pca")))[PC_signif_idx_ref], 
                         k.anchor = k.anchor,
                         k.filter = k.filter,
                         k.score = k.score,
                         max.features=max.features, #default
                         verbose=T)
      }, error = function(err) {
       message(paste0("FindTransferAnchors failed with error ", err))
        if (ncol(seuratObjQuery)<=max(k.anchor, k.filter, k.score)) {
          c(k.anchor, k.filter, k.score)[which.max(k.anchor, k.filter, k.score)] <- ncol(seuratObjQuery)
        }
        FindTransferAnchors(reference = seuratObjRef, 
                            query = seuratObjQuery, 
                            project.query = T, # this is the second difference
                            reduction = "pcaproject",
                            project.query = F, # i.e. project reference
                            npcs = NULL,
                            l2.norm=TRUE,
                            dims= (1:ncol(Embeddings(object=seuratObjRef, reduction="pca")))[PC_signif_idx_ref], # Which dimensions to use from the reduction to specify the neighbor search space. Does this refer to ref or combined?
                            k.anchor = k.anchor,
                            k.filter = k.filter,
                            k.score = k.score,
                            max.features=max.features, #default
                            verbose=T)
     })

    predictions <- TransferData(anchorset = anchors, 
                                refdata = seuratObjRef[[colLabels]][,1], 
                                reduction = "pcaproject",
                                l2.norm=TRUE,
                                dims= (1:ncol(Embeddings(object=seuratObjRef, reduction="pca")))[PC_signif_idx_ref], 
                                k.weight = k.weight,
                                sd.weight = sd.weight,
                                verbose=T)
    
    seuratObjQuery <- AddMetaData(object=seuratObjQuery, 
                                metadata=predictions)
    
    Idents(seuratObjQuery) <- seuratObjQuery$predicted.id
    ### EDITED 20190516 ###
    # Plot prediction score
    p <- FeaturePlot(object=seuratObjQuery, feature = "prediction.score.max")
    saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_", seuratObjName, "_transferlabel_predictionScoreMax.pdf"), w=12, h=12)
    ###
    
    # Do not subset out; rather include the prediction score in the label and user can subset manually
    if (F) {
    # Subset out cells with poor predictions
      if (sum(seuratObjQuery@meta.data$prediction.score.max < minPredictionScore)>0) {
        #idx <- seuratObjQuery@meta.data$prediction.score.max >= minPredictionScore
        
        # Plot data before filtering on prediction.score.max
        seuratObjQuery <- ScaleData(seuratObjQuery,
                                    #vars.to.regress=vars.to.regress,
                                    do.scale=T,
                                    do.center=T,
                                    block.size=15000,
                                    min.cells.to.block=5000,
                                    verbose=T)
        
        seuratObjQuery <- RunPCA(object= seuratObjQuery,
                                 npcs= min(n_comp, 30),
                                 weight.by.var=F,
                                 verbose=F,
                                 seed.use = randomSeed)
        
        perplexity1 = min(30, ncol(GetAssayData(object=seuratObjQuery))%/%1.5)
        
        seuratObjQuery <- RunTSNE(object = seuratObjQuery, 
                                  reduction = "pca", 
                                  dims = 1:min(n_comp, 30),#(1:n_comp)[PC_signif_idx][1:min(30,sum(PC_signif_idx))], # no need to use all PCs for t-SNE
                                  seed.use = randomSeed,
                                  check_duplicates=F,
                                  #do.fast=T,
                                  perplexity=perplexity1)
        
        # Plot 
        p <- DimPlot(object=seuratObjQuery, 
                      label = T, 
                      pt.size = 1, 
                      no.legend = F, 
                      group.by = "predicted.id",
                      plot.title = paste0(seuratObjName, " by transferred label")) # + xlab("t-SNE 1") + ylab("t-SNE 2"
        saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_", seuratObjName, "_tSNE_clust_transferlabel_preFilt.pdf"), w=12, h=12)
        
        p <- FeaturePlot(object=seuratObjQuery, feature = "prediction.score.max")
        saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_", seuratObjName, "_transferlabel_predictionScore_preFilt.pdf"), w=12, h=12)
        
        ## Now remove the cells that didn't align well
        seuratObjRm <- subset(x= seuratObjQuery, 
                                 subset = prediction.score.max<minPredictionScore)
        # save the bad cells before removing them
        saveMeta(savefnc=saveRDS,seuratObjRm, file = paste0(dir_RObjects, prefix_data, "_", prefix_run, "_transferlabel_prediction.score.max_sub_", minPredictionScore, "_seurat_obj.RDS.gz"), compress="gzip")
        rm(seuratObjRm)
        
        # Overwrite the full seurat object with the subset of cells that aligned well
        seuratObjQuery <- subset(x= seuratObjQuery, 
                                 subset = prediction.score.max>=minPredictionScore)
      
       #Idents(seuratObjQuery) <- seuratObjQuery@meta.data$predicted.id
  
        seuratObjQuery <- ScaleData(seuratObjQuery,
                                  #vars.to.regress=vars.to.regress,
                                  do.scale=T,
                                  do.center=T,
                                  block.size=15000,
                                  min.cells.to.block=5000,
                                  verbose=T)
        
        seuratObjQuery <- RunPCA(object= seuratObjQuery,
                               npcs= n_comp,
                               weight.by.var=T,
                               verbose=F,
                               seed.use = randomSeed)
        
        perplexity1 = min(30, ncol(GetAssayData(object=seuratObjQuery))%/%1.5)
        
        seuratObjQuery <- RunTSNE(object = seuratObjQuery, 
                                  reduction = "pca", 
                                  dims = (1:n_comp)[PC_signif_idx][1:min(30,sum(PC_signif_idx))], # no need to use all PCs for t-SNE
                                  seed.use = randomSeed,
                                  check_duplicates=F,
                                  #do.fast=T,
                                  perplexity=perplexity1)
      }
    }
    # TODO: check whether use of PC_signif_idx is appropriate here
    return(seuratObjQuery)
  }
  #list_iterable = list("X"=list_seurat_obj)
  list_seurat_obj <- mapply(FUN = fun, "seuratObjQuery"=list_seurat_obj, 
                            seuratObjName = names(list_seurat_obj), 
                            "PC_signif_idx"=listPC_signif_idx, SIMPLIFY = F)
  #list_seurat_obj <- safeParallel(fun=fun, list_iterable=list_iterable)
}

#######################################################################
########################## T-SNE PLOTS ################################
#######################################################################

message("Plotting t-SNE")
invisible(mapply(function(seurat_obj, name) {
  p1 <- DimPlot(object=seurat_obj,  
                reduction="tsne",
                dims=c(1,2),
                #do.return = T, 
                label=F,
                pt.size = 1, 
                group.by = "sample_ID", 
                no.legend=F, 
                plot.title=paste0(name, " by sample"))
  saveMeta(savefnc=ggsave, plot=p1, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_tSNEPlot_sample.pdf"), w=10, h=10)
  
  if (!is.null(res_primary)) {
    p2 <- DimPlot(object=seurat_obj, 
                  label = T, 
                  pt.size = 1, 
                  no.legend = F, 
                  group.by = paste0("RNA_snn_res.",res_primary),
                  plot.title = paste0(name, " by cluster")) # + xlab("t-SNE 1") + ylab("t-SNE 2"
    saveMeta(savefnc=ggsave, plot=p2, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_", name, "_tSNEPlot_clust.pdf"), w=10, h=10)
  }
  
  if (!is.null(seurat_obj@meta.data$predicted.id)) {
    p3 <- DimPlot(object=seurat_obj, 
                  label = T, 
                  pt.size = 1, 
                  no.legend = F, 
                  group.by = "predicted.id",
                  plot.title = paste0(name, " by transferred label")) # + xlab("t-SNE 1") + ylab("t-SNE 2"
    saveMeta(savefnc=ggsave, plot=p3, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_", name, "_tSNEPlot_clust_transferlabel.pdf"), w=12, h=12)
    
    p3 <- FeaturePlot(object=seurat_obj, feature = "prediction.score.max")
    saveMeta(savefnc=ggsave, plot=p3, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_featurePlot_clust_transferlabel_predictionScore.pdf"), w=12, h=12)
    
  }
  # p1 <- TSNEPlot(seurat_obj, do.return = T, pt.size = 1, group.by = "condition", no.legend=F, plot.title=paste0(name, " by condition"))
  # ggsave(p1, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_tSNEPlot_condition.pdf"), w=10, h=10)
  
}, seurat_obj = list_seurat_obj,
name = names(list_seurat_obj),
SIMPLIFY=F
))

######################################################################
####################### RESET DEFAULT ASSAY TO RNA ###################
######################################################################
# see https://satijalab.org/seurat/faq.html 4.

if (!all(sapply(c(align_group_IDs,align_specify),is.null))) {
  list_seurat_obj <- lapply(list_seurat_obj, function(seurat_obj) {DefaultAssay(seurat_obj) <- "RNA"; return(seurat_obj)})
}

######################################################################
####################### FIND CLUSTER MARKERS #########################
######################################################################

if (!is.null(res_primary)){
  
  message("Computing differentially expressed genes as cluster markers")
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_FindMarkers.txt")
  # for each final seurat object, for each cluster, a dataframe of markers
  list_list_markers <- lapply(list_seurat_obj, function(seurat_obj) {
    clusters = names(table(Idents(seurat_obj)))
    list_iterable = list("X"=clusters)
    fun = function(cluster) {tryCatch({
      FindMarkers(seurat_obj,  
                  #cells.1=colnames(seurat_obj)[Idents(seurat_obj)==cluster],
                  #cells.2=NULL,
                  ident.1 = cluster,
                  only.pos = T,
                  #ident.2 = clusters[clusters!=cluster],
                  test.use  ="MAST",
                  max.cells.per.ident=1000,
                  random.seed=randomSeed,
                  #latent.vars = if (!is.null(merge_specify) | !is.null(merge_group_IDs)) "sample_ID" else NULL,
                  verbose = T)
    }, 
    error = function(err) {
      NA_character_
    })}
    list_markers=NULL
    list_markers <- lapply(FUN=fun, "X"=list_iterable[[1]])#, simplify=F, outfile=outfile, seurat_obj=seurat_obj, res_primary=res_primary)
    #list_markers <- safeParallel(fun=fun, list_iterable=list_iterable, simplify=F, outfile=outfile, seurat_obj=seurat_obj, res_primary=res_primary)
  })
  
  # rbind the outputs from each cluster into a single df per highest-level seurat_obj, adding a 'cluster' column.
  list_markers <- mapply(function(list_markers, seurat_obj) {
    list_markers <- mapply(function(df_markers, cluster) {
      if (!all(sapply(df_markers, is.na))) {
        cbind("gene" = rownames(df_markers), "cluster"=rep(cluster, nrow(df_markers)), df_markers)
      } else {
        NA_character_
      }
    },df_markers=list_markers, cluster=names(table(Idents(seurat_obj))), SIMPLIFY=F)
    list_markers <- list_markers[!sapply(list_markers, function(markers) all(is.na(markers)))]
    markers <- Reduce(x=list_markers, f=rbind)
    rownames(markers) <- NULL
    return(markers)
  }, list_markers = list_list_markers, seurat_obj = list_seurat_obj, SIMPLIFY = F)
  
  names(list_markers) <- names(list_seurat_obj)
}

######################################################################
######################## FEATUREPLOTS ################################
######################################################################

message("Making featureplots")
if (!is.null(feats_to_plot)) {
  invisible(mapply(function(seurat_obj, name) {
    if (feats_plot_separate) {
      lapply(feats_to_plot, function(feature) {
        if (any(grepl(feature, rownames(GetAssayData(object=seurat_obj)))) | any(grepl(feature, colnames(seurat_obj@meta.data)))) {
          tryCatch({
            #pdf(file = paste0(dir_plots, prefix_data,"_", prefix_run,"_", name, "_", feature,"_featurePlot.pdf"), w=8, h=8)
            p<-FeaturePlot(object=seurat_obj, 
                        #dims=c(1,2), 
                        #reduction=c("tsne"),
                        features = feature)#, 
                        #cols=c("grey95", "blue"), 
                        #label=T,
                        #pt.size = 1)
            saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data,"_", prefix_run,"_", name, "_", feature,"_featurePlot.pdf"), w=8, h=8)
          }, error = function(err) warning(paste0(name, ": FeaturePlot failed with error ", err)))
          #try(dev.off())
        } else { 
          warning(paste0(feature, " not found in ", name, " data"))}
        #ggsave(p1, filename =  paste0(dir_plots, prefix_data,"_",prefix_run,"_", name, "_", feature,"_featurePlot.pdf"), w=8, h=8)
      })
    } else {
      feats_found_idx <- sapply(feats_to_plot, function(feature) any(grepl(feature, rownames(GetAssayData(object=seurat_obj)))))
      if (all(feats_found_idx)) {
        tryCatch({
          p<- FeaturePlot(object=seurat_obj, 
                      features = feats_to_plot[1:2], 
                      #dims=c(1,2),
                      #reduction="tsne",
                      cols=c("grey92","blue", "violetred1"), 
                      blend=T,  
                      label=T,
                      pt.size = 1)
          saveMeta(savefnc=ggsave, plot=p,filename =  paste0(dir_plots, prefix_data,"_",prefix_run,"_", name, "_", paste0(feats_to_plot, collapse = "_"), "_featurePlot.pdf"), w=8, h=8)
        }, error = function(err) warning(paste0(name, ": FeaturePlot failed with error ", err)))
      } else {
        warning(paste0("FeaturePlot failed because ", paste0(feats_to_plot[!feats_found_idx], collapse = " ") , " weren't found in ", name))
      }
    }
  }, 
  seurat_obj = list_seurat_obj,
  name = names(list_seurat_obj),
  SIMPLIFY=F))
}

############################## SOUPX PLOTS ###########################

if (run_SoupX) {
  
  fun = function(seurat_obj, gene) {
    p<-FeaturePlot(seurat_obj, features.plot = gene, cols.use=c("grey95", "blue"), pt.size = 1, no.legend = F)
    saveMeta(savefnc=ggsave, plot=p, filename = paste0(dir_plots, prefix_data, "_", prefix_run, "_", sample_ID, "_", gene, "_SoupX_featurePlot.pdf"), width = 12, height=10)
  }
  
  for (sample_ID in sample_IDs) {
    
    topgenes_path <- paste0(dir_scratch, prefix_data, "_", prefix_run, "_", sample_ID, "_SoupX_topgenes.RDS")
    
    if (file.exists(topgenes_path)) {
      
      topgenes <- readRDS(file = topgenes_path)
      
      list_iterable = list("seurat_obj" = list_seurat_obj, 
                  "gene"=topgenes)
      
      safeParallel(fun=fun, 
                   list_iterable=list_iterable, 
                   dir_plots=dir_plots, 
                   prefix_data=prefix_data, 
                   prefix_run=prefix_run, 
                   sample_ID=sample_ID)
      
      rm(topgenes)
      file.remove(topgenes_path)
    }
  }
}

########################## VIOLIN PLOTS ##############################

message("making nCount_RNA and percent.mito and percent.ribo violin plots by sample")
invisible(mapply(function(seurat_obj, name) {
  p <- VlnPlot(object = seurat_obj, 
          features = "nCount_RNA", 
          #do.sort = F, 
          same.y.lims=T, 
          #legend.position = "bottom", 
          #single.legend = F, 
          #do.return=T, 
          group.by = "sample_ID")
  saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_nCount_RNA_VlnPlot_sample.pdf"), w=10, h=10)
  p<-VlnPlot(object = seurat_obj, 
          features = "percent.mito", 
          #do.sort = F, 
          same.y.lims=T, 
          #legend.position = "bottom", 
          #single.legend = F, 
          #do.return=T, 
          group.by = "sample_ID")
  saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_percent.mito_VlnPlot_sample.pdf"), w=10, h=10)
  p<-VlnPlot(object = seurat_obj, 
          features = "percent.ribo", 
          #do.sort = F, 
          same.y.lims=T, 
          #legend.position = "bottom", 
          #single.legend = F, 
          #do.return=T, 
          group.by = "sample_ID")
  saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_percent.ribo_VlnPlot_sample.pdf"), w=10, h=10)
  
}, seurat_obj = list_seurat_obj,
name = names(list_seurat_obj),
SIMPLIFY=F
))

#### Facet plots of samples in each cluster and clusters in each sample ####

message("plotting sample proportions in each cluster and cluster proportions for each sample")
invisible(mapply(function(seurat_obj, name){
  df = data.frame("sample_ID" =as.character(seurat_obj$sample_ID), "cluster"=as.character(Idents(seurat_obj)))
  
  p<-ggplot(df, aes(sample_ID, cluster)) + geom_bar(stat = "identity") + facet_wrap(.~cluster) 
  saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_samples_per_cluster.pdf"), w=24, h=15)
  
  p<- ggplot(df, aes(cluster, sample_ID)) + geom_bar(stat = "identity") + facet_wrap(.~sample_ID) 
    theme(strip.text.x = element_text(size=8, angle=75),
          strip.text.y = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="red", fill="#CCCCFF"))
  saveMeta(savefnc=ggsave, plot=p, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_cluster_per_sample.pdf"), w=24, h=15)
}, seurat_obj = list_seurat_obj, name=names(list_seurat_obj),SIMPLIFY=F))


if (!is.null(path_transferLabelsRef)) {
  message("plotting sample proportions in each cluster and cluster proportions for each sample")
  invisible(mapply(function(seurat_obj, name){
    df = data.frame("sample_ID" =as.character(seurat_obj$sample_ID), "predicted.id"=as.character(seurat_obj@meta.data$predicted.id))
    
    p<- ggplot(df, aes(sample_ID, predicted.id)) + geom_bar(stat = "identity") + facet_wrap(.~predicted.id) 
    saveMeta(savefnc=ggsave, plot=p,filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_samples_per_predicted.id.pdf"), w=24, h=15)
    
    p<- ggplot(df, aes(predicted.id, sample_ID)) + geom_bar(stat = "identity") + facet_wrap(.~sample_ID) 
    theme(strip.text.x = element_text(size=8, angle=75),
          strip.text.y = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="red", fill="#CCCCFF"))
    saveMeta(savefnc=ggsave, plot=p,filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_predicted.id_per_sample.pdf"), w=24, h=15)
  }, seurat_obj = list_seurat_obj, name=names(list_seurat_obj),SIMPLIFY=F))
}

######################################################################
########################## OUTPUT TABLES, ROBJECTS ###################
######################################################################

if (!is.null(res_primary)) {
  ############################## CLUSTER GENE MARKERS ##################
  message("Writing out gene markers")

  tryCatch({
    saveMeta(savefnc=openxlsx::write.xlsx, x = list_markers,
                          file = paste0(dir_tables, prefix_data, "_", prefix_run, "_markers.xlsx"))
    }, error =function(err){ 
      message(paste0("openxlsx::write.xlsx failed with error", err))
      invisible(mapply(function(markers,name) saveMeta(savefnc=write.csv, x=markers, file=paste0(dir_tables, prefix_data, "_", prefix_run, "_", name, "_markers.csv"), quote = F, row.names=T), markes=list_markers, name=names(list_seurat_obj), SIMPLIFY=F))

    })
  
}
  
########################## SAVE SEURAT OBJECTS #######################

outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_outputTables_RObjects.txt")
list_iterable=list("seurat_obj" = list_seurat_obj, 
          "name" = names(list_seurat_obj))
fun = function(seurat_obj, name) {
  invisible({
    saveMeta(savefnc=saveRDS, object= seurat_obj, file = paste0(dir_RObjects, prefix_data, "_", prefix_run, "_", name, "_seurat_obj.RDS.gz"), compress = "gzip")
    })
  }
invisible(safeParallel(fun=fun, list_iterable=list_iterable, timeout = 1200, outfile=outfile))

######################################################################
######################### LOG SESSION ################################
######################################################################

setwd(dir_log)
saveMeta(doPrint=T)
