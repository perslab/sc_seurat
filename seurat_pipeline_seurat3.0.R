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
              help = "Path(s) to metadata file(s) in one of the standard (compressed) character separated formats. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g. ''c('<path1>','<path2>')''. If not given, takes any metadata stored within the path_data object(s). [default %default]"),  
  make_option("--n_cells_loaded", type="character", default = 'c(9000)',
              help = "Approximately how many cells loaded in each sample? Takes a vector, in quotes, with one value per group of samples or a single common value. Used by doubletFinder to estimate doublet rate, [default %default]"),
  make_option("--n_cells_recovered", type="character", default = NULL,
              help = "Approximately how many cells recovered per sample? Takes a vector, in quotes, with one value per group of samples or a single common value. Used to cut off cells in QC. If left as NULL, a sensible value is computed on basis on n_cells_loaded [default %default]"),
  make_option("--use_filtered_feature_bc_matrix", type="logical", default = T,
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
  make_option("--percent.ribo_max", type="numeric", default = 0.15,
              help = "Maximum proportion of ribosomal genes tolerated in a cell, [default %default]"),
  make_option("--rm_sc_multiplets", type="logical", default = T,
              help = "Use DoubletFinder to remove suspected multiplets? [default %default]"),
  make_option("--vars.to.regress", type="character", default='c("nCount_RNA", "percent.mito", "percent.ribo")',
              help="Provide arguments to Seurat's ScaleData function in the form of a vector in quotes, defaults to ''c('nCount_RNA', 'percent.mito', 'percent.ribo')'' [default %default]"),
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
  make_option("--n_PC", type="integer", default = 50L,
              help = "How many principal components? [default %default]"),
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
  make_option("--minPredictionScore", type="numeric", default = 0.7,
              help = "If path_transferLabelsRef is given, minimum prediction score (0-1), [default %default]"),
  make_option("--feats_to_plot", type="character", default = 'c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo","Malat1")',
              help = "Features to plot. Format as a quoted character with a vector of values, without whitespace. [default %default]"),
  make_option("--feats_plot_separate", type="logical", default = T,
              help = "Plot features separately or in one plot? If FALSE, only the first two features will be plotted [default %default]"),
  make_option("--RAM_Gb_max", type="integer", default=200,
              help = "Upper limit on Gb RAM available. Taken into account when setting up parallel processes. [default %default]"),
  make_option("--path_runLog", type="character", default=NULL,
              help = "Path to file to log the run and the git commit. If left as NULL, write to a file called runLog.text in the dirLog [default %default]")
  
)

######################################################################
######################### UTILITY FUNCTIONS ##########################
######################################################################
  
source(file="/projects/jonatan/tools/functions-src/utility_functions.R")

######################################################################
########################### PACKAGES #################################
######################################################################

ipak(c("optparse", "Matrix", "Seurat", "ggplot2", "scales", "dplyr", "parallel", "reshape", "reshape2", "cowplot"))#, "pSI", "loomR", "doubletFinder")
stopifnot(as.character(packageVersion("Seurat"))=='3.0.0.9000')

######################################################################
########################### SC FUNCTIONS #############################
######################################################################

source(file="/projects/jonatan/tools/functions-src/functions_sc.R")

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
n_PC <- opt$n_PC
use_jackstraw <- opt$use_jackstraw
#n_CC <- opt$n_CC
res_primary <- opt$res_primary
res_to_calc <- opt$res_to_calc ; if (!is.null(res_to_calc)) res_to_calc <- eval(parse(text = res_to_calc))
path_transferLabelsRef <- opt$path_transferLabelsRef
colLabels <- opt$colLabels
minPredictionScore <- opt$minPredictionScore
feats_to_plot <- opt$feats_to_plot ; if (!is.null(feats_to_plot)) feats_to_plot <- eval(parse(text = feats_to_plot))
feats_plot_separate <- opt$feats_plot_separate
#n_cores <- opt$n_cores
RAM_Gb_max <- opt$RAM_Gb_max
path_runLog <- opt$path_runLog

######################################################################
######################## CONDITIONED PACKAGES  #######################
######################################################################

#if (!is.null(res_primary)) ipak(c("xlsx"))
if (!is.null(res_primary)) ipak(c("openxlsx"))

if (run_SoupX) ipak(c('SoupX', 'LTMGSCA'))
#'   
#'   ## https://github.com/constantAmateur/SoupX/blob/master/R/load10X.R
#'   
#'   #' Load a collection of 10X data-sets
#'   #'
#'   #' Loads unfiltered 10X data from each data-set and identifies which droplets are cells using the cellranger defaults.
#'   #'
#'   #' @export
#'   #' @param dataDirs Vector of top level cellranger output directories (the directory that contains the "raw_feature_bc_matrix" folder).
#'   #' @param channelNames To make droplet barcodes unique across experiment, each channel needs its own unique label.  If not given, this is set numerically.
#'   #' @param ... Extra parameters passed to \code{SoupChannel} construction function.
#'   #' @return A SoupChannelList object containing the count tables for each 10X dataset.
#'   #' @seealso SoupChannel SoupChannelList estimateSoup
#'   #' @importFrom Seurat Read10X
#'   
#'   load10X = function(dataDirs, 
#'                      use_filtered_feature_bc_matrix, 
#'                      sample_ID,
#'                      n_cells,
#'                      nCount_RNA_min, 
#'                      nCount_RNA_max,
#'                      dir_plots,
#'                      prefix_data,
#'                      prefix_run,
#'                      channelNames=NULL,...){
#'     if(is.null(channelNames))
#'       channelNames = sprintf('Channel%d',seq_along(dataDirs))
#'     channels = list()
#'     for(i in seq_along(dataDirs)){
#'       message(sprintf("Loading data for 10X channel %s from %s",channelNames[i],dataDirs[i]))
#'       dataDir = dataDirs[i]
#'       #Get reference
#'       #ref = list.files(file.path(dataDir))#,'raw_feature_bc_matrix'))
#'       #Load the 10X data
#'       tod = Read10X(file.path(dataDir))#,'raw_feature_bc_matrix/'))
#'       #tod = Read10X(file.path(dataDir,'raw_feature_bc_matrix',ref))
#'       
#'       ### CHANGED ###
#'       if (use_filtered_feature_bc_matrix) {
#'         #Get the barcodes that cell ranger considers to contain cells
#'         cells = read.delim(file.path(dataDir,'filtered_feature_bc_matrix',ref,'barcodes.tsv'),sep='\t',header=FALSE)
#'         cells = gsub('-1','',cells[,1])
#'         
#'         #Get the index in the big table
#'         cellIdx = match(cells,colnames(tod))
#'         
#'       } else {
#'         
#'         tod %>% colSums -> nCount_RNA_sums
#'         nCount_RNA_sums %>% rank -> rank_tmp
#'         which(rank_tmp >= (length(rank_tmp)-n_cells)) -> cellIdx
#'         cells <- colnames(tod)[cellIdx]
#'         nCount_RNA_sums[cellIdx] %>% sort(decreasing=T) -> nCount_RNA_sums_sort
#'         
#'         # Plot nCount_RNA histogram # TODO test this
#'         ggplot(data.frame(nCount_RNA=nCount_RNA_sums_sort, barcode=1:length(nCount_RNA_sums_sort)), aes(barcode,nCount_RNA)) + geom_line() + 
#'           scale_x_continuous(trans='log10', limits=c(1,n_cells+n_cells%/%5), breaks = c(sapply(1:4, function(x) 10^x))) + 
#'           scale_y_continuous(trans='log10', limits=c(1,nCount_RNA_sums_sort[1]), breaks = c(sapply(1:4, function(x) 10^x))) + 
#'           geom_hline(yintercept=nCount_RNA_min) + geom_hline(yintercept=nCount_RNA_max) + 
#'           ggtitle(paste0(sample_ID, ": UMI counts vs. barcodes, top ", as.character(n_cells), " barcodes"))
#'         ggsave(filename=paste0(dir_plots, prefix_data,"_", prefix_run,"_", sample_ID,"_nCount_RNA_vs_barcode.pdf"), w=12, h=8)
#'         
#'       }
#'       #############
#'       channels[[channelNames[i]]] = SoupChannel(tod,
#'                                                 tod[,cellIdx,drop=FALSE],
#'                                                 channelName=channelNames[i],
#'                                                 #ref=ref,
#'                                                 path=dataDir,
#'                                                 dataType='10X')
#'     }
#'     channels = SoupChannelList(channels)
#'     return(channels)
#'   }
#' }

######################################################################
############################ OPTIONS #################################
######################################################################

options(stringsAsFactors = F)

######################################################################
############################ CONSTANTS ###############################
######################################################################

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

dir_current <- paste0(LocationOfThisScript(), "/")

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

randomSeed <- 12345
set.seed(randomSeed)

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

# Get dirs_sample and sample_IDs by searching dirs_project_10x
if (!is.null(dirs_project_10x)) {
  
  #ref_transcript_cell <- if (flag_organism == "mmusculus") "mm10" else if (flag_organism=="hsapiens") "hg19"
  #ref_transcript_nuclei <- if (flag_organism == "mmusculus") "mm10-1\\.2\\.0_premrna" else if (flag_organism=="hsapiens") "GRCh38" 

  # Do initial rough search
  fun <- function(dir_project_10x) {
    dir(path = dir_project_10x, full.names = T, recursive = T, include.dirs = T)
  }
  args <- list("X" = dirs_project_10x)
  dirs_all <- safeParallel(fun=fun, args=args)
  dirs_all <- unlist(dirs_all, use.names=F)
  
  # For some reason, giving the pattern argument to dir() returns nothing
  pattern0 = if (!use_filtered_feature_bc_matrix | run_SoupX) { 
    paste0(".*outs/raw_feature_bc_matrix$")
    #paste0(".*outs/filtered_feature_bc_matrix/", ref_transcript_cell, "$|.*outs/filtered_feature_bc_matrix, "/", ref_transcript_nuclei, "$")
  } else { 
    paste0(".*outs/filtered_feature_bc_matrix$")#/", ref_transcript_cell, "$|.*outs/filtered_feature_bc_matrix$")#, ref_transcript_nuclei, "$")
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
  sample_IDs <- gsub(paste0("/scratch|/nfsdata|/data|/sc-10x|/data-runs|-\\d{4}_cells/|outs|/raw_feature_bc_matrix|/filtered_feature_bc_matrix|raw_feature_bc_matrix|filtered_feature_bc_matrix/"), "", dirs_sample)
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
  list_data_tmp <- safeParallel(arg=list("dir_sample"=dirs_sample), fun=Read10X, outfile=outfile)

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

  args = list("dir_sample"=dirs_sample, "sample_ID"= sample_IDs)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_loadData.txt")
  list_data_tmp <- safeParallel(fun=fun, args=args, outfile=outfile)
}

# Get filtered cell idx - unless we loaded cellranger's filtered matrix and are not using SoupX

if (run_SoupX & use_filtered_feature_bc_matrix)  {
  fun = function(dir_sample, data_tmp) {
    #Get the barcodes that cell ranger considers to contain cells
    path_barcodes <- gsub("raw_feature_bc_matrix/$", "filtered_feature_bc_matrix/barcodes.tsv.gz", dir_sample)
    cells = read.delim(path_barcodes,sep='\t',header=FALSE)
    cells = gsub('-1','',cells[,1])
    
    #Get the index in the big table
    mcellIdx <- match(cells,colnames(data_tmp))
  }
  
  args = list("dir_sample"=dirs_sample, "data_tmp"=list_data_tmp)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_get_filtered_cell_idx.txt")
  
  list_cellIdx <- safeParallel(fun=fun, args=args, outfile=outfile)
  
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
    ggsave(filename=paste0(dir_plots, prefix_data,"_", prefix_run,"_", sample_ID,"_nCount_RNA_vs_barcode.pdf"), w=12, h=8)
    return(cellIdx)
  } 
  
  args = list("n_cells"=n_cells_recovered, "data_tmp"=list_data_tmp, "sample_ID"=sample_IDs)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_get_top_cell_idx.txt")
  
  list_cellIdx <- safeParallel(fun=fun, args=args, outfile=outfile)

} else {
  list_cellIdx <- NULL
}

# SoupX ambient RNA filtering
if (run_SoupX) {
  fun= function(data_tmp, cellIdx, channelName) {
    SoupChannel(tod=data_tmp,
                toc=tod[,cellIdx,drop=FALSE],
                channelName=channelName,
                #ref=ref,
                #path=dataDir,
                dataType='10X')
    }
  args = list("data_tmp" = list_data_tmp, "cellIdx"=list_cellIdx, "channelName" = sample_IDs)
  outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_SoupChannel.txt")
  list_channels <- safeParallel(fun=fun, args=args, outfile=outfile)
  
  scl = SoupChannelList(list_channels)

  if (!is.null(SoupX_genes)) {
    list_topgenes <- lapply(1:length(sample_IDs), function(i) SoupX_genes)
  }  else {
    # since no suitably specifically expressed genes provided by user
    # find bimodally distributed genes 
    scl <- inferNonExpressedGenes(scl)
    
    fun <- function(sample_ID) {
      if (any(scl$channels[[sample_ID]]$nonExpressedGenes$isUseful)) {
        nonExpressedGenes_useful <- scl$channels[[sample_ID]]$nonExpressedGenes[scl$channels[[sample_ID]]$nonExpressedGenes$isUseful,]
        nonExpressedGenes_useful <- nonExpressedGenes_useful[order(nonExpressedGenes_useful$extremity, decreasing = T),]
        
        if (!sum(nonExpressedGenes_useful$extremity>0.7)>1) return(NA) 
        
        topgenes <- rownames(nonExpressedGenes_useful)[nonExpressedGenes_useful$extremity>0.7]
        topgenes <- topgenes[1:min(10, length(topgenes))]
        toc_log <- log(as.matrix(scl$toc))
        list_fits <- lapply(topgenes, function(genes) SeparateKRpkmNew(x=toc_log[genes,, drop=F], n=100,q=0,k=2,err = 1e-10))
        idx_topgenes_ok <- sapply(list_fits, function(fit) fit[2,1]>0.4)
        topgenes <- topgenes[idx_topgenes_ok]
        
        if (length(topgenes)<=1) return(NA)
        
        return(topgenes) 
  
      } else {
        return(NA)
      }
    }
    args= list("sample_ID"=sample_IDs)
    outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_SoupX_findTopGenes.txt")
    list_topgenes <- safeParallel(fun=fun, args=args, outfile=outfile)
  }
  
  # If we have some suitable genes to estimate ambient RNA contamination, filter the expression matrix
  if (!all(sapply(list_topgenes, function(vec_topgenes) all(is.na(vec_topgenes)), simplify = T))) {
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
        ggsave(filename = paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_SoupX_plotMarkerDistribution.pdf"), w=15, h=8)
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
        saveRDS(object = topgenes, 
                file = paste0(dir_RObjects, prefix_data, "_", prefix_run, "_", sample_ID, "_SoupX_topgenes.RDS"))
        return(data_tmp)
      }, error= function(err) {
        data_tmp
      })
    }
      
    outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_SoupXAdjustcounts.txt")
    args=list("data_tmp"=list_data_tmp, "sample_ID" = sample_IDs, topgenes = list_topgenes)
    list_data_tmp <- safeParallel(fun=fun, args=args, outfile=outfile, scl=scl) 
  }
} else if (!use_filtered_feature_bc_matrix){
  fun = function(data_tmp, cellIdx) {
    data_tmp[,cellIdx, drop=F]
  }
  args=list(data_tmp=list_data_tmp, cellIdx = list_cellIdx)
  outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_filterDataTmp.txt")
  list_data_tmp <- safeParallel(fun=fun,args=args, outfile = outfile)
}

# Now compute percent.mito, percent ribo
fun = function(data_tmp, sample_ID) {
  
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
  
  ggsave(filename =  paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_nCount_RNA_nFeature_RNA.pdf"), w=12, h=8)
  ## percent.mito, percent.ribo
  ggplot(metadata, aes(percent.mito, percent.ribo)) + 
    geom_point(shape=1) + 
    geom_vline(xintercept=percent.mito_max) +
    geom_hline(yintercept = percent.ribo_max) + 
    #scale_y_continuous(breaks = seq(from=0.1, to=0.9, by=0.1)) +
    #scale_x_continuous(breaks = seq(from=0.1, to=0.9, by=0.1)) + 
    ggtitle(paste0(sample_ID, ": prop. ribo vs. prop. mito, top ", as.character(n_cells), " barcodes"))
  ggsave(filename =  paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_percent.mito_percent.ribo.pdf"), w=12, h=8)
  
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

args=list("sample_ID"=sample_IDs, data_tmp = list_data_tmp)
outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_QCandMakeSeuratObj.txt")

list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile)

# outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_load_data_do_QC.txt")
# fun <- function(sample_ID, dir_sample, n_cells) {
#   
#   idents <- NULL
#   
#   if (!is.null(dirs_project_10x)) { 
#   ########################### SOUPX ####################################
#     if (run_SoupX) {
#       
#       #dir_sample_SoupX <- gsub("/raw_feature_bc_matrix/|/filtered_feature_bc_matrix/", "", dir_sample)
#       #dir_sample_SoupX <- gsub(paste0("/raw_feature_bc_matrix|/filtered_feature_bc_matrix|",ref_transcript_cell ,"/|", ref_transcript_nuclei, "/"), "", dir_sample)
#       scl_success <-T 
#       tryCatch({
#         scl <- load10X(dataDirs =  dir_sample,#dir_sample_SoupX, 
#                      channelNames = sample_ID,
#                      use_filtered_feature_bc_matrix = use_filtered_feature_bc_matrix,
#                      n_cells= n_cells, 
#                      nCount_RNA_min=nCount_RNA_min, 
#                      nCount_RNA_max=nCount_RNA_max, 
#                      sample_ID=sample_ID, 
#                      dir_plots=dir_plots, 
#                      prefix_data=prefix_data, 
#                      prefix_run=prefix_run)
#         }, error=function(err){
#                        scl_success <-F
#                      })
#     
#       if (scl_success) {
#         if (is.null(SoupX_genes)) {
#         # infer soup specific (bimodally distributed) genes since none provided by user
#         
#         scl <- inferNonExpressedGenes(scl)
#         
#         if (any(scl$channels[[sample_ID]]$nonExpressedGenes$isUseful)) {
#           nonExpressedGenes_useful <- scl$channels[[sample_ID]]$nonExpressedGenes[scl$channels[[sample_ID]]$nonExpressedGenes$isUseful,]
#           nonExpressedGenes_useful <- nonExpressedGenes_useful[order(nonExpressedGenes_useful$extremity, decreasing = T),]
#           greenlight = sum(nonExpressedGenes_useful$extremity>0.7)>1
#   
#         } else {
#           greenlight=F
#         }
#         
#         if (greenlight) {
#           topgenes <- rownames(nonExpressedGenes_useful)[nonExpressedGenes_useful$extremity>0.7]
#           topgenes <- topgenes[1:min(10, length(topgenes))]
#           
#           # Fit a mix of two truncated Normal distribution using LTMGSCA package; filter genes by upper mode fit >0.35
#           toc_log <- log(as.matrix(scl$toc))
#           list_fits <- lapply(topgenes, function(genes) SeparateKRpkmNew(x=toc_log[genes,, drop=F], n=100,q=0,k=2,err = 1e-10))
#           idx_topgenes_ok <- sapply(list_fits, function(fit) fit[2,1]>0.4)
#           topgenes <- topgenes[idx_topgenes_ok]
#           
#           if (length(topgenes)<=1) {
#             greenlight=F
#           } 
#         }
#       } else {
#         topgenes <- SoupX_genes 
#         greenlight=T
#       }
#     } else {
#       greenlight=F
#     }
#       
#       if (greenlight) {
#         # This plot shows:
#         # 1. for each selected gene sample 100 cells
#         # 2. for each cell, compute the expected normalized expression based on soup
#         # 3. for each cell, compute the observed normalized expression
#         # 4. for each cell, compute the log ratio of observed to expected
#         # 5. where the ratio is high, call the gene as actually expressed (not just by soup)
#         #       Ideally the 'actually expressed' should look like a mode
#         # 6. 
#         
#         tryCatch({
#           plotMarkerDistribution(scl, sample_ID, topgenes)
#           ggsave(filename = paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_SoupX_plotMarkerDistribution.pdf"), w=15, h=8)
#     
#           #estimate the contamination fraction
#           scl = calculateContaminationFraction(scl, sample_ID, list(genes = topgenes), excludeMethod = 'pCut', exCut=0.5)
#           #plotChannelContamination(scl, "Channel1")
#           #ggsave(filename = paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_SoupX_plotChannelContamination.pdf"), w=15, h=8)
#           
#           # Cell level contamination fraction
#           scl = interpolateCellContamination(scl, sample_ID, interpolationMethod = c("lowess"))
#           scl = adjustCounts(scl) # adds adjusted (raw) counts matrix in scl$atoc. 
#           # Alternative is strainCells, which modifies expression fraction. Columns (cells all sum to 1) 
#           
#           # What genes were set to zero min most cells? (output to log file)
#           cntSoggy = rowSums(scl$toc > 0)
#           cntAdjusted = rowSums(scl$atoc > 0)
#           mostZeroed <- head(sort(x=(cntSoggy - cntAdjusted)/cntSoggy, decreasing=T), n= 20)
#             
#           # Continue with adjusted counts matrix
#           data_tmp <- scl$atoc
#           
#           # Save topgenes temporarily for plotting later
#           saveRDS(object = topgenes, 
#                   file = paste0(dir_scratch, prefix_data, "_", prefix_run, "_", sample_ID, "_SoupX_topgenes.RDS"))
#         }, error=function(err) {
#                   message(paste0(sample_ID, ": plotMarkerDistribution failed with error: ", err))
#                   data_tmp <- Seurat::Read10X(data.dir = dir_sample)
#                 })
#       } else {
#         data_tmp <- Seurat::Read10X(data.dir = dir_sample)
#       }
#     } else {
#       data_tmp <- Seurat::Read10X(data.dir = dir_sample)
#       
#       if (!use_filtered_feature_bc_matrix) {
#         # Take top n cells
#         data_tmp %>% colSums -> nCount_RNA_sums
#         nCount_RNA_sums %>% rank -> rank_tmp
#         which(rank_tmp >= (length(rank_tmp)-n_cells)) -> idx_cells
#         #which(rank_tmp <= (length(rank_tmp)-5000) & rank_tmp > (length(rank_tmp)-10000)) -> idx_cells
#         
#         nCount_RNA_sums[idx_cells] %>% sort(decreasing=T) -> nCount_RNA_sums_sort
#         
#         # Plot nCount_RNA histogram # TODO test this
#         ggplot(data.frame(nCount_RNA=nCount_RNA_sums_sort, barcode=1:length(nCount_RNA_sums_sort)), aes(barcode,nCount_RNA)) + geom_line() + 
#           scale_x_continuous(trans='log10', limits=c(1,n_cells+n_cells%/%5), breaks = c(sapply(1:4, function(x) 10^x))) + 
#           scale_y_continuous(trans='log10', limits=c(1,nCount_RNA_sums_sort[1]), breaks = c(sapply(1:4, function(x) 10^x))) + 
#           geom_hline(yintercept=nCount_RNA_min) + geom_hline(yintercept=nCount_RNA_max) + 
#           ggtitle(paste0(sample_ID, ": UMI counts vs. barcodes, top ", as.character(n_cells), " barcodes"))
#         ggsave(filename=paste0(dir_plots, prefix_data,"_", prefix_run,"_", sample_ID,"_nCount_RNA_vs_barcode.pdf"), w=12, h=8)
#         
#         # filter the matrix to n cells loaded
#         data_tmp <- data_tmp[, idx_cells] 
#       }
#       
#     }
#     
#     colnames(data_tmp) <- make.unique(paste0(sample_ID, "_", colnames(data_tmp)))
#     
#   } else {
#    
#     data_tmp <- load_obj(f=dir_sample)
#     # If seurat object take raw data
#     seurat_meta <- NULL
#     if (class(data_tmp)=="seurat") {
#       #seurat_metadata <- if(is.null(dim(data_tmp@meta.data))) data_tmp@meta.data else NULL  # TODO
#       idents <- Idents(data_tmp)
#       data_tmp <- GetAssayData(data_tmp, slot="counts")
#       colnames(data_tmp) <- make.unique(paste0(sample_ID, "_", colnames(data_tmp)))
#   
#     } else {
#       if (all(c(1,2,3) %in% rownames(data_tmp))) {
#         if (!is.null(data_tmp[["X"]]))  {
#           rownames(data_tmp) <- data_tmp[["X"]] 
#           data_tmp[["X"]] <- NULL
#         } else {
#             stop("cannot identify expression data gene names. Try to add as a column 'X'")
#           }
#       }
#     }
#     # Do not filter down to n cells loaded
#   }
#   
#   # Compute and add percent mito and percent ribo as metadata
#   mito.genes <- grepl(pattern = "^mt-", x = rownames(data_tmp), ignore.case=T)
#   ribo.genes <- grepl(pattern = "^Rp[sl][[:digit:]]", x = rownames(data_tmp), ignore.case=T)
#   colSums_tmp <- colSums(x = data_tmp)
#   
#   metadata <- data.frame(percent.mito=colSums(x = data_tmp[mito.genes,])/colSums_tmp, 
#                          percent.ribo = colSums(x = data_tmp[ribo.genes,])/colSums_tmp, 
#                          row.names=colnames(data_tmp))
#   # Add nCount_RNA
#   metadata[["nCount_RNA"]] <- colSums(data_tmp)
#   
#   # add idents
#   if (!is.null(idents)) metadata[["idents"]] <- idents
#     
#   # Add nFeature_RNA
#   metadata[["nFeature_RNA"]] <-  colSums(data_tmp>0)
#   metadata[["sample_ID"]] <-  rep(sample_ID, nrow(metadata))
#   
#   # Add sample IDs
#   rownames(metadata) <- colnames(data_tmp)
#   
#   # Plot QC metrics
#   ## nCount_RNA_nFeature_RNA
#   ggplot(metadata, aes(nCount_RNA, nFeature_RNA)) + 
#     geom_point(shape=1) + 
#     scale_x_continuous(trans='log10', limits=c(1,max(metadata[["nCount_RNA"]])), breaks = c(sapply(1:4, function(x) 10^x))) +
#     scale_y_continuous(trans='log2', limits=c(1, max(metadata[["nFeature_RNA"]])), breaks = c(sapply(1:4, function(x) 10^x))) + 
#     geom_vline(xintercept=nCount_RNA_min) + 
#     geom_vline(xintercept=nCount_RNA_max) +
#     geom_hline(yintercept=nFeature_RNA_min) + 
#     geom_hline(yintercept=nFeature_RNA_max) +
#     ggtitle(paste0(sample_ID, ": nFeature_RNA vs nCount_RNA counts vs. barcodes, top ", as.character(n_cells), " barcodes"))
#   
#   ggsave(filename =  paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_nCount_RNA_nFeature_RNA.pdf"), w=12, h=8)
#   ## percent.mito, percent.ribo
#   ggplot(metadata, aes(percent.mito, percent.ribo)) + 
#     geom_point(shape=1) + 
#     geom_vline(xintercept=percent.mito_max) +
#     geom_hline(yintercept = percent.ribo_max) + 
#     #scale_y_continuous(breaks = seq(from=0.1, to=0.9, by=0.1)) +
#     #scale_x_continuous(breaks = seq(from=0.1, to=0.9, by=0.1)) + 
#     ggtitle(paste0(sample_ID, ": prop. ribo vs. prop. mito, top ", as.character(n_cells), " barcodes"))
#   ggsave(filename =  paste0(dir_plots,prefix_data,"_",prefix_run,"_", sample_ID, "_percent.mito_percent.ribo.pdf"), w=12, h=8)
#   
#   # Filter data matrix on QC metrics
#   idx_QC_ok <- metadata[["nCount_RNA"]] >= nCount_RNA_min & 
#     metadata[["nCount_RNA"]] <= nCount_RNA_max & 
#     metadata[["nFeature_RNA"]] >= nFeature_RNA_min & 
#     metadata[["nFeature_RNA"]] <= nFeature_RNA_max & 
#     metadata[["percent.mito"]] <= percent.mito_max & 
#     metadata[["percent.ribo"]] <= percent.ribo_max
#   
#   if (sum(idx_QC_ok)>=50) {
#     data_tmp <- data_tmp[,idx_QC_ok]
#     metadata <- metadata[idx_QC_ok,]
#     rownames(metadata) <- colnames(data_tmp)
#     
#     seurat_obj <- CreateSeuratObject(counts = data_tmp, 
#                                      project = paste0(prefix_data, "_", prefix_run), 
#                                      #min.cells = -Inf, 
#                                      min.features = -Inf, 
#                                      #normalization.method = "LogNormalize", 
#                                      #scale.factor = 1e4, 
#                                      #do.scale=F, 
#                                      #do.center=F, 
#                                      meta.data = metadata)
#     
#     # We scale and center here as we do not wish to regress out technical confounders before doubletFinder - see preprint p. 9: https://www.biorxiv.org/content/biorxiv/early/2018/06/20/352484.full.pdf 
#     return(seurat_obj)
#   } else {
#     warning (paste0(sample_ID, " had fewer than 50 cells after QC and was discarded"))
#     return(NULL)
#   }
# }
# 
# message("Loading data")
# 
# args <- list("dir_sample" = dirs_sample, 
#              "sample_ID" = sample_IDs, 
#              "n_cells" = n_cells_recovered)
# 
# list_seurat_obj <- safeParallel(fun=fun, 
#                                 args=args,
#                                 simplify=F,
#                                 outfile=outfile, 
#                                 dir_plots=dir_plots, 
#                                 prefix_data=prefix_data,
#                                 prefix_run=prefix_run,
#                                 nCount_RNA_min=nCount_RNA_min,
#                                 nCount_RNA_max=nCount_RNA_max, 
#                                 nFeature_RNA_min=nFeature_RNA_min, 
#                                 nFeature_RNA_max=nFeature_RNA_max, 
#                                 percent.ribo_max=percent.ribo_max, 
#                                 percent.mito_max=percent.mito_max,
#                                 use_filtered_feature_bc_matrix=use_filtered_feature_bc_matrix, 
#                                 dir_scratch=dir_scratch,
#                                 n_cells = n_cells_recovered[1]) # needed otherwise not found

names(list_seurat_obj) <- sample_IDs

idx_sample_ok <- !sapply(list_seurat_obj, is.null)

samples_dropped <- paste0(sample_IDs[!idx_sample_ok], collapse=" ")

if (nchar(samples_dropped)>0) {
  samples_dropped_df <- data.frame(samples = samples_dropped, data = prefix_data, run = prefix_run, percent.mito_max=percent.mito_max, percent.ribo_max=percent.ribo_max)
  write.csv(samples_dropped_df, file=paste0(dir_log, prefix_data, "_", prefix_run, "_samples_dropped.csv"), quote = F, row.names = F)
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
  if (length(list_metadata)==length(list_seurat_obj)) {
    message("Adding metadata")
    outfile = paste0(dir_log, prefix_data,"_",prefix_run,"_load_metadata.txt")
    fun = load_obj 
    args = list("X"=paths_metadata)
    list_metadata <- safeParallel(fun=fun, args=args, outfile=outfile)
  } else {
    warning("Number of metadata files must equal number of expression data matrices under dirs_project_10x / paths_data")
  }
}

######################################################################
############################ NORMALISE ###############################
######################################################################

outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_NormalizeData.txt")
fun = function(seurat_obj) NormalizeData(object = seurat_obj, normalization.method="LogNormalize",verbose=T)
args = list("X"=list_seurat_obj)
list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile)

######################################################################
########## FIND HIGHLY VAR GENES ###########################
######################################################################

outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_log_FindVariableFeatures.txt")
fun = function(seurat_obj) tryCatch({FindVariableFeatures(object = seurat_obj,
                                                nfeatures = nAnchorFeatures,
                                                selection.method="vst")}, 
                                    error= function(err) {
                                      message("FindVariableFeatures failed with selection.method='vst', trying with 'mean.var.plot'")
                                      FindVariableFeatures(object = seurat_obj,
                                                           nfeatures = nAnchorFeatures,
                                                           selection.method="mean.var.plot")
                                    })
args = list("X"=list_seurat_obj)
list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile)

######################################################################
########### SCALE, CENTER, REGRESS OUT CONFOUNDERS ###################
######################################################################

#if (is.null(RAM_Gb_max)) RAM_Gb_max=200
#additional_Gb = max(as.numeric(sapply(list_seurat_obj, FUN = function(x) object.size(x), simplify = T)))/1024^3
#obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
#n_cores <- max(1, min(detectCores()%/%3, RAM_Gb_max %/% (obj_size_Gb + additional_Gb))-1)
#n_cores <- min(n_cores, length(list_seurat_obj))

outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_ScaleData.txt")
fun = function(seurat_obj) ScaleData(object = seurat_obj, 
                                     vars.to.regress = vars.to.regress, 
                                     features= VariableFeatures(seurat_obj),
                                     block.size=15000,
                                     min.cells.to.block=5000,
                                     verbose=T)
args = list("X"=list_seurat_obj)
#list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile)
list_seurat_obj <- lapply(FUN = fun, X =args[[1]])#
######################################################################
################################ DO PCA ##############################
######################################################################

message("Computing PCA")
outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_PCA.txt")
args = list("seurat_obj"=list_seurat_obj, "obj_name"=names(list_seurat_obj))
fun = function(seurat_obj, obj_name) {
  tryCatch({
  pcs.compute1 = min(n_PC, min(ncol(GetAssayData(seurat_obj, slot="scale.data")), length(VariableFeatures(object=seurat_obj)))%/%2)
  RunPCA(object = seurat_obj, 
         weight.by.var = F,
         pcs.compute = pcs.compute1, 
         do.print = F,
         verbose=T)
                },   error = function(err1) {
                  pcs.compute2=min(n_PC, min(ncol(GetAssayData(seurat_obj, slot = "scale.data")), length(VariableFeatures(object=seurat_obj)))%/%3)
                  warning(paste0("RunPCA failed using ", pcs.compute1, " components with error ", err1, ", maybe due to cell count: ", dim(GetAssayData(object=seurat_obj))[2], ". Trying with ", pcs.compute2, " components"))
                  tryCatch({RunPCA(object = seurat_obj, 
                                   weight.by.var = F,
                                   pcs.compute = pcs.compute2, 
                                   do.print = F, 
                                   verbose=T)}, error = function(err2) {
                                     warning(paste0(obj_name, ": PCA failed again with ", pcs.compute2, " components with error: ", err2, ". Maybe due to cell count: ", dim(GetAssayData(object= seurat_obj))[2], ". Returning seurat object without PCA"))
                                     seurat_obj})})
}
list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile, n_PC =n_PC)

######################################################################
############################## DO TSNE ###############################
######################################################################

message("Computing t-SNE")
outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_TSNE.txt")
args = list("seurat_obj"=list_seurat_obj, "obj_name" = names(list_seurat_obj))
fun = function(seurat_obj, obj_name) {
   
    perplexity1 = min(30, min(ncol(Embeddings(object=seurat_obj, reduction="pca")),ncol(GetAssayData(object=seurat_obj)))%/%1.5)
    tryCatch({
      RunTSNE(object = seurat_obj, 
                      reduction = "pca", 
                      dims = 1:min(30,ncol(Embeddings(object=seurat_obj, reduction="pca"))), # no need to use all PCs for t-SNE
                      seed.use = randomSeed,
                      check_duplicates=F,
                      #do.fast=T,
                      perplexity=perplexity1)
      }, error = function(err1) {
                        perplexity2 = perplexity1 = min(30, min(ncol(Embeddings(object=seurat_obj, reduction="pca")),ncol(GetAssayData(object=seurat_obj)))%/%2.5)
                        warning(paste0(obj_name, ": RunTSNE failed with error: ", err1, ", using perplexity ", perplexity1, ". Maybe due to cell count: ", dim(GetAssayData(object=seurat_obj))[2], ". Trying with perplexity ", perplexity2))
                        tryCatch({RunTSNE(object = seurat_obj, 
                                          reduction = "pca", 
                                          dims = 1:min(30,ncol(Embeddings(object=seurat_obj, reduction="pca"))), # no need to use all PCs for t-SNE
                                          #do.fast=T,
                                          seed.use = randomSeed,
                                          check_duplicates=F,
                                          perplexity = perplexity2)
                        }, error= function(err2) {
                          warning(paste0("RunTSNE failed again with error ", err2, ", returning original seurat object"))
                          seurat_obj})
                      })
  }
list_seurat_obj = safeParallel(fun=fun, args=args, outfile=outfile)
  
######################################################################
######################### FIND DOUBLETS ##############################
######################################################################

if (rm_sc_multiplets) {
  
  message("Detecting doublets")
  
  # doublet rate data source: 10x https://pdfs.semanticscholar.org/presentation/0c5c/ad2edbb8f0b78cdba1c78b4324213ac20ab3.pdf
  
  #train_rate <- c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6)
  #train_cells_loaded <- c(870, 1700, 3500, 5300, 7000, 8700, 10500, 12200, 14000, 15700, 17400)
  # train_rate <- c(0,0.6,0.8,5.1,7.6,10.5)
  # train_cells_loaded <- c(150, 1530,1700, 10000, 17000, 19000)
  
  #doub_rate_lm <- lm(train_rate~train_cells_loaded)
  #n_expected_doub <- round((doub_rate_lm$coefficients[1]+doub_rate_lm$coefficients[2]*n_cells_loaded)/100*n_cells_loaded, 0)
  
  # Get doublet predictions as metadata
  
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_doubletFinder.txt")
  fun = function(seurat_obj, obj_name) {
      tryCatch({
        proportion.NN = max(0.02, 10/ncol(seurat_obj))
        doubletFinder(seu = seurat_obj, 
                      #expected.doublets = n_expected_doub, 
                      proportion.artificial = 0.25, 
                      proportion.NN=proportion.NN)
        }, 
        error = function(err) {
        warning(paste0(obj_name, ": doubletFinder failed with error: ", err, ".  Maybe due to cell count: ", dim(GetAssayData(object=seurat_obj))[2], ". Returning seurat object without predicted doublets"))
        seurat_obj
        })
  }
  args = list("seurat_obj"=list_seurat_obj, "obj_name" = names(list_seurat_obj))
  
  list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile, doubletFinder=doubletFinder)#,n_expected_doub=n_expected_doub)
    
  ##### Plot the doublets #####
  
  args = list("seurat_obj"=list_seurat_obj, "obj_name" = names(list_seurat_obj))
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_doubletFinderPlots.txt")
  fun= function(seurat_obj, obj_name) {
      try({
        DimPlot(object=seurat_obj,  
        reduction="tsne",
        dims=c(1,2),
        #do.return = T, 
        #label=F,
        #pt.size = 1, 
        group.by = "pANNPredictions", 
        #no.legend=F, 
        plot.title=paste0(obj_name, " predicted doublets"))
        ggsave(filename=paste0(dir_plots, prefix_data, "_", prefix_run, "_", obj_name, "_TSNE_doublets.pdf"), w = 10, h = 10)
        })
  }
  invisible(safeParallel(fun=fun, args=args, outfile=outfile))
  #try(invisible(mapply(FUN = fun, seurat_obj=list_seurat_obj, obj_name = names(list_seurat_obj), SIMPLIFY=F)))
 
  # Save a log 
  log_prop_singlets <- sapply(X= list_seurat_obj, FUN = function(seurat_obj){
    tryCatch({
      round(sum(seurat_obj$pANNPredictions=="Singlet")/ncol(seurat_obj),2)}, error = function(err) NA)
  }, simplify = T)
  
  log_prop_singlets_df <- matrix(log_prop_singlets, nrow=1) %>% data.frame(.,stringsAsFactors = F, row.names = NULL)
  colnames(log_prop_singlets_df) <- names(list_seurat_obj)
  write.table(log_prop_singlets, file=paste0(dir_log, prefix_data, "_", prefix_run, "_prop_singlets.tab"))
  
  message("Removing suspected doublets")
  # Filter out doublets
  outfile = paste0(dir_log, prefix_data, "_", prefix_run, "_doubletFinder_subset.txt")
  fun = function(seurat_obj, obj_name) {
    tryCatch({
      subset(object = seurat_obj, subset= "pANNPredictions"=="Singlet")
      #SubsetData(object=seurat_obj, subset.name = "pANNPredictions", accept.value="Singlet", subset.raw=T)
    }, error = function(err) {
      warning(paste0(obj_name, ": filtering out doublets failed with error: ", err, ". Maybe doubletFinder failed, possibly due to an upstream PCA or t-SNE failure."))
      seurat_obj
    })
  }
  args=list(seurat_obj=list_seurat_obj, obj_name = names(list_seurat_obj))  
  list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile)
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


if (!is.null(merge_group_IDs) | !is.null(merge_specify)) {
  
  ######################################################################
  ########################## NORMALISE #################################
  ######################################################################
  
  fun = function(seurat_obj) NormalizeData(object=seurat_obj)
  args = list("X"=list_seurat_obj)
  list_seurat_obj <- lapply(FUN=fun, "X"=args[[1]])
  
  ######################################################################
  ####################### SCALE, CENTER ################################
  ######################################################################
  
  # if (is.null(RAM_Gb_max)) RAM_Gb_max=200
  # additional_Gb = max(as.numeric(sapply(list_seurat_obj, FUN = function(x) object.size(x), simplify = T)))/1024^3
  # obj_size_Gb <- as.numeric(sum(sapply(ls(envir = .GlobalEnv), function(x) object.size(x=eval(parse(text=x)))))) / 1024^3
  # n_cores <- max(1, min(detectCores()%/%3, RAM_Gb_max %/% (obj_size_Gb + additional_Gb))-1)
  # n_cores <- min(n_cores, length(list_seurat_obj))
  
  list_seurat_obj <- lapply(X=list_seurat_obj, 
                            FUN = function(seurat_obj) {
                              ScaleData(object = seurat_obj,
                                        #vars.to.regress = vars.to.regress, 
                                        block.size=15000,
                                        min.cells.to.block=5000,
                                        verbose=T) 
                                        #do.par = T, 
                                        #num.cores = n_cores)
                              })
  invisible(gc())
  
  ######################################################################
  #################### FIND HIGHLY VAR GENES ###########################
  ######################################################################
  
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_log_FindVariableFeatures.txt")
  fun = function(seurat_obj) FindVariableFeatures(object = seurat_obj, 
                                                  do.plot=F,
                                                  nfeatures = nAnchorFeatures)
  args = list("X"=list_seurat_obj)
  list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile)
  
}
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
                             dims = 1:n_PC, 
                             k.anchor=k.anchor,
                             k.filter=k.filter,
                             k.score=k.score)
      
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
                               dims = 1:n_PC, 
                               k.anchor=k.anchor,
                               k.filter=k.filter,
                               k.score=k.score)}, 
               error = function(err1) {
                  message(paste0("FindIntegrationAnchors failed again with error ", err1, ", returning NULL"))
                  return(NULL)
        })
    })
  }
  
  args = list(list_seurat_obj_tmp = group_list_seurat_obj_tmp, name = names(group_list_seurat_obj_tmp))
  
  message("Finding achors")  
  
  list_anchors <- safeParallel(fun=fun, args=args, outfile=outfile)
  
  names(list_anchors) <- names(group_list_seurat_obj_tmp)
 
  # filter out those groups that failed
  list_anchors <- list_anchors[!all(sapply(list_anchors, is.null))]
  
  message("Integrating data")
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_log_IntegrateData.txt")
  fun = function(anchors, name) {
    tryCatch({
      IntegrateData(anchorset = anchors, 
                    dims = 1:n_PC)
    }, error=function(err) {
      message(paste0(name, ": dropped because IntegrateData failed with error ", err))
      return(NULL)
    })
  }
  args = list(anchors = list_anchors, name = names(list_anchors))
  
  
  list_seurat_align <- safeParallel(fun=fun, args=args, outfile=outfile, n_PC=n_PC)
  names(list_seurat_align) <- align_group_IDs
  
  # Filter again
  list_seurat_obj <- list_seurat_align[!sapply(list_seurat_align, is.null)]
  rm(list_seurat_align, group_list_seurat_obj_tmp)
  
}

######################################################################
#################### SCALE DATA AFTER INTEGRATING ####################
######################################################################

if (!is.null(align_group_IDs) | !is.null(align_specify)) {
  
  message("scaling integrated data")
  fun = function(seurat_obj, name) {
    DefaultAssay(object = seurat_obj) <- "integrated"
    seurat_obj <- ScaleData(object = seurat_obj, 
                            block.size=15000,
                            min.cells.to.block = 10000,
                            verbose = T)
  
  }
  args=list("seurat_obj" = list_seurat_obj, 
            "name" = names(list_seurat_obj))
  outfile = paste0(dir_log, prefix_data, "_", prefix_run, "_log_scale_integrated_data.txt")
  list_seurat_obj<-safeParallel(fun=fun, args=args,outfile=outfile)
}

######################################################################
####### FIND HIGHLY VAR FEATURES (AFTER MERGING) ####################
######################################################################

if (!is.null(merge_specify) | !is.null(merge_group_IDs)) {
  
  outfile=paste0(dir_log,prefix_data,"_",prefix_run,"_FindVariableFeatures.txt")
  fun = function(seurat_obj) {
      FindVariableFeatures(object = seurat_obj, 
                           do.plot=F,
                           nfeatures = nAnchorFeatures)
    }
  args = list("X" = list_seurat_obj)
  list_seurat_obj <- safeParallel(fun=fun,args=args, outfile=outfile)

}
######################################################################
#################### PCA (AFTER MERGING OR ALIGNING) #################
######################################################################

if (!is.null(align_group_IDs) | !is.null(align_specify) !is.null(merge_group_IDs) | !is.null(merge_specify)) {
  message("Computing PCA")
  fun =  function(seurat_obj, obj_name) {
      tryCatch({
        pcs.compute1 = min(n_PC, min(ncol(seurat_obj), length(VariableFeatures(object=seurat_obj)))%/%2)
        seurat_obj<-RunPCA(object = seurat_obj, 
               weight.by.var = F,
               npcs = pcs.compute1, 
               seed.use=randomSeed,
               verbose=T)
        # project to get loadings for all genes, not just variable
        seurat_obj <- ProjectDim(object = seurat_obj)
      },   error = function(err) {
        pcs.compute2=min(n_PC, min(ncol(seurat_obj), length(VariableFeatures(object=seurat_obj)))%/%3)
        warning(paste0("RunPCA failed using ", pcs.compute1, " components with error ", err1, ", maybe due to cell count: ", dim(GetAssayData(object=seurat_obj))[2], ". Trying with ", pcs.compute2, " components"))
        tryCatch({
          seurat_obj<-RunPCA(object = seurat_obj, 
                         weight.by.var = F,
                         npcs = pcs.compute2, 
                         seed.use=randomSeed,
                         verbose=T)
          # project so we get loadings for all genes, not just variable 
          ProjectDim(object = seurat_obj)
        }, error = function(err2) {
                           warning(paste0(obj_name, ": PCA failed again with ", pcs.compute2, " components with error: ", err2, ". Maybe due to cell count: ", dim(GetAssayData(object=seurat_obj))[2], ". Returning seurat object without PCA"))
                           seurat_obj})
      })}
  args = list("seurat_obj"=list_seurat_obj, "obj_name"=names(list_seurat_obj))
  outfile=paste0(dir_log,prefix_data,"_",prefix_run,"_PCA2.txt")
  list_seurat_obj<- safeParallel(fun=fun, args=args, outfile=outfile)
}  
######################################################################
############################ JACKSTRAW ###############################
######################################################################

if (use_jackstraw) {
  message("Using Jackstraw resampling to determine significant principal components")
  
  pvalThreshold = 0.05

  list_PC_signif_idx <- mapply(function(seurat_obj, name) {
    
    PC_signif_idx <- rep(TRUE, ncol(Loadings(object=seurat_obj, reduction="pca"))) # default
      
    prop.freq <- max(0.016, round(4/length(VariableFeatures(object=seurat_obj)),3)) # to ensure we have at least 3 samples so the algorithm works well
    # see https://github.com/satijalab/seurat/issues/5
    pcs.compute = ncol(Loadings(object=seurat_obj, reduction="pca"))
    seurat_obj <- JackStraw(object = seurat_obj, 
                            reduction="pca",
                            dims = pcs.compute,
                            num.replicate = 250, 
                            #display.progress = T,
                            #do.par = T,
                            verbose=T,
                            #num.cores = n_cores,
                            prop.freq = prop.freq)

    seurat_obj <- ScoreJackStraw(object=seurat_obj, reduction="pca", dims=1:pcs.compute, do.plot=T)
    ggsave(filename = paste0(dir_plots, prefix_data, "_", prefix_run, "_", name, "_JackStrawPlot.pdf"), width=10, height=10)
    PC_signif_idx <- seurat_obj@reductions$pca@jackstraw$overall.p.values[,2]< pvalThreshold
    return(PC_signif_idx)
  }, seurat_obj=list_seurat_obj, name=names(list_seurat_obj), SIMPLIFY = F)
} else {
  list_PC_signif_idx <- lapply(1:length(list_seurat_obj), function(seurat_obj){
    rep(TRUE, n_PC)
  })
} 
######################################################################
############################ T-SNE ###################################
######################################################################

if (is.null(path_transferLabelsRef)) { # if not we'll compute later anyway
  
  message("Computing t-SNE")
  
  fun =  function(seurat_obj, obj_name, PC_signif_idx) {
    n_comp = n_PC#if (is.null(align_group_IDs) & is.null(align_specify)) n_PC else n_CC
    perplexity1 = min(30, min(n_comp,ncol(GetAssayData(object=seurat_obj)))%/%1.5)
    tryCatch({RunTSNE(object = seurat_obj, 
                      tsne.method="Rtsne",
                      reduction = "pca",#if (is.null(align_group_IDs) & is.null(align_specify)) "pca" else "cca.aligned", 
                      dims= (1:n_comp)[PC_signif_idx], # no need to use all PCs for t-SNE
                      seed.use = randomSeed,
                      #do.fast=T,
                      perplexity=perplexity1,
                      check_duplicates=F)
    }, error = function(err1) {
      perplexity2 = min(30, min(n_comp,ncol(GetAssayData(object=seurat_obj)))%/%2.5)
      warning(paste0(obj_name, ": RunTSNE failed with error: ", err1, ", using perplexity ", perplexity1, ". Maybe due to cell count: ", dim(GetAssayData(object=seurat_obj))[2], ". Trying with perplexity ", perplexity2))
      tryCatch({RunTSNE(object = seurat_obj, 
                        reduction = "pca",#if (is.null(align_group_IDs) & is.null(align_specify)) "pca" else "cca.aligned", 
                        dims = (1:n_comp)[PC_signif_idx], # no need to use all PCs for t-SNE
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
  args=list("seurat_obj"=list_seurat_obj,
            "obj_name" = names(list_seurat_obj), 
            "PC_signif_idx" = list_PC_signif_idx)
  list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile)

}
######################################################################
############################ FIND CLUSTERS ###########################
######################################################################

if (!is.null(res_primary)) {
  
  message("Finding neighbors")
  fun = function(seurat_obj) {
    FindNeighbors(object=seurat_obj, verbose=T)
  }
  args = list("X"=list_seurat_obj)
  list_seurat_obj <- safeParallel(fun=fun, args=args)
  
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
      args= list("X"=res_to_calc)
      
      list_seurat_obj_res <- safeParallel(fun=fun, args=args, outfile=outfile, seurat_obj=seurat_obj,  randomSeed=randomSeed)
                                 
      # Save all the cluster assignments as meta data in the original seurat object
      for (i in seq(1:length(res_to_calc))) {
        seurat_obj[[paste0("clust.res.", res_to_calc[i])]] <- Idents(object=list_seurat_obj_res[[i]])
      }
      # No need to have this massive list in session!
      rm(list_seurat_obj_res)
      
      Idents(object = seurat_obj) <- paste0("clust.res.", res_primary)
      
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
    args = list("seurat_obj"=list_seurat_obj, "PC_signif_idx"=list_PC_signif_idx)
    list_seurat_obj <- safeParallel(fun=fun, args=args, outfile=outfile, res_primary=res_primary, randomSeed=randomSeed)
  }
}

######################################################################
######################### TRANSFER LABELS ############################
######################################################################

if (!is.null(path_transferLabelsRef)) {
  seuratObjRef <- load_obj(path_transferLabelsRef)
  
  if (length(seuratObjRef@assays$RNA@var.features)==0) {
    seuratObjRef <- FindVariableFeatures(object=seuratObjRef,
                                         selection.method="vst",
                                         nfeatures = nAnchorFeatures,
                                         verbose=T)
  }
  
  if (any(seuratObjRef@meta.data$percent.mito, is.null)) {
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
  
  if (all(dim(seuratObjRef@assays$RNA@scale.data))==0) {
    seuratObjRef <- ScaleData(seuratObjRef,
                              #vars.to.regress=vars.to.regress,
                              block.size=15000,
                              min.cells.to.block=5000,
                              verbose=T)
  }
  
  if (is.null(seuratObjRef@reductions$pca)) { 
    seuratObjRef <- RunPCA(object= seuratObjRef,
                           npcs= n_PC,
                           weight.by.var=T,
                           verbose=F,
                           seed.use = randomSeed)
  }
  
  fun = function(seuratObjQuery) {
    anchors <- tryCatch({
      FindTransferAnchors(reference = seuratObjRef, 
                         query = seuratObjQuery,
                         reduction = "pcaproject",
                         npcs = n_PC,
                         l2.norm=TRUE,
                         dims= 1:ncol(Embeddings(object=seuratObjRef, reduction="pca")), 
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
                            npcs = n_PC,
                            l2.norm=TRUE,
                            dims= 1:ncol(Embeddings(object=seuratObjRef, reduction="pca")), 
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
                                dims= 1:ncol(Embeddings(object=seuratObjRef, reduction="pca")), 
                                k.weight = k.weight,
                                sd.weight = sd.weight,
                                verbose=T)
    
    seuratObjQuery <- AddMetaData(object=seuratObjQuery, 
                                metadata=predictions)
    
    # Subset out cells with poor predictions
    if (sum(seuratObjQuery@meta.data$prediction.score.max < minPredictionScore)>0) {
      idx <- seuratObjQuery@meta.data$prediction.score.max >= minPredictionScore
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
                             npcs= n_PC,
                             weight.by.var=T,
                             verbose=F,
                             seed.use = randomSeed)
      
      perplexity1 = min(30, min(ncol(Embeddings(object=seuratObjQuery, reduction="pca")),
                                ncol(GetAssayData(object=seuratObjQuery)))%/%1.5)

      seuratObjQuery <- RunTSNE(object = seuratObjQuery, 
                                reduction = "pca", 
                                dims = 1:min(30,ncol(Embeddings(object=seuratObjQuery, reduction="pca"))), # no need to use all PCs for t-SNE
                                seed.use = randomSeed,
                                check_duplicates=F,
                                #do.fast=T,
                                perplexity=perplexity1)
      
    }
    
    return(seuratObjQuery)
  }
  args = list("X"=list_seurat_obj)
  list_seurat_obj <- safeParallel(fun=fun, args=args)
}

######################################################################
####################### FIND CLUSTER MARKERS #########################
######################################################################

pvalThreshold=0.05

if (!is.null(res_primary)){
  
  message("Computing differentially expressed genes as cluster markers")
  outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_FindMarkers.txt")
  # for each final seurat object, for each cluster, a dataframe of markers
  list_list_markers <- lapply(list_seurat_obj, function(seurat_obj) {
      clusters = names(table(Idents(seurat_obj)))
      args = list("X"=clusters)
      fun = function(cluster) {tryCatch({
      FindMarkers(seurat_obj,  
                  #cells.1=colnames(seurat_obj)[Idents(seurat_obj)==cluster],
                  #cells.2=NULL,
                  ident.1 = cluster,
                  #ident.2 = clusters[clusters!=cluster],
                  test.use  ="MAST",
                  #random.seed=randomSeed,
                  #latent.vars = if (!is.null(merge_specify) | !is.null(merge_group_IDs)) "sample_ID" else NULL,
                  verbose = T)
        }, 
        error = function(err) {
          NA_character_
      })}
      list_markers=NULL
      list_markers <- safeParallel(fun=fun, args=args, simplify=F, outfile=outfile, seurat_obj=seurat_obj, res_primary=res_primary)
      })
  
  # rbind the outputs from each cluster into a single df per highest-level seurat_obj, adding a 'cluster' column.
  list_markers <- mapply(function(list_markers, seurat_obj) {
    list_markers <- mapply(function(df_markers, cluster) {
      if (!all(sapply(df_markers, is.na))) {
        cbind("gene" = rownames(df_markers), "cluster"=rep(cluster, nrow(df_markers)), df_markers)
      } else {
        NA_character_
      }
    },df_markers=list_markers, cluster=names(table(seurat_obj[[paste0("clust.res.", res_primary)]])), SIMPLIFY=F)
    list_markers <- list_markers[!sapply(list_markers, function(markers) all(is.na(markers)))]
    markers <- Reduce(x=list_markers, f=rbind)
    rownames(markers) <- NULL
    return(markers)
  }, list_markers = list_list_markers, seurat_obj = list_seurat_obj, SIMPLIFY = F)

  names(list_markers) <- names(list_seurat_obj)
}

######################################################################
############################# PLOTS #################################
######################################################################

########################## PLOT T-SNE ################################

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
  ggsave(p1, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_tSNEPlot_sample.pdf"), w=10, h=10)
  
  if (!is.null(res_primary)) {
    p2 <- DimPlot(object=seurat_obj, 
                  label = T, 
                  pt.size = 1, 
                  no.legend = F, 
                  group.by = paste0("clust.res.",res_primary),
                  plot.title = paste0(name, " by cluster")) # + xlab("t-SNE 1") + ylab("t-SNE 2"
    ggsave(p2, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_", name, "_tSNEPlot_clust.pdf"), w=10, h=10)
  }
  
  
  if (!is.null(seurat_obj@meta.data$predicted.id)) {
    p3 <- DimPlot(object=seurat_obj, 
                  label = T, 
                  pt.size = 1, 
                  no.legend = F, 
                  group.by = "predicted.id",
                  plot.title = paste0(name, " by transferred label")) # + xlab("t-SNE 1") + ylab("t-SNE 2"
    ggsave(p3, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_", name, "_tSNEPlot_clust_transferlabel.pdf"), w=10, h=10)
  }
  # p1 <- TSNEPlot(seurat_obj, do.return = T, pt.size = 1, group.by = "condition", no.legend=F, plot.title=paste0(name, " by condition"))
  # ggsave(p1, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_tSNEPlot_condition.pdf"), w=10, h=10)
  
}, seurat_obj = list_seurat_obj,
name = names(list_seurat_obj),
SIMPLIFY=F
))

# list_seurat_obj <- lapply(list_seurat_obj, function(seurat_obj) {
#   metadata <- ifelse(grepl("GF", seurat_obj@meta.data$sample_ID), yes="GF", no = "CR")
#   metadata <- data.frame("condition" = metadata, row.names = rownames(seurat_obj@meta.data))
#   AddMetaData(seurat_obj, metadata)
# })
############################## FEATUREPLOTS ##############################

message("Making featureplots")
if (!is.null(feats_to_plot)) {
  invisible(mapply(function(seurat_obj, name) {
    if (feats_plot_separate) {
      lapply(feats_to_plot, function(feature) {
        if (any(grepl(feature, rownames(GetAssayData(object=seurat_obj)))) | any(grepl(feature, colnames(seurat_obj@meta.data)))) {
          tryCatch({
            #pdf(file = paste0(dir_plots, prefix_data,"_", prefix_run,"_", name, "_", feature,"_featurePlot.pdf"), w=8, h=8)
            FeaturePlot(object=seurat_obj, 
                        #dims=c(1,2), 
                        #reduction=c("tsne"),
                        features = feature)#, 
                        #cols=c("grey95", "blue"), 
                        #label=T,
                        #pt.size = 1)
            ggsave(filename =  paste0(dir_plots, prefix_data,"_", prefix_run,"_", name, "_", feature,"_featurePlot.pdf"), w=8, h=8)
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
          FeaturePlot(object=seurat_obj, 
                      features = feats_to_plot[1:2], 
                      #dims=c(1,2),
                      #reduction="tsne",
                      cols=c("grey92","blue", "violetred1"), 
                      blend=T,  
                      label=T,
                      pt.size = 1)
          ggsave(filename =  paste0(dir_plots, prefix_data,"_",prefix_run,"_", name, "_", paste0(feats_to_plot, collapse = "_"), "_featurePlot.pdf"), w=8, h=8)
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
    FeaturePlot(seurat_obj, features.plot = gene, cols.use=c("grey95", "blue"), pt.size = 1, no.legend = F)
    ggsave(filename = paste0(dir_plots, prefix_data, "_", prefix_run, "_", sample_ID, "_", gene, "_SoupX_featurePlot.pdf"), width = 12, height=10)
  }
  
  for (sample_ID in sample_IDs) {
    
    topgenes_path <- paste0(dir_scratch, prefix_data, "_", prefix_run, "_", sample_ID, "_SoupX_topgenes.RDS")
    
    if (file.exists(topgenes_path)) {
      
      topgenes <- readRDS(file = topgenes_path)
      
      args = list("seurat_obj" = list_seurat_obj, 
                  "gene"=topgenes)
      
      safeParallel(fun=fun, 
                   args=args, 
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
  VlnPlot(object = seurat_obj, 
          features = "nCount_RNA", 
          #do.sort = F, 
          same.y.lims=T, 
          #legend.position = "bottom", 
          #single.legend = F, 
          #do.return=T, 
          group.by = "sample_ID")
  ggsave(filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_nCount_RNA_VlnPlot_sample.pdf"), w=10, h=10)
  VlnPlot(object = seurat_obj, 
          features = "percent.mito", 
          #do.sort = F, 
          same.y.lims=T, 
          #legend.position = "bottom", 
          #single.legend = F, 
          #do.return=T, 
          group.by = "sample_ID")
  ggsave(filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_percent.mito_VlnPlot_sample.pdf"), w=10, h=10)
  VlnPlot(object = seurat_obj, 
          features = "percent.ribo", 
          #do.sort = F, 
          same.y.lims=T, 
          #legend.position = "bottom", 
          #single.legend = F, 
          #do.return=T, 
          group.by = "sample_ID")
  ggsave(filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_percent.ribo_VlnPlot_sample.pdf"), w=10, h=10)
  
}, seurat_obj = list_seurat_obj,
name = names(list_seurat_obj),
SIMPLIFY=F
))

#### Facet plots of samples in each cluster and clusters in each sample ####

message("plotting sample proportions in each cluster and cluster proportions for each sample")
invisible(mapply(function(seurat_obj, name){
  df = data.frame("sample_ID" =as.character(seurat_obj$sample_ID), "cluster"=as.character(Idents(seurat_obj)))
  
  ggplot(df, aes(sample_ID, cluster)) + geom_bar(stat = "identity") + facet_wrap(.~cluster) 
  ggsave(filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_samples_per_cluster.pdf"), w=24, h=15)
  
  ggplot(df, aes(cluster, sample_ID)) + geom_bar(stat = "identity") + facet_wrap(.~sample_ID) 
    theme(strip.text.x = element_text(size=8, angle=75),
          strip.text.y = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="red", fill="#CCCCFF"))
  ggsave(filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_cluster_per_sample.pdf"), w=24, h=15)
}, seurat_obj = list_seurat_obj, name=names(list_seurat_obj),SIMPLIFY=F))


if (!is.null(path_transferLabelsRef)) {
  message("plotting sample proportions in each cluster and cluster proportions for each sample")
  invisible(mapply(function(seurat_obj, name){
    df = data.frame("sample_ID" =as.character(seurat_obj$sample_ID), "predicted.id"=as.character(seurat_obj@meta.data$predicted.id))
    
    ggplot(df, aes(sample_ID, predicted.id)) + geom_bar(stat = "identity") + facet_wrap(.~predicted.id) 
    ggsave(filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_samples_per_predicted.id.pdf"), w=24, h=15)
    
    ggplot(df, aes(predicted.id, sample_ID)) + geom_bar(stat = "identity") + facet_wrap(.~sample_ID) 
    theme(strip.text.x = element_text(size=8, angle=75),
          strip.text.y = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="red", fill="#CCCCFF"))
    ggsave(filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_predicted.id_per_sample.pdf"), w=24, h=15)
  }, seurat_obj = list_seurat_obj, name=names(list_seurat_obj),SIMPLIFY=F))
  
}
######################################################################
########################## OUTPUT TABLES, ROBJECTS ###################
######################################################################

if (!is.null(res_primary)) {
  ############################## CLUSTER GENE MARKERS ##################
  message("Writing out gene markers")

  # invisible(tryCatch({
  #   xlsx::write.xlsx2(x = list_markers,
  #                     file = paste0(dir_tables, prefix_data, "_", prefix_run, "_", name, "_markers.xlsx"))}, error =function(err){ 
  #   message(paste0("xlsx failed with error", err))
  #   write.csv(markers, file=paste0(dir_tables, prefix_data, "_", prefix_run, "_", name, "_markers.csv"), quote = F, row.names=T)
  # }))
  
  tryCatch({
      openxlsx::write.xlsx(x = list_markers,
                          file = paste0(dir_tables, prefix_data, "_", prefix_run, "_markers.xlsx"))
    }, error =function(err){ 
      message(paste0("openxlsx::write.xlsx failed with error", err))
      invisible(mapply(function(markers,name) write.csv(markers, file=paste0(dir_tables, prefix_data, "_", prefix_run, "_", name, "_markers.csv"), quote = F, row.names=T), markes=list_markers, name=names(list_seurat_obj), SIMPLIFY=F))

    })
  
}
  
########################## SAVE SEURAT OBJECTS #######################

outfile=paste0(dir_log, prefix_data, "_", prefix_run, "_outputTables_RObjects.txt")
args=list("seurat_obj" = list_seurat_obj, 
          "name" = names(list_seurat_obj))
fun = function(seurat_obj, name) {invisible(saveRDS(seurat_obj, file = paste0(dir_RObjects, prefix_data, "_", prefix_run, "_", name, "_seurat_obj.RDS.gz"), compress = "gzip"))}
invisible(safeParallel(fun=fun, args=args, outfile=outfile))


########################## SAVE OPTIONS ##############################

# params_run <- unlist(opt, recursive=F)# all the user options/defaults
# params_run <- c(c("data"=prefix_data, "run"=prefix_run), params_run)
# matrix(unlist(params_run), nrow=1) %>% as.data.frame(row.names=NULL, stringsAsFactors=F) -> params_run_df
# colnames(params_run_df) <- names(params_run)
# params_run_path = sprintf("%s%s_%s_params_run.tab", dir_log, prefix_data, prefix_run)
# append_params_run = F#file.exists(params_run_path) # NB: if left as NULL, they won't get included!
# write.table(params_run_df, 
#             file=params_run_path, 
#             quote = F, 
#             sep = "\t", 
#             row.names=F, 
#             append = append_params_run, 
#             col.names = !append_params_run)

######################################################################
############### LOG PARAMETERS AND FILE VERSION ######################
######################################################################

as.character(Sys.time()) %>% gsub("\\ ", "_",.) %>% gsub("\\:", ".", .) ->tStop

if (is.null(path_runLog)) path_runLog <- paste0(dirLog, "_preservation_runLog.txt")

dirCurrent = paste0(LocationOfThisScript(), "/") # need to have this function defined
setwd(dirCurrent) # this should be a git directory

# get the latest git commit
gitCommitEntry <- try(system2(command="git", args=c("log", "-n 1 --oneline"), stdout=TRUE))

# Write to text file
cat(text = "\n" , file =  path_runLog, append=T, sep = "\n")
cat(text = "##########################" , file =  path_runLog, append=T, sep = "\n")

cat(text = prefixOut , file =  path_runLog, append=T, sep = "\n")
cat(sessionInfo()[[1]]$version.string, file=path_runLog, append=T, sep="\n")
if (!"try-error" %in% class(gitCommitEntry)) cat(text = paste0("git commit: ", gitCommitEntry) , file =  path_runLog, append=T, sep = "\n")
cat(text = paste0("tStart: ", tStart) , file =  path_runLog, append=T, sep = "\n")
cat(text = paste0("tStop: ", tStop) , file =  path_runLog, append=T, sep = "\n")

# output parameters (assumping use of optparse package)
cat(text = "\nPARAMETERS: " , file =  path_runLog, append=T, sep = "\n")
for (i in 1:length(opt)) {
  cat(text = paste0(names(opt)[i], "= ", opt[[i]]) , file =  path_runLog, append=T, sep = "\n")
}

cat(text = "##########################" , file =  path_runLog, append=T, sep = "\n")

######################################################################
############################### WRAP UP ##############################
######################################################################
save.image(file = paste0(dir_out, prefix_data, "_", prefix_run,"_finalSessionImage.RData.gz"), compress="gzip")
message("Script done!")