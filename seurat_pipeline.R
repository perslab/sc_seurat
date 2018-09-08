# Seurat pipeline
# Usage: time ./seurat_pipeline.R --dir_project_10x /nfsdata//nfsdata/data/sc-10x/data-runs/180604-perslab-germfree/ --dir_out /projects/jonatan/tmp-germfree --prefix_data "perslab_scnuclei_pilot" --prefix_run "seurat_1" --merge_IDs "c('all')" --n_cores 10  --res_to_calc "c(0.6,0.8,1.0,1.2,1.5,2.0)" --feats_to_plot "c('Ghrl','Lep')"

# TODO
# in RunCCA, group.by some metadata given as argument? e.g. protocol?
# check Dylan and Mette's CCA scripts
# add doubletFinder
# add geneplots (nUMI, percent.mito, percent.ribo)
# check timshel's pipeline
# add UMAP
# in FindConservedMarkers,  do we need # grouping.var = align_ID ?
# Identify celltypes

######################################################################
########################### OptParse #################################
######################################################################

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--dir_project_10x", type="character", default = NULL,
              help = "Corresponds to cellranger's $DATARUNS_DIR/$PROJECT_ID. Contains $SAMPLE_NAME subfolders with /outs/filtered_gene_bc_matrices/mm10 subdirs. [default %default]"),  
  make_option("--paths_data", type="character", default = NULL,
              help = "Path(s) to single counts data matrix stored in standard (gz compressed) delimited text format or loom. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g. ''c('<path1>','<path2>')'',  [default %default]"),  
  make_option("--dir_out", type="character",
              help = "Project directory to which to write files. Should include subdirectories /tables, /RObjects, /plots, /log"),  
  make_option("--flag_datatype", type="character", default = "sc",
              help = "Accepts arguments sc or bulk, [default %default]"), 
  make_option("--prefix_data", type="character", default = "Seurat_out",
              help = "Dataset prefix for output files, [default %default]"), 
  make_option("--prefix_run", type="character", default=substr(gsub("-","",as.character(Sys.Date())),3,1000),
              help = "Run prefix for output files, [default %default]"), 
  make_option("--paths_metadata", type="character", default = NULL,
              help = "Path(s) to metadata file(s) in one of the standard (compressed) character separated formats. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g. ''c('<path1>','<path2>')''. If not given, takes any metadata stored within the path_data object(s). [default %default]"),  
  make_option("--merge_IDs", type="character", default = NULL,
              help = "Merge samples by some identifiers? Precedes any sample alignment. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g. ''c('GF','CR')''. The group identifier must be part of the names of the sample sub-directories. Use 'c('all')' to merge all samples, [default %default]"), 
  make_option("--align_IDs", type="character", default = NULL,
              help = "Align samples by some identifiers? Succeeds any sample merge. Takes a vector, in single (double) quotes, of characters, in double (single) quotes, without whitespace, e.g.''c('tissue.1','tissue.2')''. Use ''c('all')'' to align all samples, [default %default]"), 
  make_option("--n_cores", type="integer", default = 15L,
              help = "Script uses the parallel package using FORK cluster; how many cores? [default %default]"),
  make_option("--n_comp", type="integer", default = 75L,
              help = "How many principal and canonical components? [default %default]"),
  make_option("--res_primary", type="double", default = 0.8,
              help = "resolution parameter to pass to Seurat::FindClusters. Used for finding gene markers and in plots. [default %default]"),
  make_option("--res_to_calc", type="character", default = NULL,
              help = "multiple resolution parameters to pass to Seurat::FindClusters. Must include res_primary. Format as a quoted character with a vector of values, without whitespace. [default %default]"),
  make_option("--feats_to_plot", type="character", default = NULL,
              help = "Features to plot. Format as a quoted character with a vector of values, without whitespace. [default %default]"),
  make_option("--feats_plot_separate", type="character", default = T,
              help = "Plot features separately or in one plot? If FALSE, only the first two features will be plotted [default %default]")
)

######################################################################
########################### PACKAGES #################################
######################################################################

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(loomR))

######################################################################
########################### GET OPTIONS ##############################
######################################################################

opt <- parse_args(OptionParser(option_list=option_list))

dir_project_10x <- opt$dir_project_10x
paths_data <- opt$paths_data ; if (!is.null(paths_data)) paths_data <- eval(parse(text=paths_data))
dir_out <- opt$dir_out
flag_datatype <- opt$flag_datatype
prefix_data <- opt$prefix_data
prefix_run <- opt$prefix_run
paths_metadata <- opt$paths_metadata ; if (!is.null(paths_metadata)) paths_metadata <- eval(parse(text=paths_metadata))
merge_IDs <- opt$merge_IDs ; if (!is.null(merge_IDs)) merge_IDs <- eval(parse(text=merge_IDs))
align_IDs <- opt$align_IDs ; if (!is.null(align_IDs)) align_IDs <- eval(parse(text=align_IDs))
n_cores <- opt$n_cores
n_comp <- opt$n_comp
res_primary <- opt$res_primary
res_to_calc <- opt$res_to_calc ; if (!is.null(res_to_calc)) res_to_calc <- eval(parse(text = res_to_calc))
feats_to_plot <- opt$feats_to_plot ; if (!is.null(feats_to_plot)) feats_to_plot <- eval(parse(text = feats_to_plot))
feats_plot_separate <- opt$feats_plot_separate

######################################################################
########################### SET OPTIONS (DEV) ########################
######################################################################

if (FALSE) {
  dir_project_10x <- "/nfsdata/data/sc-10x/data-runs/180214-perslab-hyp_nuclei_optimization/"
  paths_data <- NULL 
  dir_out <- "/projects/jonatan/2018_Th_neurons/"
  flag_datatype <- "sc"
  prefix_data <- "perslab_hyp_nuclei"
  prefix_run <- "Th"
  paths_metadata <- NULL
  merge_IDs <- NULL 
  align_IDs <- "hyp_nucl"
  n_cores <- 20
  n_comp <- 75
  res_primary <- 0.8 # the primary resolution - will be used for cluster marker identification.
  res_to_calc <- c(0.4, 0.6, 0.8, 1.2, 1.6, 2, 4, 10) # Alternative resolutions to calculate for FindClusters. OBS: res_primary must be included in res2calculate
  feats_to_plot <- c("Th","Ghrl","Lep")
  feats_plot_separate <- T
}

######################################################################
########################### SET OPTIONS ##############################
######################################################################

options(stringsAsFactors = F)

######################################################################
############################ CONSTANTS ###############################
######################################################################

if (is.null(dir_out)) {
  pos <- tail(gregexpr("/", seurat_path)[[1]],2)[1]
  dir_out = substr(seurat_path, 1, pos)
}

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

LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
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

dir_current = paste0(LocationOfThisScript(), "/")

flag_date = substr(gsub("-","",as.character(Sys.Date())),3,1000)

######################################################################
########################### VERIFY INPUT #############################
######################################################################

if (!xor(is.null(dir_project_10x), is.null(paths_data))) stop("Provide dir_project_10x or path_data")

######################################################################
######################## DEFINE FUNCTIONS ############################
######################################################################

load_obj <- function(f) {
  # Utility function for loading an object stored in any of 
  # .RData, .RDS, .loom, .csv, .txt, .tab, .delim
  # File may be gzip compressed 
  
  compressed = F
  
  if (grepl(pattern = "\\.gz|\\.gzip", x=f))  {
    compressed <- T
    f = paste0("gzfile('",f,"')")
  }
  
  if (grepl(pattern = "\\.RDS|\\.rds", x = f)) {
    readRDS(file=if(compressed) eval(parse(text=f)) else f)
  } else if (grepl(pattern="\\.RData|\\.Rdata", x=f)) { 
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
  } else if (grepl(pattern="\\.loom", x=f)) {
    connect(filename=f, mode = "r+")
  } else if (grepl(pattern = "\\.csv", x=f)) {
    read.csv(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
  } else if (grepl(pattern = "\\.tab", x=f)) {
    read.table(file=if(compressed) eval(parse(text=f)) else f, sep="\t", stringsAsFactors = F, quote="", header=T) 
  } else if (grepl(pattern = "\\.txt", x=f)) {
    read.delim(file=if(compressed) eval(parse(text=f)) else f, stringsAsFactors = F, quote="", header=T)
  }
}


######################################################################
############# LOAD 10x DATA AND CREATE SEURAT OBJECT #################
######################################################################

if (!is.null(dir_project_10x)) {
  dirs_sample <- dir(path = dir_project_10x, full.names = T, recursive = T, no.. = T, include.dirs = T)
  dirs_sample <- dirs_sample[grep(pattern = "/outs/filtered_gene_bc_matrices/mm10$", x = dirs_sample)]
  dirs_sample <- paste0(dirs_sample, "/")
  
  sample_IDs <- gsub("/nfsdata/data/sc-10x/data-runs/|-5000_cells/outs/filtered_gene_bc_matrices/mm10/", "", dirs_sample)
  sample_IDs <- gsub(".*/", "", sample_IDs)
  
} else {
  dirs_sample <- paths_data 
  sample_IDs <- gsub(".*/|\\.RData|\\.RDS|\\.loom|\\.csv|\\.tab|\\.txt", "", dirs_sample)
}

cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_load_raw_data.txt"))

if (!is.null(dir_project_10x)) {
  list_raw <- parLapplyLB(cl, dirs_sample, function(dir_sample) Seurat::Read10X(data.dir = dir_sample))
} else { 
  list_raw <- parLapplyLB(cl, dirs_sample, function(dir_sample) load_obj(f=dir_sample))
}

list_seurat_obj <- parLapplyLB(cl, list_raw, function(rawData) CreateSeuratObject(raw.data = rawData, 
                                                                                  min.cells = 3, 
                                                                                  min.genes = 200, 
                                                                                  project = paste0(prefix_data, "_", prefix_run)))
stopCluster(cl)

rm(list_raw)
names(list_seurat_obj) <- sample_IDs

######################################################################
####################### ADD META DATA  ###############################
######################################################################

cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, prefix_data,"_",prefix_run,"_add_metadata.txt"))

# For each data object, add sample name as a column in metadata
list_seurat_obj <- clusterMap(cl, function(seurat_obj, sample_ID)
{ AddMetaData(object = seurat_obj, metadata = data.frame(sample_ID = rep(sample_ID, times=nrow(seurat_obj@meta.data)), row.names = rownames(seurat_obj@meta.data)))
}, 
seurat_obj = list_seurat_obj,
sample_ID = sample_IDs,
SIMPLIFY=F,
.scheduling=c("dynamic"))

if (!is.null(paths_metadata)) {
  if (length(paths_metadata) != length(list_seurat_obj)) stop("Length of paths_metadata must equal number of data matrices under dir_project_10x / paths_data") 
  list_metadata <- parLapplyLB(cl, paths_metadata, load_obj)
  list_seurat_obj <- clusterMap(cl, function(seurat_obj, metadata)
  { AddMetaData(object = seurat_obj, metadata = data.frame(metadata, row.names = rownames(seurat_obj@meta.data)))
  }, 
  seurat_obj = list_seurat_obj,
  metadata = list_metadata,
  SIMPLIFY=F,
  .scheduling=c("dynamic"))
}

stopCluster(cl)

######################################################################
############################# MERGE ##################################
######################################################################

if (!is.null(merge_IDs)) {
  if (length(merge_IDs)==1) {
    for (i in 1:length(list_seurat_obj)) {
      if (i == 1) {
        # list_seurat_merged is of length 1, just for consistency with other code
        list_seurat_merged <- list_seurat_obj[i] # single brackets intended
      } else {
        list_seurat_merged <- list(MergeSeurat(object1 = list_seurat_merged[[1]], 
                                     object2 = list_seurat_obj[[i]],
                                     project = paste0(prefix_data, "_", prefix_run),
                                     min.cells = -Inf,
                                     min.genes = -Inf,
                                     do.normalize = F,
                                     add.cell.id2 = names(list_seurat_obj)[i])) 
      }
    }
  } else {
    
    cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, "MergeSeurat.txt"))
    
    list_vec_merge_idx <- parLapplyLB(cl, merge_IDs, function(id) grep(pattern=id, 
                                                     x = names(list_seurat_obj),
                                                     ignore.case=T))
    list_list_seurat_obj_tmp <- list()
    
    list_list_seurat_obj_tmp[[merge_IDs[i]]] <- parLapplyLB(cl, list_vec_merge_idx, function(vec_merge_idx) list_seurat_obj[vec_merge_idx]) 
    
    list_seurat_merged <- clusterMap(cl, list_list_seurat_obj_tmp, function(list_seurat_obj_tmp) {
      
      for (i in 1:length(list_seurat_obj_tmp)) {
        if (i == 1) {
          seurat_merged <- list_seurat_obj_tmp[[i]]
        } else {
          seurat_merged <- MergeSeurat(object1 = seurat_merged, 
                                       object2 = list_seurat_obj_tmp[[i]],
                                       project = paste0(prefix_data, "_", prefix_run),
                                       min.cells = -Inf,
                                       min.genes = -Inf,
                                       do.normalize = F,
                                       add.cell.id2 = names(list_seurat_obj_tmp)[i])
        }
      }
      
    })
  
    
    rm(list_list_seurat_obj_tmp)
  }
  
  names(list_seurat_merged) <- merge_IDs
  list_seurat_obj <- list_seurat_merged
  rm(list_seurat_merged)
  
  ## Add merge_IDs as metadata
  list_seurat_obj <- clusterMap(cl, function(seurat_obj, merge_ID)
  { AddMetaData(object = seurat_obj, metadata = data.frame(merge_ID = rep(merge_ID, times=nrow(seurat_obj@meta.data)), row.names = rownames(seurat_obj@meta.data)))
  }, 
  seurat_obj = list_seurat_obj,
  merge_ID = merge_IDs, SIMPLIFY=F,
  .scheduling=c("dynamic"))
  stopCluster(cl)
}
######################################################################
####################### DO QUALITY CONTROL ###########################
######################################################################

if (flag_datatype == "sc") {
  
  cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, "ribo_mito_filter.txt"))
  
  list_seurat_obj <- parLapplyLB(cl, list_seurat_obj, function(seurat_obj) {
    mito.genes <- grep(pattern = "^mt-", x = rownames(x = seurat_obj@raw.data), value = TRUE, ignore.case=T)
    ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = seurat_obj@raw.data), value = TRUE)
    percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito.genes, ])/Matrix::colSums(seurat_obj@raw.data)
    percent.ribo <- Matrix::colSums(seurat_obj@raw.data[ribo.genes, ])/Matrix::colSums(seurat_obj@raw.data)
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.ribo, col.name = "percent.ribo")
    seurat_obj
  }) 
  
  clusterMap(cl, function(seurat_obj, sample) {
    
    f1 <- VlnPlot(object = seurat_obj, features.plot = c("nGene", "nUMI", "percent.mito","percent.ribo"), nCol = 4, do.return = T)
    par(mfrow = c(1, 2))
    ggsave(plot = f1, filename =  paste0(dir_plots,"_",prefix_data,"_", prefix_run,"_", sample,"_violinplot.pdf"), w=12, h=8)
    
    pdf(file =  paste0(dir_plots,flag_date,"_",prefix_data,"_",prefix_run,"_",  sample,"_geneplot_nUMI_percent.mito.pdf"), w=12, h=8)
    GenePlot(object = seurat_obj, gene1 = "nUMI", gene2 = "percent.mito")
    dev.off()
    
    pdf(file =  paste0(dir_plots,flag_date,"_",prefix_data,"_",prefix_run,"_",  sample,"_geneplot_nUMI_percent.ribo.pdf"), w=12, h=8)
    GenePlot(object = seurat_obj, gene1 = "nUMI", gene2 = "percent.ribo")
    dev.off()
    
    pdf(file=paste0(dir_plots,flag_date,"_",prefix_data,"_", prefix_run,"_", sample,"_geneplot_percent.mito_percent.ribo.pdf"), w=12, h=8)
    GenePlot(object = seurat_obj, gene1 = "percent.mito", gene2 = "percent.ribo")
    dev.off()
    
    pdf(file=paste0(dir_plots,flag_date,"_",prefix_data,"_", prefix_run,"_", sample,"_geneplot_nUMI_nGene.pdf"), w=12, h=8)
    GenePlot(object = seurat_obj, gene1 = "nUMI", gene2 = "nGene")
    dev.off()
    
  }, seurat_obj = list_seurat_obj,  sample = names(list_seurat_obj), SIMPLIFY=F)
  
  # Filter out cells with high mito / ribo content
  list_seurat_obj_cellnames <- lapply(list_seurat_obj, function(seurat_obj) as.character(seurat_obj@ident))
  
  list_seurat_obj <- clusterMap(cl, function(seurat_obj, sample) {
    
    seurat_obj_copy <- seurat_obj
    
    seurat_obj <- try(FilterCells(object = seurat_obj, 
                                  subset.names = c("nGene", "percent.mito", "percent.ribo"), 
                                  low.thresholds = c(100, -Inf, -Inf), 
                                  high.thresholds = c(5000,0.15,0.15)))
    if (class(seurat_obj)=="try-error") {
      warning(paste0(sample, ": FilterCells failed, returning unfiltered object"))
      return(seurat_obj_copy) 
    } else return(seurat_obj)
  },  seurat_obj=list_seurat_obj, sample=names(list_seurat_obj), SIMPLIFY=F)

  
  df_filtered <- clusterMap(cl, function(old, new) {
    vec_filtered <- setdiff(x = old, y = as.character(new@ident))
    vec_filtered
  },
  old = lapply(list_seurat_obj, function(seurat_obj) seurat_obj@ident), 
  new = list_seurat_obj, 
  SIMPLIFY=T, .scheduling = c("dynamic"))
                        
  write.table(x = df_filtered, file = paste0(dir_log, prefix_data, "_", prefix_run, "_mito_ribo_cells_filtered_out.csv"), quote = F, col.names = T)
  rm(df_filtered)
  
  stopCluster(cl)
  
}

######################################################################
########## NORMALIZE, SCALE AND FIND HIGHLY VAR GENES ###########################
######################################################################

cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, "Normalize_Scale_HVG.txt"))

list_seurat_obj <- parLapply(cl, list_seurat_obj, function(seurat_obj) NormalizeData(object = seurat_obj, 
                                                                                     normalization.method = "LogNormalize",  
                                                                                     scale.factor = 1e4, 
                                                                                     display.progress=T))

list_seurat_obj <- parLapply(cl, list_seurat_obj, function(seurat_obj) {
  
  FindVariableGenes(object = seurat_obj, 
                    mean.function = ExpMean, 
                    dispersion.function = LogVMR, 
                    x.low.cutoff = 0.0125, 
                    x.high.cutoff = 3, 
                    y.cutoff = 0.5,
                    do.plot=F, display.progress =F)
})

stopCluster(cl)

list_seurat_obj <- lapply(list_seurat_obj, function(seurat_obj) ScaleData(object = seurat_obj, 
                             genes.use = NULL, 
                             vars.to.regress = c("nUMI", "percent.mito", "percent.ribo"), 
                             model.use="linear", 
                             do.scale=T, 
                             do.center=T, 
                             do.par=T, 
                             num.cores=n_cores))  

######################################################################
################## FIND HVG UNION AS FEATURES FOR ALIGNMENT ##########
######################################################################
## Select genes as features for Canonical Correlation Analysis-based alignment
# we will take the union of the top 2k variable genes in each dataset for alignment 

if (!is.null(align_IDs)) {
  
  cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, "GetHighlyVariableGenes.txt"))
  
  list_hvg <- parLapply(cl, list_seurat_obj, function(seurat_obj) {
    hvg <- rownames(x= head(x=seurat_obj@hvg.info, n = 2000))
    return(hvg)
  })
  stopCluster(cl)
  
  hvg.union <- Reduce(f = dplyr::union, x = list_hvg) #nb: dplyr::union removes duplicates
  
  list_rownames <- lapply(list_seurat_obj, function(seurat_obj) rownames(seurat_obj@scale.data))
  
  for (rnames in list_rownames) {
    hvg.union <- intersect(hvg.union, rnames)
  }
}
######################################################################
######################### DO ALIGNMENT ###############################
######################################################################

if (!is.null(align_IDs)) {
  
  cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, "AlignCCA.txt"))
  
  if (length(align_IDs)==1) {
    for (i in 1:length(list_seurat_obj)) {
      list_seurat_obj <- list(parLapply(cl, list_seurat_obj, function(seurat_obj) {
        RunMultiCCA(object.list=list_seurat_obj, 
                    genes.use=hvg.union, 
                    add.cell.ids = names(list_seurat_obj), 
                    niter = 25,
                    num.ccs = n_comp, 
                    standardize = TRUE) }))
      }
    
  } else if (length(align_IDs==2)) {
    
    list_seurat_obj <- parLapply(cl, list_seurat_obj, function(seurat_obj) {
      
      RunCCA(object = list_seurat_obj[[1]],
             object2 = list_seurat_obj[[2]], 
             group1 = names(list_seurat_obj)[[1]], 
             group2 = names(list_seurat_obj)[[2]]) })
    
    
  } else {

    list_vec_align_idx <- parLapplyLB(cl, align_IDs, function(id) grep(pattern=id, 
                                                                       x = names(list_seurat_obj),
                                                                       ignore.case=T))
    
    list_list_seurat_obj_tmp[[align_IDs[i]]] <- parLapplyLB(cl, list_vec_align_idx, function(vec_align_idx) list_seurat_obj[vec_align_idx]) 
    
    list_seurat_align <- clusterMap(cl, list_list_seurat_obj_tmp, function(list_seurat_obj_tmp) {
      
    for (i in 1:length(list_seurat_obj_tmp)) {
      RunMultiCCA(object.list=list_seurat_obj, 
                  genes.use=hvg.union, 
                  add.cell.ids = names(list_seurat_obj_tmp), 
                  niter = 25,
                  num.ccs = n_comp, 
                  standardize = TRUE) 
      }
    })
  }
    
  stopCluster(cl)
  
  rm(list_list_seurat_obj_tmp)

  names(list_seurat_align) <- align_IDs
  list_seurat_obj <- list_seurat_align
  rm(list_seurat_align)
  
}

######################################################################
############################# PLOT CCAs ###############################
######################################################################
## Visualize results of CCA plot CC1 versus CC2 and look at a violin plot

if (!is.null(align_IDs)) {
  
  cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, "plotCCA.txt"))
  
  invisible(clusterMap(cl, function(seurat_obj, name) {
    par(mfrow = c(1, 2))
    p1 <- DimPlot(object = seurat_obj,
                  reduction.use = "cca",
                  group.by = "sample_ID",
                  pt.size = 0.5,
                  do.return = TRUE)
    
    p2 <- VlnPlot(object = seurat_obj,
                  features.plot = "CC1",
                  group.by = "sample_ID",
                  do.return = TRUE,
                  point.size.use = 0.1)
    plot_grid(p1, p2)
    
    ggsave(filename =  paste0(dir_plots,prefix_data,"_",prefix_run,"_", name, "_CC1_vs_CC2_and_violin.pdf"), w=12, h=8)
    
    # MetageneBicorPlot shows shared correlation strength vs CCs
    p3 <- MetageneBicorPlot(seurat_obj, 
                           grouping.var = "sample_ID", 
                           dims.eval = 1:n_comp, 
                           display.progress = F)
    ggsave(plot = p3, filename = paste0(dir_plots,prefix_data,"_",prefix_run,"_", name, "_metageneBicorPlot.pdf"), w=12, h=8)
    
    pdf(file =  paste0(dir_plots,prefix_data,"_",prefix_run,"_", name, "_CCA_DimHeatMap.pdf"), w=18, h=12)
    
    DimHeatmap(object = seurat_obj, 
               reduction.type = "cca", 
               cells.use = 500, 
               dim.use = 1:n_comp, 
               do.balanced = TRUE)
    dev.off()
    
    
  }, seurat_obj = list_seurat_obj, 
  name = names(list_seurat_obj), 
  SIMPLIFY=F, 
  .scheduling = c("dynamic")))
  
  stopCluster(cl)
}

######################################################################
###################### KEEP WHAT CAN BE ALIGNED ######################
######################################################################

if (!is.null(align_IDs)) {
  
  cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, "filter_CCA.txt"))
  
  clusterMap(cl, function(seurat_obj, name) {
    
  seurat_obj <- CalcVarExpRatio(object = seurat_obj, 
                                     reduction.type = "pca", 
                                     grouping.var = "sample_ID", # TODO: what is the correct parameter?
                                     dims.use = 1:n_comp)
  
  # seurat_discard <- SubsetData(object = seurat_obj, 
  #                              subset.name = "var.ratio.pca", 
  #                              accept.high = 0.5)
  
  seurat_obj <- SubsetData(object = seurat_obj, 
                                subset.name = "var.ratio.pca", 
                                accept.low = 0.50)
  
  seurat_obj <- AlignSubspace(object=seurat_obj, 
                                   reduction.type = "cca", 
                                   grouping.var = "sample_ID", 
                                   dims.align = 1:n_comp)
  
  ## Visualize the aligned CCAs
  p1 <- VlnPlot(object = seurat_obj, 
                features.plot = "ACC1", 
                group.by = "sample_ID", #TODO
                do.return = TRUE, 
                point.size.use = 0.1)
  
  p2 <- VlnPlot(object = seurat_obj, 
                features.plot = "ACC2", 
                group.by = "sample_ID", #TODO
                do.return = TRUE, 
                point.size.use = 0.1)
  
  plot_grid(p1, p2)
  
  ggsave(filename =  paste0(dir_plots,data_prefix,"_",run_prefix,"_", name, "_aligned_CCs_vlnPlots.pdf"), w=12, h=8)
  
  return(seurat_obj)
  
  }
  , seurat_obj = list_seurat_obj, 
  name = names(list_seurat_obj), 
  SIMPLIFY=F, 
  .scheduling = c("dynamic"))
}

######################################################################
################### DO DIMENSIONALITY REDUCTION ######################
######################################################################

cl <- makeCluster(spec=n_cores, type = "FORK", outfile=paste0(dir_log, "dimReduction.txt"))

  list_seurat_obj <- parLapply(cl, list_seurat_obj, function(seurat_obj) {
    RunPCA(object = seurat_obj, 
           pcs.compute = min(n_comp, min(dim(seurat_obj@data))), 
           do.print = F)
  })
  
  list_seurat_obj <- parLapply(cl, list_seurat_obj, function(seurat_obj) {
    RunTSNE(object = seurat_obj, 
            reduction.use = if (is.null(align_IDs)) "pca" else "cca.aligned", 
            dims.use = 1:min(50,min(n_comp, min(dim(seurat_obj@data)))),
            do.fast=T)
  })

stopCluster(cl)

######################################################################
###################### FIND CLUSTERS, MARKERS ########################
######################################################################

if (!is.null(res_to_calc)) {
  
  cl <- makeCluster(min(length(res_to_calc),n_cores), type = "FORK", outfile=paste0(dir_log, "FindClusters_multiple_res_FindMarkers.txt"))
  
  list_seurat_obj <- lapply(list_seurat_obj, function(seurat_obj) {
    # Run FindClusters in parallel using different resolution parameter values
    list_seurat_obj_res <- parLapply(cl,
                                     res_to_calc,
                                     function(x) FindClusters(object = seurat_obj,
                                                              reduction.type = if (is.null(align_IDs)) "pca" else "cca.aligned",
                                                              dims.use = 1:n_comp,
                                                              print.output = 0,
                                                              save.SNN = T,
                                                              resolution = x))
    # Save all the cluster assignments as meta data in the original seurat object
    for (i in seq(1:length(res_to_calc))) {
      seurat_obj <- AddMetaData(seurat_obj, list_seurat_obj_res[[i]]@ident, col.name = paste0("clust.res.", res_to_calc[i]))
    }
    # No need to have this massive list in session!
    rm(list_seurat_obj_res)
    
    seurat_obj <- SetAllIdent(object = seurat_obj, id = paste0("clust.res.", res_primary))
  
    return(seurat_obj)
  })
  
} else {
  
  cl <- makeCluster(n_cores, type = "FORK", outfile=paste0(dir_log, "FindClusters_markers.txt"))
  
  list_seurat_obj <- parLapplyLB(cl, list_seurat_obj, function(seurat_obj) FindClusters(object = seurat_obj,
                                              reduction.type = if (is.null(align_IDs)) "pca" else "cca.aligned",
                                              dims.use = 1:n_comp,
                                              print.output = 0,
                                              save.SNN = T,
                                              resolution = res_primary))
}

####################### FIND CLUSTER MARKERS #########################

if (is.null(align_IDs)) {
  
  list_list_markers <- lapply(list_seurat_obj, function(seurat_obj) {
    clusters = names(table(seurat_obj@meta.data[[paste0("res.", res_primary)]]))
    list_markers <- parLapply(cl, clusters, function(cluster) tryCatch({FindMarkers(seurat_obj,  
                                                                                    do.print = T,
                                                                                    ident.1 = cluster, 
                                                                                    print.bar = T)}, 
                                                                       error=function(err){
                                                                         warning(paste0(cluster, ": FindMarkers failed with error "))
                                                                         as.character(err)}))  
  })
  
} else if (!is.null(align_IDs)) {
  list_list_markers <- lapply(list_seurat_obj, function(seurat_obj) {
    clusters = names(table(seurat_obj@meta.data[[paste0("res.", res_primary)]]))
    list_markers <- parLapply(cl, clusters, function(cluster) tryCatch({FindConservedMarkers(seurat_obj,  
                                                                                    do.print = T,
                                                                                    ident.1 = cluster, 
                                                                                    print.bar = T)}, # TODO:  do we need # grouping.var = sample_ID or align_ID ?
                                                                       error=function(err){
                                                                         warning(paste0(cluster, ": FindConservedMarkers failed with error "))
                                                                         as.character(err)}))  
  })
  
  
}

# Name the nested list entries

list_list_markers <- clusterMap(cl, function(list_markers, seurat_obj) {
  names(list_markers) <- names(table(seurat_obj@meta.data[[paste0("res.", res_primary)]]))
  list_markers
}, list_markers = list_list_markers, seurat_obj = list_seurat_obj, SIMPLIFY=F, .scheduling=c("dynamic"))


stopCluster(cl)

# name the outer list entries
names(list_list_markers) <- names(list_seurat_obj)

######################################################################
######################### IDENTIFY CELLTYPES #########################
######################################################################

# TODO

######################################################################
############################# PLOTS #################################
######################################################################

cl <- makeCluster(n_cores, type = "FORK", outfile=paste0(dir_log, "DrawPlots.txt"))

########################## PLOT T-SNE ################################

invisible(clusterMap(cl, function(seurat_obj, name) {
  p1 <- TSNEPlot(seurat_obj, do.return = T, pt.size = 1, group.by = "sample_ID", no.legend=F, plot.title=paste0(name, " by sample"))
  ggsave(p1, filename =  paste0(dir_plots, prefix_data,"_",prefix_run, "_", name, "_tSNEPlot_sample.pdf"), w=8, h=8)
  p2 <- TSNEPlot(seurat_obj, do.return = T, do.label = T, pt.size = 1, no.legend = F, plot.title = paste0(name, " by cluster")) # + xlab("t-SNE 1") + ylab("t-SNE 2"
  ggsave(p2, filename =  paste0(dir_plots, prefix_data, "_", prefix_run, "_", name, "_tSNEPlot_clust.pdf"), w=8, h=8)
}, seurat_obj = list_seurat_obj,
name = names(list_seurat_obj),
SIMPLIFY=F,
.scheduling = c("dynamic")
))

############################## FEATUREPLOTS ##############################

if (!is.null(feats_to_plot)) {

  invisible(clusterMap(cl, function(seurat_obj, name) {
    if (feats_plot_separate) {
      lapply(feats_to_plot, function(feature) {
        tryCatch({
        pdf(file = paste0(dir_plots, prefix_data,"_", prefix_run,"_", name, "_", feature,"_featurePlot.pdf"), w=8, h=8)
        FeaturePlot(seurat_obj, features.plot = feature, pt.size = 1, no.legend = F)
        }, error = function(err) warning(paste0(name, ": FeaturePlot failed for ", feature, ". Maybe it wasn't found? Error message: ", err)))
        try(dev.off())
        #ggsave(p1, filename =  paste0(dir_plots, prefix_data,"_",prefix_run,"_", name, "_", feature,"_featurePlot.pdf"), w=8, h=8)
      })
    } else {
      p1 <- FeaturePlot(seurat_obj, features.plot = feats_to_plot[1:2], overlay=T,  pt.size = 1, no.legend = F,do.return = T)
      ggsave(p1, filename =  paste0(dir_plots, prefix_data,"_",prefix_run,"_", name, "_", paste0(feats_to_plot, collapse = "_"), "_featurePlot.pdf"), w=8, h=8)
    }
  }, 
  seurat_obj = list_seurat_obj,
  name = names(list_seurat_obj),
  SIMPLIFY=F,
  .scheduling = c("dynamic")
  ))

}

stopCluster(cl)

######################################################################
########################## OUTPUT TABLES, ROBJECTS ###################
######################################################################


############################## CLUSTER GENE MARKERS ##################

cl <- makeCluster(n_cores, type = "FORK", outfile=paste0(dir_log, "outputTables_RObjects.txt"))

invisible(mapply(function(list_markers, name) {
                     clusterMap(cl, function(markers, cluster) {
                       write.csv(markers, file=paste0(dir_tables, prefix_data, "_", prefix_run, "_", name, "_", cluster, "_markers.csv"), quote = F, row.names=T)
                       }, 
                       markers = list_markers, 
                       cluster=names(list_markers), 
                       SIMPLIFY=F, 
                       .scheduling= c("dynamic"))
  }, 
 list_markers = list_list_markers, 
 name = names(list_seurat_obj), SIMPLIFY=F))

########################## SAVE SEURAT OBJECTS #######################

clusterMap(cl, function(seurat_obj, name) {
  invisible(saveRDS(seurat_obj, file = paste0(dir_RObjects, prefix_data, "_", prefix_run, "_", name, "_seurat_obj.RDS"), compress = "gzip"))
}, 
seurat_obj = list_seurat_obj, 
name = names(list_seurat_obj), 
SIMPLIFY=F, 
.scheduling= c("dynamic"))

stopCluster(cl)

######################################################################
############################### WRAP UP ##############################
######################################################################

save.image(paste0(dir_scratch, prefix_data, "_", prefix_run, "_session_image", flag_date), compress = "gzip")

message("Script done!")