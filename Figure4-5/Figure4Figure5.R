suppressMessages(library(SingleCellExperiment))
suppressMessages(library(ggplot2))
suppressMessages(library(igraph))
suppressMessages(library(rhdf5))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(stringr))
suppressMessages(library(ggrepel))
suppressMessages(library(ggthemes))
suppressMessages(library(gghalves))
suppressMessages(library(tidytext))
suppressMessages(library(reshape2))
suppressMessages(library(pheatmap))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(mixtools))

options(stringsAsFactors = FALSE)

suppressMessages(source("../scover_helper_functions.R"))
heatmap_colors <- colorRampPalette(c("magenta", "black", "yellow"))(100)

outdir <- "output/"
dir.create(outdir)

# Cell type labels for dataset. Row names are pool names, $cell_type1 contains cell type labels.
# Also contains information on which dataset the different pools come from, as this dataset is composed of two datasets.
curr_colData <- read.csv("data/20200912_tm_all_organs_pool80_5u5d_colData.tsv", sep="\t")
curr_colData <- as.data.frame(curr_colData)

# Adding category labels for the included tabula muris (TM) cell types
category_annot_tm <- c("epithelial", "epithelial", # bladder
                       "connective", "immune",
                       "macroglial", "pericyte",
                       "endothelial", "neuronal", 
                       "macroglial", "macroglial", # latter = precursor
                       "epithelial", "epithelial", # colon
                       "epithelial", "immune", # colon, fat
                       "endothelial", "immune", 
                       "connective", "immune", 
                       "immune", "immune", 
                       "immune", "muscle", # fat, heart
                       "endothelial", "endothelial", 
                       "adipose", "connective", 
                       "immune", "muscle",
                       "epithelial", "endothelial", # kidney, liver 
                       "epithelial", "endothelial", # liver, lung
                       "connective", "epithelial",
                       "epithelial", "epithelial", # mammary
                       "connective", "immune", # mammary, marrow
                       "immune", "immune",
                       "immune", "immune",
                       "immune", "immune",
                       "immune", "immune", # mammary, muscle
                       "endothelial", "connective", # latter = stem
                       "muscle", "muscle", # latter = stem
                       "endocrine", "exocrine", # pancreas
                       "endocrine", "epithelial",
                       "endocrine", "endocrine",
                       "epithelial", "epithelial", # skin
                       "epithelial", "immune", # skin, spleen (latter = stem)
                       "immune", "immune", # spleen, thymus
                       "epithelial", "epithelial", # tongue
                       "epithelial", "immune", # trachea
                       "connective")  
names(category_annot_tm) <- unique(curr_colData$cell_type1)
curr_colData$Category <- category_annot_tm[curr_colData$cell_type1]
curr_colData <- curr_colData[,c(1,2,5)]
colnames(curr_colData) <- c("Celltype", "Origin", "Category")


# Original dataset. Contains pooled expression values that will be used for 
# the influence score VS TF expression correlations
curr_data <- read.csv("data/20200912_tm_all_organs_pool80_5u5d.tsv.gz", 
                      sep="\t")
# Get values by separating comma-separated values in the $ind column
lc <- matrix(nrow=nrow(curr_data), ncol=length(str_split(curr_data$ind[1], ",")[[1]]))
for(i in 1:nrow(lc)){
  lc[i,] <- str_split(curr_data$ind[i], ",")[[1]]
}
lc <- t(apply(lc, 1, as.numeric))
rownames(lc) <- curr_data$gene
colnames(lc) <- rownames(curr_colData)
curr_colData <- as.data.frame(curr_colData)
curr_rowData <- data.frame(row.names=rownames(lc), feature_symbol=rownames(lc))
curr_sce <- SingleCellExperiment(assays=list(logcounts=as.matrix(lc)), colData=as.data.frame(curr_colData), 
                                 rowData=as.data.frame(curr_rowData))
curr_sce <- curr_sce[,!is.na(curr_sce$Celltype)]
num_celltypes <- length(unique(curr_sce$Celltype))

# Network output: "leave-one-out" (LOO) influence scores for 10 different models
influence_scores_hdf5 <- H5Fopen("scover_output/20201013_All_leave_change_HDF5_tm_all_7_organs_pool80_5u5d_600mot_3.h5")
q <- h5ls(influence_scores_hdf5)
prefixes <- unique(sapply(q$name, function(x) {str_split(x, "_")[[1]][2] } ))
d <- 600
all_LOO <- list()
for(i in 1:length(prefixes)){
  curr_prefix <- prefixes[i]
  
  curr_LOO <- paste0("LOO_",curr_prefix)
  curr_LOO <- rhdf5::h5read(influence_scores_hdf5, curr_LOO)
  
  colnames(curr_LOO) <- paste0(curr_prefix, "_", 0:(d-1))
  
  all_LOO[[curr_prefix]] <- curr_LOO
}
h5closeAll()
all_LOO_mat <- do.call(cbind, all_LOO)
all_LOO_mat <- t(all_LOO_mat)
dim(all_LOO_mat)
colnames(all_LOO_mat) <- rownames(curr_colData)

# Loading tomtom database
db <- readLines("data/Mus_musculus.meme")
db <- db[str_detect(db, "^MOTIF")]
motif_codes <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][2] })
motif_names <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][3] })
names(motif_names) <- motif_codes
# Loading tomtom alignments
all_tom_table <- read.delim("scover_output/all_tomtom/tomtom.tsv", sep="\t",quote = "", header = T, 
                            comment.char = "#")
all_tom_table$Target_ID <- motif_names[all_tom_table$Target_ID]
all_tom_table$Target_ID <- sapply(all_tom_table$Target_ID, function(x){
  str_remove(str_split(x, "\\)")[[1]][1], "\\(")
})
# Remove TFs that are not expressed in SCE
expressed_genes <- rownames(curr_sce)[rowSums(logcounts(curr_sce))>0]
all_tom_table <- all_tom_table[all_tom_table$Target_ID %in% expressed_genes,]
all_tom_table <- all_tom_table[all_tom_table$q.value < 0.05, ] # retain significant only
all_targets <- unique(all_tom_table$Target_ID)

# Construct graph using Tomtom alignments
# Get edge connections from table
edge_vect <- c()
for(i in 1:nrow(all_tom_table)){
  edge_vect <- c(edge_vect, all_tom_table[i,]$Query_ID, all_tom_table[i,]$Target_ID)
}
edge_weights <- all_tom_table$q.value
g <- graph(edge_vect, directed = F)
V(g)$type <- bipartite_mapping(g)$type # Make it a bipartite graph
# Add log-transformed q-values as edge weights
E(g)$weight <- -log2(edge_weights) 
# Cluster graph
g <- simplify(g)
fgc <- igraph::cluster_walktrap(g)
fgc_groups <- igraph::groups(fgc)



# Get cluster reproducibilities from graph
n_candidates <- 10
group_reproducibilities <- sapply(fgc_groups,FUN=function(x, ncan=n_candidates){
  prefixes <- stringr::str_extract(x, "^best|(^[0-9]+)")
  prefixes <- prefixes[!is.na(prefixes)]
  num_exp_found <- length(unique(prefixes))
  return( (num_exp_found/ncan) )
})
cluster_motif_detectors <- sapply(fgc_groups, FUN=function(x, query_options=all_tom_table$Query_ID){
  return(x[x %in% query_options])
})
cluster_not_motif_detectors <- sapply(fgc_groups, FUN=function(x, query_options=all_tom_table$Query_ID){
  return(x[!(x %in% query_options)])
})
cluster_motif_detectors <- sapply(X = cluster_motif_detectors, FUN=function(X){
  return(paste0(X, collapse=" "))
})
cluster_alignments <- sapply(X = cluster_not_motif_detectors, FUN=function(X){
  return(paste0(X, collapse=" "))
})
cluster_df <- data.frame(cluster_motif_detectors, cluster_alignments, cluster_reproducibility=group_reproducibilities)
cluster_df <- cluster_df[order(cluster_df$cluster_reproducibility, decreasing = T),]


celltypes <- unique(str_remove(colnames(all_LOO_mat), "_pool[0-9]+$"))
celltypes_mean_LOO <- matrix(nrow=nrow(all_LOO_mat), ncol=length(celltypes))
for(i in 1:length(celltypes)){
  curr_LOO_mat <- all_LOO_mat[,str_detect(colnames(all_LOO_mat), fixed(celltypes[i])), drop=FALSE]
  celltypes_mean_LOO[,i] <- apply(curr_LOO_mat, 1, mean)
}
colnames(celltypes_mean_LOO) <- celltypes
rownames(celltypes_mean_LOO) <- rownames(all_LOO_mat)

# Convert to z-scores 
celltypes_mean_LOO_z <- to_z(celltypes_mean_LOO)

# Get normalised "cell type entropy" scores across cell types (shannon ent)
celltype_mean_LOO_entropies <- apply(celltypes_mean_LOO, 1, shan_ent)
num_groups <- ncol(celltypes_mean_LOO)
theoretical_max <- -log2(1/num_groups)
celltype_mean_LOO_entropies <- celltype_mean_LOO_entropies / theoretical_max
motif_impacts <- apply(celltypes_mean_LOO, 1, function(x) sum(x))
motif_summary <- data.frame(row.names = names(celltype_mean_LOO_entropies), entropy=celltype_mean_LOO_entropies, impact=motif_impacts)

hist(log10(0.1+motif_summary$impact+abs(min(motif_summary$impact))))
mix_model <- normalmixEM(log10(5+log10(0.1+motif_summary$impact+abs(min(motif_summary$impact)))), k=2)

motif_summary$post_1 <- mix_model$posterior[,1]
motif_summary$post_2 <- mix_model$posterior[,2]

if((sum(motif_summary$post_1 > 0.6) + sum(motif_summary$post_2 > 0.6)) < nrow(motif_summary)){
  cat("No clear bimodal distribution: likely not a lot of 'dead motifs'\n")
  
  good_motifs <- rownames(motif_summary)
} else {
  cat("Some 'dead motifs' found\n")
  
  good_cluster <- which(mix_model$mu == max(mix_model$mu))
  good_motifs <- rownames(motif_summary)[mix_model$posterior[,good_cluster] > 0.9]
}

motif_summary$is_good_motif <- ifelse(rownames(motif_summary)%in%good_motifs, "yes", "no")
motif_summary$model <- str_remove(rownames(motif_summary), "_[0-9]+$")
motif_summary$ent_imp <- motif_summary$entropy*motif_summary$impact

celltypes_mean_LOO <- celltypes_mean_LOO[good_motifs,,drop=FALSE]
celltype_mean_LOO_entropies <- celltype_mean_LOO_entropies[good_motifs]
motif_summary <- motif_summary[good_motifs,,drop=FALSE]
all_LOO_mat <- all_LOO_mat[good_motifs,,drop=FALSE]


# Generate annotation for heatmaps
row_annot <- data.frame(entropy=celltype_mean_LOO_entropies, row.names = rownames(motif_summary))
motif_cluster <- c()
cluster_df_motifs <- stringr::str_split(cluster_df$cluster_motif_detectors, " ")
for(i in 1:nrow(row_annot)){
  found <- F
  for(j in 1:length(cluster_df_motifs)){
    if(rownames(row_annot)[i] %in% cluster_df_motifs[[j]]){
      motif_cluster <- c(motif_cluster, rownames(cluster_df)[j]) 
      found <- T
      break
    }
  }
  if(!found){
    motif_cluster <- c(motif_cluster, "not aligned")
  }
}
motif_cluster <- factor(motif_cluster, levels = c("not aligned", rownames(cluster_df)))
row_annot$motif_cluster <- motif_cluster

cluster_reprod <- c()
for(i in 1:nrow(row_annot)){
  if(row_annot$motif_cluster[i] == "not aligned") {
    cluster_reprod <- c(cluster_reprod, "NA")
    next
  }
  cluster_reprod <- c(cluster_reprod, cluster_df[as.character(row_annot$motif_cluster[i]),]$cluster_reproducibility)
}
row_annot$cluster_reprod <- cluster_reprod

sub_selection <- row_annot$motif_cluster != "not aligned" & 
  row_annot$cluster_reprod >= 0.5

celltypes_mean_LOO_z <- to_z(celltypes_mean_LOO)

celltypes_mean_LOO_z <- celltypes_mean_LOO_z[sub_selection,]
celltypes_mean_LOO <- celltypes_mean_LOO[sub_selection,]
row_annot <- row_annot[sub_selection,]
all_LOO_mat_selection <- all_LOO_mat[sub_selection,]

# Set default motif cluster annotations
cluster_names <- sapply(cluster_df$cluster_alignments, FUN=function(x){
  if(str_detect(x, "\\(")){
    if(length(str_split(x, " ")[[1]]) > 1){
      if(length(str_split(x, " ")[[1]]) == sum(str_detect(str_split(x, " ")[[1]], "\\("))){
        curr_names <- str_remove_all(sapply(str_split(str_split(x, " ")[[1]], "_"), FUN=function(x) { return(x[1]) }), "\\(|\\)")
        if(length(curr_names) <= 5){
          names <- curr_names
        } else {
          names <- c(curr_names[1:5], "...")
        }
      } else {
        curr_names <- str_split(x, " ")[[1]]
        curr_names <- str_remove_all(sapply(str_split(str_split(x, " ")[[1]], "_"), FUN=function(x) { return(x[1]) }), "\\(|\\)")
        if(length(curr_names) <= 5){
          names <- curr_names
        } else {
          names <- c(curr_names[1:5], "...")
        }
      }
    } else {
      curr_names <- str_split(x, " ")[[1]]
      names <- curr_names
    }
  } else {
    curr_names <- str_split(x, " ")[[1]]
    if(length(curr_names) <= 5){
      names <- curr_names
    } else {
      names <- c(curr_names[1:5], "...")
    }
  }
  return(paste(names, collapse = "/"))
})

cluster_df$cluster_annot <- cluster_names

cluster_annot <- c()
for(i in 1:nrow(row_annot)){
  if(row_annot$motif_cluster[i] == "not aligned") {
    cluster_annot <- c(cluster_annot, "NA")
    next
  }
  cluster_annot <- c(cluster_annot, cluster_df[as.character(row_annot$motif_cluster[i]),]$cluster_annot)
}
row_annot$cluster_annot <- cluster_annot

row_annot$cluster_annot[row_annot$motif_cluster == 12] <- "ETS motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 2] <- "Egr/KLF motif families"
# row_annot$cluster_annot[row_annot$motif_cluster == 10] <- "ETS motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 5] <- "bHLH/bZIP family"
row_annot$cluster_annot[row_annot$motif_cluster == 16] <- "Yy1"
row_annot$cluster_annot[row_annot$motif_cluster == 6] <- "Zfp143/Tbx2/Six5"

colnames(row_annot) <- c("Motif cell type entropy", "Motif cluster", "Motif cluster reproducibility", 
                         "Motif cluster annotation")

ann_colors = list(
  `Motif cluster annotation` = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33", "#040B99")[1:length(unique(row_annot$`Motif cluster annotation`))]
)
names(ann_colors$`Motif cluster annotation`) <- names(table(row_annot$`Motif cluster annotation`)[order(table(row_annot$`Motif cluster annotation`), decreasing = T)])


# Subcluster motif families 

k_max <- 10
mot_fams <- unique(row_annot$`Motif cluster annotation`)[
  order(unique(row_annot$`Motif cluster annotation`))]
print(mot_fams)

# Subcluster bHLH family
graphics.off()
curr_fam <- c("bHLH/bZIP family")
curr_annot <- row_annot[row_annot$`Motif cluster annotation` %in% curr_fam,]
wss <- sapply(1:k_max, 
              function(k){kmeans(celltypes_mean_LOO_z[row_annot$`Motif cluster annotation` %in% curr_fam,], k, nstart=50,iter.max = 15 )$tot.withinss})
plot(wss)
d1 <- diff(wss); k <- which.max(abs(diff(d1) / diff(wss[-1]))) 
k <- 2
curr_annot$subcluster <- cutree(hclust(dist(celltypes_mean_LOO_z[row_annot$`Motif cluster annotation` %in% curr_fam,])), k = k)[rownames(curr_annot)]
row_annot$`Motif cluster annotation 2` <- row_annot$`Motif cluster annotation`
for(i in 1:nrow(row_annot)){
  if(rownames(row_annot)[i] %in% rownames(curr_annot)){
    curr_mot <- rownames(row_annot)[i]
    row_annot$`Motif cluster annotation 2`[i] <- paste(curr_fam, 1:k)[curr_annot[curr_mot,]$subcluster]
  }
}
ann_colors = list(
  `Motif cluster annotation 2` = stata_pal()(length(unique(row_annot$`Motif cluster annotation 2`)))
)
names(ann_colors$`Motif cluster annotation 2`) <- names(table(row_annot$`Motif cluster annotation 2`))[order(names(table(row_annot$`Motif cluster annotation 2`)), decreasing = F)]

# Subcluster ETS family
graphics.off()
curr_fam <- "ETS motif family"
curr_annot <- row_annot[row_annot$`Motif cluster annotation` %in% curr_fam,]
wss <- sapply(1:k_max, 
              function(k){kmeans(celltypes_mean_LOO_z[row_annot$`Motif cluster annotation` %in% curr_fam,], k, nstart=50,iter.max = 15 )$tot.withinss})
plot(wss)
d1 <- diff(wss); k <- which.max(abs(diff(d1) / diff(wss[-1]))) 
print(k)
curr_annot$subcluster <- cutree(hclust(dist(celltypes_mean_LOO_z[row_annot$`Motif cluster annotation` %in% curr_fam,])), k = k)[rownames(curr_annot)]
for(i in 1:nrow(row_annot)){
  if(rownames(row_annot)[i] %in% rownames(curr_annot)){
    curr_mot <- rownames(row_annot)[i]
    row_annot$`Motif cluster annotation 2`[i] <- paste(curr_fam, 1:k)[curr_annot[curr_mot,]$subcluster]
  }
}
ann_colors = list(
  `Motif cluster annotation 2` = stata_pal()(length(unique(row_annot$`Motif cluster annotation 2`)))
)
names(ann_colors$`Motif cluster annotation 2`) <- names(table(row_annot$`Motif cluster annotation 2`))[order(names(table(row_annot$`Motif cluster annotation 2`)), decreasing = F)]


# Subclustered
row_annot_backup <- row_annot
row_annot$`Motif cluster annotation` <- row_annot$`Motif cluster annotation 2`
repr_cl_names <- unique(row_annot$`Motif cluster annotation 2`)

pool_exps <- list()
cluster_TFs <- list()
all_clusters_corr_dfs <- list()
all_cluster_corrs <- list()
all_cluster_pvals <- list()
all_cluster_mean_exps <- list()
for(i in 1:length(repr_cl_names)){
  curr_cluster <- repr_cl_names[i]
  pool_exps[[curr_cluster]] <- colMeans(all_LOO_mat_selection[row_annot$`Motif cluster annotation 2` == curr_cluster,])
  
  curr_motifs <- rownames(row_annot)[row_annot$`Motif cluster annotation 2` == curr_cluster]
  curr_cluster_tfs <- unique(all_tom_table[all_tom_table$Query_ID %in% curr_motifs,]$Target_ID)
  curr_cluster_tfs <- curr_cluster_tfs[curr_cluster_tfs %in% rownames(curr_sce)]
  curr_cluster_annot <- curr_cluster
  
  curr_pool_exps <- pool_exps[[curr_cluster]]
  curr_tom_table <- all_tom_table[all_tom_table$Query_ID %in% curr_motifs,]
  
  curr_cluster_df <- logcounts(curr_sce[curr_cluster_tfs,])
  if(nrow(curr_cluster_df) != 0){
    curr_cluster_df <- curr_cluster_df[,names(curr_pool_exps), drop=FALSE]
    
    curr_corrs <- c()
    curr_pvals <- c()
    curr_clusters_corr_dfs <- list()
    if(!is.null(nrow(curr_cluster_df))){
      for(j in 1:nrow(curr_cluster_df)){
        curr_corrs <- c(curr_corrs, cor(curr_pool_exps, curr_cluster_df[j,], method="spearman"))
        curr_pvals <- c(curr_pvals, cor.test(curr_pool_exps, curr_cluster_df[j,], method="spearman", exact = FALSE)$p.value)
        curr_corr_df <- data.frame(x=curr_pool_exps, y=curr_cluster_df[j,])
        curr_corr_df$cell_type <- curr_sce$Celltype
        curr_corr_df$Category <- category_annot_tm[curr_sce$Celltype]
        curr_corr_df$cell_type_label <- curr_corr_df$Category
        curr_corr_df$cell_type_label[duplicated(curr_corr_df$cell_type_label)] <- NA
        curr_clusters_corr_dfs[[rownames(curr_cluster_df)[j]]] <- curr_corr_df
      }
    } else {
      curr_corrs <- c(curr_corrs, cor(curr_pool_exps, curr_cluster_df, method="spearman"))
      curr_pvals <- c(curr_pvals, cor.test(curr_pool_exps, curr_cluster_df, method="spearman", exact = FALSE)$p.value)
      curr_clusters_corr_dfs[[rownames(curr_cluster_df)]] <- data.frame(x=curr_pool_exps, y=curr_cluster_df)
    }
    
    names(curr_corrs) <- rownames(curr_cluster_df)
    names(curr_pvals) <- rownames(curr_cluster_df)
    
    all_cluster_corrs[[curr_cluster]] <- curr_corrs
    all_cluster_pvals[[curr_cluster]] <- curr_pvals
    all_clusters_corr_dfs[[curr_cluster]] <- curr_clusters_corr_dfs
    all_cluster_mean_exps[[curr_cluster]] <- rowMeans(curr_cluster_df)
  } else {
    cat("nrow(curr_cluster_df) == 0")
  }
}
all_cluster_mean_exps_df <- melt(all_cluster_mean_exps) %>% magrittr::set_colnames(c("expression", "cluster"))
all_cluster_mean_exps_df$cluster_annot <- all_cluster_mean_exps_df$cluster
all_corrs_df <- melt(all_cluster_corrs) %>% magrittr::set_colnames(c("corr", "cluster"))
all_cluster_mean_exps_df$corr <- all_corrs_df$corr
all_pvals_df <- melt(all_cluster_pvals) %>% magrittr::set_colnames(c("pval", "cluster"))
all_cluster_mean_exps_df$pval <- all_pvals_df$pval
all_cluster_pvals_fdr <- lapply(all_cluster_pvals, FUN=p.adjust, method="fdr")
all_pvals_fdr_df <- melt(all_cluster_pvals_fdr) %>% magrittr::set_colnames(c("pval_fdr", "cluster"))
all_cluster_mean_exps_df$pval_fdr <- all_pvals_fdr_df$pval_fdr
all_cluster_corrs_df <- stack(all_cluster_corrs)
all_names <- c()
for(i in all_cluster_corrs){
  all_names <- c(all_names, names(i))
}
all_cluster_mean_exps_df$TF <- all_names
all_cluster_mean_exps_df$expression_2 <- all_cluster_mean_exps_df$expression
all_cluster_mean_exps_df$expression_2[all_cluster_mean_exps_df$pval_fdr > 0.05] <- NA
all_cluster_mean_exps_df$tf_labels <- NA
top_n <- 5
for(i in 1:length(unique(all_cluster_mean_exps_df$cluster_annot))){
  curr_cluster <- unique(all_cluster_mean_exps_df$cluster)[i]
  curr_corr_tf_df <- all_cluster_mean_exps_df[all_cluster_mean_exps_df$cluster == curr_cluster,,drop=FALSE]
  curr_corr_tf_df <- curr_corr_tf_df[order(abs(curr_corr_tf_df$corr), decreasing = TRUE),,drop=FALSE]
  annot_tfs <- curr_corr_tf_df$TF[1:(ifelse(nrow(curr_corr_tf_df) < top_n, nrow(curr_corr_tf_df), top_n))]
  all_cluster_mean_exps_df$tf_labels[all_cluster_mean_exps_df$cluster == curr_cluster &
                                       all_cluster_mean_exps_df$TF %in% annot_tfs] <- 
    all_cluster_mean_exps_df$TF[all_cluster_mean_exps_df$cluster == curr_cluster &
                                  all_cluster_mean_exps_df$TF %in% annot_tfs]
}
all_cluster_mean_exps_df$tf_labels[all_cluster_mean_exps_df$pval_fdr > 0.05] <- ""

ann_colors$`Motif cluster annotation` <- ann_colors$`Motif cluster annotation 2`

# Fig. S4a =====
pheatmap::pheatmap(celltypes_mean_LOO_z[order(row_annot$`Motif cluster annotation`),,drop=FALSE], 
                   cluster_rows = F, cluster_cols = T, show_rownames = F,
                   annotation_row = row_annot[,c("Motif cell type entropy", "Motif cluster annotation")], 
                   border_color=NA, annotation_colors=ann_colors, 
                   color = heatmap_colors,
                   angle_col = 45, cellwidth = 12, cellheight = .4, width=17, height=20, 
                   filename=paste0(outdir, "/FigS4a.pdf"))
graphics.off()


# Fig. 5c =====
ggplot(all_cluster_mean_exps_df, aes(x=cluster_annot, y=corr, color=expression_2)) +
  geom_jitter(width = 0) + 
  geom_label_repel(label=all_cluster_mean_exps_df$tf_labels, size=4) +
  theme_bw(base_size=14) + 
  theme_Nice() + theme(legend.position = "right") +
  labs(x="Motif cluster name", 
       y="Spearman R", color="Mean TF expression across pools")
ggsave(filename=paste0(outdir, "/Fig5c.pdf"), 
       width = 11, height=6, useDingbats=FALSE)


# Aggregate for pools
pool_LOO_aggregates <- matrix(nrow=length(unique(row_annot$`Motif cluster annotation`)), ncol=ncol(all_LOO_mat))
for(i in 1:nrow(pool_LOO_aggregates)){
  curr_cluster_annot <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))][i]
  pool_LOO_aggregates[i,] <- colSums(all_LOO_mat[rownames(row_annot),][row_annot$`Motif cluster annotation` == curr_cluster_annot,])
}
rownames(pool_LOO_aggregates) <- unique(all_cluster_mean_exps_df$cluster_annot)[order(unique(all_cluster_mean_exps_df$cluster_annot))]
colnames(pool_LOO_aggregates) <- colnames(all_LOO_mat)



# Read activations scores

motif_hdf5 <- H5Fopen("scover_output/All_motif_activations_HDF5_tm_all_7_organs_pool80_5u5d_600mot_3.h5")
q <- h5ls(motif_hdf5)
prefixes <- unique(sapply(q$name, function(x) {str_split(x, "_")[[1]][2] } ))

d <- 600

all_motif_activations <- list()
for(i in 1:length(prefixes)){
  curr_prefix <- prefixes[i]
  
  curr_motif_act <- paste0("ACTIVATIONS_",curr_prefix)
  curr_motif_act <- rhdf5::h5read(motif_hdf5, curr_motif_act)
  
  rownames(curr_motif_act) <- paste0(curr_prefix, "_", 0:(d-1))
  
  all_motif_activations[[curr_prefix]] <- curr_motif_act
}

h5closeAll()
sapply(all_motif_activations, dim)
all_motif_activations[["best"]] <- all_motif_activations[["best"]][1:600,1:1918]
sapply(all_motif_activations, dim)
all_motif_activations_mat <- do.call(rbind, all_motif_activations)

# Also add motif annotation to motif_summary
motif_summary$annotation <- sapply(rownames(motif_summary), FUN=function(x){
  if(x %in% rownames(row_annot)){
    return(row_annot[x,]$`Motif cluster annotation`)
  } else {
    return("Not aligned")
  }
})

motif_impacts_2 <- apply(all_LOO_mat, 1, function(x) sum(x))
motif_summary$impact_2 <- motif_impacts_2[rownames(motif_summary)]

# Now add information to motif summary that tells how many of the motifs are above 0.5 * max(activation)
motifs_above_half_max <- apply(all_motif_activations_mat, 1, FUN=function(x){
  sum(x > (0.5 * max(x))) / length(x)
})
motif_summary$motifs_above_half_max <- motifs_above_half_max[rownames(motif_summary)]

row_annot$`Percentage of promoters activated` <- motifs_above_half_max[rownames(row_annot)]*100

celltypes_mean_LOO_aggregates <- matrix(nrow=length(unique(row_annot$`Motif cluster annotation`)), ncol=ncol(celltypes_mean_LOO))
for(i in 1:nrow(celltypes_mean_LOO_aggregates)){
  curr_cluster_annot <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))][i]
  celltypes_mean_LOO_aggregates[i,] <- colSums(celltypes_mean_LOO[row_annot$`Motif cluster annotation` == curr_cluster_annot,])
}
rownames(celltypes_mean_LOO_aggregates) <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))]
colnames(celltypes_mean_LOO_aggregates) <- colnames(celltypes_mean_LOO)

celltypes_mean_LOO_top100 <- celltypes_mean_LOO[order(abs(rowSums(celltypes_mean_LOO)), decreasing = TRUE)[1:100],]
row_annot_top100 <- row_annot[order(rowSums(celltypes_mean_LOO), decreasing = TRUE)[1:100],]
celltypes_mean_LOO_aggregates_z <- to_z(celltypes_mean_LOO_aggregates)
curr_colData <- as.data.frame(curr_colData)
curr_colData$Celltype <- curr_colData$Celltype
curr_colData$Category <- category_annot_tm[curr_colData$Celltype]

curr_annot_color <- list(
  Category = stata_pal()(length(unique(curr_colData$Category))) #,
  #Dataset = stata_pal()(13)[c(12,13)]
)
names(curr_annot_color$Category) <- unique(curr_colData$Category)
colnames(curr_colData) <- c("Celltype", "Dataset", "Category")

curr_annot_col <- data.frame(row.names = unique(curr_colData$Celltype),
                             Dataset = sapply(unique(curr_colData$Celltype), FUN=function(x){str_split(x, " ")[[1]][1]}))
curr_annot_col$Category <- category_annot_tm[rownames(curr_annot_col)]
curr_colData$Category <- category_annot_tm[curr_colData[,1]]

# Fig S4b =====
pheatmap(celltypes_mean_LOO_aggregates, color=heatmap_colors, 
         angle_col=45, cellwidth=12, cellheight=12, 
         annotation_col = curr_annot_col[,c("Category"), drop=FALSE], 
         annotation_colors = curr_annot_color, border_color=NA, 
         filename=paste0(outdir, "/FigS4b.pdf"))
graphics.off()

category_means <- matrix(nrow=nrow(celltypes_mean_LOO), ncol=length(unique(curr_colData$Category)))
for(i in 1:ncol(category_means)){
  curr_cat <- unique(curr_colData$Category)[i]
  category_means[,i] <- rowMeans(celltypes_mean_LOO[,category_annot_tm[colnames(celltypes_mean_LOO)] == curr_cat,drop=FALSE])
}
colnames(category_means) <- unique(curr_colData$Category)
rownames(category_means) <- rownames(celltypes_mean_LOO)
category_means_z <- to_z(category_means)

category_means_aggregates <- matrix(nrow=length(unique(row_annot$`Motif cluster annotation`)), ncol=ncol(category_means))
for(i in 1:nrow(category_means_aggregates)){
  curr_cluster_annot <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))][i]
  category_means_aggregates[i,] <- colSums(category_means[row_annot$`Motif cluster annotation` == curr_cluster_annot,])
}
rownames(category_means_aggregates) <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))]
colnames(category_means_aggregates) <- colnames(category_means)
category_means_aggregates_z <- to_z(category_means_aggregates)

# Add heatmap annotation information: 
# 1) amount of motifs in family 
# 2) mean percent of promoters above 1/2 max activation 
# 3) summed mean influence scores for motif family motifs
curr_category_annot <- data.frame(row.names=rownames(category_means_aggregates_z), # 1)
                                  amount_motifs = rep("", nrow(category_means_aggregates_z)))
curr_category_annot$amount_motifs <- sapply(rownames(category_means_aggregates_z), 
                                            FUN=function(x){
                                              sum(row_annot$`Motif cluster annotation` == x)
                                            })
mean_perc_act <- row_annot %>% group_by(`Motif cluster annotation`) %>% # 2)
  summarise(mean=mean(`Percentage of promoters activated`)) %>% data.frame(row.names=1)
curr_category_annot$mean_perc_act <- mean_perc_act$mean
sum_loo <- c() # 3)
for(i in 1:length(unique(row_annot$`Motif cluster annotation`))){
  curr_motif_fam <- unique(row_annot$`Motif cluster annotation`)[i]
  sum_loo <- c(sum_loo, 
               sum(rowMeans(all_LOO_mat[rownames(row_annot[row_annot$`Motif cluster annotation` == curr_motif_fam,]),])))  
}
names(sum_loo) <- unique(row_annot$`Motif cluster annotation`)
curr_category_annot$mean_weight <- sum_loo[rownames(curr_category_annot)]
colnames(curr_category_annot) <- c("Motifs in cluster", "Activated promoters (%)", "Summed mean influence")

# Fig 4a =====
pheatmap(category_means_aggregates_z, color=heatmap_colors, 
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, 
         useDingbats=FALSE,
         annotation_row = curr_category_annot,
         filename = paste0(outdir, "Fig4a.pdf"))
graphics.off()
# then, annotation was added using values in "curr_category_annot"

# Do PCA:
prcomp_mat <- prcomp(t(all_LOO_mat_selection))
eigs <- prcomp_mat$sdev^2
eigs[1]/sum(eigs)
eigs[2]/sum(eigs)
variances_expl <- eigs/sum(sum(eigs))
prcomp_mat <- prcomp_mat$x
prcomp_mat <- as.data.frame(prcomp_mat)
prcomp_mat$Celltype <- curr_colData[rownames(prcomp_mat),]$Celltype
prcomp_mat$Category <- curr_colData[rownames(prcomp_mat),]$Category
prcomp_mat$Origin <- curr_colData[rownames(prcomp_mat),]$Origin
# Add labels randomly:
prcomp_mat <- prcomp_mat[sample(1:nrow(prcomp_mat), replace = FALSE),]
prcomp_mat$ggrepel_labels <- prcomp_mat$Category
prcomp_mat$ggrepel_labels[duplicated(prcomp_mat$ggrepel_labels)] <- ""

# Fig 4c =====
ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=Category)) + geom_point() +
  scale_color_stata() + theme_Nice(angled = FALSE) + theme(legend.position = "right") +
  labs(x=paste0("PC1 (", round(variances_expl[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl[2]*100, 2), "%)")) +
  geom_label_repel(label=prcomp_mat$ggrepel_labels, size=4, show.legend = FALSE) + 
  coord_fixed()
ggsave(paste0(outdir, "Fig4c.pdf"), width=7, height=7,
       useDingbats=FALSE)


# Aggregate scores
LOO_mat_melted <- melt(all_LOO_mat_selection) %>% magrittr::set_colnames(c("Motif", "Pool", "Weight"))
LOO_mat_melted$Celltype <- curr_colData[LOO_mat_melted$Pool,1]
LOO_mat_melted$`Motif cluster annotation` <- row_annot[as.character(LOO_mat_melted$Motif),]$`Motif cluster annotation`

all_LOO_mat_selection_aggregates <- matrix(nrow=length(unique(LOO_mat_melted$`Motif cluster annotation`)), 
                                           ncol=ncol(all_LOO_mat_selection))
for(i in 1:nrow(all_LOO_mat_selection_aggregates)){
  curr_motif_family <- unique(LOO_mat_melted$`Motif cluster annotation`)[i]
  all_LOO_mat_selection_aggregates[i,] <- colSums(all_LOO_mat_selection[rownames(all_LOO_mat_selection) 
                                                                        %in% rownames(row_annot[row_annot$`Motif cluster annotation` ==
                                                                                                  curr_motif_family,]),])
}
rownames(all_LOO_mat_selection_aggregates) <- unique(LOO_mat_melted$`Motif cluster annotation`)
colnames(all_LOO_mat_selection_aggregates) <- colnames(all_LOO_mat_selection)

all_LOO_mat_selection_aggregates_melt <- melt(all_LOO_mat_selection_aggregates) %>% 
  magrittr::set_colnames(c("Motif cluster annotation", "Pool", "Aggregate of motif influence scores"))
all_LOO_mat_selection_aggregates_melt$Celltype <- curr_colData[all_LOO_mat_selection_aggregates_melt$Pool,1]
all_LOO_mat_selection_aggregates_melt$Category <- category_annot_tm[all_LOO_mat_selection_aggregates_melt$Celltype]

all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` <- factor(x=as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`),
                                                                           levels=unique(as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))[order(unique(as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)))])


# Fig 5a =====
ggplot(all_LOO_mat_selection_aggregates_melt, 
       aes(x=reorder_within(Category,`Aggregate of motif influence scores`,`Motif cluster annotation`), y=`Aggregate of motif influence scores`, 
           fill=`Category`)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_x_reordered() +
  scale_fill_stata() + theme_bw(base_size=14) + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x = "Category", y="Aggregate of motif influence scores") +
  theme_Nice() +
  theme(legend.position = "none") +
  facet_wrap(~`Motif cluster annotation`, scales="free")
ggsave(paste0(outdir, "/Fig5a.pdf"),
       width=11,height=8, useDingbats=FALSE)

# Fig 4d =====
pheatmap(cor(t(all_LOO_mat_selection_aggregates), method="spearman"), 
         color=heatmap_colors,
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, useDingbats=FALSE,
         filename = paste0(outdir, "Fig4d.pdf"))
graphics.off()


# GO term analysis: correlations of motif influence scores with expression of GO-term related genes:
GO_terms <- readLines("data/go_scfind.tsv")
names(GO_terms) <- sapply(GO_terms, FUN=function(x){return(str_split(x, "\t")[[1]][1])})
GO_terms <- sapply(GO_terms, FUN=function(x){
  str_split(str_split(x, "\t")[[1]][2], ",")[[1]]
})
curr_GO_lengths <- c()
pb <- txtProgressBar(0, length(GO_terms), style=3)
for(i in 1:length(GO_terms)){
  curr_GO_name <- names(GO_terms)[i]
  curr_GO_genes <- GO_terms[[i]]
  curr_GO_genes <- unique(curr_GO_genes[curr_GO_genes %in% rownames(curr_sce)])
  curr_GO_length <- length(curr_GO_genes)
  curr_GO_lengths <- c(curr_GO_lengths, curr_GO_length)
  setTxtProgressBar(pb, i)
}

# Select GO terms with at least 1 gene and fewer than 50 genes in set (that are found
# in the current experiment)
GO_selection <- GO_terms[which(curr_GO_lengths < 50 & curr_GO_lengths > 0)]
curr_GO_lengths_selection <- curr_GO_lengths[curr_GO_lengths < 50 & curr_GO_lengths > 0]
GO_selection_corrs <- list()
GO_selection_pvals <- list()
pb <- txtProgressBar(0, length(GO_selection), style=3)
for(j in 1:length(GO_selection)){ # Correlate expression of GO term genes to motif influence scores. This can take a while
  curr_GO_name <- names(GO_selection)[j]
  curr_GO_genes <- GO_selection[[j]]
  curr_GO_genes <- unique(curr_GO_genes[curr_GO_genes %in% rownames(curr_sce)])
  curr_GO_expression <- colMeans(logcounts(curr_sce[curr_GO_genes,]))
  
  curr_GO_corrs <- c()
  curr_GO_pvals <- c()
  for(i in 1:length(unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))){
    curr_cluster_annot <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)[i]
    curr_melt_aggregates <- all_LOO_mat_selection_aggregates_melt[all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` == curr_cluster_annot,]
    curr_GO_corrs <- c(curr_GO_corrs, cor(curr_melt_aggregates$`Aggregate of motif influence scores`, curr_GO_expression, 
                                          method="spearman"))
    curr_GO_pvals <- c(curr_GO_pvals, p.adjust(cor.test(curr_melt_aggregates$`Aggregate of motif influence scores`, 
                                                        curr_GO_expression, 
                                                        method="spearman", 
                                                        exact = FALSE)$p.value, method="fdr"))
  }
  names(curr_GO_corrs) <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)
  names(curr_GO_pvals) <- names(curr_GO_corrs)
  
  GO_selection_corrs[[j]] <- curr_GO_corrs
  GO_selection_pvals[[j]] <- curr_GO_pvals
  
  setTxtProgressBar(pb, j)
}
names(GO_selection_corrs) <- names(GO_selection)
names(GO_selection_pvals) <- names(GO_selection)

GO_melt <- melt(GO_selection_corrs)
GO_melt$motif_family <- names(GO_selection_corrs[[1]])
GO_melt_p <- melt(GO_selection_pvals)
GO_melt$p_corr <- GO_melt_p$value
colnames(GO_melt) <- c("Spearman R", "GO term", "Motif cluster annotation", "Corrected p-value")
GO_melt_cast <- dcast(GO_melt, formula=`GO term`~`Motif cluster annotation`, value.var="Spearman R")
rownames(GO_melt_cast) <- GO_melt_cast[,1]
GO_melt_cast <- GO_melt_cast[,-c(1)]

GO_melt_cast_pval <- dcast(GO_melt, formula=`GO term`~`Motif cluster annotation`, value.var="Corrected p-value")
GO_melt_cast_pval <- as.data.frame(GO_melt_cast_pval)
rownames(GO_melt_cast_pval) <- GO_melt_cast_pval[,1]
GO_melt_cast_pval <- GO_melt_cast_pval[,-c(1)]

bottom_top_quantiles <- quantile(as.numeric(unlist(GO_melt_cast)), c(0.01, 0.99))
has_no_significant <- apply(GO_melt_cast, 1, FUN=function(x){
  return((sum(x < bottom_top_quantiles[1]) + 
            sum(x > bottom_top_quantiles[2])) 
         == 0)
})
sum(!has_no_significant)

# Cluster annotation
clust_annot <- cutree(hclust(dist(GO_melt_cast[!has_no_significant,])), k=2)
GO_row_annot <- data.frame(row.names=names(clust_annot),
                           "GO term cluster"=as.character(clust_annot))
colnames(GO_row_annot) <- "GO term cluster"
GO_row_annot_col <- list(
  "GO term cluster" = stata_pal()(2)
)
names(GO_row_annot_col$`GO term cluster`) <- c("1", "2")


# Fig 5b =====
pheatmap::pheatmap(GO_melt_cast[!has_no_significant,],
                   show_rownames = FALSE, border_color = NA,
                   color=heatmap_colors, cellwidth = 12,
                   angle_col = 45, cellheight = .2, 
                   annotation_row = GO_row_annot, 
                   annotation_colors = GO_row_annot_col,
                   filename = paste0(outdir, "/Fig5b.pdf"),
                   useDingbats=FALSE)
graphics.off()
significant_cluster <- cutree(hclust(dist(GO_melt_cast[!has_no_significant,])), k=2)
# The annotation for GO terms was added in Illustrator, and it was found by sampling from these two:
cat(names(significant_cluster[significant_cluster == 1])[sample(length(significant_cluster[significant_cluster == 1]), size=10)], sep="\n")
cat(names(significant_cluster[significant_cluster == 2])[sample(length(significant_cluster[significant_cluster == 2]), size=10)], sep="\n")


# Supplementary figs ==================
# Tbx Zfp143 Bach2

# Motif weight / FUBP1 expression plot
Zfp143Tbx2Six5_agg_across_pools <- colSums(all_LOO_mat_selection[row_annot$`Motif cluster annotation` == "Zfp143/Tbx2/Six5",])
Tbx2_exp_across_pools <- logcounts(curr_sce)["Tbx2",]
Tbx2_corr_df <- data.frame(agg_LOO = Zfp143Tbx2Six5_agg_across_pools, 
                           exp=Tbx2_exp_across_pools, 
                           cat=curr_colData$Category)

Tbx2_plot <- ggplot(Tbx2_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled = FALSE) +
  labs(color="Category", 
       x="Zfp143/Tbx2/Six5", 
       y="Tbx2 expression") + 
  scale_color_stata() + theme(aspect.ratio=1) + scale_x_continuous(breaks=c(0.1,0.2,0.3))

Zfp143_exp_across_pools <- logcounts(curr_sce)["Zfp143",]
Zfp143_corr_df <- data.frame(agg_LOO = Zfp143Tbx2Six5_agg_across_pools, 
                             exp=Zfp143_exp_across_pools, 
                             cat=curr_colData$Category)
Zfp143_plot <- ggplot(Zfp143_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled = FALSE) +
  labs(color="Category", 
       x="Zfp143/Tbx2/Six5", 
       y="Zfp143 expression") + 
  scale_color_stata() + theme(aspect.ratio=1) + scale_x_continuous(breaks=c(0.2,0.3))

bHLHbZIPfamily1_agg_across_pools <- colSums(all_LOO_mat_selection[row_annot$`Motif cluster annotation` == "bHLH/bZIP family 1",])
Bach2_exp_across_pools <- logcounts(curr_sce)["Bach2",]
Bach2_corr_df <- data.frame(agg_LOO = bHLHbZIPfamily1_agg_across_pools, 
                            exp=Bach2_exp_across_pools, 
                            cat=curr_colData$Category)

Bach2_plot <- ggplot(Bach2_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled = FALSE) + theme(legend.position="right") +
  labs(color="Category", 
       x="bHLH/bZIP family 1", 
       y="Bach2 expression") + 
  scale_color_stata() + theme(aspect.ratio=1) + scale_x_continuous(breaks = round(seq(0.5, 1, by = 0.25),2))

# Fig 5d =====
Tbx2_plot | Zfp143_plot | Bach2_plot # requires "patchwork" package
ggsave(paste0(outdir, "/Fig5d.pdf"), width=8, height=4, useDingbats=FALSE)

# Correlations were added using 
cor(Tbx2_corr_df$agg_LOO, Tbx2_corr_df$exp, method="spearman")
cor(Zfp143_corr_df$agg_LOO, Zfp143_corr_df$exp, method="spearman")
cor(Bach2_corr_df$agg_LOO, Bach2_corr_df$exp, method="spearman")
