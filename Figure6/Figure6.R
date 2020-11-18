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

setwd("~/scanem_pytorch/data_gen/FigureGithub/Figure6/")

suppressMessages(source("../scanem_helper_functions.R"))
heatmap_colors <- colorRampPalette(c("magenta", "black", "yellow"))(100)

outdir <- "output/"
dir.create(outdir)

# Starting with plots for P0

# Cell type labels for dataset. Row names are pool names, $cell_type1 contains cell type labels.
# Also contains information on which dataset the different pools come from, as this dataset is composed of two datasets.
curr_colData <- read.csv("data/20200710_ss_p0_exinonly_pool100_noprom_8u5d_colData.tsv", sep="\t")
curr_colData <- as.data.frame(curr_colData)
P0_category_table <- c("2/3", "2/3", "3/4/5", "4", "4/5", 
                           "5", "5/6", "6", 
                           "Inhibitory", "Inhibitory")
names(P0_category_table) <- unique(curr_colData$cell_type1)
curr_colData$Category <- P0_category_table[curr_colData$cell_type1]
curr_colData$Celltype <- curr_colData$cell_type1

# This contains the SingleCellExperiment object with pooled expression information for P0
curr_sce_p0 <- readRDS("data/P0_gene_expression_pooled.RDS")
num_celltypes <- length(unique(curr_sce_p0$cell_type1))

threshold <- 0.05

# Network output: "leave-one-out" (LOO) influence scores for 10 different models
influence_scores_hdf5 <- H5Fopen("scanem_output/20201015_All_leave_change_HDF5_ss_p0_exin_pool100_noprom_8u5d_newscanem_2.h5")
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

# Read db
db <- readLines("data/Mus_musculus.meme")
db <- db[str_detect(db, "^MOTIF")]
codes <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][2] })
names <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][3] })
names(names) <- codes

all_tom_table <- read.delim("scanem_output/p0_all_tomtom/tomtom.tsv", sep="\t",quote = "", header = T, 
                            comment.char = "#")
all_tom_table$Target_ID <- names[all_tom_table$Target_ID]
all_tom_table$Target_ID <- sapply(all_tom_table$Target_ID, function(x){
  str_remove(str_split(x, "\\)")[[1]][1], "\\(")
})
# Remove motifs that are expressed in SCE
expressed_genes <- rownames(curr_sce_p0)[rowSums(counts(curr_sce_p0))>0]
all_tom_table <- all_tom_table[all_tom_table$Target_ID %in% expressed_genes,]
# Only use alignments below alignment threshold
all_tom_table <- all_tom_table[all_tom_table$q.value < threshold, ]
all_targets <- unique(all_tom_table$Target_ID)

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

celltypes_mean_LOO <- celltypes_mean_LOO[good_motifs,,drop=FALSE]
celltype_mean_LOO_entropies <- celltype_mean_LOO_entropies[good_motifs]
motif_summary <- motif_summary[good_motifs,,drop=FALSE]
all_LOO_mat <- all_LOO_mat[good_motifs,,drop=FALSE]

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

row_annot$cluster_annot[row_annot$motif_cluster == 6] <- "Ctcf"
row_annot$cluster_annot[row_annot$motif_cluster == 10] <- "Meis2"
row_annot$cluster_annot[row_annot$motif_cluster == 13] <- "bHLH motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 3] <- "Egr/KLF motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 1] <- "Rfx motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 7] <- "Mef2c"
row_annot$cluster_annot[row_annot$motif_cluster == 5] <- "Nfia/Nfib"

colnames(row_annot) <- c("Motif cell type entropy", "Motif cluster", "Motif cluster reproducibility", 
                         "Motif cluster annotation")
ann_colors = list(
  `Motif cluster annotation` = stata_pal()(length(unique(row_annot$`Motif cluster annotation`))) # c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33", "#040B99")
)
names(ann_colors$`Motif cluster annotation`) <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))]

# Load motif activation scores
all_motif_activations_mat <- readRDS("data/P0_motif_activations_mat.RDS")
motifs_above_half_max <- apply(all_motif_activations_mat, 1, FUN=function(x){
  sum(x > (0.5 * max(x))) / length(x)
})
motif_summary$motifs_above_half_max <- motifs_above_half_max[rownames(motif_summary)]
row_annot$`Percentage of regions activated` <- motifs_above_half_max[rownames(row_annot)] * 100

celltypes_mean_LOO_aggregates <- matrix(nrow=length(unique(row_annot$`Motif cluster annotation`)), ncol=ncol(celltypes_mean_LOO))
for(i in 1:nrow(celltypes_mean_LOO_aggregates)){
  curr_cluster_annot <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))][i]
  celltypes_mean_LOO_aggregates[i,] <- colSums(celltypes_mean_LOO[row_annot$`Motif cluster annotation` == curr_cluster_annot,])
}
rownames(celltypes_mean_LOO_aggregates) <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))]
colnames(celltypes_mean_LOO_aggregates) <- colnames(celltypes_mean_LOO)
celltypes_mean_LOO_aggregates_z <- to_z(celltypes_mean_LOO_aggregates)

curr_colData$Category <- P0_category_table[curr_colData$cell_type1]
colnames(curr_colData) <- c("Celltype", "Category")

curr_annot_color <- list(
  Category = stata_pal()(length(unique(curr_colData$Category)))
)
names(curr_annot_color$Category) <- unique(curr_colData$Category)

curr_annot_col <- data.frame(row.names = unique(curr_colData$Celltype))
curr_annot_col$Category <- P0_category_table[rownames(curr_annot_col)]

curr_annot_row <- data.frame(row.names = rownames(celltypes_mean_LOO_aggregates))
curr_annot_row$amount_motifs <- sapply(rownames(celltypes_mean_LOO_aggregates), 
                                       FUN=function(x){
                                         sum(row_annot$`Motif cluster annotation` == x)
                                       })

perc_activated_per_family <- row_annot %>% group_by(`Motif cluster annotation`) %>%
  summarise(mean_perc_activ = mean(`Percentage of regions activated`)) %>%
  as.data.frame()
mean_perc_activ <- perc_activated_per_family$mean_perc_activ
names(mean_perc_activ) <- perc_activated_per_family$`Motif cluster annotation`
curr_annot_row$`Activated regions (%)` <- mean_perc_activ[rownames(curr_annot_row)] * 100

sum_loo <- c()
for(i in 1:length(unique(row_annot$`Motif cluster annotation`))){
  curr_motif_fam <- unique(row_annot$`Motif cluster annotation`)[i]
  sum_loo <- c(sum_loo, 
               sum(rowMeans(all_LOO_mat[rownames(row_annot[row_annot$`Motif cluster annotation` == curr_motif_fam,]),])))  
}
names(sum_loo) <- unique(row_annot$`Motif cluster annotation`)
curr_annot_row$mean_weight <- sum_loo[rownames(curr_annot_row)]

colnames(curr_annot_row) <- c("Motifs in cluster", "Activated regions (%)", "Summed mean influence")

curr_annot_col_2 <- curr_annot_col
colnames(curr_annot_col_2) <- c("Excitatory layer")
curr_annot_color$`Excitatory layer` <- curr_annot_color$Category

# Fig 6a, first half =====
pheatmap(celltypes_mean_LOO_aggregates_z[,P0_category_table != "Inhibitory"], color=heatmap_colors,
         angle_col=45, cellwidth=12, cellheight=12,
         annotation_col = curr_annot_col_2,
         annotation_colors = curr_annot_color, cluster_cols = FALSE, 
         border_color=NA, width = 10, height=10, 
         annotation_row = curr_annot_row,
         filename=paste0(outdir, "/Fig6a_1.pdf"))
graphics.off()
# values are added from curr_annot_row

Hlf_melt <- melt(logcounts(curr_sce_p0)["Hlf",])
Hlf_melt$Category <- curr_sce_p0$Category
Hlf_melt_ex <- Hlf_melt[curr_sce_p0$Category != "Inhibitory",]

# Fig 6c, first half =====
ggplot(Hlf_melt_ex, aes(x=Category, y=value, fill=Category)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_fill_stata() + theme_bw(base_size=14) +
  labs(x = "Excitatory layer", y="Log-transformed expression") +
  theme_Nice()
ggsave(paste0(outdir, "/Fig6c_1.pdf"),
       width=3.5,height=3.5, useDingbats=FALSE)






# Now for adult

# Cell type labels for dataset. Row names are pool names, $cell_type1 contains cell type labels.
# Also contains information on which dataset the different pools come from, as this dataset is composed of two datasets.
curr_colData <- read.csv("data/20200710_ss_adult_exinonly_pool100_noprom_8u5d_colData.tsv", sep="\t")
curr_colData <- as.data.frame(curr_colData)
Adult_category_table <- c("2", "3", "3", "4", "4", 
                          "5", "5", "5", "6", 
                          "Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory")
names(Adult_category_table) <- unique(curr_colData$cell_type1)
curr_colData$Category <- Adult_category_table[curr_colData$cell_type1]
curr_colData$Celltype <- curr_colData$cell_type1

# This contains the SingleCellExperiment object with pooled expression information for P0
curr_sce_adult <- readRDS("data/ADULT_gene_expression_pooled.RDS")
num_celltypes <- length(unique(curr_sce_adult$cell_type1))

threshold <- 0.05

# Network output: "leave-one-out" (LOO) influence scores for 10 different models
influence_scores_hdf5 <- H5Fopen("scanem_output/20201015_All_leave_change_HDF5_ss_adult_exinonly_pool100_noprom_8u5d_2.h5")
q <- h5ls(influence_scores_hdf5)
prefixes <- unique(sapply(q$name, function(x) {str_split(x, "_")[[1]][2] } ))
d <- 300
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

# Read db
db <- readLines("data/Mus_musculus.meme")
db <- db[str_detect(db, "^MOTIF")]
codes <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][2] })
names <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][3] })
names(names) <- codes

all_tom_table <- read.delim("scanem_output/adult_all_tomtom/tomtom.tsv", sep="\t",quote = "", header = T, 
                            comment.char = "#")
all_tom_table$Target_ID <- names[all_tom_table$Target_ID]
all_tom_table$Target_ID <- sapply(all_tom_table$Target_ID, function(x){
  str_remove(str_split(x, "\\)")[[1]][1], "\\(")
})
# Remove motifs that are expressed in SCE
expressed_genes <- rownames(curr_sce_adult)[rowSums(counts(curr_sce_adult))>0]
all_tom_table <- all_tom_table[all_tom_table$Target_ID %in% expressed_genes,]
# Only use alignments below alignment threshold
all_tom_table <- all_tom_table[all_tom_table$q.value < threshold, ]
all_targets <- unique(all_tom_table$Target_ID)

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

celltypes_mean_LOO <- celltypes_mean_LOO[good_motifs,,drop=FALSE]
celltype_mean_LOO_entropies <- celltype_mean_LOO_entropies[good_motifs]
motif_summary <- motif_summary[good_motifs,,drop=FALSE]
all_LOO_mat <- all_LOO_mat[good_motifs,,drop=FALSE]

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

row_annot$cluster_annot[row_annot$motif_cluster == 7] <- "bHLH motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 12] <- "T-box motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 1] <- "Egr/KLF motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 3] <- "Mef2c"
row_annot$cluster_annot[row_annot$motif_cluster == 6] <- "Nfia/Nfib"
row_annot$cluster_annot[row_annot$motif_cluster == 2] <- "Rfx motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 16] <- "Ctcf"

row_annot$cluster_annot[row_annot$motif_cluster == 15] <- "bZIP motif family 1"
row_annot$cluster_annot[row_annot$motif_cluster == 5] <- "bZIP motif family 2"

colnames(row_annot) <- c("Motif cell type entropy", "Motif cluster", "Motif cluster reproducibility", 
                         "Motif cluster annotation")
ann_colors = list(
  `Motif cluster annotation` = stata_pal()(length(unique(row_annot$`Motif cluster annotation`))) # c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33", "#040B99")
)
names(ann_colors$`Motif cluster annotation`) <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))]

# Load motif activation scores
all_motif_activations_mat <- readRDS("data/ADULT_motif_activations_mat.RDS")
motifs_above_half_max <- apply(all_motif_activations_mat, 1, FUN=function(x){
  sum(x > (0.5 * max(x))) / length(x)
})
motif_summary$motifs_above_half_max <- motifs_above_half_max[rownames(motif_summary)]
row_annot$`Percentage of regions activated` <- motifs_above_half_max[rownames(row_annot)] * 100

celltypes_mean_LOO_aggregates <- matrix(nrow=length(unique(row_annot$`Motif cluster annotation`)), ncol=ncol(celltypes_mean_LOO))
for(i in 1:nrow(celltypes_mean_LOO_aggregates)){
  curr_cluster_annot <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))][i]
  celltypes_mean_LOO_aggregates[i,] <- colSums(celltypes_mean_LOO[row_annot$`Motif cluster annotation` == curr_cluster_annot,])
}
rownames(celltypes_mean_LOO_aggregates) <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))]
colnames(celltypes_mean_LOO_aggregates) <- colnames(celltypes_mean_LOO)
celltypes_mean_LOO_aggregates_z <- to_z(celltypes_mean_LOO_aggregates)

curr_colData$Category <- Adult_category_table[curr_colData$cell_type1]
colnames(curr_colData) <- c("Celltype", "Category")

curr_annot_color <- list(
  Category = stata_pal()(length(unique(curr_colData$Category)))
)
names(curr_annot_color$Category) <- unique(curr_colData$Category)

curr_annot_col <- data.frame(row.names = unique(curr_colData$Celltype))
curr_annot_col$Category <- Adult_category_table[rownames(curr_annot_col)]

curr_annot_row <- data.frame(row.names = rownames(celltypes_mean_LOO_aggregates))
curr_annot_row$amount_motifs <- sapply(rownames(celltypes_mean_LOO_aggregates), 
                                       FUN=function(x){
                                         sum(row_annot$`Motif cluster annotation` == x)
                                       })

perc_activated_per_family <- row_annot %>% group_by(`Motif cluster annotation`) %>%
  summarise(mean_perc_activ = mean(`Percentage of regions activated`)) %>%
  as.data.frame()
mean_perc_activ <- perc_activated_per_family$mean_perc_activ
names(mean_perc_activ) <- perc_activated_per_family$`Motif cluster annotation`
curr_annot_row$`Activated regions (%)` <- mean_perc_activ[rownames(curr_annot_row)] * 100

sum_loo <- c()
for(i in 1:length(unique(row_annot$`Motif cluster annotation`))){
  curr_motif_fam <- unique(row_annot$`Motif cluster annotation`)[i]
  sum_loo <- c(sum_loo, 
               sum(rowMeans(all_LOO_mat[rownames(row_annot[row_annot$`Motif cluster annotation` == curr_motif_fam,]),])))  
}
names(sum_loo) <- unique(row_annot$`Motif cluster annotation`)
curr_annot_row$mean_weight <- sum_loo[rownames(curr_annot_row)]

colnames(curr_annot_row) <- c("Motifs in cluster", "Activated regions (%)", "Summed mean influence")

curr_annot_col_2 <- curr_annot_col
colnames(curr_annot_col_2) <- c("Excitatory layer")
curr_annot_color$`Excitatory layer` <- curr_annot_color$Category

# Fig 6a, second half =====
pheatmap(celltypes_mean_LOO_aggregates_z[,Adult_category_table != "Inhibitory"], color=heatmap_colors,
         angle_col=45, cellwidth=12, cellheight=12,
         annotation_col = curr_annot_col_2,
         annotation_colors = curr_annot_color, cluster_cols = FALSE, 
         border_color=NA, width = 10, height=10, 
         annotation_row = curr_annot_row,
         filename=paste0(outdir, "/Fig6a_2.pdf"))
graphics.off()
# values are added from curr_annot_row

curr_sce_adult$Category <- Adult_category_table[curr_sce_adult$cell_type1]
Hlf_melt <- melt(logcounts(curr_sce_adult)["Hlf",])
Hlf_melt$Category <- curr_sce_adult$Category
Hlf_melt_ex <- Hlf_melt[curr_sce_adult$Category != "Inhibitory",]

# Fig 6c, first half =====
ggplot(Hlf_melt_ex, aes(x=Category, y=value, fill=Category)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_fill_stata() + theme_bw(base_size=14) +
  labs(x = "Excitatory layer", y="Log-transformed expression") +
  theme_Nice(angled = FALSE)
ggsave(paste0(outdir, "/Fig6c_2.pdf"),
       width=3.5,height=3.5, useDingbats=FALSE)





all_LOO_mat_selection_aggregates_melt_p0 <- readRDS("data/20201018_all_LOO_mat_selection_aggregates_melt_P0.RDS")
all_LOO_mat_selection_aggregates_melt_adult <- readRDS("data/20201018_all_LOO_mat_selection_aggregates_melt_adult.RDS")

all_LOO_mat_selection_aggregates_melt_adult$Category <- Adult_category_table[all_LOO_mat_selection_aggregates_melt_adult$Celltype]

# Fig 6b =====
ggplot(all_LOO_mat_selection_aggregates_melt_adult[all_LOO_mat_selection_aggregates_melt_adult$Category != "Inhibitory" &
                                                     all_LOO_mat_selection_aggregates_melt_adult$`Motif cluster annotation` == "Hlf",], 
       aes(x=Category, y=`Aggregate of motif weights`, 
           fill=`Category`)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_x_reordered() +
  scale_fill_stata() + theme_bw(base_size=14) +
  labs(x = "Excitatory layer", y="Aggregate of motif influence scores") +
  theme_Nice(angled = FALSE)
ggsave(paste0(outdir, "/Fig6b.pdf"),
       width=3.5,height=3.5, useDingbats=FALSE)





curr_sce_p0 <- curr_sce_p0[,curr_sce_p0$Category != "Inhibitory"]
curr_sce_adult <- curr_sce_adult[,curr_sce_adult$Category != "Inhibitory"]

get_melted_exp <- function(genes, sce, layer_annot=NULL){
  curr_mat <- melt(colSums(logcounts(sce[genes,])))
  colnames(curr_mat) <- "Average expression"
  curr_mat$Celltype <- colData(sce)$cell_type1
  if(!is.null(layer_annot)){
    curr_mat$Layer <- layer_annot[curr_mat$Celltype]
  }
  return(curr_mat)
}
Nfia_Nfib_excit_p0 <- get_melted_exp(c("Nfia", "Nfib"), curr_sce_p0, P0_category_table)
Nfia_Nfib_excit_adult <- get_melted_exp(c("Nfia", "Nfib"), curr_sce_adult, Adult_category_table)

all_LOO_mat_selection_aggregates_melt_p0_excit <- all_LOO_mat_selection_aggregates_melt_p0[!str_detect(all_LOO_mat_selection_aggregates_melt_p0$Celltype, "^In"),]
all_LOO_mat_selection_aggregates_melt_adult_excit <- all_LOO_mat_selection_aggregates_melt_adult[!str_detect(all_LOO_mat_selection_aggregates_melt_adult$Celltype, "^In"),]
all_LOO_mat_selection_aggregates_melt_p0_excit$Pool <- as.character(all_LOO_mat_selection_aggregates_melt_p0_excit$Pool)
all_LOO_mat_selection_aggregates_melt_p0_excit$Layer <- sapply(all_LOO_mat_selection_aggregates_melt_p0_excit$Category, FUN=function(x){str_remove(x,"Layer ")})
all_LOO_mat_selection_aggregates_melt_adult_excit$Pool <- as.character(all_LOO_mat_selection_aggregates_melt_adult_excit$Pool)
all_LOO_mat_selection_aggregates_melt_adult_excit$Layer <- all_LOO_mat_selection_aggregates_melt_adult_excit$Category

Nfia_Nfib_LOO_p0 <- all_LOO_mat_selection_aggregates_melt_p0_excit[all_LOO_mat_selection_aggregates_melt_p0_excit$`Motif cluster annotation` == "Nfia/Nfib",]
Nfia_Nfib_LOO_adult <- all_LOO_mat_selection_aggregates_melt_adult_excit[all_LOO_mat_selection_aggregates_melt_adult_excit$`Motif cluster annotation` == "Nfia/Nfib",]

Nfia_Nfib_df_p0 <- data.frame(expression = Nfia_Nfib_excit_p0$`Average expression`, 
                              agg_LOO = Nfia_Nfib_LOO_p0$`Aggregate of motif weights`, 
                              pool = Nfia_Nfib_LOO_p0$Pool, 
                              layer = Nfia_Nfib_LOO_p0$Layer)
Nfia_Nfib_df_adult <- data.frame(expression = Nfia_Nfib_excit_adult$`Average expression`, 
                                 agg_LOO = Nfia_Nfib_LOO_adult$`Aggregate of motif weights`, 
                                 pool = Nfia_Nfib_LOO_adult$Pool, 
                                 layer = Nfia_Nfib_LOO_adult$Layer)
p1 <- ggplot(Nfia_Nfib_df_p0, aes(x=expression, y=agg_LOO, color=layer)) + 
  geom_point() + theme_bw(base_size=12) + 
  theme_Nice() +
  labs(x="Nfia/Nfib", y="", 
       color="Layer") + scale_color_stata() 
p2 <- ggplot(Nfia_Nfib_df_adult, aes(x=expression, y=agg_LOO, color=layer)) + 
  geom_point() + theme_bw(base_size=12) + 
  theme_Nice() +
  labs(x="Nfia/Nfib", y="", 
       color="Layer") + scale_color_stata() 
p1/p2

Mef2c_excit_p0 <- get_melted_exp(c("Mef2c"), curr_sce_p0, P0_category_table)
Mef2c_excit_adult <- get_melted_exp(c("Mef2c"), curr_sce_adult, Adult_category_table)
Mef2c_LOO_p0 <- all_LOO_mat_selection_aggregates_melt_p0_excit[all_LOO_mat_selection_aggregates_melt_p0_excit$`Motif cluster annotation` == "Mef2c",]
Mef2c_LOO_adult <- all_LOO_mat_selection_aggregates_melt_adult_excit[all_LOO_mat_selection_aggregates_melt_adult_excit$`Motif cluster annotation` == "Mef2c",]
Mef2c_df_p0 <- data.frame(expression = Mef2c_excit_p0$`Average expression`, 
                          agg_LOO = Mef2c_LOO_p0$`Aggregate of motif weights`, 
                          pool = Mef2c_LOO_p0$Pool, 
                          layer = Mef2c_LOO_p0$Layer)
Mef2c_df_adult <- data.frame(expression = Mef2c_excit_adult$`Average expression`, 
                             agg_LOO = Mef2c_LOO_adult$`Aggregate of motif weights`, 
                             pool = Mef2c_LOO_adult$Pool, 
                             layer = Mef2c_LOO_adult$Layer)
p3 <- ggplot(Mef2c_df_p0, aes(x=expression, y=agg_LOO, color=layer)) + 
  geom_point() + theme_bw(base_size=12) + 
  theme_Nice() +
  labs(x="Mef2c", y="", 
       color="Layer") + scale_color_stata() 
p4 <- ggplot(Mef2c_df_adult, aes(x=expression, y=agg_LOO, color=layer)) + 
  geom_point() + theme_bw(base_size=12) + 
  theme_Nice() +
  labs(x="Mef2c", y="", 
       color="Layer") + scale_color_stata() 

Rfx3_excit_p0 <- get_melted_exp(c("Rfx3"), curr_sce_p0, P0_category_table)
Rfx3_excit_adult <- get_melted_exp(c("Rfx3"), curr_sce_adult, Adult_category_table)
Rfx3_LOO_p0 <- all_LOO_mat_selection_aggregates_melt_p0_excit[all_LOO_mat_selection_aggregates_melt_p0_excit$`Motif cluster annotation` == "Rfx motif family",]
Rfx3_LOO_adult <- all_LOO_mat_selection_aggregates_melt_adult_excit[all_LOO_mat_selection_aggregates_melt_adult_excit$`Motif cluster annotation` == "Rfx motif family",]
Rfx3_df_p0 <- data.frame(expression = Rfx3_excit_p0$`Average expression`, 
                         agg_LOO = Rfx3_LOO_p0$`Aggregate of motif weights`, 
                         pool = Rfx3_LOO_p0$Pool, 
                         layer = Rfx3_LOO_p0$Layer)
Rfx3_df_adult <- data.frame(expression = Rfx3_excit_adult$`Average expression`, 
                            agg_LOO = Rfx3_LOO_adult$`Aggregate of motif weights`, 
                            pool = Rfx3_LOO_adult$Pool, 
                            layer = Rfx3_LOO_adult$Layer)
p5 <- ggplot(Rfx3_df_p0, aes(x=expression, y=agg_LOO, color=layer)) + 
  geom_point() + theme_bw(base_size=12) + 
  theme_Nice() +
  labs(x="Rfx3", y="", 
       color="Layer") + scale_color_stata() 
p6 <- ggplot(Rfx3_df_adult, aes(x=expression, y=agg_LOO, color=layer)) + 
  geom_point() + theme_bw(base_size=12) + 
  theme_Nice() +
  labs(x="Rfx3", y="", 
       color="Layer") + scale_color_stata() 


Tcf4_excit_p0 <- get_melted_exp(c("Tcf4"), curr_sce_p0, P0_category_table)
Tcf4_excit_adult <- get_melted_exp(c("Tcf4"), curr_sce_adult, Adult_category_table)
Tcf4_LOO_p0 <- all_LOO_mat_selection_aggregates_melt_p0_excit[all_LOO_mat_selection_aggregates_melt_p0_excit$`Motif cluster annotation` == "bHLH motif family",]
Tcf4_LOO_adult <- all_LOO_mat_selection_aggregates_melt_adult_excit[all_LOO_mat_selection_aggregates_melt_adult_excit$`Motif cluster annotation` == "bHLH motif family",]
Tcf4_df_p0 <- data.frame(expression = Tcf4_excit_p0$`Average expression`, 
                         agg_LOO = Tcf4_LOO_p0$`Aggregate of motif weights`, 
                         pool = Tcf4_LOO_p0$Pool, 
                         layer = Tcf4_LOO_p0$Layer)
Tcf4_df_adult <- data.frame(expression = Tcf4_excit_adult$`Average expression`, 
                            agg_LOO = Tcf4_LOO_adult$`Aggregate of motif weights`, 
                            pool = Tcf4_LOO_adult$Pool, 
                            layer = Tcf4_LOO_adult$Layer)
p7 <- ggplot(Tcf4_df_p0, aes(x=expression, y=agg_LOO, color=layer)) + 
  geom_point() + theme_bw(base_size=12) + 
  theme_Nice() + theme(legend.position="right") +
  labs(x="Tcf4", y="", 
       color="Layer") + scale_color_stata() 
p8 <- ggplot(Tcf4_df_adult, aes(x=expression, y=agg_LOO, color=layer)) + 
  geom_point() + theme_bw(base_size=12) + 
  theme_Nice() + theme(legend.position="right") +
  labs(x="Tcf4", y="", 
       color="Layer") + scale_color_stata()

# Fig 6f =====

(p1/p2)|(p3/p4)|(p5/p6)|(p7/p8)
ggsave(filename=paste0(outdir, "/Fig6f.pdf"), 
       width=8,height=4, useDingbats=FALSE)
# All correlations found using e.g. 
cor(Tcf4_df_p0$expression,  Tcf4_df_p0$agg_LOO, method="spearman")  



# Coefficient of variation
cv <- function(x){
  sd(x)/mean(x)
}

row_annot_P0 <- readRDS("data/20201016_rownames_LOO_P0.RDS")
row_annot_ADULT <- readRDS("data/20201016_rownames_LOO_ADULT.RDS")
celltypes_mean_LOO_P0 <- readRDS("data/20201016_celltypes_mean_LOO_P0.RDS")
celltypes_mean_LOO_ADULT <- readRDS("data/20201016_celltypes_mean_LOO_ADULT.RDS")

# Limit to excitatory celltypes
celltypes_mean_LOO_P0 <- celltypes_mean_LOO_P0[,!str_detect(colnames(celltypes_mean_LOO_P0), "^In")]
celltypes_mean_LOO_ADULT <- celltypes_mean_LOO_ADULT[,!str_detect(colnames(celltypes_mean_LOO_ADULT), "^In")]

res1_mean_cv <- sapply(unique(row_annot_P0$`Motif cluster annotation`), FUN=function(x){
  return(mean(apply(celltypes_mean_LOO_P0[row_annot_P0[
    rownames(celltypes_mean_LOO_P0),]$`Motif cluster annotation` == x,],
    1, 
    FUN=function(q) { cv(q[!str_detect(names(q), "^In")])} )))
})

res2_mean_cv <- sapply(unique(row_annot_ADULT$`Motif cluster annotation`), FUN=function(x){
  return(mean(apply(celltypes_mean_LOO_ADULT[row_annot_ADULT[
    rownames(celltypes_mean_LOO_ADULT),]$`Motif cluster annotation` == x,],
    1, 
    FUN=function(q) { cv(q[!str_detect(names(q), "^In")])} )))
})

CV_df <- data.frame(meancv=c(res1_mean_cv, res2_mean_cv),
                    experiment = factor(c(rep("P0", length(res1_mean_cv)), rep("Adult", length(res2_mean_cv))), levels=c("P0", "Adult")),
                    motif_family = c(unique(row_annot_P0$`Motif cluster annotation`), unique(row_annot_ADULT$`Motif cluster annotation`)))


# Fig 6e =====
CV_df[CV_df$motif_family %in% intersect(names(res1_mean_cv), names(res2_mean_cv)),] %>%
  ggplot(aes(x=motif_family, y=meancv, fill=experiment)) +
  geom_col(position="dodge") + theme_bw(base_size=14) + 
  theme_Nice() + theme(legend.position = "right") +
  labs(x="Motif family", y="Mean CV", fill="Origin") +
  scale_fill_stata()
ggsave(filename=paste0(outdir, "Fig6e.pdf"),
       width=4.5, height=4.5)


# Fig 6d was plotted using Tomtom from the MEME suite:
# Shobhit Gupta, JA Stamatoyannopolous, Timothy Bailey and William Stafford Noble, 
# "Quantifying similarity between motifs", Genome Biology, 8(2):R24, 2007. 
# Afterwards, it was formatted using Illustrator 

