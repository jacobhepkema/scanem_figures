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

options(stringsAsFactors = FALSE)

suppressMessages(source("../scanem_helper_functions.R"))

outdir <- "output/"
dir.create(outdir)

# Cell type labels for kidney dataset. Row names are pool names, $cell_type1 contains cell type labels.
# Also contains information on which dataset the different pools come from, as this dataset is composed of two datasets.
curr_colData <- read.csv("data/20200825_combined_mature_fetal_kidney_pool120_5u5d_colData.tsv", sep="\t")
curr_colData <- as.data.frame(curr_colData)

# Adding category labels for the kidney cell types
Category_table <- c("Immune", "Nephron progenitor",
                    "Immune", "Immune", 
                    "Nephron epithelium", "Nephron progenitor", 
                    "Nephron progenitor", "Endothelium",
                    "Stroma", "Stroma",
                    "Immune", "Nephron epithelium", 
                    "Immune", "Immune", 
                    "Nephron progenitor", "Immune", 
                    "Immune", "Stroma", 
                    "Stroma", "Immune", 
                    "Immune", "Nephron progenitor",
                    "Nephron epithelium", "Immune", 
                    "Nephron progenitor", "Immune", 
                    "Nephron progenitor", "Immune", 
                    "Immune", "Stroma", 
                    "Stroma", "Nephron progenitor", 
                    "Nephron progenitor", "Nephron epithelium", 
                    "Nephron progenitor", "Stroma",
                    
                    "Endothelium", "Immune", 
                    "Immune", "Immune", 
                    "Nephron epithelium", "Endothelium", 
                    "Nephron epithelium", "Nephron progenitor",
                    "Endothelium", "Nephron epithelium", 
                    "Immune", "Immune", 
                    "Immune", "Stroma", 
                    "Immune", "Immune", 
                    "Nephron epithelium", "Endothelium", 
                    "Endothelium", "Nephron epithelium", 
                    "Nephron epithelium", "Nephron epithelium", 
                    "Nephron epithelium", "Nephron epithelium", 
                    "Nephron epithelium")
names(Category_table) <- unique(curr_colData$cell_type1)
curr_colData$Category <- Category_table[curr_colData$cell_type1]
colnames(curr_colData) <- c("Celltype", "Origin", "Category")


# Original dataset. Contains pooled expression values that will be used for 
# the influence score VS TF expression correlations
curr_data <- read.csv("data/20200825_combined_mature_fetal_kidney_pool120_5u5d.tsv.gz", 
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
influence_scores_hdf5 <- H5Fopen("scanem_output/All_leave_one_out_change_HDF5_combined_matfet_pool120_5u5d_600mot.h5")
prefixes <- unique(sapply(h5ls(influence_scores_hdf5)$name, function(x) {str_split(x, "_")[[1]][2] } ))
d <- 600 # Number of motifs in model
all_LOO <- list()
for(i in 1:length(prefixes)){ # Cycle through 10 candidate models
  curr_prefix <- prefixes[i]
  curr_LOO <- paste0("LOO_",curr_prefix)
  curr_LOO <- rhdf5::h5read(influence_scores_hdf5, curr_LOO)
  colnames(curr_LOO) <- paste0(curr_prefix, "_", 0:(d-1))
  all_LOO[[curr_prefix]] <- curr_LOO
}
h5closeAll()
all_LOO_mat <- do.call(cbind, all_LOO) # Combine LOO scores for different models
all_LOO_mat <- t(all_LOO_mat)
dim(all_LOO_mat) # 6000 motifs by 529 pools
colnames(all_LOO_mat) <- rownames(curr_colData)

# Loading tomtom database
db <- readLines("data/Homo_sapiens.meme")
db <- db[str_detect(db, "^MOTIF")]
motif_codes <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][2] })
motif_names <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][3] })
names(motif_names) <- motif_codes
# Loading tomtom alignments
all_tom_table <- read.delim("scanem_output/all_tomtom/tomtom.tsv", sep="\t",quote = "", header = T, 
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

# Get mean influence scores in cell types
celltypes_mean_LOO <- matrix(nrow=nrow(all_LOO_mat), ncol=length(unique(curr_colData$Celltype)))
for(i in 1:ncol(celltypes_mean_LOO)){
  curr_LOO_mat <- all_LOO_mat[,curr_colData$Celltype == unique(curr_colData$Celltype)[i], drop=FALSE]
  celltypes_mean_LOO[,i] <- apply(curr_LOO_mat, 1, mean)
}
colnames(celltypes_mean_LOO) <- unique(curr_colData$Celltype)
rownames(celltypes_mean_LOO) <- rownames(all_LOO_mat)
celltypes_mean_LOO <- as.data.frame(celltypes_mean_LOO)

# Convert to z-scores 
celltypes_mean_LOO_z <- to_z(celltypes_mean_LOO)

# Get normalised "cell type entropy" scores across cell types (shannon ent)
celltype_mean_LOO_entropies <- apply(celltypes_mean_LOO, 1, shan_ent)
num_groups <- ncol(celltypes_mean_LOO)
theoretical_max <- -log2(1/num_groups)
celltype_mean_LOO_entropies <- celltype_mean_LOO_entropies / theoretical_max
motif_impacts <- apply(celltypes_mean_LOO, 1, function(x) sum(abs(x)))
motif_summary <- data.frame(row.names = names(celltype_mean_LOO_entropies), entropy=celltype_mean_LOO_entropies, impact=motif_impacts)
motif_summary$model <- str_remove(rownames(motif_summary), "_[0-9]+$")

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

# Select reproducible motifs that aligned to actual motifs
sub_selection <- row_annot$motif_cluster != "not aligned" & 
  row_annot$cluster_reprod >= 0.5
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
# Adjust to better names
row_annot$cluster_annot[row_annot$motif_cluster == 7] <- "ETS motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 11] <- "ZBTB33/BRCA1"
row_annot$cluster_annot[row_annot$motif_cluster == 2] <- "EGR/KLF motif families"
row_annot$cluster_annot[row_annot$motif_cluster == 9] <- "YY1"
row_annot$cluster_annot[row_annot$motif_cluster == 3] <- "bZIP motif family"
row_annot$cluster_annot[row_annot$motif_cluster == 1] <- "ZFX"

# Better colnames for heatmap annotation
colnames(row_annot) <- c("Motif cell type entropy", "Motif cluster", "Motif cluster reproducibility", 
                         "Motif cluster annotation")
# Generate colors for heatmap annotation using stata_pal() from ggthemes
ann_colors = list(
  `Motif cluster annotation` = stata_pal()(length(unique(row_annot$`Motif cluster annotation`)))
)
names(ann_colors$`Motif cluster annotation`) <- names(table(row_annot$`Motif cluster annotation`)[order(table(row_annot$`Motif cluster annotation`), decreasing = T)])

# Subcluster ETS and update annotations

k_max <- 10
mot_fams <- unique(row_annot$`Motif cluster annotation`)[
  order(unique(row_annot$`Motif cluster annotation`))]
graphics.off()

print(mot_fams)
curr_fam <- c("ETS motif family")

# Elbow method for determining number of clusters
curr_annot <- row_annot[row_annot$`Motif cluster annotation` %in% curr_fam,]
wss <- sapply(1:k_max, 
              function(k){kmeans(celltypes_mean_LOO_z[row_annot$`Motif cluster annotation` %in% curr_fam,], k, nstart=50,iter.max = 15 )$tot.withinss})

d1 <- diff(wss); k <- which.max(abs(diff(d1) / diff(wss[-1])))
print(k)

curr_annot$subcluster <- cutree(hclust(dist(celltypes_mean_LOO_z[row_annot$`Motif cluster annotation` %in% curr_fam,])), k = k)[rownames(curr_annot)]
row_annot$`Motif cluster annotation 2` <- row_annot$`Motif cluster annotation`
for(i in 1:nrow(row_annot)){
  if(rownames(row_annot)[i] %in% rownames(curr_annot)){
    curr_mot <- rownames(row_annot)[i]
    row_annot$`Motif cluster annotation 2`[i] <- paste("ETS motif family", 1:k)[curr_annot[curr_mot,]$subcluster]
  }
}
row_annot$`Motif cluster annotation` <- row_annot$`Motif cluster annotation 2`
row_annot$`Motif cluster annotation 2` <- NULL
ann_colors = list(
  `Motif cluster annotation` = stata_pal()(length(unique(row_annot$`Motif cluster annotation`)))
)
names(ann_colors$`Motif cluster annotation`) <- names(table(row_annot$`Motif cluster annotation`))[order(names(table(row_annot$`Motif cluster annotation`)), decreasing = F)]


# Fig. S2a =====
heatmap_colors <- colorRampPalette(c("magenta", "black", "yellow"))(100)
pheatmap::pheatmap(celltypes_mean_LOO_z[order(row_annot$`Motif cluster annotation`),,drop=FALSE], 
                   cluster_rows = F, cluster_cols = T, show_rownames = F,
                   annotation_row = row_annot[,c("Motif cell type entropy", "Motif cluster annotation")], 
                   color = heatmap_colors, 
                   border_color=NA, annotation_colors=ann_colors, 
                   angle_col = 45, cellwidth = 12, cellheight = .4, width=17, height=20,
                   filename=paste0(outdir, "/FigS2a.pdf")) 




# Read activations scores

motif_hdf5 <- H5Fopen("scanem_output/20201020_All_motif_activations_HDF5_combined_matfet_pool120_5u5d_600mot.h5")
q <- h5ls(motif_hdf5)

all_motif_activations <- list()
for(i in 1:length(prefixes)){
  curr_prefix <- prefixes[i]
  
  curr_motif_act <- paste0("ACTIVATIONS_",curr_prefix)
  curr_motif_act <- rhdf5::h5read(motif_hdf5, curr_motif_act)
  
  rownames(curr_motif_act) <- paste0(curr_prefix, "_", 0:(d-1))
  
  all_motif_activations[[curr_prefix]] <- curr_motif_act
}
h5closeAll()

dimensions <- sapply(all_motif_activations, dim)
# Get rid of the odd sequence that prevents concatenation
all_motif_activations <- lapply(all_motif_activations, FUN=function(x){
  x <- x[1:d, 1:min(dimensions[2,])]
  return(x)
})
sapply(all_motif_activations, dim) # now the second dimensions are the same
# Make into one matric:
all_motif_activations_mat <- do.call(rbind, all_motif_activations)

# Also add motif annotation to motif_summary
motif_summary$annotation <- sapply(rownames(motif_summary), FUN=function(x){
  if(x %in% rownames(row_annot)){
    return(row_annot[x,]$`Motif cluster annotation`)
  } else {
    return("Not aligned")
  }
})

# Now add information to motif summary that tells how many of the motifs are above 0.5 * max(activation)
motifs_above_half_max <- apply(all_motif_activations_mat, 1, FUN=function(x){
  sum(x > (0.5 * max(x))) / length(x)
})
motif_summary$motifs_above_half_max <- motifs_above_half_max[rownames(motif_summary)]


# Aggregate scores
graphics.off()
celltypes_mean_LOO_aggregates <- matrix(nrow=length(unique(row_annot$`Motif cluster annotation`)), ncol=ncol(celltypes_mean_LOO))
for(i in 1:nrow(celltypes_mean_LOO_aggregates)){
  curr_cluster_annot <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))][i]
  celltypes_mean_LOO_aggregates[i,] <- colSums(celltypes_mean_LOO[row_annot$`Motif cluster annotation` == curr_cluster_annot,])
}
rownames(celltypes_mean_LOO_aggregates) <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))]
colnames(celltypes_mean_LOO_aggregates) <- colnames(celltypes_mean_LOO)
# Z transform:
celltypes_mean_LOO_aggregates_z <- to_z(celltypes_mean_LOO_aggregates)

# Column annotations
curr_annot_color <- list(
  Category = stata_pal()(length(unique(curr_colData$Category))),
  Dataset = stata_pal()(7)[c(6,7)]
)
names(curr_annot_color$Category) <- unique(curr_colData$Category)
names(curr_annot_color$Dataset) <- unique(curr_colData$Origin)

curr_annot_col <- data.frame(row.names = unique(curr_colData$Celltype),
                             Dataset = sapply(unique(curr_colData$Celltype), FUN=function(x){str_split(x, " ")[[1]][1]}))
curr_annot_col$Category <- Category_table[rownames(curr_annot_col)]

# Fig. S2b =====
pheatmap(celltypes_mean_LOO_aggregates, color=heatmap_colors,
         angle_col=45, cellwidth=12, cellheight=12,
         annotation_col = curr_annot_col,
         border_color = NA,
         annotation_colors = curr_annot_color, 
         filename = paste0(outdir,"/FigS2b.pdf"))

# Calculate means in categories
category_means <- matrix(nrow=nrow(celltypes_mean_LOO), ncol=length(unique(curr_colData$Category)))
for(i in 1:ncol(category_means)){
  curr_cat <- unique(curr_colData$Category)[i]
  category_means[,i] <- rowMeans(celltypes_mean_LOO[,Category_table[colnames(celltypes_mean_LOO)] == curr_cat,drop=FALSE])
}
colnames(category_means) <- unique(curr_colData$Category)
rownames(category_means) <- rownames(celltypes_mean_LOO)
# Z-transform:
category_means_z <- to_z(category_means)

# Calculate aggregate means for the categories (aggregates of motifs, means across categories)
category_means_aggregates <- matrix(nrow=length(unique(row_annot$`Motif cluster annotation`)), ncol=ncol(category_means))
for(i in 1:nrow(category_means_aggregates)){
  curr_cluster_annot <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))][i]
  category_means_aggregates[i,] <- colSums(category_means[row_annot$`Motif cluster annotation` == curr_cluster_annot,])
}
rownames(category_means_aggregates) <- unique(row_annot$`Motif cluster annotation`)[order(unique(row_annot$`Motif cluster annotation`))]
colnames(category_means_aggregates) <- colnames(category_means)
# Z-transform:
category_means_aggregates_z <- to_z(category_means_aggregates)

# df for motif category annotation
curr_category_annot <- data.frame(row.names=rownames(category_means_aggregates_z), 
                                  amount_motifs = rep("", nrow(category_means_aggregates_z)))
curr_category_annot$amount_motifs <- sapply(rownames(category_means_aggregates_z), 
                                            FUN=function(x){
                                              sum(row_annot$`Motif cluster annotation` == x)
                                            })
colnames(curr_category_annot) <- "Motifs in cluster"

# Get motif activation scores
row_annot$halfmax <- motifs_above_half_max[rownames(row_annot)]
mean_halfmax_df <- row_annot %>% group_by(`Motif cluster annotation`) %>%
  summarize(mean_halfmax=mean(halfmax)) %>% as.data.frame()
mean_halfmax_df$mean_halfmax <- mean_halfmax_df$mean_halfmax * 100
colnames(mean_halfmax_df) <- c("Motif cluster annotation", "Activated promoters (%)")
# Add to motif category annot df
curr_category_annot$`Activated promoters (%)` <- mean_halfmax_df$`Activated promoters (%)`

# Add mean influence score to motif category annot df
sum_loo <- c()
for(i in 1:length(unique(row_annot$`Motif cluster annotation`))){
  curr_motif_fam <- unique(row_annot$`Motif cluster annotation`)[i]
  sum_loo <- c(sum_loo, 
               sum(rowMeans(all_LOO_mat[rownames(row_annot[row_annot$`Motif cluster annotation` == curr_motif_fam,]),])))  
}
names(sum_loo) <- unique(row_annot$`Motif cluster annotation`)
curr_category_annot$`Summed mean influence` <- sum_loo[rownames(curr_category_annot)]

# Fig 2a =====
pheatmap(category_means_aggregates_z, color=heatmap_colors, 
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, 
         filename = paste0(outdir, "Fig2a.pdf"),
         useDingbats=FALSE,
         annotation_row = curr_category_annot)
# Note: the numbers were later added in Illustrator using curr_category_annot:
curr_category_annot


# Do PCA:
graphics.off()
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

# Fig 2c =====
ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=Category, shape=Origin)) + 
  geom_point(size=2) +
  scale_color_stata() + theme_bw(base_size=14) +
  labs(x=paste0("PC1 (", round(variances_expl[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl[2]*100, 2), "%)")) +
  coord_fixed() +
  theme_Nice(angled=FALSE) + theme(legend.position="right")
ggsave(paste0(outdir, "Fig2c.pdf"), width=7, height=7,
       useDingbats=FALSE)


# Load pseudotime values for individual cells

pseudotime_values <- read.csv("data/Kidney_pseudotime_data.tsv", sep="\t")
colnames(pseudotime_values) <- c("UMAP1", "UMAP2", "Pool", "Cell type", "Proximal tubule pseudotime",
                                 "prtub")

# Fig 2d =====
ggplot(pseudotime_values, aes(x=UMAP1, y=UMAP2, color=`Proximal tubule pseudotime`)) + 
  geom_point() + theme_bw(base_size=14) +
  scale_color_gradient(high="#31e03d", low="black", breaks=c(0,.5,1)) +
  theme_Nice(angled = FALSE) + theme(legend.position = "right") + 
  coord_equal()
ggsave(paste0(outdir, "Fig2d.pdf"), height = 6, width=6, useDingbats = FALSE)


# Melt motif leave-one-out scores 
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
  magrittr::set_colnames(c("Motif cluster annotation", "Pool", "Aggregate of motif weights"))
all_LOO_mat_selection_aggregates_melt$Celltype <- curr_colData[all_LOO_mat_selection_aggregates_melt$Pool,1]
all_LOO_mat_selection_aggregates_melt$Category <- Category_table[all_LOO_mat_selection_aggregates_melt$Celltype]
# Order:
all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` <- factor(x=as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`),
                                                                           levels=unique(as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))[order(unique(as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)))])
graphics.off()

# Fig 3a =====

ggplot(all_LOO_mat_selection_aggregates_melt, 
       aes(x=reorder_within(Category,`Aggregate of motif weights`,`Motif cluster annotation`), y=`Aggregate of motif weights`, 
           fill=`Category`)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_x_reordered() +
  scale_fill_stata() + theme_bw(base_size=14) + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x = "Category", y="Aggregate of motif influence scores") +
  theme_Nice() +
  facet_wrap(~`Motif cluster annotation`, scales="free", ncol = 5)

ggsave(paste0(outdir, "/Fig3a.pdf"),
       width=16,height=6, useDingbats=FALSE)



# Load pseudotime values for pools
pool_pseudotime_values <- readRDS("data/20200924_pseudotime_values_fet_kid.RDS")
# Limit to pools where at least 20 cells have pseudotime values
pool_pseudotime_values <- pool_pseudotime_values[pool_pseudotime_values$cells_in_pool >= 20,] 
pseudotime_analysis_pools <- rownames(pool_pseudotime_values)

options(stringsAsFactors = FALSE)
pseudotime_corrs <- c()
for(i in 1:length(unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))){
  curr_cluster_annot <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)[i]
  print(curr_cluster_annot)
  
  curr_melt_aggregates <- all_LOO_mat_selection_aggregates_melt[all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` == curr_cluster_annot,]
  
  curr_melt_aggregates$Pool <- as.character(curr_melt_aggregates$Pool)
  rownames(curr_melt_aggregates) <- curr_melt_aggregates$Pool
  curr_melt_aggregates <- curr_melt_aggregates[pseudotime_analysis_pools,]
  
  pseudotime_corrs <- c(pseudotime_corrs, cor(curr_melt_aggregates$`Aggregate of motif weights`, 
                                              pool_pseudotime_values$pseudotime_values, 
                                              method="spearman"))
  
  cor_df <- data.frame(row.names = rownames(pool_pseudotime_values), 
                       motif_w_agg = curr_melt_aggregates$`Aggregate of motif weights`,
                       pool_pseudotime_values = pool_pseudotime_values$pseudotime_values)
  cor_df$Celltype <- curr_colData[rownames(cor_df),1]
  cor_df$Category <- Category_table[cor_df$Celltype]
  cor_df$cell_type_label <- cor_df$Category
  cor_df$cell_type_label[duplicated(cor_df$cell_type_label)] <- ""
  # q <- ggplot(cor_df, aes(x=motif_w_agg, y=pool_pseudotime_values, color=Celltype)) +
  #   geom_point() +
  #   theme_bw(base_size=14) + xlab("Aggregate of motif weights in pool") + 
  #   ylab(paste0("Mean pseudotime value in pools")) +
  #   labs(title = paste0(curr_cluster_annot,
  #                       " (Spearman R = ", round(pseudotime_corrs[i],
  #                                                digits=2), ")"), 
  #        color="Cell type") +
  #   geom_label_repel(label=cor_df$cell_type_label, show.legend = FALSE) +
  #   scale_color_stata()
  # plot(q)
}
names(pseudotime_corrs) <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)

# Get correlation with proliferation marker gene expression 
proliferation_marker <- colMeans(logcounts(curr_sce[c("MKI67", 
                                                      "PLK1", 
                                                      "E2F1", 
                                                      "FOXM1", 
                                                      "MKI67", 
                                                      "MCM2", 
                                                      "MCM7",
                                                      "BUB1",
                                                      "CCNE1",
                                                      "CCND1",
                                                      "CCNB1",
                                                      "TOP2A"),pseudotime_analysis_pools]))
prolif_corrs <- c()
for(i in 1:length(unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))){
  curr_cluster_annot <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)[i]
  
  curr_melt_aggregates <- all_LOO_mat_selection_aggregates_melt[all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` == curr_cluster_annot,]
  curr_melt_aggregates$Pool <- as.character(curr_melt_aggregates$Pool)
  rownames(curr_melt_aggregates) <- curr_melt_aggregates$Pool
  curr_melt_aggregates <- curr_melt_aggregates[pseudotime_analysis_pools,]
  
  prolif_corrs <- c(prolif_corrs, cor(curr_melt_aggregates$`Aggregate of motif weights`, 
                                      proliferation_marker, method="spearman"))
  
  cor_df <- data.frame(row.names = pseudotime_analysis_pools, 
                       motif_w_agg = curr_melt_aggregates$`Aggregate of motif weights`,
                       proliferation_marker = proliferation_marker)
  cor_df$Celltype <- curr_colData[rownames(cor_df),1]
  cor_df$Category <- Category_table[cor_df$Celltype]
  cor_df$cell_type_label <- cor_df$Category
  cor_df$cell_type_label[duplicated(cor_df$cell_type_label)] <- ""
  # q <- ggplot(cor_df, aes(x=motif_w_agg, y=proliferation_marker, color=Celltype)) +
  #   geom_point() +
  #   theme_bw(base_size=14) + xlab("Aggregate of motif weight across pools") + 
  #   ylab(paste0("Average proliferation marker expression across pools")) +
  #   labs(title = paste0(curr_cluster_annot,
  #                       " (Spearman R = ", round(prolif_corrs[i],
  #                                                digits=2), ")"), 
  #        color="Cell type") +
  #   geom_label_repel(label=cor_df$cell_type_label) +
  #   scale_color_stata()
  # plot(q)
}
names(prolif_corrs) <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)

# Get GO term related to proximal tubule development
GO_terms_2 <- read.csv("data/QuickGO-annotations-1600958410345-20200924.tsv", sep="\t")
GO_terms_2 <- GO_terms_2[GO_terms_2$TAXON.ID == 9606,]
GO_terms_2 <- GO_terms_2$SYMBOL
GO_terms_2 <- unique(GO_terms_2)
GO_terms_2 <- GO_terms_2[GO_terms_2 %in% rownames(curr_sce)]
GO_marker_2 <- colMeans(logcounts(curr_sce[GO_terms_2,pseudotime_analysis_pools]))
GO_corrs_2 <- c()
for(i in 1:length(unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))){
  curr_cluster_annot <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)[i]
  print(curr_cluster_annot)
  
  curr_melt_aggregates <- all_LOO_mat_selection_aggregates_melt[all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` == curr_cluster_annot,]
  
  curr_melt_aggregates$Pool <- as.character(curr_melt_aggregates$Pool)
  rownames(curr_melt_aggregates) <- curr_melt_aggregates$Pool
  curr_melt_aggregates <- curr_melt_aggregates[pseudotime_analysis_pools,]
  
  GO_corrs_2 <- c(GO_corrs_2, cor(curr_melt_aggregates$`Aggregate of motif weights`, 
                                  GO_marker_2, method="spearman"))
  
  cor_df <- data.frame(row.names = pseudotime_analysis_pools, 
                       motif_w_agg = curr_melt_aggregates$`Aggregate of motif weights`,
                       GO_marker_2 = GO_marker_2)
  cor_df$Celltype <- curr_colData[rownames(cor_df),1]
  cor_df$Category <- Category_table[cor_df$Celltype]
  cor_df$cell_type_label <- cor_df$Category
  cor_df$cell_type_label[duplicated(cor_df$cell_type_label)] <- ""
  # q <- ggplot(cor_df, aes(x=motif_w_agg, y=proliferation_marker, color=Celltype)) +
  #   geom_point() +
  #   theme_bw(base_size=14) + xlab("Aggregate of motif weight across pools") + 
  #   ylab(paste0("Average proliferation marker expression across pools")) +
  #   labs(title = paste0(curr_cluster_annot,
  #                       " (Spearman R = ", round(prolif_corrs[i],
  #                                                digits=2), ")"), 
  #        color="Cell type") +
  #   geom_label_repel(label=cor_df$cell_type_label) +
  #   scale_color_stata()
  # plot(q)
}
names(GO_corrs_2) <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)

pseudotime_proliferation_corr_df <- rbind(pseudotime_corrs, prolif_corrs, GO_corrs_2)
rownames(pseudotime_proliferation_corr_df) <- c("Correlation with proximal tubule pseudotime", 
                                                "Correlation with proliferation",
                                                "Correlation with GO:proximal tubule development")

# Fig 2e =====
pheatmap::pheatmap(pseudotime_proliferation_corr_df, cluster_rows = FALSE, 
                   color=heatmap_colors, border_color = NA, cellwidth = 12, cellheight = 12,
                   angle_col=45, filename = paste0(outdir, "/Fig2e.pdf"), width=7, height=3)
graphics.off()


# Correlations of TF family influence scores with individual genes
FUBP1_agg_across_pools <- colSums(all_LOO_mat_selection[row_annot$`Motif cluster annotation` == "FUBP1",])
FUBP1_exp_across_pools <- logcounts(curr_sce)["FUBP1",]
FUBP1_corr_df <- data.frame(agg_LOO = FUBP1_agg_across_pools, 
                            exp=FUBP1_exp_across_pools, 
                            cat=curr_colData$Category)

ETS2_agg_across_pools <- colSums(all_LOO_mat_selection[row_annot$`Motif cluster annotation` == "ETS motif family 2",])
ELF1_exp_across_pools <- logcounts(curr_sce)["ELF1",]
OSR1_exp_across_pools <- logcounts(curr_sce)["OSR1",]
OSR2_exp_across_pools <- logcounts(curr_sce)["OSR2",]
ETV4_exp_across_pools <- logcounts(curr_sce)["ETV4",]
ETV5_exp_across_pools <- logcounts(curr_sce)["ETV5",]
ETS2_corr_df <- data.frame(agg_LOO = ETS2_agg_across_pools, 
                           exp=ELF1_exp_across_pools, 
                           cat=curr_colData$Category)
OSR1_corr_df <- data.frame(agg_LOO = ETS2_agg_across_pools, 
                           exp=OSR1_exp_across_pools, 
                           cat=curr_colData$Category)
OSR2_corr_df <- data.frame(agg_LOO = ETS2_agg_across_pools, 
                           exp=OSR2_exp_across_pools, 
                           cat=curr_colData$Category)
ETV4_corr_df <- data.frame(agg_LOO = ETS2_agg_across_pools, 
                           exp=ETV4_exp_across_pools, 
                           cat=curr_colData$Category)
ETV5_corr_df <- data.frame(agg_LOO = ETS2_agg_across_pools, 
                           exp=ETV5_exp_across_pools, 
                           cat=curr_colData$Category)

FUBP1_plot <- ggplot(FUBP1_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="FUBP1", 
       y="FUBP1 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0, 0.02, 0.04)) + theme(aspect.ratio=1)
ETS2_plot <- ggplot(ETS2_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE)+
  labs(color="Category", 
       x="ETS motif family 2", 
       y="ELF1 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0.25,0.35,0.45)) + theme(aspect.ratio=1)
OSR1_plot <- ggplot(OSR1_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) + 
  labs(color="Category", 
       x="ETS motif family 2", 
       y="OSR1 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0.25,0.35,0.45)) + theme(aspect.ratio=1)
OSR2_plot <- ggplot(OSR2_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="ETS motif family 2", 
       y="OSR2 expression") + 
  scale_color_stata()+ scale_x_continuous(breaks=c(0.25,0.35,0.45)) + theme(aspect.ratio=1)
ETV4_plot <- ggplot(ETV4_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="ETS motif family 2", 
       y="ETV4 expression") + 
  scale_color_stata()+ scale_x_continuous(breaks=c(0.25,0.35,0.45)) + theme(aspect.ratio=1)
ETV5_plot <- ggplot(ETV5_corr_df, aes(x=agg_LOO, y=exp, color=cat)) +
  geom_point() + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="ETS motif family 2", 
       y="ETV5 expression") + 
  scale_color_stata()+ scale_x_continuous(breaks=c(0.25,0.35,0.45)) + theme(aspect.ratio=1) + theme(legend.position="right")

# Fig 3b =====
(FUBP1_plot | ETS2_plot | OSR1_plot) / (OSR2_plot | ETV4_plot | ETV5_plot)
ggsave(paste0(outdir, "/Fig3b.pdf"), width=8, height=4, useDingbats=FALSE)
# Correlations were obtained using e.g. cor(ETV5_corr_df$agg_LOO, ETV5_corr_df$exp, method="spearman")



# Get all correlations between TFs and motif family weights
repr_cl_names <- unique(row_annot$`Motif cluster annotation`)
pool_exps <- list()
cluster_TFs <- list()
all_clusters_corr_dfs <- list()
all_cluster_corrs <- list()
all_cluster_pvals <- list()
all_cluster_mean_exps <- list()
for(i in 1:length(repr_cl_names)){
  curr_cluster <- repr_cl_names[i]
  pool_exps[[curr_cluster]] <- colMeans(all_LOO_mat_selection[row_annot$`Motif cluster annotation` == curr_cluster,])
  curr_motifs <- rownames(row_annot)[row_annot$`Motif cluster annotation` == curr_cluster]
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
        curr_corr_df$Category <- Category_table[curr_sce$Celltype]
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
all_cluster_corrs_df <- stack(all_cluster_corrs)
all_cluster_pvals_fdr <- lapply(all_cluster_pvals, FUN=p.adjust, method="fdr")
all_pvals_fdr_df <- melt(all_cluster_pvals_fdr) %>% magrittr::set_colnames(c("pval_fdr", "cluster"))
all_cluster_mean_exps_df$pval_fdr <- all_pvals_fdr_df$pval_fd
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


# Fig 3c =====
ggplot(all_cluster_mean_exps_df, aes(x=cluster_annot, y=corr, color=expression_2)) +
  geom_jitter(width = 0) + 
  geom_label_repel(label=all_cluster_mean_exps_df$tf_labels, size=4) +
  theme_bw(base_size=14) + 
  theme_Nice() + theme(legend.position = "right") + 
  labs(x="Motif cluster name", 
       y="Spearman R", color="Mean TF expression across pools") 
ggsave(filename=paste0(outdir, "/Fig3c.pdf"), 
       width = 12, height=6, useDingbats=FALSE)
