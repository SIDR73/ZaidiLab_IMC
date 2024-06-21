rm(list = ls())

####REQUIRED LIBRARIES####
library(reshape2)
library(pals)
library(ggplot2)
library(Hmisc)
library(devtools)
library(ComplexHeatmap)
library(ggiraphExtra)
library(diffcyt)
# library(ggvoronoi)
library(ggpubr)
library(multcomp)
library(sf)
library(clusterSim)
library(circlize)
library(RColorBrewer)
library(stringr)
library(igraph)
library(readxl)
library(dplyr)
library(packcircles)
library(gridExtra)
library(limma)
library(qgraph)
library(flowCore)
library(pheatmap)
library(matrixStats)
library(ggbreak)
library(readxl)

setwd("C:/Users/sidha/OneDrive/Sid Stuff/PROJECTS/Jaffee Lab/tiKPC_IMC_Analysis/Phenograph/Phenograph")

work<-getwd()

metaDataFile = paste0(work,"/Config/metadata_tikpc.xlsx")
panelDataFile = paste0(work,"/Config/tikpc_cleanpanel.xlsx")
dataDirectory = paste0(work,"/Data")

require(scales);require(readxl);require(plyr);require(dplyr);require(DataEditR); 
require(Rphenograph);require(Hmisc); require(ComplexHeatmap); require(pals); require(matrixStats);
require(reshape2); require(ggplot2)

# Run this only when initial RDS is generated
output<-readRDS("global_data_tikpc.RDS")
data_full <- data.frame(output[1])
data <- data.matrix(output[2])
data01 <- output[3]
csv_full <- output[4]

## Read-in metadata and clean =======
ifelse(grepl(metaDataFile,pattern='.xlsx'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
md$file_name <- factor(md$file_name)
md$image_name <- factor(md$image_name)
md$sample_id <- factor(md$sample_id)
md$File_order <- factor(md$File_order)
md$Type <- factor(md$Type)
md$Mouse <- factor(md$Mouse)

# ##input image id into metadata
# image_id<-c()
# for (i in 1:length(md$file_name)){
#   tempfile <- read.csv(paste0(dataDirectory,"/",md$file_name[i]))
#   df<- as.data.frame(cbind(paste0(md$file_name[i]), unique(tempfile$ImageId)))
#   image_id<-rbind(image_id,df)
# }
# # From Melissa - I changed initial image_id$V2 to image_id$V1 because 
# # there was otherwise an empty vector.. even doing this didn't get it working
# md$ImageId <- image_id$V1[match(image_id$V1,md$file_name)]


# Error Handling: Check that all files in metadata present in dataDirectory
if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])){
  print(paste('ERR: not all filenames in metadata present in data folder - missing',
              md$file_name[!which(data$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),
                                                                                   pattern = '.csv')])],'Subsetting...'))
  md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])),]
}

#Set up levels ======
samplevels=c("7986_panin_1",
             "7986_norm_1",
             "8084_panin_1",
             "8084_norm_1",
             "8084_norm_2",
             "8084_norm_3",
             "8084_norm_4",
             "8084_norm_5",
             "8084_norm_6",
             "8084_norm_7",
             "7999_tc_1",
             "7999_te_1",
             "7999_te_2",
             "7999_tc_2",
             "7999_la_1",
             "7999_la_2",
             "7999_tc_3",
             "7999_te_3",
             "7999_tc_4",
             "7999_tc_5",
             "7999_te_4",
             "7999_tc_6",
             "7811_tc_1",
             "7811_te_1",
             "7811_tc_2",
             "7811_te_2",
             "7811_tc_3",
             "7811_norm_1",
             "8019_norm_1",
             "8019_panin_1",
             "8019_la_1",
             "8019_panin_2",
             "8016_panin_1",
             "8016_panin_2",
             "7731_la_1",
             "7731_la_2",
             "7731_panin_1",
             "7813_te_1",
             "7813_tc_1",
             "7813_te_2",
             "7813_te_3",
             "7813_tc_2",
             "7785_la_1",
             "7785_panin_1",
             "7785_norm_1",
             "8082_te_1",
             "8082_la_1",
             "8082_la_2",
             "8082_te_2",
             "8082_tc_1",
             "8082_te_3",
             "8180_la_1",
             "8180_tc_1",
             "8180_te_1",
             "8180_te_2",
             "8180_te_3",
             "8180_tc_2",
             "8180_te_4",
             "8180_tc_3",
             "8180_te_5",
             "8180_tc_4")

tcellimagelevels=c("7986_ROI_001_P",
                   "8084_ROI_001_P",
                   "8084_ROI_002_N1",
                   "8084_ROI_002_N2",
                   "8084_ROI_003_N1",
                   "8084_ROI_003_N4",
                   "7999_ROI_001_TC",
                   "7999_ROI_001_TE",
                   "7999_ROI_002_TE",
                   "7999_ROI_003_TC",
                   "7999_ROI_005_LA1",
                   "7999_ROI_005_LA2",
                   "7999_ROI_005_TC",
                   "7999_ROI_005_TE",
                   "7999_ROI_006_TC",
                   "7999_ROI_007_TC",
                   "7999_ROI_008_TE",
                   "7999_ROI_009_TC",
                   "7811_ROI_001_split_TC",
                   "7811_ROI_001_split_TE",
                   "7811_ROI_001_TC",
                   "7811_ROI_001_TE",
                   "7811_ROI_002_TC",
                   "7811_ROI_003_N",
                   "8019_ROI_001_N",
                   "8019_ROI_003_P",
                   "8019_ROI_004_LA",
                   "8019_ROI_004_P",
                   "8016_ROI_001_P",
                   "8016_ROI_002_P",
                   "7731_ROI_001_LA1",
                   "7731_ROI_001_LA2",
                   "7731_ROI_001_P",
                   "7813_ROI_001_TE",
                   "7813_ROI_002_TC",
                   "7813_ROI_002_TE",
                   "7813_ROI_003_TE",
                   "7813_ROI_005_TC",
                   "7785_ROI_001_LA",
                   "7785_ROI_001_P",
                   "7785_ROI_004_N",
                   "8082_ROI_001_TE",
                   "8082_ROI_002_LA1",
                   "8082_ROI_002_LA2",
                   "8082_ROI_002_TE",
                   "8082_ROI_003_TC",
                   "8082_ROI_003_TE",
                   "8180_ROI_001_LA",
                   "8180_ROI_001_TC",
                   "8180_ROI_001_TE",
                   "8180_ROI_002_TE",
                   "8180_ROI_003_TE",
                   "8180_ROI_004_TC",
                   "8180_ROI_004_TE",
                   "8180_ROI_005_TC",
                   "8180_ROI_005_TE",
                   "8180_ROI_006_TC")

namelevels=c("7986_ROI_001_P",
             "7986_ROI_4_N",
             "8084_ROI_001_P",
             "8084_ROI_002_N1",
             "8084_ROI_002_N2",
             "8084_ROI_003_N1",
             "8084_ROI_003_N2",
             "8084_ROI_003_N3",
             "8084_ROI_003_N4",
             "8084_ROI_003_N5",
             "7999_ROI_001_TC",
             "7999_ROI_001_TE",
             "7999_ROI_002_TE",
             "7999_ROI_003_TC",
             "7999_ROI_005_LA1",
             "7999_ROI_005_LA2",
             "7999_ROI_005_TC",
             "7999_ROI_005_TE",
             "7999_ROI_006_TC",
             "7999_ROI_007_TC",
             "7999_ROI_008_TE",
             "7999_ROI_009_TC",
             "7811_ROI_001_split_TC",
             "7811_ROI_001_split_TE",
             "7811_ROI_001_TC",
             "7811_ROI_001_TE",
             "7811_ROI_002_TC",
             "7811_ROI_003_N",
             "8019_ROI_001_N",
             "8019_ROI_003_P",
             "8019_ROI_004_LA",
             "8019_ROI_004_P",
             "8016_ROI_001_P",
             "8016_ROI_002_P",
             "7731_ROI_001_LA1",
             "7731_ROI_001_LA2",
             "7731_ROI_001_P",
             "7813_ROI_001_TE",
             "7813_ROI_002_TC",
             "7813_ROI_002_TE",
             "7813_ROI_003_TE",
             "7813_ROI_005_TC",
             "7785_ROI_001_LA",
             "7785_ROI_001_P",
             "7785_ROI_004_N",
             "8082_ROI_001_TE",
             "8082_ROI_002_LA1",
             "8082_ROI_002_LA2",
             "8082_ROI_002_TE",
             "8082_ROI_003_TC",
             "8082_ROI_003_TE",
             "8180_ROI_001_LA",
             "8180_ROI_001_TC",
             "8180_ROI_001_TE",
             "8180_ROI_002_TE",
             "8180_ROI_003_TE",
             "8180_ROI_004_TC",
             "8180_ROI_004_TE",
             "8180_ROI_005_TC",
             "8180_ROI_005_TE",
             "8180_ROI_006_TC")

typelevels=c("Normal",
              "PanIN",
              "Tumor Edge",
             "Tumor Core",
             "Lymphoid Aggregate")

mouselevels=c("1",
              "2",
              "3",
              "4",
              "5",
              "6",
              "7",
              "8",
              "9",
              "10",
              "11")

## Read csv into csv_raw =========
csv_raw <- lapply(paste0(dataDirectory,"/",md$file_name),read.csv)
csv_raw_full <- plyr::ldply(csv_raw, rbind)

#export to clean panel + rename Name as sample_id
rawcolnames <- c()
rawcolnames$name <- colnames(csv_raw_full)
write.csv(rawcolnames, 'rawpanel.csv') # MADE AN EDIT

#clean csv_raw_full dataframe to csv_full containing analysis markers only
cleanpanel <- read_xlsx('tikpc_cleanpanel.xlsx')
colnames(csv_raw_full) <- cleanpanel$clean_names
panel <- cleanpanel$clean_names[cleanpanel$analysis > 0]
csv_full <- csv_raw_full[,colnames(csv_raw_full) %in% panel]

#sort panels into different categories BINARY CLASSIFICATION
subtype_markers <- cleanpanel$clean_names[cleanpanel$subtype == 1]
functional_markers <- cleanpanel$clean_names[cleanpanel$functional == 1]
otherparameters <- cleanpanel$clean_names[cleanpanel$other ==1]
cluster_by <- cleanpanel$clean_names[cleanpanel$cluster_by == 1]
tcellmarkers <- cleanpanel$clean_names[cleanpanel$tcell == 1]
macmarkers <- cleanpanel$clean_names[cleanpanel$mac == 1]
dendmarkers <- cleanpanel$clean_names[cleanpanel$dendritic == 1]
stromamarkers <- cleanpanel$clean_names[cleanpanel$stroma == 1]
bcellmarkers <- cleanpanel$clean_names[cleanpanel$bcell == 1]
neutmarkers <- cleanpanel$clean_names[cleanpanel$neutrophil == 1]

#can call other phenotypes 


#Cluster heatmap for unannotated clusters======
data_full <- csv_full
data <- data.matrix(csv_full[,-1])
data <- asinh(data[, union(subtype_markers,functional_markers)] / 0.8)

#phenograph clustering of data
rng <- colQuantiles(data, probs = c(0.05, 0.95))
data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01[data01 < 0] <- 0; data01[data01 > 1] <- 1;data01 <-data01[,union(subtype_markers,functional_markers)]
data01[is.na(data01)] <- 0

set.seed(1234)
phenographout<-Rphenograph(data01)
data_full$cluster<-factor(membership(phenographout[[2]])) # Store clusters here

cluster_mean <- data.frame(data01, cluster = data_full$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_mat<-as.matrix(cluster_mean[,union(subtype_markers,functional_markers)])

rownames(cluster_mean_mat)<-1:nrow(cluster_mean_mat)

cluster_scaled<-t(scale(t(cluster_mean_mat)))

rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#save as RDS file
global_data <- list(data_full, data, data01, csv_full)
saveRDS(global_data, "global_data_tikpc.RDS")

## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_full$cluster)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))


pdf("clusterheatmap_tikpc_SID1.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="tiKPC Phenograph Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 


#cluster heatmap for merged annotations ===========
clusterMergeFile = paste0(work,"/Config/merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("B",
                "Epith",
                "Mac",
                "CD57+",
                "Stromal",
                "T")


colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(data_full$cluster, cluster_merging$original_cluster)
data_full$cluster1m <- cluster_merging$new_cluster[mm1]

cluster_mean_merged <- data.frame(data01, cluster = data_full$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,union(subtype_markers,functional_markers)])

cluster_scaled_merged<-t(scale(t(cluster_mean_merged_mat)))

rownames(cluster_scaled_merged)<-1:nrow(cluster_scaled_merged)

## Annotation for the merged clusters

## ERROR

if(!is.null(clusterMergeFile)){
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
  annotation_row$Merged <- cluster_merging$new_cluster
  color_clusters2 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Merged))
  names(color_clusters2) <- levels(cluster_merging$new_cluster)
  annotation_colors$Merged <- color_clusters2
}

## Colors for the heatmap

legend_breaks = seq(from = 0, to = 1, by = 0.2)

clusternames<-clusterlevels

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusternames))

names(colorassigned)<-clusternames

rownames(cluster_scaled_merged) <- cluster_mean_merged$cluster

color_list = list(clusters=colorassigned)

color_list_byoriginal = colorassigned[match(unique(cluster_merging$new_cluster),names(colorassigned))]

cp<-rowAnnotation(#clusters=clusternames,
                  col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_full$cluster1m)),
                    gp = gpar(fill=colorassigned),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

pdf("clusterheatmap_IP01merged_SID1.pdf",width=10,height=4)
Heatmap(cluster_scaled_merged[clusterlevels,],
        column_title="IP01 Phenograph Merged Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled_merged)), to = round(max(cluster_scaled_merged)))),
        width = ncol(cluster_scaled_merged)*unit(4, "mm"), 
        height = nrow(cluster_scaled_merged)*unit(4, "mm"))
dev.off()



#T cell sub clusters=======
# data_immune <- data_full[data_full$cluster %in% c(4, 10, 19),]
data_immune <- data_full[data_full$cluster %in% c(4, 19),]
data_immune_expr <- data.matrix(asinh(data_immune[,tcellmarkers]/0.8)) # T to t

rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,tcellmarkers]

# Sid's Implementation for handling NaN values:
data_immune_expr01[is.na(data_immune_expr01)] <- 0    # Replace NA values with 0
data_immune_expr01[is.nan(data_immune_expr01)] <- 0   # Replace NaN values with 0
data_immune_expr01[is.infinite(data_immune_expr01)] <- 0  # Replace Inf values with 0


set.seed(1234)
phenographout_immune<-Rphenograph(data_immune_expr01)
data_immune$cluster<-factor(membership(phenographout_immune[[2]]))
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,tcellmarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#save as RDS file
tcell_data <- list(data_immune, data_immune_expr, data_immune_expr01)
saveRDS(tcell_data, "tcell_data_SID.RDS")
readRDS("tcell_data_SID.RDS")

## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune$cluster)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75,
                       width = unit(2,"cm")))


pdf("clusterheatmap_TRegs_SID2.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="Phenograph T Regulatory Cell Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"),
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off()

### MERGED T Cell Plots

clusterMergeFile = paste0(work,"/Config/merge_tcells.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("CD3+ CD4+ TCF1+",
                "CD8+ TCF1+ GZMB+",
                "CD3+ CD8+ GZMB+",
                "CD3+ TCF1+",
                "GZMB+ PD1+",
                "CD3+ FOXP3+ TOX+",
                "CD3+ BCL6+ PD1+",
                "CD3+ GZMB+ PD1+",
                "CD3+ TOX+ CD11c+",
                "CD11c+ GZMB+",
                "CD3+ CD4+ FOXP3+ TOX+ PD1+")

clusterlevels <- factor(clusterlevels, levels = unique(data_immune$cluster1m))

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(data_immune$cluster, cluster_merging$original_cluster)
data_immune$cluster1m <- cluster_merging$new_cluster[mm1]

cluster_mean_merged <- data.frame(data_immune_expr01, cluster = data_immune$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

# Create matrix w T cell markers for heatmap
cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,tcellmarkers])

cluster_scaled_merged<-t(scale(t(cluster_mean_merged_mat)))

rownames(cluster_scaled_merged)<-1:nrow(cluster_scaled_merged)

## Annotation for the merged clusters

if(!is.null(clusterMergeFile)){
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
  
  # Sid Edit
  # annotation_row <- data.frame(Cluster = data_immune$cluster1m)
  # annotation_row <- annotation_row[!duplicated(annotation_row$Cluster), ]
  
  annotation_row$Merged <- cluster_merging$new_cluster
  color_clusters2 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Merged))
  names(color_clusters2) <- levels(cluster_merging$new_cluster)
  annotation_colors$Merged <- color_clusters2
}

## Colors for the heatmap

legend_breaks = seq(from = 0, to = 1, by = 0.2)
####
clusternames<-clusterlevels
####
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusternames))

names(colorassigned)<-clusternames
####
rownames(cluster_scaled_merged)<-cluster_mean_merged$cluster

color_list = list(clusters=colorassigned)

color_list_byoriginal = colorassigned[match(unique(cluster_merging$new_cluster),names(colorassigned))]

cp<-rowAnnotation(clusters=clusternames,
  col=color_list,
  gp = gpar(col = "white", lwd = .5),
  counts= anno_barplot(
    as.vector(table(data_immune$cluster1m)),
    gp = gpar(fill=colorassigned),
    border = F,
    bar_width = 0.75, 
    width = unit(2,"cm")))


rownames(cluster_scaled_merged) <- cluster_mean_merged$cluster

pdf("clusterheatmap_merged_TReg_SID2.pdf",width=10,height=4)
Heatmap(cluster_scaled_merged[clusterlevels,],
        column_title="tiKPC Phenograph Merged T_Reg Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled_merged)), to = round(max(cluster_scaled_merged)))),
        width = ncol(cluster_scaled_merged)*unit(4, "mm"),
        height = nrow(cluster_scaled_merged)*unit(4, "mm"))
dev.off()




####Abundance plots======
tumorlevels=c("PanIN",
                "Normal",
                "Tumor Core",
                "Tumor Edge",
                "Lymphoid Aggregate")
counts_table<-table(data_full$cluster, data_full$sample_id)
# counts_table<-table(annotation_row$cluster, data_full$sample_id)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <-as.data.frame.matrix(props_table)

#Densities (make sure area csv follows the order in counts table)
areas <- read_xlsx(paste0(work,'/Config/Areas.xlsx')) ## Modified from IP01_area.xlsx
densities <- t(t(counts)/areas$TotalArea) ## ISSUE: THE LENGTHS ARE DIFFERENT SO SHORTER LIST CYCLES BACK

write.csv(counts,'counts.csv')
write.csv(props,'props.csv')
write.csv(densities, 'densities.csv')

ggdf <- melt(data.frame(cluster = rownames(props), props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$tumor <- factor(md$Type[match(ggdf$sample_id,md$sample_id)], levels=tumorlevels)

#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities), densities, check.names = FALSE),
             id.vars = "cluster", value.name = "densities", 
             variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdfd$tumor <- factor(md$Type[match(ggdf$sample_id,md$sample_id)], levels=tumorlevels)

#plot box plots
#% CELLS
ggp2<-ggplot(ggdfd,aes(x=tumor,y=densities,fill=tumor))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=6,scales="free")+
  ylab("Density (# of cells)")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11)
  )
pdf('Density_box_unannotated_SID1.pdf',width=12,height=14)
ggp2
dev.off()

####Abundance plots======
counts_table_m<-table(data_full$cluster1m, data_full$sample_id)
props_table_m <- t(t(counts_table_m) / colSums(counts_table_m)) * 100
counts_m <- as.data.frame.matrix(counts_table_m)
props_m <-as.data.frame.matrix(props_table_m)

densities_m <- t(t(counts_m)/areas$TotalArea)


ggdf <- melt(data.frame(cluster = rownames(props_m), props_m, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$tumor <- factor(md$Type[match(ggdf$sample_id,md$sample_id)], levels=tumorlevels)
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)

#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities_m), densities_m, check.names = FALSE),
              id.vars = "cluster", value.name = "densities_m", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdfd$tumor <- factor(md$Type[match(ggdf$sample_id,md$sample_id)], levels=tumorlevels)
ggdfd$cluster <- factor(ggdfd$cluster, levels=clusterlevels)

#plot box plots
#% CELLS
ggp2<-ggplot(ggdf,aes(x=tumor,y=proportion,fill=tumor))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=4,scales="free")+
  ylab("Abundance (% of total cells)")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11)
  )
pdf('Abundance_box_merged_SID1.pdf',width=6,height=5)
ggp2
dev.off()

ggp2<-ggplot(ggdfd,aes(x=tumor,y=densities_m,fill=tumor))+
  geom_boxplot(outlier.shape=NA, lwd=0.25)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=4,scales="free")+
  ylab("Density (# of cells)")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11)
  )
pdf('Density_box_merged_SID1.pdf',width=6,height=5)
ggp2
dev.off()
