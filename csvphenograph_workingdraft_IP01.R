rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()

metaDataFile = paste0(work,"/Config/metadata_IP01.xlsx")
panelDataFile = paste0(work,"/Config/panel.xlsx")
dataDirectory = paste0(work,"/Data")

require(scales);require(readxl);require(plyr);require(dplyr);require(DataEditR); 
require(Rphenograph);require(Hmisc); require(ComplexHeatmap); require(pals); require(matrixStats);
require(reshape2); require(ggplot2)

output<-readRDS("global_data_IP01.RDS")
data_full <- data.frame(output[1])
data <- data.matrix(output[2])
data01 <- output[3]
csv_full <- output[4]

## Read-in metadata and clean =======
ifelse(grepl(metaDataFile,pattern='.xlsx'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
md$file_name <- factor(md$file_name)
md$File_order <- factor(md$File_order)
md$Tumor <- factor(md$Tumor)

##input image id into metadata
image_id<-c()
for (i in 1:length(md$file_name)){
  tempfile <- read.csv(paste0(dataDirectory,"/",md$file_name[i]))
  df<- as.data.frame(cbind(paste0(md$file_name[i]), unique(tempfile$ImageId)))
  image_id<-rbind(image_id,df)
}
md$ImageId <- image_id$V2[match(image_id$V1,md$file_name)]


## Make sure all files in metadata present in datadirectory
if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])){
  print(paste('ERR: not all filenames in metadata present in data folder - missing',
              md$file_name[!which(data$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),
                                                                                   pattern = '.csv')])],'Subsetting...'))
  md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])),]
}

#Set up levels ======
samplevels=c("01_C1","01_C2","01_C3","01_I3","01_I2","01_I1")


tumorlevels=c("IPMN",
              "Transition",
              "Tumor")

## Read csv into csv_raw =========
csv_raw <- lapply(paste0(dataDirectory,"/",md$file_name),read.csv)
csv_raw_full <- plyr::ldply(csv_raw, rbind)
csv_raw_full$ImageId <- md$sample_id[match(csv_raw_full$ImageId,md$ImageId)]

#export to clean panel + rename ImageId as sample_id
rawcolnames <- c()
rawcolnames$name <- colnames(csv_raw_full)
rawcolnames$sum <- colSums(csv_raw_full)
write.csv(rawcolnames, 'rawpanel.csv')

#clean csv_raw_full dataframe to csv_full containing analysis markers only
cleanpanel <- read_xlsx('cleanpanel.xlsx')
colnames(csv_raw_full) <- cleanpanel$clean_names
panel <- cleanpanel$clean_names[cleanpanel$analysis > 0]
csv_full <- csv_raw_full[,colnames(csv_raw_full) %in% panel]



#sort panels into different categories
subtype_markers <- cleanpanel$clean_names[cleanpanel$subtype == 1]
functional_markers <- cleanpanel$clean_names[cleanpanel$functional == 1]
otherparameters <- cleanpanel$clean_names[cleanpanel$other ==1]
cluster_by <- cleanpanel$clean_names[cleanpanel$cluster_by == 1]
Tcellmarkers <- cleanpanel$clean_names[cleanpanel$tcell == 1]
Macmarkers <- cleanpanel$clean_names[cleanpanel$mac == 1]
NKmarkers <- cleanpanel$clean_names[cleanpanel$nk == 1]
Stromamarkers <- cleanpanel$clean_names[cleanpanel$stroma == 1]


#Cluster heatmap for unannotated clusters======
data_full <- csv_full
data <- data.matrix(csv_full[,-1])
data <- asinh(data[, union(subtype_markers,functional_markers)] / 0.8)

#phenograph clustering of data
rng <- colQuantiles(data, probs = c(0.05, 0.95))
data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01[data01 < 0] <- 0; data01[data01 > 1] <- 1;data01 <-data01[,union(subtype_markers,functional_markers)]

set.seed(1234)
phenographout<-Rphenograph(data01)
data_full$cluster<-factor(membership(phenographout[[2]]))

cluster_mean <- data.frame(data01, cluster = data_full$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_mat<-as.matrix(cluster_mean[,union(subtype_markers,functional_markers)])

rownames(cluster_mean_mat)<-1:nrow(cluster_mean_mat)

cluster_scaled<-t(scale(t(cluster_mean_mat)))

rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#save as RDS file
global_data <- list(data_full, data, data01, csv_full)
saveRDS(global_data, "global_data_IP01.RDS")

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


pdf("clusterheatmap_IP01.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="IP01 Phenograph Clusters",
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

clusternames<-sort(unique(cluster_merging$new_cluster))

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))

names(colorassigned)<-clusternames

color_list = list(clusters=colorassigned)

color_list_byoriginal = colorassigned[match(unique(cluster_merging$new_cluster),names(colorassigned))]

cp<-rowAnnotation(clusters=clusternames,
                  col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_full$cluster1m)),
                    gp = gpar(fill=colorassigned),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

pdf("clusterheatmap_IP01merged.pdf",width=10,height=4)
Heatmap(cluster_scaled_merged,
        column_title="IP01 Phenograph Merged Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = F,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled_merged)*unit(4, "mm"), 
        height = nrow(cluster_scaled_merged)*unit(4, "mm"))
dev.off()



#T cell sub clusters=======
data_immune <- data_full[data_full$cluster %in% c(4, 10, 19),]
data_immune_expr <- data.matrix(asinh(data_immune[,Tcellmarkers]/0.8))

rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,Tcellmarkers]


set.seed(1234)
phenographout_immune<-Rphenograph(data_immune_expr01)
data_immune$cluster<-factor(membership(phenographout_immune[[2]]))
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,Tcellmarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#save as RDS file
tcell_data <- list(data_immune, data_immune_expr, data_immune_expr01)
saveRDS(tcell_data, "tcell_data.RDS")
readRDS("tcell_data.RDS")

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


pdf("clusterheatmap_IP01Tcells.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="IP01 Phenograph T Clusters",
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



####Abundance plots======
counts_table<-table(data_full$cluster, data_full$sample_id)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <-as.data.frame.matrix(props_table)

#Densities
areas <- read_xlsx(paste0(work,'/Config/IP01_area.xlsx'))
densities <- t(t(counts)/areas$TotalArea)

write.csv(counts,'counts.csv')
write.csv(props,'props.csv')
write.csv(densities, 'densities.csv')

ggdf <- melt(data.frame(cluster = rownames(props), props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$tumor <- factor(md$Tumor[match(ggdf$sample_id,md$sample_id)], levels=tumorlevels)

#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities), densities, check.names = FALSE),
             id.vars = "cluster", value.name = "densities", 
             variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdfd$tumor <- factor(md$Tumor[match(ggdf$sample_id,md$sample_id)], levels=tumorlevels)

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
pdf('Density_box_unannotated.pdf',width=12,height=14)
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
ggdf$tumor <- factor(md$Tumor[match(ggdf$sample_id,md$sample_id)], levels=tumorlevels)

#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities_m), densities_m, check.names = FALSE),
              id.vars = "cluster", value.name = "densities_m", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdfd$tumor <- factor(md$Tumor[match(ggdf$sample_id,md$sample_id)], levels=tumorlevels)


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
pdf('Abundance_box_merged.pdf',width=6,height=5)
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
pdf('Density_box_merged.pdf',width=6,height=5)
ggp2
dev.off()

