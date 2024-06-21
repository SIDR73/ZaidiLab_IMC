#T cell sub clusters=======
data_immune <- data_full[data_full$cluster %in% c(19,4,16,8,30,5),]
data_immune_expr <- data.matrix(asinh(data_immune[,tcellmarkers]/0.8))

rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,tcellmarkers]
data_immune_expr01[is.na(data_immune_expr01)] <- 0


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
saveRDS(tcell_data, "tcell_data.RDS")
tcell_data<-readRDS("tcell_data.RDS")
data_immune <- data.frame(tcell_data[1])
data_immune_expr <- data.frame(tcell_data[2])
data_immune_expr01 <- data.frame(tcell_data[3])


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


pdf("clusterheatmap_tikpc_tcellsSID2.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="tiKPC Phenograph T Clusters",
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
                "CD3+ CD4+ FOXP3+ TOX+ PD1+",
                "CD8+ CD11c+",
                "CD3+")



colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(data_immune$cluster, cluster_merging$original_cluster)
data_immune$cluster1m <- cluster_merging$new_cluster[mm1]

cluster_mean_merged <- data.frame(data_immune_expr01, cluster = data_immune$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,tcellmarkers])

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
####
clusternames<-clusterlevels
####
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusternames))

names(colorassigned)<-clusternames
####
rownames(cluster_scaled_merged)<-cluster_mean_merged$cluster

color_list = list(clusters=colorassigned)

color_list_byoriginal = colorassigned[match(unique(cluster_merging$new_cluster),names(colorassigned))]

cp<-rowAnnotation(#clusters=clusternames,
                  col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_immune$cluster1m)),
                    gp = gpar(fill=colorassigned),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

pdf("clusterheatmap_merged_tcellSID2.pdf",width=10,height=4)
Heatmap(cluster_scaled_merged[clusterlevels,],
        column_title="tiKPC Phenograph Merged T Clusters",
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
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled_merged)*unit(4, "mm"), 
        height = nrow(cluster_scaled_merged)*unit(4, "mm"))
dev.off()


####Abundance plots======
counts_table<-table(data_immune$cluster, data_immune$sample_id)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <-as.data.frame.matrix(props_table)

#Densities (manually re-order Areas table to match the counts table order)
# I had to make a separate Tcell area file because not every sample in the original table has T cells
areas <- read_xlsx(paste0(work,'/Config/tcellareas.xlsx'))
densities <- t(t(counts)/areas$TotalArea)

mdtcells <- md[md$image_name %in% unique(colnames(props)),]

ggdf <- melt(data.frame(cluster = rownames(props), props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=tcellimagelevels)
ggdf$type <- factor(mdtcells$Type[match(ggdf$sample_id,mdtcells$image_name)], levels=typelevels)
ggdf$mouse <- factor(mdtcells$Mouse[match(ggdf$sample_id,mdtcells$image_name)], levels=mouselevels)



#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities), densities, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdf$sample_id, levels=tcellimagelevels)
ggdfd$type <- factor(mdtcells$Type[match(ggdf$sample_id,mdtcells$image_name)], levels=typelevels)
ggdfd$mouse <- factor(mdtcells$Mouse[match(ggdf$sample_id,mdtcells$image_name)], levels=mouselevels)

#plot box plots
#% CELLS
ggp2<-ggplot(ggdf,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=6,scales="free")+
  ylab("Abundance (% of cells)")+
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
pdf('Abundance_box_tcellSID2.pdf',width=8,height=4)
ggp2
dev.off()

#Densities
ggp2<-ggplot(ggdfd,aes(x=type,y=densities,fill=type))+
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
pdf('Density_box_tcell_unannotatedSID2.pdf',width=8,height=4)
ggp2
dev.off()


####Abundance plots======
counts_table_m<-table(data_immune$cluster1m, data_immune$sample_id)
props_table_m <- t(t(counts_table_m) / colSums(counts_table_m)) * 100
counts_m <- as.data.frame.matrix(counts_table_m)
props_m <-as.data.frame.matrix(props_table_m)


areas <- read_xlsx(paste0(work,'/Config/tcellareas.xlsx'))
densities_m <- t(t(counts_m)/areas$TotalArea)

mdtcells <- md[md$image_name %in% unique(colnames(props_m)),]

ggdf <- melt(data.frame(cluster = rownames(props_m), props_m, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=tcellimagelevels)
ggdf$type <- factor(mdtcells$Type[match(ggdf$sample_id,mdtcells$image_name)], levels=typelevels)
ggdf$mouse <- factor(mdtcells$Mouse[match(ggdf$sample_id,mdtcells$image_name)], levels=mouselevels)

#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities_m), densities_m, check.names = FALSE),
              id.vars = "cluster", value.name = "densities_m", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdf$sample_id, levels=tcellimagelevels)
ggdfd$type <- factor(mdtcells$Type[match(ggdf$sample_id,mdtcells$image_name)], levels=typelevels)
ggdfd$mouse <- factor(mdtcells$Mouse[match(ggdf$sample_id,mdtcells$image_name)], levels=mouselevels)


#plot box plots
#% CELLS
ggp2<-ggplot(ggdf,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=4,scales="free")+
  ylab("Abundance (% of total cells)")+
  ggtitle("tiKPC Proportion T Cell Merged")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=4, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11)
  )
pdf('Abundance_box_merged_tcell1SID2.pdf',width=13,height=7)
ggp2
dev.off()

# Remove lymphoid aggregate samples
ggdfd_subset <- ggdfd[ggdf$type != "Lymphoid Aggregate",]

# First specify the comparisons for stats
type_comp <- list(c("Normal", "PanIN"), c("Normal", "Tumor Edge"), c("Normal", "Tumor Core"),
                  c("PanIn", "Tumor Edge"), c("PanIn", "Tumor Core"), c("Tumor Edge", "Tumor Core"))

# Make the stats table yourself
wc_type_stats <- compare_means(densities_m ~ type, ggdfd, method = "wilcox.test", group.by = "cluster")
write.csv(wc_type_stats, file = "wc_type_stats_tcells_subset.csv")

ggp2<-ggplot(ggdfd_subset,aes(x=type,y=densities_m,fill=type))+
  geom_boxplot(outlier.shape=NA, lwd=0.25)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=4,scales="free")+
  ylab("Density (# of cells)")+
  ggtitle("tiKPC Density T Cell Merged")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=4, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 8))+
          stat_compare_means(comparisons = type_comp, method = "wilcox.test", tip.length = 0,
                             label = "p.signif",
                             hide.ns = F,
                             bracket.size = 0.5,
                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                symbols = c("****", "***", "**", "*", "ns")), size = 2)
  
pdf('Density_box_merged_tcell1_subsetSID2.pdf',width=9,height=7)
ggp2
dev.off()

