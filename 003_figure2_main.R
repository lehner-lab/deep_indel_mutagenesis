

############################################################################
################ script to re-produce the Figure 2  plots #################
############################################################################

#########################################################################################################
#########################################################################################################
##### heatmaps aPCA Figure 2

############ GRB2-SH3: run function to plot the heatmap for substitutions
result<-plot_heatmap_aPCA_grb2_subs()
grb2_fold_subs_hm <- result[[1]]
grb2_seq_code <- result[[2]]

## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_fold_subs_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:52,
           y = min(grb2_seq_code$Pos),
           label = 1:52,
           color=grb2_seq_code[c(1:52),]$color_seq,
           vjust = 2.8,
           size=4) +
  annotate(geom = "text",
           x = 1:52,
           y = min(grb2_seq_code$Pos),
           label = WT.Vectorised,
           vjust = 4.2,
           size=4)

############ PSD95-PDZ3: run function to plot the heatmap for substitutions
result<-plot_heatmap_aPCA_pdz3_subs()
pdz3_fold_subs_hm <- result[[1]]
pdz3_seq_code <- as.data.frame(result[[2]])

## pdz3 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_fold_subs_hm +
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:52,
           y = min(pdz3_seq_code$Pos),
           label = 1:52,
           color=pdz3_seq_code[c(1:52),]$color_seq,
           vjust = 2.8,
           size=4) +
  annotate(geom = "text",
           x = 1:52,
           y = min(pdz3_seq_code$Pos),
           label = WT.Vectorised,
           vjust = 4.2,
           size=4)

############ GRB2-SH3: run function to plot the heatmap for all insertions
result<-plot_heatmap_aPCA_grb2_allins()
grb2_fold_allins_hm <- result[[1]]

## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV "
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_fold_allins_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:53,
           y = min(grb2_seq_code$Pos),
           label = 1:53,
           color=grb2_seq_code[c(1:53),]$color_seq,
           vjust = 2.8,
           size=4) +
  annotate(geom = "text",
           x = 1:53,
           y = min(grb2_seq_code$Pos),
           label = WT.Vectorised,
           vjust = 4.2,
           size=4)

############ PSD95-PDZ3: run function to plot the heatmap for all insertions
result<-plot_heatmap_aPCA_pdz3_allins()
pdz3_fold_allins_hm <- result[[1]]

## grb2 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV "
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_fold_allins_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:53,
           y = min(pdz3_seq_code$Pos),
           label = 1:53,
           color=pdz3_seq_code[c(1:53),]$color_seq,
           vjust = 2.8,
           size=4) +
  annotate(geom = "text",
           x = 1:53,
           y = min(pdz3_seq_code$Pos),
           label = WT.Vectorised,
           vjust = 4.2,
           size=4)

############ GRB2-SH3: run function to plot the heatmap for insertion repeats
result<-plot_heatmap_aPCA_grb2_ins()
grb2_fold_ins_hm <- result[[1]]

## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_fold_ins_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:53,
           y = min(grb2_seq_code$Pos),
           label = 1:53,
           color=grb2_seq_code[c(1:53),]$color_seq,
           vjust = 2.8,
           size=4) 

############ GRB2-SH3: run function to plot the heatmap for deletions
result<-plot_heatmap_aPCA_grb2_dels()
grb2_fold_dels_hm <- result[[1]]

## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_fold_dels_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:52,
           y = min(grb2_seq_code$Pos),
           label = 1:52,
           color=grb2_seq_code[c(1:52),]$color_seq,
           vjust = 2.8,
           size=4) 

############ PSD95-PDZ3: run function to plot the heatmap for insertion repeats
result<-plot_heatmap_aPCA_pdz3_ins()
pdz3_fold_ins_hm <- result[[1]]

## grb2 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_fold_ins_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:53,
           y = min(pdz3_seq_code$Pos),
           label = 1:53,
           color=pdz3_seq_code[c(1:53),]$color_seq,
           vjust = 2.8,
           size=4) 

############ PSD95-PDZ3: run function to plot the heatmap for deletions
result<-plot_heatmap_aPCA_pdz3_dels()
pdz3_fold_dels_hm <- result[[1]]

## grb2 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_fold_dels_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:52,
           y = min(pdz3_seq_code$Pos),
           label = 1:52,
           color=pdz3_seq_code[c(1:52),]$color_seq,
           vjust = 2.8,
           size=4) 

#########################################################################################################
#########################################################################################################
##### scatterplots aPCA Figure 2b

############ GRB2-SH3
## calculate mean sub/position for GRB2-SH3 aPCA
grb2_mean_sub <- calculate_meansub_pos("GRB2-SH3")


## merge with deletion data 
grb2_mut_type_cor <- merge(scaled_variants_aPCA[scaled_variants_aPCA$type == "singleDEL" & 
                                                  scaled_variants_aPCA$domain == "GRB2-SH3",
                                                c("domain", "Pos", "scaled_fitness", "scaled_sigma")],
                           grb2_mean_sub,
                           by = c("Pos")) 

## now add ins_before and ins_afer: 1x CCC
grb2_ins_before <- scaled_variants_aPCA[scaled_variants_aPCA$type == "singleINS" & 
                                          scaled_variants_aPCA$domain == "GRB2-SH3",
                                        c("Pos", "scaled_fitness", "scaled_sigma")]
grb2_ins_after <- grb2_ins_before
grb2_ins_after$Pos <- grb2_ins_after$Pos -1

## merge to the mut_type df: grb2_mut_type_cor
grb2_mut_type_cor <- merge(grb2_mut_type_cor,
                           grb2_ins_before,
                           by= c("Pos"),
                           all.x = T)

grb2_mut_type_cor <- merge(grb2_mut_type_cor,
                           grb2_ins_after,
                           by= c("Pos"),
                           all.x = T)

colnames(grb2_mut_type_cor)[3:4] <- c("scaled_fitness_del", "scaled_sigma_del")
colnames(grb2_mut_type_cor)[7:8] <- c("scaled_fitness_ins_before", "scaled_sigma_ins_before")
colnames(grb2_mut_type_cor)[9:10] <- c("scaled_fitness_ins_after", "scaled_sigma_ins_after")

## plot del vs subs
scatter_for_fig2(grb2_mut_type_cor, 
                 grb2_mut_type_cor$scaled_fitness_mean,
                 grb2_mut_type_cor$scaled_fitness_del,
                 grb2_mut_type_cor$scaled_sigma_del)
cor.test(grb2_mut_type_cor$scaled_fitness_mean, grb2_mut_type_cor$scaled_fitness_del, method="pearson")

## plot ins_after vs subs
scatter_for_fig2(grb2_mut_type_cor, 
                 grb2_mut_type_cor$scaled_fitness_mean,
                 grb2_mut_type_cor$scaled_fitness_ins_after,
                 grb2_mut_type_cor$scaled_sigma_ins_after)
cor.test(grb2_mut_type_cor$scaled_fitness_mean, grb2_mut_type_cor$scaled_fitness_ins_after, method="pearson")

## plot ins_after vs dels
scatter_for_fig2(grb2_mut_type_cor, 
                 grb2_mut_type_cor$scaled_fitness_del,
                 grb2_mut_type_cor$scaled_fitness_ins_after,
                 grb2_mut_type_cor$scaled_sigma_ins_after)+
  geom_errorbar(data=grb2_mut_type_cor,
                aes(x=scaled_fitness_del,
                    y=scaled_fitness_ins_after,
                    xmin=scaled_fitness_del-scaled_sigma_del,
                    xmax=scaled_fitness_del+scaled_sigma_del, color="grey60"))

cor.test(grb2_mut_type_cor$scaled_fitness_del, grb2_mut_type_cor$scaled_fitness_ins_after, method="pearson")

############ PSD95-PDZ3
## calculate mean sub/position for PSD95-PDZ3 aPCA
pdz3_mean_sub <- calculate_meansub_pos("PSD95-PDZ3")

## merge with deletion data 
pdz3_mut_type_cor <- merge(scaled_variants_aPCA[scaled_variants_aPCA$type == "singleDEL" & 
                                                  scaled_variants_aPCA$domain == "PSD95-PDZ3",
                                                c("domain", "Pos", "scaled_fitness", "scaled_sigma")],
                           pdz3_mean_sub,
                           by = c("Pos")) 

## now add ins_before and ins_afer: 1x CCC
pdz3_ins_before <- scaled_variants_aPCA[scaled_variants_aPCA$type == "singleINS" & 
                                          scaled_variants_aPCA$domain == "PSD95-PDZ3",
                                        c("Pos", "scaled_fitness", "scaled_sigma")]
pdz3_ins_after <- pdz3_ins_before
pdz3_ins_after$Pos <- pdz3_ins_after$Pos -1

## merge to the mut_type df: pdz3_mut_type_cor
pdz3_mut_type_cor <- merge(pdz3_mut_type_cor,
                           pdz3_ins_before,
                           by= c("Pos"),
                           all.x = T)

pdz3_mut_type_cor <- merge(pdz3_mut_type_cor,
                           pdz3_ins_after,
                           by= c("Pos"),
                           all.x = T)

colnames(pdz3_mut_type_cor)[3:4] <- c("scaled_fitness_del", "scaled_sigma_del")
colnames(pdz3_mut_type_cor)[7:8] <- c("scaled_fitness_ins_before", "scaled_sigma_ins_before")
colnames(pdz3_mut_type_cor)[9:10] <- c("scaled_fitness_ins_after", "scaled_sigma_ins_after")

## plot del vs subs
scatter_for_fig2(pdz3_mut_type_cor, 
                 pdz3_mut_type_cor$scaled_fitness_mean,
                 pdz3_mut_type_cor$scaled_fitness_del,
                 pdz3_mut_type_cor$scaled_sigma_del)
cor.test(pdz3_mut_type_cor$scaled_fitness_mean, pdz3_mut_type_cor$scaled_fitness_del, method="pearson")

## plot ins_after vs subs
scatter_for_fig2(pdz3_mut_type_cor, 
                 pdz3_mut_type_cor$scaled_fitness_mean,
                 pdz3_mut_type_cor$scaled_fitness_ins_after,
                 pdz3_mut_type_cor$scaled_sigma_ins_after)
cor.test(pdz3_mut_type_cor$scaled_fitness_mean, pdz3_mut_type_cor$scaled_fitness_ins_after, method="pearson")

## plot ins_after vs dels
scatter_for_fig2(pdz3_mut_type_cor, 
                 pdz3_mut_type_cor$scaled_fitness_del,
                 pdz3_mut_type_cor$scaled_fitness_ins_after,
                 pdz3_mut_type_cor$scaled_sigma_ins_after)+
  geom_errorbar(data=pdz3_mut_type_cor,
                aes(x=scaled_fitness_del,
                    y=scaled_fitness_ins_after,
                    xmin=scaled_fitness_del-scaled_sigma_del,
                    xmax=scaled_fitness_del+scaled_sigma_del, color="grey60"))
cor.test(pdz3_mut_type_cor$scaled_fitness_del, pdz3_mut_type_cor$scaled_fitness_ins_after, method="pearson")

#########################################################################################################
#########################################################################################################
##### correlation histograms tsuboyama Figure 2c: mutation types
## subs vs deletions
tsuboyama_subVSdel<-c()
for (i in unique(tsuboyama_nat_doms_all$pdb_name)){
  if (nrow(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,])>3){
    temp<-cor.test(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$ddG_ML_subs,
                   tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$ddG_ML_dels,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsuboyama_subVSdel <- tsuboyama_subVSdel %>% rbind(temp2)
  }
}

## subs vs ins_after
ins_after <- tsuboyama_nat_doms_all[,c("Pos", "pdb_name", "ddG_ML_ins")]
ins_after$Pos <- ins_after$Pos -1
subs <- tsuboyama_nat_doms_all[,c("Pos", "pdb_name", "ddG_ML_subs")]

ins_afterVSsubs <- merge(ins_after,
                         subs, by=c("pdb_name", "Pos"))

tsuboyama_subVSins_after<-c()
for (i in unique(ins_afterVSsubs$pdb_name)){
  if (nrow(ins_afterVSsubs[ins_afterVSsubs$pdb_name==i,])>3){
    temp<-cor.test(ins_afterVSsubs[ins_afterVSsubs$pdb_name==i,]$ddG_ML_subs,
                   ins_afterVSsubs[ins_afterVSsubs$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_afterVSsubs[ins_afterVSsubs$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsuboyama_subVSins_after <- tsuboyama_subVSins_after %>% rbind(temp2)
  }
}



## dels vs ins_after
ins_afterVSdel <- merge(ins_after,
                        tsuboyama_nat_doms_all[,c("pdb_name", "Pos", "ddG_ML_dels")],
                        by = c("pdb_name", "Pos"))

tsuboyama_delsVSins_after<-c()
for (i in unique(ins_afterVSdel$pdb_name)){
  if (nrow(ins_afterVSdel[ins_afterVSdel$pdb_name==i,])>3){
    temp<-cor.test(ins_afterVSdel[ins_afterVSdel$pdb_name==i,]$ddG_ML_dels,
                   ins_afterVSdel[ins_afterVSdel$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_afterVSdel[ins_afterVSdel$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsuboyama_delsVSins_after <- tsuboyama_delsVSins_after %>% rbind(temp2)
  }
}

############################################
## joined histogram for main figure. 
color1 <- adjustcolor("darkorange3",alpha.f = 0.8)
color2 <- adjustcolor("antiquewhite4",alpha.f = 0.8)
color3 <- adjustcolor("wheat",alpha.f = 0.4)


hist(tsuboyama_delsVSins_after$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

## make sure its the same subset of doms
hist(tsuboyama_subVSdel[which(tsuboyama_delsVSins_after$domain %in% tsuboyama_subVSdel$domain),]$cor,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)

hist(tsuboyama_subVSins_after$cor,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color3,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)

# Add a legend to the left side of the plot
legend("left",
       legend = c("ins_afterVSdel", "subVSdel", "subVSins_after"),
       fill = c(color1, color2, color3),  # Transparent colors
       border = "black",
       bty = "n"
)

mlv(na.omit(tsuboyama_delsVSins_after$cor), method="naive")
mlv(na.omit(tsuboyama_subVSdel[which(tsuboyama_delsVSins_after$domain %in% tsuboyama_subVSdel$domain),]$cor), method="naive")
mlv(na.omit(tsuboyama_subVSins_after$cor), method="naive")



#########################################################################################################
#########################################################################################################
##### correlation histograms tsuboyama Figure 2c: mutation types but structure vs loop

## define loop and structure
struc<-c("Strand", "AlphaHelix", "310Helix")
loop<-c("ntermini", "loop","ctermini")

## ins_after vs dels
ins_afterVSdel <- merge(ins_after,
                        tsuboyama_nat_doms_all[,c("pdb_name", "Pos", "ddG_ML_dels", "simple_struc")],
                        by=c("pdb_name", "Pos"))

ins_afterVSdel_loop <- ins_afterVSdel[ins_afterVSdel$simple_struc %in% loop,]
ins_afterVSdel_struc <- ins_afterVSdel[ins_afterVSdel$simple_struc %in% struc,]

ins_afterVSdel_loop_cor<-c()
for (i in unique(ins_afterVSdel_loop$pdb_name)){
  if (nrow(ins_afterVSdel_loop[ins_afterVSdel_loop$pdb_name==i ,])>3){
    temp<-cor.test(ins_afterVSdel_loop[ins_afterVSdel_loop$pdb_name==i,]$ddG_ML_dels,
                   ins_afterVSdel_loop[ins_afterVSdel_loop$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_afterVSdel_loop[ins_afterVSdel_loop$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    ins_afterVSdel_loop_cor <- ins_afterVSdel_loop_cor %>% rbind(temp2)
  }
}

ins_afterVSdel_struc_cor<-c()
for (i in unique(ins_afterVSdel_struc$pdb_name)){
  if (nrow(ins_afterVSdel_struc[ins_afterVSdel_struc$pdb_name==i ,])>3){
    temp<-cor.test(ins_afterVSdel_struc[ins_afterVSdel_struc$pdb_name==i,]$ddG_ML_dels,
                   ins_afterVSdel_struc[ins_afterVSdel_struc$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_afterVSdel_struc[ins_afterVSdel_struc$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    ins_afterVSdel_struc_cor <- ins_afterVSdel_struc_cor %>% rbind(temp2)
  }
}

### plot
color1 <- adjustcolor("black",alpha.f = 0.8)
color2 <- adjustcolor("white",alpha.f = 0.4)

########
hist(ins_afterVSdel_struc_cor$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     ylim = c(0,30),
     main = NULL,
     xlim = c(-1, 1),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

hist(ins_afterVSdel_loop_cor$cor,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)

legend("left",
       legend = c("structure", "loop"),
       fill = c(color1, color2),  # Transparent colors
       border = "black",
       bty = "n",
       cex = 1.5
)

mlv(na.omit(ins_afterVSdel_struc_cor$cor), method="naive")
mlv(na.omit(ins_afterVSdel_loop_cor$cor), method="naive")


## ins_after vs subs
ins_afterVSsub <- merge(ins_after,
                        tsuboyama_nat_doms_all[,c("pdb_name", "Pos", "ddG_ML_subs", "simple_struc")],
                        by=c("pdb_name", "Pos"))

ins_afterVSsub_loop <- ins_afterVSsub[ins_afterVSsub$simple_struc %in% loop,]
ins_afterVSsub_struc <- ins_afterVSsub[ins_afterVSsub$simple_struc %in% struc,]

ins_afterVSsub_loop_cor<-c()
for (i in unique(ins_afterVSsub_loop$pdb_name)){
  if (nrow(ins_afterVSsub_loop[ins_afterVSsub_loop$pdb_name==i ,])>3){
    temp<-cor.test(ins_afterVSsub_loop[ins_afterVSsub_loop$pdb_name==i,]$ddG_ML_subs,
                   ins_afterVSsub_loop[ins_afterVSsub_loop$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_afterVSdel_loop[ins_afterVSdel_loop$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    ins_afterVSsub_loop_cor <- ins_afterVSsub_loop_cor %>% rbind(temp2)
  }
}

ins_afterVSsub_struc_cor<-c()
for (i in unique(ins_afterVSsub_struc$pdb_name)){
  if (nrow(ins_afterVSsub_struc[ins_afterVSsub_struc$pdb_name==i ,])>3){
    temp<-cor.test(ins_afterVSsub_struc[ins_afterVSsub_struc$pdb_name==i,]$ddG_ML_subs,
                   ins_afterVSsub_struc[ins_afterVSsub_struc$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_afterVSsub_struc[ins_afterVSsub_struc$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    ins_afterVSsub_struc_cor <- ins_afterVSsub_struc_cor %>% rbind(temp2)
  }
}

### plot
color1 <- adjustcolor("black",alpha.f = 0.8)
color2 <- adjustcolor("white",alpha.f = 0.4)

########
hist(ins_afterVSsub_struc_cor$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

hist(ins_afterVSsub_loop_cor$cor,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)

legend("left",
       legend = c("structure", "loop"),
       fill = c(color1, color2),  # Transparent colors
       border = "black",
       bty = "n",
       cex = 1.5
)

mlv(na.omit(ins_afterVSsub_struc_cor$cor), method="naive")
mlv(na.omit(ins_afterVSsub_loop_cor$cor), method="naive")


## dels vs subs
subVSdel_loop <- tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% loop,]
subVSdel_struc <- tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% struc,]

subVSdel_loop_cor<-c()
for (i in unique(subVSdel_loop$pdb_name)){
  if (nrow(subVSdel_loop[subVSdel_loop$pdb_name==i ,])>3){
    temp<-cor.test(subVSdel_loop[subVSdel_loop$pdb_name==i,]$ddG_ML_dels,
                   subVSdel_loop[subVSdel_loop$pdb_name==i,]$ddG_ML_subs,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(subVSdel_loop[subVSdel_loop$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    subVSdel_loop_cor <- subVSdel_loop_cor %>% rbind(temp2)
  }
}

subVSdel_struc_cor<-c()
for (i in unique(subVSdel_struc$pdb_name)){
  if (nrow(subVSdel_struc[subVSdel_struc$pdb_name==i ,])>3){
    temp<-cor.test(subVSdel_struc[subVSdel_struc$pdb_name==i,]$ddG_ML_dels,
                   subVSdel_struc[subVSdel_struc$pdb_name==i,]$ddG_ML_subs,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(subVSdel_struc[subVSdel_struc$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    subVSdel_struc_cor <- subVSdel_struc_cor %>% rbind(temp2)
  }
}

### plot
color1 <- adjustcolor("black",alpha.f = 0.8)
color2 <- adjustcolor("white",alpha.f = 0.4)

########
hist(subVSdel_struc_cor$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

hist(subVSdel_loop_cor$cor,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)

legend("left",
       legend = c("structure", "loop"),
       fill = c(color1, color2),  # Transparent colors
       border = "black",
       bty = "n",
       cex = 1.5
)

mlv(na.omit(subVSdel_struc_cor$cor), method="naive")
mlv(na.omit(subVSdel_loop_cor$cor), method="naive")


#########################################################################################################
#########################################################################################################
##### standard deviation per residue plots for all insertions and all substitutions for GRB2-SH3 and PSD95-PDZ3

## for GRB2-SH3
grb2_subs<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="GRB2-SH3" &
                                scaled_variants_aPCA$type=="substitution",c("Pos", "scaled_fitness", "scaled_sigma", "Mut")]

grb2_allins<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="GRB2-SH3" &
                                  scaled_variants_aPCA$type=="allins",c("Pos", "scaled_fitness", "scaled_sigma", "insID")]


###  measure of dispersion: standard deviation
variation_data_pos_grb2_subs <- grb2_subs %>%
  group_by(Pos) %>%
  summarize(
    sd_sub_pos = sd(scaled_fitness)
  )

variation_data_pos_grb2_ins <- grb2_allins %>%
  group_by(Pos) %>%
  summarize(
    sd_ins_pos = sd(scaled_fitness)
  )

## plot
color1 <- adjustcolor("#414487FF", alpha.f = 0.8) # deletion color
color2 <- adjustcolor("#2A788EFF", alpha.f = 0.8) # insertion before color

ggplot() +
  geom_bar(data=variation_data_pos_grb2_ins, aes(x=Pos,y = sd_ins_pos, fill = "Insertion"), stat = "identity", position = "dodge", alpha = 0.5) +
  geom_bar(data=variation_data_pos_grb2_subs, aes(x=Pos,y = sd_sub_pos, fill = "Substitution"), stat = "identity", position = "dodge", alpha = 0.5) +
  labs(x = "Position", y = "Standard deviation") +
  scale_fill_manual(values = c("Insertion" = color2, "Substitution" = color1)) +
  scale_x_continuous(breaks = seq(1, 53, by = 1)) +  # Set the x-axis tick positions
  theme_minimal()+
  theme(legend.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        axis.text = element_text(size=11),
        axis.title.y = element_blank())


## for PSD95-PDZ3
pdz3_subs<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="PSD95-PDZ3" &
                                  scaled_variants_aPCA$type=="substitution",c("Pos", "scaled_fitness", "scaled_sigma", "Mut")]

pdz3_allins<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="PSD95-PDZ3" &
                                    scaled_variants_aPCA$type=="allins",c("Pos", "scaled_fitness", "scaled_sigma", "insID")]


###  measure of dispersion: standard deviation
variation_data_pos_pdz3_subs <- pdz3_subs %>%
  group_by(Pos) %>%
  summarize(
    sd_sub_pos = sd(scaled_fitness)
  )

variation_data_pos_pdz3_ins <- pdz3_allins %>%
  group_by(Pos) %>%
  summarize(
    sd_ins_pos = sd(scaled_fitness)
  )

## plot
color1 <- adjustcolor("#414487FF", alpha.f = 0.8) # deletion color
color2 <- adjustcolor("#2A788EFF", alpha.f = 0.8) # insertion before color

ggplot() +
  geom_bar(data=variation_data_pos_pdz3_ins, aes(x=Pos,y = sd_ins_pos, fill = "Insertion"), stat = "identity", position = "dodge", alpha = 0.5) +
  geom_bar(data=variation_data_pos_pdz3_subs, aes(x=Pos,y = sd_sub_pos, fill = "Substitution"), stat = "identity", position = "dodge", alpha = 0.5) +
  labs(x = "Position", y = "Standard deviation") +
  scale_fill_manual(values = c("Insertion" = color2, "Substitution" = color1)) +
  scale_x_continuous(breaks = seq(1, 53, by = 1)) +  # Set the x-axis tick positions
  theme_minimal()+
  theme(legend.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        axis.text = element_text(size=11),
        axis.title.y = element_blank())

#########################################################################################################
#########################################################################################################
##### tolerance  per Mut plots for all insertions and all substitutions for GRB2-SH3 and PSD95-PDZ3

## GRB2-SH3
scatter_for_fig2_2(grb2_subs, grb2_allins)

## PSD95-PDZ3
colnames(pdz3_allins)[4]<-"Mut"
scatter_for_fig2_2(pdz3_subs, pdz3_allins)


#########################################################################################################
#########################################################################################################
##### tolerance  pos plots for all insertions for GRB2-SH3 and PSD95-PDZ3

# for GRB2-SH3 insertions 
df<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="GRB2-SH3" &
                           scaled_variants_aPCA$type=="allins",c("Pos", "scaled_fitness", "scaled_sigma", "insID")]

# Calculate the mean and standard deviation of scaled_fitness for each Pos
summary_df <- df %>%
  group_by(Pos) %>%
  summarise(min_fitness = min(scaled_fitness),
            max_fitness = max(scaled_fitness))

### find top label
max_scaled_fitness_rows <- df %>%
  group_by(Pos) %>%
  summarise(min_fitness = min(scaled_fitness),
            max_fitness = max(scaled_fitness),
            min_insID = insID[which.min(scaled_fitness)],
            max_insID = insID[which.max(scaled_fitness)])

label_df <- max_scaled_fitness_rows[,c("Pos", "max_fitness", "max_insID")]
label_df_min <- max_scaled_fitness_rows[,c("Pos", "min_fitness", "min_insID")]

## plot
insID_effect_var_plot()


# for PSD95-PDZ3 insertions 
df<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="PSD95-PDZ3" &
                           scaled_variants_aPCA$type=="allins",c("Pos", "scaled_fitness", "scaled_sigma", "insID")]

# Calculate the mean and standard deviation of scaled_fitness for each Pos
summary_df <- df %>%
  group_by(Pos) %>%
  summarise(min_fitness = min(scaled_fitness),
            max_fitness = max(scaled_fitness))

### find top label
max_scaled_fitness_rows <- df %>%
  group_by(Pos) %>%
  summarise(min_fitness = min(scaled_fitness),
            max_fitness = max(scaled_fitness),
            min_insID = insID[which.min(scaled_fitness)],
            max_insID = insID[which.max(scaled_fitness)])

label_df <- max_scaled_fitness_rows[,c("Pos", "max_fitness", "max_insID")]
label_df_min <- max_scaled_fitness_rows[,c("Pos", "min_fitness", "min_insID")]

## plot
insID_effect_var_plot()


#########################################################################################################
#########################################################################################################
##### del sub plots. 

delsubs <- scaled_variants_aPCA[scaled_variants_aPCA$mut_type  == "delSub", c("Pos", "domain", "scaled_fitness", "scaled_sigma")]
singledels <- scaled_variants_aPCA[scaled_variants_aPCA$type  == "singleDEL", c("Pos", "domain", "scaled_fitness", "scaled_sigma")]

cor_delsub_del <- merge(singledels,
                        delsubs,
                        by=c("Pos", "domain"))

colnames(cor_delsub_del)[3:6] <- c("scaled_fitness_del", 
                                   "scaled_sigma_del", 
                                   "scaled_fitness_delsub", 
                                   "scaled_sigma_delsub")

## plot
ggplot(data=cor_delsub_del,
       aes(x=scaled_fitness_delsub,
           y=scaled_fitness_del))+
  geom_errorbar(aes(ymax = scaled_fitness_del + scaled_sigma_del,
                    ymin = scaled_fitness_del - scaled_sigma_del), color = "grey30", alpha = 0.2)+
  geom_errorbar(aes(xmax = scaled_fitness_delsub + scaled_sigma_delsub,
                    xmin = scaled_fitness_delsub - scaled_sigma_delsub), color = "grey30", alpha = 0.2)+
  geom_point()+
  scale_color_manual(values = c("grey30","black"))+
  theme_classic()+
  xlab("abundance: delSub") +
  ylab("abundance: single deletion") +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 18, face = "bold")
  )

cor.test(cor_delsub_del$scaled_fitness_delsub,
         cor_delsub_del$scaled_fitness_del,
         method="pearson")



