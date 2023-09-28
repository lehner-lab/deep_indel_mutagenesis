############################################################################
################ script to re-produce the Figure 5  plots #################
############################################################################

#########################################################################################################
#########################################################################################################
##### density plots of indel and substitution effects from aPCA and bPCA b)

# isolate data from aPCA and bPCA for GRB2-SH3 and plot using custom functions
aPCA_and_bPCA_density_maps(scaled_variants_bPCA[scaled_variants_bPCA$domain == "GRB2-SH3" & scaled_variants_bPCA$mut_type %in% c("deletions","insertions", "substitutions"),],
                           scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3" & scaled_variants_aPCA$mut_type %in% c("deletions","insertions", "substitutions"),])

# isolate data from aPCA and bPCA for PSD95-PDZ3 and plot using custom functions
aPCA_and_bPCA_density_maps(scaled_variants_bPCA[scaled_variants_bPCA$domain == "PSD95-PDZ3" & scaled_variants_bPCA$mut_type %in% c("deletions","insertions", "substitutions"),],
                           scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3" & scaled_variants_aPCA$mut_type %in% c("deletions","insertions", "substitutions"),])

#########################################################################################################
#########################################################################################################
##### heatmaps of substitution, insertion and deletion effect on protein-protein binding in c) 

############ GRB2-SH3: run function to plot the heatmap for substitutions
result<-plot_heatmap_bPCA_grb2_subs()
grb2_bind_subs_hm <- result[[1]]
grb2_seq_code <- result[[2]]

## remove pos 0 (=wt)
grb2_seq_code<-grb2_seq_code[c(2:54),]

## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_bind_subs_hm+
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
result<-plot_heatmap_bPCA_pdz3_subs()
pdz3_bind_subs_hm <- result[[1]]
pdz3_seq_code <- as.data.frame(result[[2]])

## remove pos 0 (=wt)
pdz3_seq_code<-pdz3_seq_code[c(2:54),]

## pdz3 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_bind_subs_hm +
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
result<-plot_heatmap_bPCA_grb2_allins()
grb2_bind_allins_hm <- result[[1]]

## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV "
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_bind_allins_hm+
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

############ PSD95-PDZ3: run function to plot the heatmap for all insertions
result<-plot_heatmap_bPCA_pdz3_allins()
pdz3_bind_allins_hm <- result[[1]]

## grb2 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV "
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_bind_allins_hm+
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

############ GRB2-SH3: run function to plot the heatmap for deletions
result<-plot_heatmap_bPCA_grb2_dels()
grb2_bind_dels_hm <- result[[1]]

## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_bind_dels_hm+
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

############ PSD95-PDZ3: run function to plot the heatmap for deletions
result<-plot_heatmap_bPCA_pdz3_dels()
pdz3_bind_dels_hm <- result[[1]]

## pdz3 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_bind_dels_hm+
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
##### aPCA and bPCA comparison per residue for both domains d)

####### GRB2-SH3:
## calculate mean and population variation (sd) for aPCA
# isolate data
all_ins_fold <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3" &
                                       scaled_variants_aPCA$type == "allins",]

all_ins_fold_mean <- aggregate(all_ins_fold$scaled_fitness, by=list(all_ins_fold$Pos), FUN=mean) 
colnames(all_ins_fold_mean)[1:2] <- c("Pos", "mean_fold")

all_ins_fold_sd <- aggregate(all_ins_fold$scaled_fitness, by=list(all_ins_fold$Pos), FUN=sd) 
colnames(all_ins_fold_sd)[1:2] <- c("Pos", "sd_fold")

all_ins_fold <- merge(all_ins_fold_mean,
                      all_ins_fold_sd, by = "Pos")

## calculate mean and population variation (sd) for bPCA
all_ins_bind <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "GRB2-SH3" &
                                       scaled_variants_bPCA$type == "allins",]

all_ins_bind_mean <- aggregate(all_ins_bind$scaled_fitness, by=list(all_ins_bind$Pos), FUN=mean) 
colnames(all_ins_bind_mean)[1:2] <- c("Pos", "mean_bind")

all_ins_bind_sd <- aggregate(all_ins_bind$scaled_fitness, by=list(all_ins_bind$Pos), FUN=sd) 
colnames(all_ins_bind_sd)[1:2] <- c("Pos", "sd_bind")

all_ins_bind <- merge(all_ins_bind_mean,
                      all_ins_bind_sd, by = "Pos")

## merge aPCA and bPCA
grb2_aPCAvsbPCA <- merge(all_ins_fold, all_ins_bind, by = "Pos")

## plot
  ggplot(grb2_aPCAvsbPCA, aes(x = Pos)) +
    geom_bar(aes(y = (mean_bind-1)), stat = "identity", fill = "#2A788EFF") +
    geom_errorbar(aes(ymin = (mean_bind-1) - sd_bind, ymax = (mean_bind-1) + sd_bind), width = 0.2, color = "#2A788EFF")+
    geom_bar(aes(y = (mean_fold-1)), stat = "identity", fill = "gray", alpha = 0.7) +
    geom_errorbar(aes(ymin = (mean_fold-1) - sd_fold, ymax = (mean_fold-1) + sd_fold), width = 0.2, color = "gray") +
  # Set axis labels
    labs(x = "position", y = "mean ins/pos") +
  # Customize the appearance
    theme_classic()+
    scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51", "52", "53"))+
    theme(axis.text.y = element_text(size=16),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = grb2_seq_code$color_seq, size=12))+
    scale_y_continuous(breaks = seq(-1.1, 0.1, by = 0.4), 
                       labels = c("-0.1", "0.3", "0.7", "1.1")) ## need to manipulate the ylabs manually


####### PSD95-PDZ3
## calculate mean and population variation (sd) for aPCA
# isolate data
all_ins_fold <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3" &
                                         scaled_variants_aPCA$type == "allins",]
  
all_ins_fold_mean <- aggregate(all_ins_fold$scaled_fitness, by=list(all_ins_fold$Pos), FUN=mean) 
colnames(all_ins_fold_mean)[1:2] <- c("Pos", "mean_fold")
  
all_ins_fold_sd <- aggregate(all_ins_fold$scaled_fitness, by=list(all_ins_fold$Pos), FUN=sd) 
colnames(all_ins_fold_sd)[1:2] <- c("Pos", "sd_fold")
  
all_ins_fold <- merge(all_ins_fold_mean,
                      all_ins_fold_sd, by = "Pos")
  
## calculate mean and population variation (sd) for bPCA
all_ins_bind <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "PSD95-PDZ3" &
                                         scaled_variants_bPCA$type == "allins",]
  
all_ins_bind_mean <- aggregate(all_ins_bind$scaled_fitness, by=list(all_ins_bind$Pos), FUN=mean) 
colnames(all_ins_bind_mean)[1:2] <- c("Pos", "mean_bind")
  
all_ins_bind_sd <- aggregate(all_ins_bind$scaled_fitness, by=list(all_ins_bind$Pos), FUN=sd) 
colnames(all_ins_bind_sd)[1:2] <- c("Pos", "sd_bind")
  
all_ins_bind <- merge(all_ins_bind_mean,
                      all_ins_bind_sd, by = "Pos")
  
## merge
pdz3_aPCAvsbPCA <- merge(all_ins_fold, all_ins_bind, by = "Pos")
  
## plot
ggplot(pdz3_aPCAvsbPCA, aes(x = Pos)) +
    geom_bar(aes(y = (mean_bind-1)), stat = "identity", fill = "#2A788EFF") +
    geom_errorbar(aes(ymin = (mean_bind-1) - sd_bind, ymax = (mean_bind-1) + sd_bind), width = 0.2, color = "#2A788EFF")+
    geom_bar(aes(y = (mean_fold-1)), stat = "identity", fill = "gray", alpha = 0.7) +
    geom_errorbar(aes(ymin = (mean_fold-1) - sd_fold, ymax = (mean_fold-1) + sd_fold), width = 0.2, color = "gray") +
    # Set axis labels
    labs(x = "position", y = "mean ins/pos") +
    # Customize the appearance
    theme_classic()+
    scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51", "52", "53"))+
    theme(axis.text.y = element_text(size=16),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = pdz3_seq_code$color_seq, size=12))+
    scale_y_continuous(breaks = seq(-1.1, 0.4, by = 0.5),
                       labels = c("-0.1", "0.4", "0.9", "1.4")) ## need to manipulate the ylabs manually



#########################################################################################################
#########################################################################################################
##### compare aPCA and bPCA scores for insertions and substitutions in e) 

###### GRB2-SH3 
# isolate all insertions for aPCA and bPCA
all_ins_fold <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3" &
                                         scaled_variants_aPCA$type == "allins",]
  
all_ins_bind <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "GRB2-SH3" &
                                         scaled_variants_bPCA$type == "allins",]
  
grb2_foldVSbind_ins <- merge(all_ins_fold[,c("Pos", "scaled_fitness", "scaled_sigma", "insID")],
                           all_ins_bind[,c("Pos", "scaled_fitness", "scaled_sigma", "insID", "rSASA")],
                           by = c("Pos", "insID"))
  
colnames(grb2_foldVSbind_ins)[3:6]<-c("scaled_fitness_fold",
                                    "scaled_sigma_fold", 
                                    "scaled_fitness_bind",
                                    "scaled_sigma_bind")
  
grb2_foldVSbind_ins <- grb2_foldVSbind_ins[!grb2_foldVSbind_ins$Pos == "53",]
  
## mark core, surface and binding_interface residues based on rSASA
grb2_foldVSbind_ins$core_surface<-""
for (i in 1:nrow(grb2_foldVSbind_ins)){
  if (grb2_foldVSbind_ins[i,]$rSASA <30){
    grb2_foldVSbind_ins[i,]$core_surface <- "core"
    }else if (grb2_foldVSbind_ins[i,]$rSASA >30){
      grb2_foldVSbind_ins[i,]$core_surface <- "surface"
    }
}

for (i in 1:nrow(grb2_foldVSbind_ins)){
  if (grb2_foldVSbind_ins[i,]$Pos %in% c(7,9,12,13,16, 32,34:35, 37, 46, 48, 50:51)){
    grb2_foldVSbind_ins[i,]$core_surface <- "binding_surface"
  }
}
  
# isolate all substitutions for aPCA and bPCA
all_ins_fold <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3" &
                                         scaled_variants_aPCA$type == "substitution",]
  
all_ins_bind <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "GRB2-SH3" &
                                         scaled_variants_bPCA$type == "substitution",]
  
grb2_foldVSbind_sub <- merge(all_ins_fold[,c("Pos", "scaled_fitness", "scaled_sigma", "Mut")],
                           all_ins_bind[,c("Pos", "scaled_fitness", "scaled_sigma", "Mut", "rSASA")],
                           by = c("Pos", "Mut"))
  
colnames(grb2_foldVSbind_sub)[3:6]<-c("scaled_fitness_fold",
                                    "scaled_sigma_fold", 
                                    "scaled_fitness_bind",
                                    "scaled_sigma_bind")

## mark core, surface and binding_interface residues based on rSASA
grb2_foldVSbind_sub$core_surface<-""
  for (i in 1:nrow(grb2_foldVSbind_sub)){
    if (grb2_foldVSbind_sub[i,]$rSASA <30){
      grb2_foldVSbind_sub[i,]$core_surface <- "core"
    }else if (grb2_foldVSbind_sub[i,]$rSASA >30){
      grb2_foldVSbind_sub[i,]$core_surface <- "surface"
    }
}
  
for (i in 1:nrow(grb2_foldVSbind_sub)){
    if (grb2_foldVSbind_sub[i,]$Pos %in% c(7,9,12,13,16, 32,34:35, 37, 46, 48, 50:51)){
      grb2_foldVSbind_sub[i,]$core_surface <- "binding_surface"
    }
}
  
# set scale for x and y axis.  
min_val <- min(c(grb2_foldVSbind_sub$scaled_fitness_bind,
                   grb2_foldVSbind_sub$scaled_fitness_fold,
                   grb2_foldVSbind_ins$scaled_fitness_bind,
                   grb2_foldVSbind_ins$scaled_fitness_fold))
max_val <- max(c(grb2_foldVSbind_sub$scaled_fitness_bind,
                   grb2_foldVSbind_sub$scaled_fitness_fold,
                   grb2_foldVSbind_ins$scaled_fitness_bind,
                   grb2_foldVSbind_ins$scaled_fitness_fold))
    
## plot: all insertions for GRB2-SH3
ggplot(data = grb2_foldVSbind_ins,
         aes(x = scaled_fitness_fold,
             y = scaled_fitness_bind,
             color = core_surface)) +
    geom_vline(xintercept = 1)+
    geom_hline(yintercept = 1)+
    geom_point(alpha = 0.6, fill = "white", size = 3) +  # Adjust alpha here (0.6 in this example)
    scale_color_manual(values = c("red3", "blue3", "green3")) +
    xlab("protein solubility") +
    ylab("protein binding") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          legend.position = "none")+
    xlim(min_val, max_val)+
    ylim(min_val, max_val)
  
  
## plot: all substitutions for GRB2-SH3
ggplot(data = grb2_foldVSbind_sub,
         aes(x = scaled_fitness_fold,
             y = scaled_fitness_bind,
             color = core_surface)) +
    geom_vline(xintercept = 1)+
    geom_hline(yintercept = 1)+
    geom_point(alpha = 0.6, fill = "white", size = 3) +  # Adjust alpha here (0.6 in this example)
    scale_color_manual(values = c("red3", "blue3", "green3")) +
    xlab("protein solubility") +
    ylab("protein binding") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          legend.position = "none")+
    xlim(min_val, max_val)+
    ylim(min_val, max_val)


###### PSD95-PDZ3
# isolate all insertions for aPCA and bPCA
all_ins_fold <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3" &
                                         scaled_variants_aPCA$type == "allins",]
  
all_ins_bind <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "PSD95-PDZ3" &
                                         scaled_variants_bPCA$type == "allins",]
  
pdz3_foldVSbind_ins <- merge(all_ins_fold[,c("Pos", "scaled_fitness", "scaled_sigma", "insID")],
                               all_ins_bind[,c("Pos", "scaled_fitness", "scaled_sigma", "insID", "rSASA")],
                               by = c("Pos", "insID"))
  
colnames(pdz3_foldVSbind_ins)[3:6]<-c("scaled_fitness_fold",
                                        "scaled_sigma_fold", 
                                        "scaled_fitness_bind",
                                        "scaled_sigma_bind")
  
pdz3_foldVSbind_ins <- pdz3_foldVSbind_ins[!pdz3_foldVSbind_ins$Pos == "53",]
  
## mark core, surface and binding_interface residues based on rSASA
pdz3_foldVSbind_ins$core_surface<-""
  for (i in 1:nrow(pdz3_foldVSbind_ins)){
    if (pdz3_foldVSbind_ins[i,]$rSASA <30){
      pdz3_foldVSbind_ins[i,]$core_surface <- "core"
    }else if (pdz3_foldVSbind_ins[i,]$rSASA >30){
      pdz3_foldVSbind_ins[i,]$core_surface <- "surface"
    }
}
  
for (i in 1:nrow(pdz3_foldVSbind_ins)){
    if (pdz3_foldVSbind_ins[i,]$Pos %in% c(12:18, 21, 29)){
      pdz3_foldVSbind_ins[i,]$core_surface <- "binding_surface"
    }
}
  
# isolate all substitutions for aPCA and bPCA
all_ins_fold <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3" &
                                         scaled_variants_aPCA$type == "substitution",]
  
all_ins_bind <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "PSD95-PDZ3" &
                                         scaled_variants_bPCA$type == "substitution",]
  
pdz3_foldVSbind_sub <- merge(all_ins_fold[,c("Pos", "scaled_fitness", "scaled_sigma", "Mut")],
                               all_ins_bind[,c("Pos", "scaled_fitness", "scaled_sigma", "Mut", "rSASA")],
                               by = c("Pos", "Mut"))
  
colnames(pdz3_foldVSbind_sub)[3:6]<-c("scaled_fitness_fold",
                                        "scaled_sigma_fold", 
                                        "scaled_fitness_bind",
                                        "scaled_sigma_bind")
  
## mark core, surface and binding_interface residues based on rSASA
pdz3_foldVSbind_sub$core_surface<-""
for (i in 1:nrow(pdz3_foldVSbind_sub)){
    if (pdz3_foldVSbind_sub[i,]$rSASA <30){
      pdz3_foldVSbind_sub[i,]$core_surface <- "core"
    }else if (pdz3_foldVSbind_sub[i,]$rSASA >30){
      pdz3_foldVSbind_sub[i,]$core_surface <- "surface"
    }
}
  
for (i in 1:nrow(pdz3_foldVSbind_sub)){
    if (pdz3_foldVSbind_sub[i,]$Pos %in% c(12:18, 21, 29)){
      pdz3_foldVSbind_sub[i,]$core_surface <- "binding_surface"
    }
}
  
# set scale for x and y axis.  
min_val <- min(c(pdz3_foldVSbind_sub$scaled_fitness_bind,
                   pdz3_foldVSbind_sub$scaled_fitness_fold,
                   pdz3_foldVSbind_ins$scaled_fitness_bind,
                   pdz3_foldVSbind_ins$scaled_fitness_fold))
max_val <- max(c(pdz3_foldVSbind_sub$scaled_fitness_bind,
                   pdz3_foldVSbind_sub$scaled_fitness_fold,
                   pdz3_foldVSbind_ins$scaled_fitness_bind,
                   pdz3_foldVSbind_ins$scaled_fitness_fold))
  
## plot: all insertions for PSD95-PDZ3
ggplot(data = pdz3_foldVSbind_ins,
         aes(x = scaled_fitness_fold,
             y = scaled_fitness_bind,
             color = core_surface)) +
    geom_vline(xintercept = 1)+
    geom_hline(yintercept = 1)+
    geom_point(alpha = 0.6, fill = "white", size = 3) +  # Adjust alpha here (0.6 in this example)
    scale_color_manual(values = c("red3", "blue3", "green3")) +
    xlab("protein solubility") +
    ylab("protein binding") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          legend.position = "none")+
    xlim(min_val, max_val)+
    ylim(min_val, max_val)
  
  
## plot: all substitutions for PSD95-PDZ3
ggplot(data = pdz3_foldVSbind_sub,
         aes(x = scaled_fitness_fold,
             y = scaled_fitness_bind,
             color = core_surface)) +
    geom_vline(xintercept = 1)+
    geom_hline(yintercept = 1)+
    geom_point(alpha = 0.6, fill = "white", size = 3) +  # Adjust alpha here (0.6 in this example)
    scale_color_manual(values = c("red3", "blue3", "green3")) +
    xlab("protein solubility") +
    ylab("protein binding") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          legend.position = "none")+
    xlim(min_val, max_val)+
    ylim(min_val, max_val)







