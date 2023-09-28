############################################################################
########## script to re-produce the Extended figure 5  plots ##############
############################################################################


#########################################################################################################
#########################################################################################################
##### insertion repeats and delSub  heatmaps for a) and b)

############ GRB2-SH3: run function to plot the heatmap for insertion repeats
result<-plot_heatmap_bPCA_grb2_ins()
grb2_bind_ins_hm <- result

## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_bind_ins_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:53,
           y = min(grb2_seq_code$Pos),
           label = 1:53,
           color=grb2_seq_code[c(1:53),]$color_seq,
           vjust = 2.9,
           size=4) 

############ PSD95-PDZ3: run function to plot the heatmap for insertion repeats
result<-plot_heatmap_bPCA_pdz3_ins()
pdz3_bind_ins_hm <- result

## grb2 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_bind_ins_hm+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:53,
           y = min(pdz3_seq_code$Pos),
           label = 1:53,
           color=pdz3_seq_code[c(1:53),]$color_seq,
           vjust = 2.9,
           size=4) 

############ plot Delsubs for GRB2-SH3
result <- plot_heatmap_bPCA_grb2_delsubs()
grb2_bind_hm_delsub <- result

## isolate the delsub substitutions
result <- find_subID_delsubs_grb2()

grb2_delsubs_1 <- result[[1]]
grb2_delsubs_2 <- result[[2]]

grb2_bind_hm_delsub+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:52,
           y = min(grb2_seq_code$Pos),
           label = 1:52,
           color=grb2_seq_code[c(1:52),]$color_seq,
           vjust = 3.5,
           size=4) +
  annotate(geom = "text",
           x = grb2_delsubs_1$Pos,
           y = min(grb2_seq_code$Pos),
           label = grb2_delsubs_1$subID,
           vjust = 5,
           size=4) +
  annotate(geom = "text",
           x = grb2_delsubs_2$Pos,
           y = min(grb2_seq_code$Pos),
           label = grb2_delsubs_2$subID,
           vjust = 6.5,
           size=4) 

############ plot Delsubs for PSD95-PDZ3
result <- plot_heatmap_bPCA_pdz3_delsubs()
pdz3_bind_hm_delsub <- result

## isolate the delsub substitutions
result <- find_subID_delsubs_pdz3()

pdz3_delsubs_1 <- result[[1]]
pdz3_delsubs_2 <- result[[2]]

pdz3_bind_hm_delsub+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off")+
  annotate(geom = "text",
           x = 1:52,
           y = min(pdz3_seq_code$Pos),
           label = 1:52,
           color=pdz3_seq_code[c(1:52),]$color_seq,
           vjust = 3.5,
           size=4) +
  annotate(geom = "text",
           x = pdz3_delsubs_1$Pos,
           y = min(pdz3_seq_code$Pos),
           label = pdz3_delsubs_1$subID,
           vjust = 5,
           size=4) +
  annotate(geom = "text",
           x = pdz3_delsubs_2$Pos,
           y = min(pdz3_seq_code$Pos),
           label = pdz3_delsubs_2$subID,
           vjust = 6.5,
           size=4) 









