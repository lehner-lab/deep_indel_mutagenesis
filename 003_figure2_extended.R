############################################################################
################ script to re-produce the Extended Figure 2  plots #########
############################################################################

##### scatterplots aPCA a)
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


## plot ins_before vs subs
scatter_for_fig2(grb2_mut_type_cor, 
                 grb2_mut_type_cor$scaled_fitness_mean,
                 grb2_mut_type_cor$scaled_fitness_ins_before,
                 grb2_mut_type_cor$scaled_sigma_ins_before)
cor.test(grb2_mut_type_cor$scaled_fitness_mean, grb2_mut_type_cor$scaled_fitness_ins_before, method="pearson")

## plot ins_before vs dels
scatter_for_fig2(grb2_mut_type_cor, 
                 grb2_mut_type_cor$scaled_fitness_del,
                 grb2_mut_type_cor$scaled_fitness_ins_before,
                 grb2_mut_type_cor$scaled_sigma_ins_before)+
  geom_errorbar(data=grb2_mut_type_cor,
                aes(x=scaled_fitness_del,
                    y=scaled_fitness_ins_before,
                    xmin=scaled_fitness_del-scaled_sigma_del,
                    xmax=scaled_fitness_del+scaled_sigma_del, color="grey60"))
cor.test(grb2_mut_type_cor$scaled_fitness_del, grb2_mut_type_cor$scaled_fitness_ins_before, method="pearson")


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



## plot ins_before vs subs
scatter_for_fig2(pdz3_mut_type_cor, 
                 pdz3_mut_type_cor$scaled_fitness_mean,
                 pdz3_mut_type_cor$scaled_fitness_ins_before,
                 pdz3_mut_type_cor$scaled_sigma_ins_before)
cor.test(pdz3_mut_type_cor$scaled_fitness_mean, pdz3_mut_type_cor$scaled_fitness_ins_before, method="pearson")


## plot ins_before vs dels
scatter_for_fig2(pdz3_mut_type_cor, 
                 pdz3_mut_type_cor$scaled_fitness_del,
                 pdz3_mut_type_cor$scaled_fitness_ins_before,
                 pdz3_mut_type_cor$scaled_sigma_ins_before)+
  geom_errorbar(data=pdz3_mut_type_cor,
                aes(x=scaled_fitness_del,
                    y=scaled_fitness_ins_before,
                    xmin=scaled_fitness_del-scaled_sigma_del,
                    xmax=scaled_fitness_del+scaled_sigma_del, color="grey60"))
cor.test(pdz3_mut_type_cor$scaled_fitness_del, pdz3_mut_type_cor$scaled_fitness_ins_before, method="pearson")


#########################################################################################################
#########################################################################################################
##### correlation histograms tsuboyama b) and c)
## subs vs ins_before
tsuboyama_subVSins_before<-c()
for (i in unique(tsuboyama_nat_doms_all$pdb_name)){
  if (nrow(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,])>3){
    temp<-cor.test(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$ddG_ML_subs,
                   tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsuboyama_subVSins_before <- tsuboyama_subVSins_before %>% rbind(temp2)
  }
}

## dels vs ins_before
tsuboyama_delsVSins_before<-c()
for (i in unique(tsuboyama_nat_doms_all$pdb_name)){
  if (nrow(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,])>3){
    temp<-cor.test(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$ddG_ML_dels,
                   tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsuboyama_delsVSins_before <- tsuboyama_delsVSins_before %>% rbind(temp2)
  }
}

############################################
## joined histogram for supplementary figure 2. 
color1 <- adjustcolor("black",alpha.f = 0.8)
color2 <- adjustcolor("white",alpha.f = 0.4)

########
## for substitutions vs ins before/after
hist(tsuboyama_subVSins_after$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     ylim = c(0,25),
     main = NULL,
     xlim = c(-1, 1),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

hist(tsuboyama_subVSins_before$cor,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)

legend("left",
       legend = c("subVSins_after", "subVSins_before"),
       fill = c(color1, color2),  # Transparent colors
       border = "black",
       bty = "n"
)

mlv(na.omit(tsuboyama_subVSins_before$cor), method="naive")

################################
## for deletions vs ins before/after
hist(tsuboyama_delsVSins_after$cor,
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

hist(tsuboyama_delsVSins_before$cor,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)

legend("left",
       legend = c("delVSins_after", "delVSins_before"),
       fill = c(color1, color2),  # Transparent colors
       border = "black",
       bty = "n"
)

mlv(na.omit(tsuboyama_delsVSins_after$cor), method="naive")
mlv(na.omit(tsuboyama_delsVSins_before$cor), method="naive")



#########################################################################################################
#########################################################################################################
##### correlation histograms tsuboyama: mutation types but structure vs loop

## define loop and structure
struc<-c("Strand", "AlphaHelix", "310Helix")
loop<-c("ntermini", "loop","ctermini")

## ins_before vs dels
ins_beforeVSdel_loop <- tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% loop,]
ins_beforeVSdel_struc <- tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% struc,]

ins_beforeVSdel_loop_cor<-c()
for (i in unique(ins_beforeVSdel_loop$pdb_name)){
  if (nrow(ins_beforeVSdel_loop[ins_beforeVSdel_loop$pdb_name==i ,])>3){
    temp<-cor.test(ins_beforeVSdel_loop[ins_beforeVSdel_loop$pdb_name==i,]$ddG_ML_dels,
                   ins_beforeVSdel_loop[ins_beforeVSdel_loop$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_beforeVSdel_loop[ins_beforeVSdel_loop$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    ins_beforeVSdel_loop_cor <- ins_beforeVSdel_loop_cor %>% rbind(temp2)
  }
}

ins_beforeVSdel_struc_cor<-c()
for (i in unique(ins_beforeVSdel_struc$pdb_name)){
  if (nrow(ins_beforeVSdel_struc[ins_beforeVSdel_struc$pdb_name==i ,])>3){
    temp<-cor.test(ins_beforeVSdel_struc[ins_beforeVSdel_struc$pdb_name==i,]$ddG_ML_dels,
                   ins_beforeVSdel_struc[ins_beforeVSdel_struc$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_beforeVSdel_struc[ins_beforeVSdel_struc$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    ins_beforeVSdel_struc_cor <- ins_beforeVSdel_struc_cor %>% rbind(temp2)
  }
}

### plot
color1 <- adjustcolor("black",alpha.f = 0.8)
color2 <- adjustcolor("white",alpha.f = 0.4)

########
hist(ins_beforeVSdel_struc_cor$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     ylim = c(0,30),
     xlim = c(-1, 1),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

hist(ins_beforeVSdel_loop_cor$cor,
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

## ins_before vs subs
ins_beforeVSsub_loop <- tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% loop,]
ins_beforeVSsub_struc <- tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% struc,]

ins_beforeVSsub_loop_cor<-c()
for (i in unique(ins_beforeVSsub_loop$pdb_name)){
  if (nrow(ins_beforeVSsub_loop[ins_beforeVSsub_loop$pdb_name==i ,])>3){
    temp<-cor.test(ins_beforeVSsub_loop[ins_beforeVSsub_loop$pdb_name==i,]$ddG_ML_subs,
                   ins_beforeVSsub_loop[ins_beforeVSsub_loop$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_beforeVSsub_loop[ins_beforeVSsub_loop$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    ins_beforeVSsub_loop_cor <- ins_beforeVSsub_loop_cor %>% rbind(temp2)
  }
}

ins_beforeVSsub_struc_cor<-c()
for (i in unique(ins_beforeVSsub_struc$pdb_name)){
  if (nrow(ins_beforeVSsub_struc[ins_beforeVSsub_struc$pdb_name==i ,])>3){
    temp<-cor.test(ins_beforeVSsub_struc[ins_beforeVSsub_struc$pdb_name==i,]$ddG_ML_subs,
                   ins_beforeVSsub_struc[ins_beforeVSsub_struc$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    
    temp2<-data.frame(domain=i,
                      rows=nrow(ins_beforeVSsub_struc[ins_beforeVSsub_struc$pdb_name==i,]),
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    ins_beforeVSsub_struc_cor <- ins_beforeVSsub_struc_cor %>% rbind(temp2)
  }
}

### plot
color1 <- adjustcolor("black",alpha.f = 0.8)
color2 <- adjustcolor("white",alpha.f = 0.4)

########
hist(ins_beforeVSsub_struc_cor$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     ylim = c(0,20),
     xlim = c(-1, 1),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

hist(ins_beforeVSsub_loop_cor$cor,
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

mlv(na.omit(ins_beforeVSsub_struc_cor$cor), method="naive")
mlv(na.omit(ins_beforeVSsub_loop_cor$cor), method="naive")


#########################################################################################################
#########################################################################################################
##### tolerance  per Mut plots for all insertions and all substitutions for GRB2-SH3 and PSD95-PDZ3 but divided by secondary structure element. 

## subset data for domain + structure and use function below to plot correlation between ins and sub (/mutation)
grb2_subs<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="GRB2-SH3" &
                                  scaled_variants_aPCA$type=="substitution",c("Pos", "scaled_fitness", "scaled_sigma", "Mut", "Structure_name")]

grb2_allins<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="GRB2-SH3" &
                                    scaled_variants_aPCA$type=="allins",c("Pos", "scaled_fitness", "scaled_sigma", "insID", "Structure_name")]

colnames(grb2_allins)[4]<-"Mut"
grb2_subs <- grb2_subs %>% mutate(simple_struc=
                                    ifelse(Structure_name == "Coil" | Structure_name == "Turn" | Structure_name == "Bridge", "loop",
                                           ifelse(Structure_name == "310Helix", "310Helix",
                                                  ifelse(Structure_name == "AlphaHelix", "AlphaHelix", 
                                                         ifelse(Structure_name == "Strand", "Strand", NA)))))
grb2_allins <- grb2_allins %>% mutate(simple_struc=
                                        ifelse(Structure_name == "Coil" | Structure_name == "Turn" | Structure_name == "Bridge", "loop",
                                               ifelse(Structure_name == "310Helix", "310Helix",
                                                      ifelse(Structure_name == "AlphaHelix", "AlphaHelix", 
                                                             ifelse(Structure_name == "Strand", "Strand", NA)))))

# plot: loop
scatter_for_fig2_2(grb2_subs[grb2_subs$simple_struc == "loop",],
                   grb2_allins[grb2_allins$simple_struc == "loop",])

# plot: strand
scatter_for_fig2_2(grb2_subs[grb2_subs$simple_struc == "Strand",],
                   grb2_allins[grb2_allins$simple_struc == "Strand",])

# plot: 310Helix
scatter_for_fig2_2(grb2_subs[grb2_subs$simple_struc == "310Helix",],
                   grb2_allins[grb2_allins$simple_struc == "310Helix",])

## now for pdz3
pdz3_subs<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="PSD95-PDZ3" &
                                  scaled_variants_aPCA$type=="substitution",c("Pos", "scaled_fitness", "scaled_sigma", "Mut", "Structure_name")]

pdz3_allins<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="PSD95-PDZ3" &
                                    scaled_variants_aPCA$type=="allins",c("Pos", "scaled_fitness", "scaled_sigma", "insID", "Structure_name")]

colnames(pdz3_allins)[4]<-"Mut"
pdz3_subs <- pdz3_subs %>% mutate(simple_struc=
                                    ifelse(Structure_name == "Coil" | Structure_name == "Turn" | Structure_name == "Bridge", "loop",
                                           ifelse(Structure_name == "310Helix", "310Helix",
                                                  ifelse(Structure_name == "AlphaHelix", "AlphaHelix", 
                                                         ifelse(Structure_name == "Strand", "Strand", NA)))))
pdz3_allins <- pdz3_allins %>% mutate(simple_struc=
                                        ifelse(Structure_name == "Coil" | Structure_name == "Turn" | Structure_name == "Bridge", "loop",
                                               ifelse(Structure_name == "310Helix", "310Helix",
                                                      ifelse(Structure_name == "AlphaHelix", "AlphaHelix", 
                                                             ifelse(Structure_name == "Strand", "Strand", NA)))))

# plot: loop
scatter_for_fig2_2(pdz3_subs[pdz3_subs$simple_struc == "loop",],
                   pdz3_allins[pdz3_allins$simple_struc == "loop",])

# plot: strand
scatter_for_fig2_2(pdz3_subs[pdz3_subs$simple_struc == "Strand",],
                   pdz3_allins[pdz3_allins$simple_struc == "Strand",])

# plot: 310Helix
scatter_for_fig2_2(pdz3_subs[pdz3_subs$simple_struc == "AlphaHelix",],
                   pdz3_allins[pdz3_allins$simple_struc == "AlphaHelix",])


#########################################################################################################
#########################################################################################################
##### tolerance  pos plots for all insertions and all substitutions for GRB2-SH3 and PSD95-PDZ3

# for GRB2-SH3 substititons. 
df<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="GRB2-SH3" &
                           scaled_variants_aPCA$type=="substitution",c("Pos", "scaled_fitness", "scaled_sigma", "Mut")]

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
            min_insID = Mut[which.min(scaled_fitness)],
            max_insID = Mut[which.max(scaled_fitness)])

label_df <- max_scaled_fitness_rows[,c("Pos", "max_fitness", "max_insID")]
label_df_min <- max_scaled_fitness_rows[,c("Pos", "min_fitness", "min_insID")]

## plot
mut_effect_var_plot()

# for PSD95-PDZ3 substititons. 
df<-scaled_variants_aPCA[scaled_variants_aPCA$domain=="PSD95-PDZ3" &
                           scaled_variants_aPCA$type=="substitution",c("Pos", "scaled_fitness", "scaled_sigma", "Mut")]

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
            min_insID = Mut[which.min(scaled_fitness)],
            max_insID = Mut[which.max(scaled_fitness)])

label_df <- max_scaled_fitness_rows[,c("Pos", "max_fitness", "max_insID")]
label_df_min <- max_scaled_fitness_rows[,c("Pos", "min_fitness", "min_insID")]

## plot
mut_effect_var_plot()

#########################################################################################################
#########################################################################################################
##### correlation all ins_before vs all subs for Supplemental Figure 2
grb2_df_ins <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3" &
                                      scaled_variants_aPCA$type == c("allins"),c("Pos", "scaled_fitness", "scaled_sigma", "insID")]

colnames(grb2_df_ins)[4] <- "Mut"

grb2_df_subs <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3" &
                                       scaled_variants_aPCA$type %in% c("substitution"), c("Pos", "scaled_fitness", "scaled_sigma", "Mut")]

grb2_df <- merge(grb2_df_ins,
                 grb2_df_subs,
                 by = c("Pos", "Mut"))

colnames(grb2_df)[3:6] <- c("scaled_fitness_ins", "scaled_sigma_ins", "scaled_fitness_subs", "scaled_sigma_subs")

ggplot(data=grb2_df,
       aes(x = scaled_fitness_ins,
           y = scaled_fitness_subs))+
  geom_point()+
  scale_color_manual(values = c("grey30","black"))+
  theme_classic()+
  xlab("abundance: insertions before") +
  ylab("abundance: substitutions") +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 18, face = "bold")
  )

## pdz3
pdz3_df_ins <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3" &
                                      scaled_variants_aPCA$type == c("allins"),c("Pos", "scaled_fitness", "scaled_sigma", "insID")]

colnames(pdz3_df_ins)[4] <- "Mut"

pdz3_df_subs <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3" &
                                       scaled_variants_aPCA$type %in% c("substitution"), c("Pos", "scaled_fitness", "scaled_sigma", "Mut")]

pdz3_df <- merge(pdz3_df_ins,
                 pdz3_df_subs,
                 by = c("Pos", "Mut"))

colnames(pdz3_df)[3:6] <- c("scaled_fitness_ins", "scaled_sigma_ins", "scaled_fitness_subs", "scaled_sigma_subs")

ggplot(data=pdz3_df,
       aes(x = scaled_fitness_ins,
           y = scaled_fitness_subs))+
  geom_point()+
  scale_color_manual(values = c("grey30","black"))+
  theme_classic()+
  xlab("abundance: insertions before") +
  ylab("abundance: substitutions") +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 18, face = "bold")
  )


#########################################################################################################
#########################################################################################################
##### correlation all ins_after vs all subs for Supplemental Figure 2
grb2_df_ins <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3" &
                                      scaled_variants_aPCA$type == c("allins"),c("Pos", "scaled_fitness", "scaled_sigma", "insID")]

colnames(grb2_df_ins)[4] <- "Mut"
grb2_df_ins$Pos <- grb2_df_ins$Pos-1

grb2_df_subs <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3" &
                                       scaled_variants_aPCA$type %in% c("substitution"), c("Pos", "scaled_fitness", "scaled_sigma", "Mut")]

grb2_df <- merge(grb2_df_ins,
                 grb2_df_subs,
                 by = c("Pos", "Mut"))

colnames(grb2_df)[3:6] <- c("scaled_fitness_ins", "scaled_sigma_ins", "scaled_fitness_subs", "scaled_sigma_subs")

ggplot(data=grb2_df,
       aes(x = scaled_fitness_ins,
           y = scaled_fitness_subs))+
  geom_point()+
  scale_color_manual(values = c("grey30","black"))+
  theme_classic()+
  xlab("abundance: insertions after") +
  ylab("abundance: substitutions") +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 18, face = "bold")
  )

## pdz3
pdz3_df_ins <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3" &
                                      scaled_variants_aPCA$type == c("allins"),c("Pos", "scaled_fitness", "scaled_sigma", "insID")]

colnames(pdz3_df_ins)[4] <- "Mut"
pdz3_df_ins$Pos <- pdz3_df_ins$Pos-1

pdz3_df_subs <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3" &
                                       scaled_variants_aPCA$type %in% c("substitution"), c("Pos", "scaled_fitness", "scaled_sigma", "Mut")]

pdz3_df <- merge(pdz3_df_ins,
                 pdz3_df_subs,
                 by = c("Pos", "Mut"))

colnames(pdz3_df)[3:6] <- c("scaled_fitness_ins", "scaled_sigma_ins", "scaled_fitness_subs", "scaled_sigma_subs")

ggplot(data=pdz3_df,
       aes(x = scaled_fitness_ins,
           y = scaled_fitness_subs))+
  geom_point()+
  scale_color_manual(values = c("grey30","black"))+
  theme_classic()+
  xlab("abundance: insertions after") +
  ylab("abundance: substitutions") +
  theme(
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 18, face = "bold")
  )



















