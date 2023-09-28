############################################################################
################ script to re-produce the Figure 3  plots #################
############################################################################

#########################################################################################################
#########################################################################################################
##### scatter plots a)

## calculate mean_sub/position
grb2_mean_sub <- calculate_meansub_pos("GRB2-SH3")
pdz3_mean_sub <- calculate_meansub_pos("PSD95-PDZ3")

## add rSASA
grb2_mean_sub <- merge(grb2_mean_sub,
                       scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3",c("Pos", "rSASA", "domain")],
                       by=c("Pos"))

grb2_mean_sub <- grb2_mean_sub[!duplicated(grb2_mean_sub),]


pdz3_mean_sub <- merge(pdz3_mean_sub,
                       scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3",c("Pos", "rSASA", "domain")],
                       by=c("Pos"))

pdz3_mean_sub <- pdz3_mean_sub[!duplicated(pdz3_mean_sub),]

grb2_pdz3_mean_sub<-rbind(grb2_mean_sub,
                          pdz3_mean_sub)

## plot
scatter_figure_3_rSASA_subs(grb2_pdz3_mean_sub)

cor.test(grb2_pdz3_mean_sub[grb2_pdz3_mean_sub$domain=="GRB2-SH3",]$scaled_fitness_mean,
         grb2_pdz3_mean_sub[grb2_pdz3_mean_sub$domain=="GRB2-SH3",]$rSASA,
         method="pearson")

cor.test(grb2_pdz3_mean_sub[grb2_pdz3_mean_sub$domain=="PSD95-PDZ3",]$scaled_fitness_mean,
         grb2_pdz3_mean_sub[grb2_pdz3_mean_sub$domain=="PSD95-PDZ3",]$rSASA,
         method="pearson")


### now for CCC ins vs rSASA. 
df <- scaled_variants_aPCA[scaled_variants_aPCA$domain %in% c("GRB2-SH3", "PSD95-PDZ3") & scaled_variants_aPCA$type == "singleINS", c("rSASA", "Pos", "domain", "scaled_fitness", "scaled_sigma")]

scatter_figure_3_rSASA_indels(df)+
  ylab("abundance (CCC ins)")
  
cor.test(df[df$domain=="GRB2-SH3",]$rSASA,
         df[df$domain=="GRB2-SH3",]$scaled_fitness,
         method="pearson")

cor.test(df[df$domain=="PSD95-PDZ3",]$rSASA,
         df[df$domain=="PSD95-PDZ3",]$scaled_fitness,
         method="pearson")

### now for single del vs rSASA. 
df <- scaled_variants_aPCA[scaled_variants_aPCA$domain %in% c("GRB2-SH3", "PSD95-PDZ3") & scaled_variants_aPCA$type == "singleDEL", c("rSASA", "Pos", "domain", "scaled_fitness", "scaled_sigma")]

scatter_figure_3_rSASA_indels(df)+
  ylab("abundance (single del)")

cor.test(df[df$domain=="GRB2-SH3",]$rSASA,
         df[df$domain=="GRB2-SH3",]$scaled_fitness,
         method="pearson")

cor.test(df[df$domain=="PSD95-PDZ3",]$rSASA,
         df[df$domain=="PSD95-PDZ3",]$scaled_fitness,
         method="pearson")


#########################################################################################################
#########################################################################################################
##### scatter plots b)

#subs
tsyboyama_rsasa_sub<-c()
for (i in unique(tsuboyama_nat_doms_all$pdb_name)){
  if (nrow(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,])>3){
    temp<-cor.test(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$rSASA,
                   tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$ddG_ML_subs,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_sub <- tsyboyama_rsasa_sub %>% rbind(temp2)
  }
}

# del
tsyboyama_rsasa_del<-c()
for (i in unique(tsuboyama_nat_doms_all$pdb_name)){
  if (nrow(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,])>3){
    temp<-cor.test(tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$rSASA,
                   tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$pdb_name==i,]$ddG_ML_dels,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_del <- tsyboyama_rsasa_del %>% rbind(temp2)
  }
}

## ins_after
ins_after <- tsuboyama_nat_doms_all[,c("Pos", "pdb_name", "ddG_ML_ins", "rSASA")]
ins_after$Pos <- ins_after$Pos -1

tsyboyama_rsasa_ins_after<-c()
for (i in unique(ins_after$pdb_name)){
  if (nrow(ins_after[ins_after$pdb_name==i,])>3){
    temp<-cor.test(ins_after[ins_after$pdb_name==i,]$rSASA,
                   ins_after[ins_after$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_ins_after <- tsyboyama_rsasa_ins_after %>% rbind(temp2)
  }
}

## plot joined hist
color1 <- adjustcolor("#414487FF", alpha.f = 0.8) 
color2 <- adjustcolor("#2A788EFF", alpha.f = 0.6) 
color3 <- adjustcolor("#7AD151FF", alpha.f = 0.6)  

hist(-tsyboyama_rsasa_sub$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = color1,  
     border = "black"  
)

hist(-tsyboyama_rsasa_del$cor,
     breaks = 30,
     add = TRUE, 
     col = color3,  
     border = "black"
)

hist(-tsyboyama_rsasa_ins_after$cor,
     breaks = 30,
     add = TRUE,  
     col = color2,  
     border = "black" 
)

legend("left",
       legend = c("subs", "ins_after", "dels"),
       fill = c(color1, color2, color3),  # Transparent colors
       border = "black",
       bty = "n"
)

mlv(na.omit(tsyboyama_rsasa_sub$cor), method="naive")
mlv(na.omit(tsyboyama_rsasa_del$cor), method="naive")
mlv(na.omit(tsyboyama_rsasa_ins_after$cor), method="naive")


#########################################################################################################
#########################################################################################################
##### helix plots c)
alphahelix_ddG<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc=="AlphaHelix",]
# count the elements. 
unique_counts_alphahelix <- aggregate(element_no_simple ~ align_to_center, data = alphahelix_ddG, FUN = function(x) length(unique(x)))

## plot: subs
plot_periodity_helix(alphahelix_ddG$ddG_ML_subs,
               alphahelix_ddG$align_to_center,
               alphahelix_ddG,
               unique_counts_alphahelix) +
  ylab("ddG substitutions")+
  xlab("realigned position")


## version with collapsed residues before pos -6 and after pos 6
## create a new df. 
result<-plot_periodity_df_collapsed(alphahelix_ddG,
                            c(min(alphahelix_ddG$align_to_center):-6),
                            c(paste(min(alphahelix_ddG$align_to_center), "-6",sep=":")),
                            c(6:max(alphahelix_ddG$align_to_center)),
                            c(paste("6",max(alphahelix_ddG$align_to_center),sep=":")),
                            c("-15:-6", "-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4", "5", "6:15"))

alphahelix_ddG_collapsed <- result[[1]]

# count the elements. 
unique_counts_alphahelix <- aggregate(element_no_simple ~ align_to_center, data = alphahelix_ddG_collapsed, FUN = function(x) length(unique(x)))


## plot: ddG subs
plot_periodity_collapsed(alphahelix_ddG_collapsed$ddG_ML_subs, 
                         alphahelix_ddG_collapsed$align_to_center, 
                         alphahelix_ddG_collapsed, 
                         "red3", 
                         unique_counts_alphahelix)+
  ylab("ddG substitutions")

## plot: ddG ins
plot_periodity_collapsed(alphahelix_ddG_collapsed$ddG_ML_ins, 
                         alphahelix_ddG_collapsed$align_to_center, 
                         alphahelix_ddG_collapsed, 
                         "red3", 
                         unique_counts_alphahelix)+
  ylab("ddG insertions")

## plot: ddG dels
plot_periodity_collapsed(alphahelix_ddG_collapsed$ddG_ML_dels, 
                         alphahelix_ddG_collapsed$align_to_center, 
                         alphahelix_ddG_collapsed, 
                         "red3", 
                         unique_counts_alphahelix)+
  ylab("ddG deletions")+
  xlab("position in helix")


##### strand plots c)
strand_ddG<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc=="Strand",]
# count the elements. 
unique_counts_strand <- aggregate(element_no_simple ~ align_to_center, data = strand_ddG, FUN = function(x) length(unique(x)))

## plot: subs
plot_periodity_strand(strand_ddG$ddG_ML_subs,
                      strand_ddG$align_to_center,
                      strand_ddG,
                      unique_counts_strand) +
  ylab("ddG substitutions")+
  xlab("position in strand")

## plot: ins
plot_periodity_strand(strand_ddG$ddG_ML_ins,
                      strand_ddG$align_to_center,
                      strand_ddG,
                      unique_counts_strand) +
  ylab("ddG insertions")+
  xlab("position in strand")

## plot: dels
plot_periodity_strand(strand_ddG$ddG_ML_dels,
                      strand_ddG$align_to_center,
                      strand_ddG,
                      unique_counts_strand) +
  ylab("ddG deletions")+
  xlab("position in strand")

#########################################################################################################
#########################################################################################################
##### violin plots e)

## group data for n/termini 
groups_ddG_ml<-rbind(data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("ntermini", "ctermini"),]$ddG_ML_subs,
                                mut="subs"),
                     data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("ntermini", "ctermini"),]$ddG_ML_ins,
                                mut="ins"),
                     data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("ntermini", "ctermini"),]$ddG_ML_dels,
                                mut="dels"))

x_axis_order <- c("subs", "ins","dels")
groups_ddG_ml$mut <- factor(groups_ddG_ml$mut, levels = x_axis_order)

#plot
plot_violins_fig3(groups_ddG_ml)
p_values <- pairwise.wilcox.test(groups_ddG_ml$ddG, groups_ddG_ml$mut, p.adjust.method = "bonferroni")

## group data for loop
groups_ddG_ml<-rbind(data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("loop"),]$ddG_ML_subs,
                                mut="subs"),
                     data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("loop"),]$ddG_ML_ins,
                                mut="ins"),
                     data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("loop"),]$ddG_ML_dels,
                                mut="dels"))

x_axis_order <- c("subs", "ins","dels")
groups_ddG_ml$mut <- factor(groups_ddG_ml$mut, levels = x_axis_order)

#plot
plot_violins_fig3(groups_ddG_ml)
p_values <- pairwise.wilcox.test(groups_ddG_ml$ddG, groups_ddG_ml$mut, p.adjust.method = "bonferroni")

  
## group data for helix
groups_ddG_ml<-rbind(data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("AlphaHelix", "310Helix"),]$ddG_ML_subs,
                                mut="subs"),
                     data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("AlphaHelix", "310Helix"),]$ddG_ML_ins,
                                mut="ins"),
                     data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("AlphaHelix", "310Helix"),]$ddG_ML_dels,
                                mut="dels"))

x_axis_order <- c("subs", "ins","dels")
groups_ddG_ml$mut <- factor(groups_ddG_ml$mut, levels = x_axis_order)

#plot
plot_violins_fig3(groups_ddG_ml)
p_values <- pairwise.wilcox.test(groups_ddG_ml$ddG, groups_ddG_ml$mut, p.adjust.method = "bonferroni")

## group data for strand
groups_ddG_ml<-rbind(data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("Strand"),]$ddG_ML_subs,
                                mut="subs"),
                     data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("Strand"),]$ddG_ML_ins,
                                mut="ins"),
                     data.frame(ddG=tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("Strand"),]$ddG_ML_dels,
                                mut="dels"))

x_axis_order <- c("subs", "ins","dels")
groups_ddG_ml$mut <- factor(groups_ddG_ml$mut, levels = x_axis_order)

#plot
plot_violins_fig3(groups_ddG_ml)
p_values <- pairwise.wilcox.test(groups_ddG_ml$ddG, groups_ddG_ml$mut, p.adjust.method = "bonferroni")



  


