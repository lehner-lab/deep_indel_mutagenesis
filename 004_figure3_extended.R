############################################################################
########   script to re-produce the Extended Figure 3  plots  ##############
############################################################################

## load packages
library(gridExtra)

#########################################################################################################
#########################################################################################################
##### hist plots a), correlation to rSASA but divided into loop and structure. 

#subs
struc<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("310Helix", "AlphaHelix", "Strand"),]
loop<-tsuboyama_nat_doms_all[!tsuboyama_nat_doms_all$simple_struc %in% c("310Helix", "AlphaHelix", "Strand"),]

tsyboyama_rsasa_sub_struc<-c()
for (i in unique(struc$pdb_name)){
  if (nrow(struc[struc$pdb_name==i,])>3){
    temp<-cor.test(struc[struc$pdb_name==i,]$rSASA,
                   struc[struc$pdb_name==i,]$ddG_ML_subs,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_sub_struc <- tsyboyama_rsasa_sub_struc %>% rbind(temp2)
  }
}

tsyboyama_rsasa_sub_loop<-c()
for (i in unique(loop$pdb_name)){
  if (nrow(loop[loop$pdb_name==i,])>3){
    temp<-cor.test(loop[loop$pdb_name==i,]$rSASA,
                   loop[loop$pdb_name==i,]$ddG_ML_subs,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_sub_loop <- tsyboyama_rsasa_sub_loop %>% rbind(temp2)
  }
}

color1 <- adjustcolor("#414487FF", alpha.f = 0.8) 
color2 <- adjustcolor("white", alpha.f = 0.8) 
hist(-tsyboyama_rsasa_sub_loop$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     ylab = "substitution",
     xlim = c(-1, 1),
     ylim = c(0,25),
     col = color2,  
     border = "black"  
)

hist(-tsyboyama_rsasa_sub_struc$cor,
     breaks = 30,
     add = TRUE, 
     col = color1,  
     border = "black"
)

legend("left",
       legend = c("loop", "structure"),
       fill = c(color2, color1),  
       border = "black",
       bty = "n"
)

## dels
tsyboyama_rsasa_del_struc<-c()
for (i in unique(struc$pdb_name)){
  if (nrow(struc[struc$pdb_name==i,])>3){
    temp<-cor.test(struc[struc$pdb_name==i,]$rSASA,
                   struc[struc$pdb_name==i,]$ddG_ML_dels,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_del_struc <- tsyboyama_rsasa_del_struc %>% rbind(temp2)
  }
}

tsyboyama_rsasa_del_loop<-c()
for (i in unique(loop$pdb_name)){
  if (nrow(loop[loop$pdb_name==i,])>3){
    temp<-cor.test(loop[loop$pdb_name==i,]$rSASA,
                   loop[loop$pdb_name==i,]$ddG_ML_dels,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_del_loop <- tsyboyama_rsasa_del_loop %>% rbind(temp2)
  }
}

color1 <- adjustcolor("#7AD151FF", alpha.f = 0.8) 
color2 <- adjustcolor("white", alpha.f = 0.8) 
hist(-tsyboyama_rsasa_del_loop$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     ylab = "deletions",
     xlim = c(-1, 1),
     ylim = c(0,25),
     col = color2,  
     border = "black"  
)

hist(-tsyboyama_rsasa_del_struc$cor,
     breaks = 30,
     add = TRUE, 
     col = color1,  
     border = "black"
)

legend("left",
       legend = c("loop", "structure"),
       fill = c(color2, color1),  
       border = "black",
       bty = "n"
)

## ins
tsyboyama_rsasa_ins_struc<-c()
for (i in unique(struc$pdb_name)){
  if (nrow(struc[struc$pdb_name==i,])>3){
    temp<-cor.test(struc[struc$pdb_name==i,]$rSASA,
                   struc[struc$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_ins_struc <- tsyboyama_rsasa_ins_struc %>% rbind(temp2)
  }
}

tsyboyama_rsasa_ins_loop<-c()
for (i in unique(loop$pdb_name)){
  if (nrow(loop[loop$pdb_name==i,])>3){
    temp<-cor.test(loop[loop$pdb_name==i,]$rSASA,
                   loop[loop$pdb_name==i,]$ddG_ML_ins,
                   method="pearson")
    temp2<-data.frame(domain=i,
                      cor=temp$estimate,
                      pvalue=temp$p.value)
    
    tsyboyama_rsasa_ins_loop <- tsyboyama_rsasa_ins_loop %>% rbind(temp2)
  }
}

color1 <- adjustcolor("#2A788EFF", alpha.f = 0.8) 
color2 <- adjustcolor("white", alpha.f = 0.8) 
hist(-tsyboyama_rsasa_ins_loop$cor,
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     ylab = "insertions",
     xlim = c(-1, 1),
     ylim = c(0,25),
     col = color2,  
     border = "black"  
)

hist(-tsyboyama_rsasa_ins_struc$cor,
     breaks = 30,
     add = TRUE, 
     col = color1,  
     border = "black"
)

legend("left",
       legend = c("loop", "structure"),
       fill = c(color2, color1),  
       border = "black",
       bty = "n"
)



#########################################################################################################
#########################################################################################################
##### b) effects of indels/subs on ntermini: divide into short and long termini
ntermini_ddG<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc=="ntermini",]
## boxplot:
unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),], FUN = function(x) length(unique(x)))

result <- termini_plot(ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),],
              ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),]$ddG_ML_subs,
              unique_counts)

short <- result[[1]]


##### collapse the last residues
# Collapse data for align_to_center from -7 to -14
ntermini_ddG$align_to_center <- ifelse(
  ntermini_ddG$align_to_center %in% c(-6:-18),
  "-18:-6",
  as.character(ntermini_ddG$align_to_center)
)
# Create a new factor variable for align_to_center with desired order
ntermini_ddG$align_to_center <- factor(
  ntermini_ddG$align_to_center,
  levels = c("-18:-6","-5","-4","-3","-2","-1")
)

unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ntermini_ddG[ntermini_ddG$element_lenght %in% 4:14, ], FUN = function(x) length(unique(x)))

result <- termini_plot(ntermini_ddG[ntermini_ddG$element_lenght%in% c(4:14),],
                        ntermini_ddG[ntermini_ddG$element_lenght%in% c(4:14),]$ddG_ML_subs,
                        unique_counts)
long <- result[[1]]
long <- long+ylab("ddG substitutions")+theme(axis.text.y = element_text(size=16),
                                             axis.title.y = element_text(size=18,angle=90))

## combine
combined_plot <- grid.arrange(long, short, ncol = 2, widths = c(3, 2))

## insertions
unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),], FUN = function(x) length(unique(x)))
result <- termini_plot(ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),],
                        ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),]$ddG_ML_ins,
                        unique_counts)

short <- result[[1]]

unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ntermini_ddG[ntermini_ddG$element_lenght %in% 4:14, ], FUN = function(x) length(unique(x)))

result <- termini_plot(ntermini_ddG[ntermini_ddG$element_lenght%in% c(4:14),],
                        ntermini_ddG[ntermini_ddG$element_lenght%in% c(4:14),]$ddG_ML_ins,
                        unique_counts)
long <- result[[1]]
long <- long+ylab("ddG insertions")+theme(axis.text.y = element_text(size=16),
                                             axis.title.y = element_text(size=18,angle=90))

## combine
combined_plot <- grid.arrange(long, short, ncol = 2, widths = c(3, 2))

## deletions
unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),], FUN = function(x) length(unique(x)))
result <- termini_plot(ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),],
                        ntermini_ddG[ntermini_ddG$element_lenght%in% c(1:3),]$ddG_ML_dels,
                        unique_counts)

short <- result[[1]]

unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ntermini_ddG[ntermini_ddG$element_lenght %in% 4:14, ], FUN = function(x) length(unique(x)))

result <- termini_plot(ntermini_ddG[ntermini_ddG$element_lenght%in% c(4:14),],
                        ntermini_ddG[ntermini_ddG$element_lenght%in% c(4:14),]$ddG_ML_dels,
                        unique_counts)
long <- result[[1]]
long <- long+ylab("ddG deletions")+theme(axis.text.y = element_text(size=16),
                                          axis.title.y = element_text(size=18,angle=90))

## combine
combined_plot <- grid.arrange(long, short, ncol = 2, widths = c(3, 2))

#########################################################################################################
#########################################################################################################
##### b) effects of indels/subs on ctermini: divide into short and long termini
ctermini_ddG<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc=="ctermini",]
## boxplot:
unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),], FUN = function(x) length(unique(x)))

result <- termini_plot(ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),],
                        ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),]$ddG_ML_subs,
                        unique_counts)

short <- result[[1]]


##### collapse the last residues
# Collapse data for align_to_center from -7 to -14
ctermini_ddG$align_to_center <- ifelse(
  ctermini_ddG$align_to_center %in% c(6:8),
  "6:8",
  as.character(ctermini_ddG$align_to_center)
)
# Create a new factor variable for align_to_center with desired order
ctermini_ddG$align_to_center <- factor(
  ctermini_ddG$align_to_center,
  levels = c("1","2","3","4","5","6:8")
)

unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ctermini_ddG[ctermini_ddG$element_lenght %in% 4:14, ], FUN = function(x) length(unique(x)))

result <- termini_plot(ctermini_ddG[ctermini_ddG$element_lenght%in% c(4:14),],
                        ctermini_ddG[ctermini_ddG$element_lenght%in% c(4:14),]$ddG_ML_subs,
                        unique_counts)
long <- result[[1]]
long <- long+theme(axis.text.y = element_text(size=16))

## combine
combined_plot <- grid.arrange(long, short, ncol = 2, widths = c(3, 2))

## insertions
unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),], FUN = function(x) length(unique(x)))
result <- termini_plot(ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),],
                        ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),]$ddG_ML_ins,
                        unique_counts)

short <- result[[1]]

unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ctermini_ddG[ctermini_ddG$element_lenght %in% 4:14, ], FUN = function(x) length(unique(x)))

result <- termini_plot(ctermini_ddG[ctermini_ddG$element_lenght%in% c(4:14),],
                        ctermini_ddG[ctermini_ddG$element_lenght%in% c(4:14),]$ddG_ML_ins,
                        unique_counts)
long <- result[[1]]
long <- long+theme(axis.text.y = element_text(size=16))

## combine
combined_plot <- grid.arrange(long, short, ncol = 2, widths = c(3, 2))

## deletions
unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),], FUN = function(x) length(unique(x)))
result <- termini_plot(ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),],
                        ctermini_ddG[ctermini_ddG$element_lenght%in% c(1:3),]$ddG_ML_dels,
                        unique_counts)

short <- result[[1]]

unique_counts <- aggregate(element_no_simple ~ align_to_center, data = ctermini_ddG[ctermini_ddG$element_lenght %in% 4:14, ], FUN = function(x) length(unique(x)))

result <- termini_plot(ctermini_ddG[ctermini_ddG$element_lenght%in% c(4:14),],
                        ctermini_ddG[ctermini_ddG$element_lenght%in% c(4:14),]$ddG_ML_dels,
                        unique_counts)
long <- result[[1]]
long <- long+theme(axis.text.y = element_text(size=16))

## combine
combined_plot <- grid.arrange(long, short, ncol = 2, widths = c(3, 2))

#########################################################################################################
#########################################################################################################
##### c) indel/sub effect on lenght, divided into different elements. 

## helix
helix_df<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("310Helix", "AlphaHelix"),]
unique_counts_helix <- aggregate(element_no_simple ~ element_lenght, data = helix_df, FUN = function(x) length(unique(x)))

plot_lenght_fig3(helix_df,
                 unique_counts_helix)+
  xlab("helix lenght")

## strand
strand_df<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc=="Strand",]
unique_counts_strand <- aggregate(element_no_simple ~ element_lenght, data = strand_df, FUN = function(x) length(unique(x)))

plot_lenght_fig3(strand_df,
                 unique_counts_strand)+
  xlab("strand lenght")

## loop
loop_df<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc=="loop",]
unique_counts_loop <- aggregate(element_no_simple ~ element_lenght, data = loop_df, FUN = function(x) length(unique(x)))

plot_lenght_fig3(loop_df,
                 unique_counts_loop)+
  xlab("loop lenght")+
  theme(legend.position = "right",
        legend.text = element_text(size=18))

#########################################################################################################
#########################################################################################################
##### d) indel/sub effect on structure before or after, divided into structure elements and loop. 
## divide into structure and loop
structure_df<-tsuboyama_nat_doms_all[tsuboyama_nat_doms_all$simple_struc %in% c("AlphaHelix", "Strand", "310Helix"),]

# structure before: structured elements. 
unique_counts <- aggregate(element_no_simple ~ structure_before, data = structure_df, FUN = function(x) length(unique(x)))
custom_order <- c("start", "ntermini", "loop", "310Helix", "AlphaHelix", "Strand")


plot_neighbour_fig3(structure_df,
                    structure_df$structure_before,
                    unique_counts,
                    unique_counts$structure_before)+
  xlab("structure element before")

# structure after: structured elements. 
unique_counts <- aggregate(element_no_simple ~ structure_after, data = structure_df, FUN = function(x) length(unique(x)))
custom_order <- c("ctermini", "loop", "310Helix", "AlphaHelix", "Strand")


plot_neighbour_fig3(structure_df,
                    structure_df$structure_after,
                    unique_counts,
                    unique_counts$structure_after)+
  xlab("structure element after")

## loop
loop_df<-tsuboyama_nat_doms_all[!tsuboyama_nat_doms_all$simple_struc %in% c("AlphaHelix", "Strand", "310Helix"),]

# structure before: structured elements. 
unique_counts <- aggregate(element_no_simple ~ structure_before, data = loop_df, FUN = function(x) length(unique(x)))
custom_order <- c("start", "310Helix", "AlphaHelix", "Strand")


plot_neighbour_fig3(loop_df,
                    loop_df$structure_before,
                    unique_counts,
                    unique_counts$structure_before)+
  xlab("structure element before")

# structure after: structured elements. 
unique_counts <- aggregate(element_no_simple ~ structure_after, data = loop_df, FUN = function(x) length(unique(x)))
custom_order <- c("end", "310Helix", "AlphaHelix", "Strand")


plot_neighbour_fig3(loop_df,
                    loop_df$structure_after,
                    unique_counts,
                    unique_counts$structure_after)+
  xlab("structure element after")







