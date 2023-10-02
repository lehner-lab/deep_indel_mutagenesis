############################################################################
######## script to re-produce the Extended Figure 4  plots #################
############################################################################

## load packages
library(viridis)

############################################################################
### a) and b): comparison of CADD and provean preformance on substitutions and indels.

# load the CADD and provean "fresh" from additional_dfs.rds; replace ".." with you folder path
tmp=readRDS(file="~/../additional_dfs.rds")

CADD_subs = tmp[[1]]
CADD_insA = tmp[[2]]
CADD_dels = tmp[[3]]
provean_dels = tmp[[11]]
provean_insA = tmp[[12]]
provean_subs = tmp[[13]]

## compare substitution predictions first: merge PROVEAN and CADD
colnames(CADD_subs)[6] <- "mut"
cadd_provean_subs <- merge(CADD_subs[,c("CADD_raw", "pdb_name", "mut")],
                           provean_subs[,c("mut", "score", "pdb_name")],
                           by = c("mut", "pdb_name"))

colnames(cadd_provean_subs)[4]<-"provean_raw"

## I ran cadd on all possible mutations for the human domains from Tsuboyama et al. regardless if the Tsuboyama et al. data frame had missing ddG for those mutations
## The same for PROVEAN. Therefore we now have more provean vs cadd comparisons than the n=42 presented in main Fig 4. Hence, we need to isolate the domains.
rows_del <- which(!cadd_provean_subs$pdb_name %in% cadd_pdbs_in_tsuboyama)
cadd_provean_subs <- cadd_provean_subs[-rows_del,]

## plot CADD vs PROVEAN substitutions
ggplot(cadd_provean_subs, aes(x = CADD_raw, y = provean_raw)) +
  geom_bin2d(aes(fill = ..count..), bins = 50) +
  scale_fill_viridis(name = "Count", guide = "legend") +
  theme_classic() +
  xlab("CADD score") +
  ylab("PROVEAN score") +
  scale_x_reverse() +
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
        axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

cor.test(-cadd_provean_subs$CADD_raw,
         cadd_provean_subs$provean_raw,
         method="pearson")


### now calculate the per domain correlations between substitution predictions for CADD and PROVEAN 
## add ddG
cadd_provean_subs <- merge(cadd_provean_subs,
                           tsuboyama_nat_doms_structure_subs[,c("mut", "pdb_name","ddG_ML")],
                           by=c("mut", "pdb_name"))

## for CADD
cor_subs_cadd<-data.frame()
for (i in unique(cadd_provean_subs$pdb_name)){
  if (nrow(cadd_provean_subs[cadd_provean_subs$pdb_name == i,])>3){
    temp<-cor.test(cadd_provean_subs[cadd_provean_subs$pdb_name == i,]$ddG_ML,
                   cadd_provean_subs[cadd_provean_subs$pdb_name == i,]$CADD_raw,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_subs_cadd<- cor_subs_cadd %>% 
      rbind(temp2)
  }
}

## for PROVEAN
cor_subs_prov<-data.frame()
for (i in unique(cadd_provean_subs$pdb_name)){
  if (nrow(cadd_provean_subs[cadd_provean_subs$pdb_name == i,])>3){
    temp<-cor.test(cadd_provean_subs[cadd_provean_subs$pdb_name == i,]$ddG_ML,
                   cadd_provean_subs[cadd_provean_subs$pdb_name == i,]$provean_raw,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_subs_prov<- cor_subs_prov %>% 
      rbind(temp2)
  }
}

## plot
breaks <- seq(-1, 1, length.out = 21)
hist(cor_subs_cadd$R, 
     breaks = breaks,
     ylab = NULL,
     xlab = NULL,
     ylim = c(0,14),
     xlim= c(-1,1),
     cex.lab = 1.5, cex.axis = 1.5,
     main = NULL,
     col = "#414487FF",  
     border = "black"  
)

# Add the third histogram with a different color
hist(-cor_subs_prov$R, 
     breaks = breaks,
     ylim = c(0,10),
     xlim= c(-1,1),
     add = TRUE, # Add the histogram to the existing plot
     col = adjustcolor("white", alpha.f = 0.4), 
     border = "black" 
)

legend("topleft", 
       legend = c("CADD", "PROVEAN"),
       fill = c("#414487FF", adjustcolor("white", alpha.f = 0.4)),
       border = "black",
       bty = "n"
)


## compare insertion predictions for CADD and PROVEAN next
cadd_provean_insA <- merge(CADD_insA[,c("CADD_raw", "Pos", "pdb_name")],
                           provean_insA[,c("score", "Pos", "pdb_name")],
                           by = c("Pos", "pdb_name"))

colnames(cadd_provean_insA)[4] <- "provean_raw"

## I ran cadd on all possible mutations for the human domains from Tsuboyama et al. regardless if the Tsuboyama et al. data frame had missing ddG for those mutations
## The same for PROVEAN. Therefore we now have more provean vs cadd comparisons than the n=42 presented in main Fig 4. Hence, we need to isolate the domains.
rows_del <- which(!cadd_provean_insA$pdb_name %in% cadd_pdbs_in_tsuboyama)
cadd_provean_insA <- cadd_provean_insA[-rows_del,]

## plot.
ggplot(cadd_provean_insA, aes(x = CADD_raw, y = provean_raw)) +
  geom_bin2d(aes(fill = ..count..), bins = 50) +
  scale_fill_viridis(name = "Count", guide = "legend") +
  theme_classic() +
  xlab("CADD score") +
  ylab("PROVEAN score") +
  scale_x_reverse() +
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
        axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

cor.test(-cadd_provean_insA$CADD_raw,
         cadd_provean_insA$provean_raw,
         method="pearson")


## compare deletions predictions for CADD and PROVEAN next
cadd_provean_dels <- merge(CADD_dels[,c("CADD_raw", "Pos", "pdb_name")],
                           provean_dels[,c("score", "Pos", "pdb_name")],
                           by = c("Pos", "pdb_name"))
colnames(cadd_provean_dels)[4] <- "provean_raw"

## I ran cadd on all possible mutations for the human domains from Tsuboyama et al. regardless if the Tsuboyama et al. data frame had missing ddG for those mutations
## The same for PROVEAN. Therefore we now have more provean vs cadd comparisons than the n=42 presented in main Fig 4. Hence, we need to isolate the domains.
rows_del <- which(!cadd_provean_dels$pdb_name %in% cadd_pdbs_in_tsuboyama)
cadd_provean_dels <- cadd_provean_dels[-rows_del,]


## plot. 
ggplot(cadd_provean_dels, aes(x = CADD_raw, y = provean_raw)) +
  geom_bin2d(aes(fill = ..count..), bins = 50) +
  scale_fill_viridis(name = "Count", guide = "legend") +
  theme_classic() +
  xlab("CADD score") +
  ylab("PROVEAN score") +
  scale_x_reverse() +
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
        axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 0),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

cor.test(-cadd_provean_dels$CADD_raw,
         cadd_provean_dels$provean_raw,
         method="pearson")

############################################################################
### plot coefficents for deletion models 5 and 5p c)

# If you ran the model yourself you can use the script below to read your .txt files 
#directory_path = "~/.."
#result <- load_coefficents_data_dels(directory_path)
#list_coef <- result
  
# If you want to plot the data for the models we ran, use this file loaded in 000: coefficients_deletions
list_coef <- coefficients_deletions


####### plot top 5 and 5p coefficiens (significant)
## isolate the data for the models
top_m5<-list_coef[["m5"]][list_coef[["m5"]]$s1 !=0,]
top_m5p<-list_coef[["m5p"]][list_coef[["m5p"]]$s1 !=0,]

# plot to see the encoded features
ggplot(data=top_m5,
       aes(x=s1,
           y= reorder(as.character(coef_name), as.numeric(s1))))+
  geom_boxplot()+
  xlab("coefficients")+
  ylab("features")+
  theme_classic()

## use the plot below to rename the features for better understanding of the coefficients 
plot_m5_coeff()


## now the same for 5p
ggplot(data=top_m5p,
       aes(x=s1,
           y= reorder(as.character(coef_name), as.numeric(s1))))+
  geom_boxplot()+
  xlab("coefficients")+
  ylab("features")+
  theme_classic()

## use the plot below to rename the features for better understanding of the coefficients 
plot_m5p_coeff()

############################################################################
### plot coefficents for insertion models 5 and 5p c)

# If you ran the model yourself you can use the script below to read your .txt files 
#directory_path = "~/.."
#result <- load_coefficents_data_dels(directory_path) ## the function for deletions can be used
#list_coef <- result

# If you want to plot the data for the models we ran, use this file loaded in 000: coefficients_insertions
list_coef <- coefficients_insertions


####### plot top 5 and 5p coefficiens (significant)
## isolate the data for the models
top_m5<-list_coef[["m5"]][list_coef[["m5"]]$s1 !=0,]
top_m5p<-list_coef[["m5p"]][list_coef[["m5p"]]$s1 !=0,]


# plot to see the encoded features
ggplot(data=top_m5,
       aes(x=s1,
           y= reorder(as.character(coef_name), as.numeric(s1))))+
  geom_boxplot()+
  xlab("coefficients")+
  ylab("features")+
  theme_classic()

## use the plot below to rename the features for better understanding of the coefficients 
plot_m5_coeff_ins()


## now the same for 5p
ggplot(data=top_m5p,
       aes(x=s1,
           y= reorder(as.character(coef_name), as.numeric(s1))))+
  geom_boxplot()+
  xlab("coefficients")+
  ylab("features")+
  theme_classic()

## rename the features for the final plot. 
plot_m5p_coeff_ins()



