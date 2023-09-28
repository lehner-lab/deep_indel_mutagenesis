############################################################################
################ script to re-produce the Figure 4  plots #################
############################################################################

############################################################################
### CADD, PROVEAN, ESMv1, DDMut.

###### CADD
###  merge subs
colnames(CADD_subs)[6] <- "mut"
CADD_subs <- merge(CADD_subs,
                         tsuboyama_nat_doms_structure_subs[,c("pdb_name", "mut", "ddG_ML")],
                         by=c("pdb_name", "mut"))

###  merge indels
CADD_insA <- merge(CADD_insA,
                           tsuboyama_nat_doms_all[,c("Pos", "pdb_name", "ddG_ML_dels", "ddG_ML_ins")],
                           by = c("Pos", "pdb_name"))

colnames(CADD_insA)[5:6] <- c("CADD_raw_ins", "PHRED_ins")

CADD_insA <- merge(CADD_insA,
                           CADD_dels[,c("CADD_raw", "PHRED", "pdb_name", "Pos")],
                           by = c("Pos", "pdb_name"))

colnames(CADD_insA)[9:10] <- c("CADD_raw_dels", "PHRED_dels")
CADD_indels <- CADD_insA

## make sure we are looking at the same subset of domains in indel and sub df. 
to_keep <- which(unique(CADD_subs$pdb_name) %in% unique(CADD_indels$pdb_name))
names_to_keep <- unique(CADD_subs$pdb_name)[to_keep]
rows_to_keep <- which(CADD_subs$pdb_name %in% names_to_keep)
CADD_subs <- CADD_subs[rows_to_keep,]

rows_to_keep <- which(CADD_indels$pdb_name %in% unique(CADD_subs$pdb_name))
CADD_indels <- CADD_indels[rows_to_keep,]

## caluclate correlation scores for the combined histograms. 
cor_subs<-data.frame()
for (i in unique(CADD_subs$pdb_name)){
  if (nrow(CADD_subs[CADD_subs$pdb_name == i,])>3){
    temp<-cor.test(CADD_subs[CADD_subs$pdb_name == i,]$ddG_ML,
                   CADD_subs[CADD_subs$pdb_name == i,]$CADD_raw,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_subs<- cor_subs %>% 
      rbind(temp2)
  }
}

cor_dels<-data.frame()
for (i in unique(CADD_indels$pdb_name)){
  if (nrow(CADD_indels[CADD_indels$pdb_name == i,])>3){
    temp<-cor.test(CADD_indels[CADD_indels$pdb_name == i,]$ddG_ML_dels,
                   CADD_indels[CADD_indels$pdb_name == i,]$CADD_raw_dels,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_dels<- cor_dels %>% 
      rbind(temp2)
  }
}

cor_ins<-data.frame()
for (i in unique(CADD_indels$pdb_name)){
  if (nrow(CADD_indels[CADD_indels$pdb_name == i,])>3){
    temp<-cor.test(CADD_indels[CADD_indels$pdb_name == i,]$ddG_ML_ins,
                   CADD_indels[CADD_indels$pdb_name == i,]$CADD_raw_ins,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_ins<- cor_ins %>% 
      rbind(temp2)
  }
}
  

## plot histogram
color1 <- adjustcolor("#414487FF", alpha.f = 0.8) 
color2 <- adjustcolor("#2A788EFF", alpha.f = 0.8) 
color3 <- adjustcolor("#7AD151FF", alpha.f = 0.8)  

hist(cor_subs$R, 
     breaks = 20,
     xlab = "Pearson's R",
     cex.lab = 1.5, cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = color1, # Set the color for the first histogram
     border = "black" # Set the border color for the first histogram
)

# Add the second histogram with a different color
hist(cor_dels$R, 
     breaks = 20,
     add = TRUE, # Add the histogram to the existing plot
     col = color3, # Set the color for the second histogram
     border = "black" # Set the border color for the second histogram
)

# Add the third histogram with a different color
hist(cor_ins$R, 
     breaks = 20,
     add = TRUE, # Add the histogram to the existing plot
     col = color2, # Set the color for the third histogram
     border = "black" # Set the border color for the third histogram
)


legend("topleft", 
       legend = c("Subs", "Dels", "Ins"),
       fill = c(color1, color3, color2),
       border = "black",
       bty = "n"
)

mlv(na.omit(cor_subs$R),method="naive")
mlv(na.omit(cor_ins$R),method="naive")
mlv(na.omit(cor_dels$R),method="naive")

##########################################
###### PROVEAN
###  merge subs
provean_subs <- merge(provean_subs,
                      tsuboyama_nat_doms_structure_subs[,c("pdb_name", "mut", "ddG_ML")],
                      by=c("pdb_name", "mut"))

###  merge indels
provean_insA <- merge(provean_insA,
                      tsuboyama_nat_doms_all[,c("Pos", "pdb_name", "ddG_ML_dels", "ddG_ML_ins")],
                      by = c("Pos", "pdb_name"))

colnames(provean_insA)[4] <- c("score_ins")

provean_insA <- merge(provean_insA,
                      provean_dels[,c("variant", "score", "pdb_name", "Pos")],
                      by = c("Pos", "pdb_name"))

colnames(provean_insA)[8] <- c("score_dels")
provean_indels <- provean_insA[,c("Pos", "pdb_name", "score_ins", "score_dels", "ddG_ML_dels", "ddG_ML_ins")]

## keep only the matching ones in indel_df
rows_to_keep <- which(provean_subs$pdb_name %in% unique(provean_indels$pdb_name))
provean_subs <- provean_subs[rows_to_keep,]

## calculate correlation scores for the combined histograms. 
cor_subs<-data.frame()
for (i in unique(provean_subs$pdb_name)){
  if (nrow(provean_subs[provean_subs$pdb_name == i,])>3){
    temp<-cor.test(provean_subs[provean_subs$pdb_name == i,]$ddG_ML,
                   provean_subs[provean_subs$pdb_name == i,]$score,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_subs<- cor_subs %>% 
      rbind(temp2)
  }
}

cor_dels<-data.frame()
for (i in unique(provean_indels$pdb_name)){
  if (nrow(provean_indels[provean_indels$pdb_name == i,])>3){
    temp<-cor.test(provean_indels[provean_indels$pdb_name == i,]$ddG_ML_dels,
                   provean_indels[provean_indels$pdb_name == i,]$score_dels,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_dels<- cor_dels %>% 
      rbind(temp2)
  }
}

cor_ins<-data.frame()
for (i in unique(provean_indels$pdb_name)){
  if (nrow(provean_indels[provean_indels$pdb_name == i,])>3){
    temp<-cor.test(provean_indels[provean_indels$pdb_name == i,]$ddG_ML_ins,
                   provean_indels[provean_indels$pdb_name == i,]$score_ins,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_ins<- cor_ins %>% 
      rbind(temp2)
  }
}

## make sure the same domains are displayed across all mut types
rows_to_keep <- which(cor_subs$name %in% unique(cor_ins$name))
cor_subs <- cor_subs[rows_to_keep,]

## plot histogram
color1 <- adjustcolor("#414487FF", alpha.f = 0.8) 
color2 <- adjustcolor("#2A788EFF", alpha.f = 0.8) 
color3 <- adjustcolor("#7AD151FF", alpha.f = 0.8)  

hist(-cor_subs$R, 
     breaks = 20,
     xlab = "Pearson's R",
     cex.lab = 1.5, cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = color1, # Set the color for the first histogram
     border = "black" # Set the border color for the first histogram
)

# Add the second histogram with a different color
hist(-cor_dels$R, 
     breaks = 30,
     add = TRUE, # Add the histogram to the existing plot
     col = color3, # Set the color for the second histogram
     border = "black" # Set the border color for the second histogram
)

# Add the third histogram with a different color
hist(-cor_ins$R, 
     breaks = 30,
     add = TRUE, # Add the histogram to the existing plot
     col = color2, # Set the color for the third histogram
     border = "black" # Set the border color for the third histogram
)


legend("topleft", 
       legend = c("Subs", "Dels", "Ins"),
       fill = c(color1, color3, color2),
       border = "black",
       bty = "n"
)

mlv(na.omit(cor_subs$R),method="naive")
mlv(na.omit(cor_ins$R),method="naive")
mlv(na.omit(cor_dels$R),method="naive")

## range of the distribution of correlations for ins
q1 <- quantile(na.omit(-cor_ins$R), 0.25)
q3 <- quantile(na.omit(-cor_ins$R), 0.75)

## range of the distribution of correlations for dels
q1 <- quantile(na.omit(-cor_dels$R), 0.25)
q3 <- quantile(na.omit(-cor_dels$R), 0.75)


##########################################
###### PROVEAN
###  merge subs
esm1v_subs <- esm1v_predictions

esm1v_subs <- merge(esm1v_subs,
                    tsuboyama_nat_doms_structure_subs[,c("pdb_name", "mut", "ddG_ML")],
                    by=c("pdb_name", "mut"))

rows_to_del <- which(is.na(esm1v_subs$ddG_ML))
esm1v_subs <- esm1v_subs[-rows_to_del,]

## caluclate correlation scores for the combined histograms. 
cor_subs<-data.frame()
for (i in unique(esm1v_subs$pdb_name)){
  if (nrow(esm1v_subs[esm1v_subs$pdb_name == i,])>3){
    temp<-cor.test(esm1v_subs[esm1v_subs$pdb_name == i,]$ddG_ML,
                   esm1v_subs[esm1v_subs$pdb_name == i,]$avg_esmv1,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_subs<- cor_subs %>% 
      rbind(temp2)
  }
}

## plot
color1 <- adjustcolor("#414487FF", alpha.f = 0.8) 
hist(-cor_subs$R, 
     breaks = 20,
     xlab = "Pearson's R",
     cex.lab = 1.5, cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = color1, # Set the color for the first histogram
     border = "black" # Set the border color for the first histogram
)

legend("topleft", 
       legend = c("Subs"),
       fill = c(color1),
       border = "black",
       bty = "n"
)

mlv(na.omit(cor_subs$R),method="naive")


##########################################
#### DDMut
###  merge subs
ddmut_subs <- ddmut_predictions
ddmut_subs <- merge(ddmut_subs,
                    tsuboyama_nat_doms_structure_subs[,c("pdb_name", "mut", "ddG_ML")],
                    by=c("pdb_name", "mut"))

## caluclate correlation scores for the combined histograms. 
cor_subs<-data.frame()
for (i in unique(ddmut_subs$pdb_name)){
  if (nrow(ddmut_subs[ddmut_subs$pdb_name == i,])>3){
    temp<-cor.test(ddmut_subs[ddmut_subs$pdb_name == i,]$ddG_ML,
                   ddmut_subs[ddmut_subs$pdb_name == i,]$score_ddmut,
                   method="pearson")
    temp2<-data.frame(name=i,
                      R=temp$estimate,
                      pvalue=temp$p.value)
    cor_subs<- cor_subs %>% 
      rbind(temp2)
  }
}

## plot
color1 <- adjustcolor("#414487FF", alpha.f = 0.8) 
hist(-cor_subs$R, 
     breaks = 20,
     xlab = "Pearson's R",
     cex.lab = 1.5, cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = color1, # Set the color for the first histogram
     border = "black" # Set the border color for the first histogram
)

legend("topleft", 
       legend = c("Subs"),
       fill = c(color1),
       border = "black",
       bty = "n"
)

mlv(na.omit(cor_subs$R),method="naive")

############################################################################
### figure c): the cross-validation for models 1-5 (deletions)

### first I need to process the data to run the cross-validation and different models.
result <- encoded_data_for_predictor()
ddG_ml_encoded <- result[[1]]

## define location where you want to save this .rds
output_location<-output_location_ddGencoded
output_name<-"ddG_ml_encoded.rds"

saveRDS(list(ddG_ml_encoded), 
        file=paste(output_location,
                   output_name, sep="/"))


############################################################################
## for the next step, you also need the ddMut file with predictions for all tsyboyama substitutions: ddmut_prediction_mean.rds
## already loaded in step 00_x

############################################################################
## the next step was run on the CRG cluster using paralisation. 
# the script to run the cross-validation for deletion effect prediction is called: ddG_deletions_published.R

## IF you want to run the predictor yourself load the data for deletions predictor
#output_folder = ""
#result <- load_deletion_prediction_data(output_folder)
#file_list<-result

# IF you want to use the data without running the predictor yourself, use this: file_list_deletions_models
file_list <- file_list_deletions_models

## models use the scripts below to get correlation coefficient for predicted vs observed for all the models
result <- cor_predictor_m1()
model1_del <- result[[1]]

result <- cor_predictor_m2()
model2_del <- result[[1]]

result <- cor_predictor_m3()
model3_del <- result[[1]]

result <- cor_predictor_m4()
model4_del <- result[[1]]

result <- cor_predictor_m5()
model5_del <- result[[1]]

result <- cor_predictor_m1p()
model1p_del <- result[[1]]

result <- cor_predictor_m5p()
model5p_del <- result[[1]]

### plot model 1, 2 and 3
####overlap 3
hist(model1_del$Pearson, 
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5, cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     col = "black", # Set the color for the first histogram
     border = "black" # Set the border color for the first histogram
)


# Add the second histogram with a different color
hist(model2_del$Pearson, 
     breaks = 20,
     add = TRUE, # Add the histogram to the existing plot
     col = adjustcolor("grey40", alpha.f=0.6), # Set the color for the second histogram
     border = "black" # Set the border color for the second histogram
)

# Add the third histogram with a different color
hist(model3_del$Pearson, 
     breaks = 30,
     add = TRUE, # Add the histogram to the existing plot
     col = adjustcolor("white", alpha.f=0.6), # Set the color for the third histogram
     border = "black" # Set the border color for the third histogram
)


legend("left", 
       legend = c("model1", "model2", "model3"),
       fill = c("black", "grey40", "white"),
       border = "black",
       bty = "n", cex = 1.5
)


mlv(na.omit(model1_del$Pearson),method="naive")
mlv(na.omit(model2_del$Pearson),method="naive")
mlv(na.omit(model3_del$Pearson),method="naive")

# plot model 4 and 5
color1 <- "#C0E8A8"  # Transparent light gray
color2 <- adjustcolor("#F5F5F5", alpha.f = 0.4)  # More transparent dark gray

# Plot the histogram for model5
hist(model5_del$Pearson,
     breaks = 30,
     xlab = NULL,
     ylab = NULL,
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     ylim = c(0, 35),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

# Add the histogram for model4 with a different color
hist(model4_del$Pearson,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)

# Add a legend to the left side of the plot
legend("left",
       legend = c("model5", "model4"),
       fill = c(color1, color2),  # Transparent colors
       border = "black",
       bty = "n", cex = 1.5
)

mlv(model4_del$Pearson, method="naive")
mlv(model5_del$Pearson, method="naive")

## model 1p and model5p
color3 <- adjustcolor("yellow3", alpha.f = 0.6)  # More transparent dark gray
hist(model5p_del$Pearson,
     breaks = 30,
     xlab = NULL,
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     ylim = c(0,25),
     ylab = NULL,
     col = color3,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

hist(model1p_del$Pearson,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)


legend("left",
       legend = c("model5p", "model1p"),
       fill = c(color3, color2),  # Transparent colors
       border = "black",
       bty = "n", cex=1.5
)

mlv(na.omit(model1p_del$Pearson),method="naive")
mlv(na.omit(model5p_del$Pearson),method="naive")

############################################################################
### figure c): the cross-validation for models 1-5 (insertions)

############################################################################
## the next step was run on the CRG cluster using paralisation. 
# the script to run the cross-validation for insetion effect prediction is called: ddG_insertions_published.R

## IF you want to run the predictor yourself load the data for insertions predictor
#output_folder = ""
#result <- load_deletion_prediction_data(output_folder)
#file_list<-result

# IF you want to use the data without running the predictor yourself, use this: file_list_insertions_models
file_list <- file_list_insertions_models

## model use the scripts below to get correlation coefficient for predicted vs observed for all the model
result <- cor_predictor_m1()
model1_ins <- result[[1]]

result <- cor_predictor_m2()
model2_ins <- result[[1]]

result <- cor_predictor_m3()
model3_ins <- result[[1]]

result <- cor_predictor_m4()
model4_ins <- result[[1]]

result <- cor_predictor_m5()
model5_ins <- result[[1]]

result <- cor_predictor_m1p()
model1p_ins <- result[[1]]

result <- cor_predictor_m5p()
model5p_ins <- result[[1]]

### plot model 1, 2 and 3
####overlap 3
hist(model1_ins$Pearson, 
     breaks = 30,
     xlab = "Pearson's R",
     cex.lab = 1.5, cex.axis = 1.5,
     main = NULL,
     ylim = c(0,25),
     xlim = c(-1, 1),
     col = "black", # Set the color for the first histogram
     border = "black" # Set the border color for the first histogram
)

# Add the second histogram with a different color
hist(model2_ins$Pearson, 
     breaks = 20,
     add = TRUE, # Add the histogram to the existing plot
     col = adjustcolor("grey40", alpha.f=0.6), # Set the color for the second histogram
     border = "black" # Set the border color for the second histogram
)

# Add the third histogram with a different color
hist(model3_ins$Pearson, 
     breaks = 30,
     add = TRUE, # Add the histogram to the existing plot
     col = adjustcolor("white", alpha.f=0.6), # Set the color for the third histogram
     border = "black" # Set the border color for the third histogram
)


legend("left", 
       legend = c("model1", "model2", "model3"),
       fill = c("black", "grey40", "white"),
       border = "black",
       bty = "n", cex=1.5
)

mlv(na.omit(model1_ins$Pearson),method="naive")
mlv(na.omit(model2_ins$Pearson),method="naive")
mlv(na.omit(model3_ins$Pearson),method="naive")


# plot model 4 and 5
color1 <- "#2A788ECC"  
color2 <- adjustcolor("#F5F5F5", alpha.f = 0.4)  

# Plot the histogram for model5
hist(model5_ins$Pearson,
     breaks = 30,
     xlab = NULL,
     ylab = NULL,
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     ylim = c(0, 25),
     col = color1,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

# Add the histogram for model4 with a different color
hist(model4_ins$Pearson,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)


# Add a legend to the left side of the plot
legend("left",
       legend = c("model4", "model5"),
       fill = c(color2, color1),  # Transparent colors
       border = "black",
       bty = "n", cex=1.5
)

mlv(na.omit(model4_ins$Pearson),method="naive")
mlv(na.omit(model5_ins$Pearson),method="naive")


## model 1p and model5p
color3 <- adjustcolor("yellow3", alpha.f = 0.6)  # More transparent dark gray
hist(model5p_ins$Pearson,
     breaks = 30,
     xlab = NULL,
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = NULL,
     xlim = c(-1, 1),
     ylim = c(0,25),
     ylab = NULL,
     col = color3,  # Set the color for model5 histogram (transparent light gray)
     border = "black"  # Set the border color for model5 histogram
)

hist(model1p_ins$Pearson,
     breaks = 30,
     add = TRUE,  # Add the histogram to the existing plot
     col = color2,  # Set the color for model4 histogram (more transparent dark gray)
     border = "black"  # Set the border color for model4 histogram
)


legend("left",
       legend = c("model5p", "model1p"),
       fill = c(color3, color2),  # Transparent colors
       border = "black",
       bty = "n", cex=1.5
)

mlv(na.omit(model1p_ins$Pearson),method="naive")
mlv(na.omit(model5p_ins$Pearson),method="naive")




