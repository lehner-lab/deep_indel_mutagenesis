############################################################################
####### script to re-produce the Extended Figure 1 plots #############
############################################################################

## load packages
library(GGally)

############################################################################
######### figure with replicate fitness correlations for aPCA

## make a list of sequences that I need to find in the fitness replicate Dimsum file (including indels+substitions+delsubs)
rep_seq <- c(scaled_variants_aPCA$aa_seq)

## load the data from the fitness replicate  DimSum file that contains fitness estimates for all 3 replicates
# load raw data data for 7 doms
input_file_path <- paste(dimsum_file_path, "aPCA_domains_variant_data_merge.RData",sep="")
load(input_file_path)

replicates<-all_variants[,c(2,17:22)]

# load raw data for grb2
input_file_path <- paste(dimsum_file_path, "grb2_fold_fitness_replicates.RData",sep="")
load(input_file_path)

## merge with 7 domains data
replicates<-rbind(replicates,
                  all_variants[,c(2,17:22)])

# load data for pdz3
input_file_path <- paste(dimsum_file_path, "pdz3_fold_fitness_replicates.RData",sep="")
load(input_file_path)

## merge with 7 domains data + GRB2-SH3 data
replicates<-rbind(replicates,
                  all_variants[,c(2,17:22)])

## keep only the final variants corresponding to scaled_variants_aPCA$aa_seq
replicates_abundance<-replicates

rows_to_keep<-which(replicates_abundance$aa_seq %in% rep_seq)
replicates_abundance<-replicates_abundance[rows_to_keep,]

## rename columns, naming them after each biological replicate
cor_matrix <- replicates_abundance[,c(1:4)]
colnames(cor_matrix)[1]<-"ID"
colnames(cor_matrix)[2]<-"S1.1"
colnames(cor_matrix)[3]<-"S1.2"
colnames(cor_matrix)[4]<-"S1.3"
cor_matrix<-as.data.frame(cor_matrix)

## plot replicate correlations
ggpairs(cor_matrix[,c(2:4)],
        columnLabels = c("Rep1", "Rep2", "Rep3"),
        upper = list(continuous = wrap('cor', size = 4)),
        xlab = "protein abundance (9 domains)",
        ylab = "protein abundance (9 domains)")


############################################################################
######### correlations to Tsuboyama et al. data for 5 proteins. 

## load the Tsuboyama data set with insA and del scores for the 5 overlapping domains and make a megred df using the make_df_aPCA_versus_Tsuboyama() function
result <- make_df_aPCA_versus_Tsuboyama()
cor_tsuboyamaVStopolska_realigned <- result[[1]]

## plot the 5 different domains. 
plot_validation("FBP11-FF1")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "FBP11-FF1",]$deltaG_ins,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "FBP11-FF1",]$scaled_fitness_ins,
         method = "pearson")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "FBP11-FF1",]$deltaG_dels,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "FBP11-FF1",]$scaled_fitness_dels,
         method = "pearson")


plot_validation("VIL1-HP")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "VIL1-HP",]$deltaG_ins,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "VIL1-HP",]$scaled_fitness_ins,
         method = "pearson")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "VIL1-HP",]$deltaG_dels,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "VIL1-HP",]$scaled_fitness_dels,
         method = "pearson")


plot_validation("CSPA-CSD")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "CSPA-CSD",]$deltaG_ins,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "CSPA-CSD",]$scaled_fitness_ins,
         method = "pearson")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "CSPA-CSD",]$deltaG_dels,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "CSPA-CSD",]$scaled_fitness_dels,
         method = "pearson")


plot_validation("CSPB-CSD")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "CSPB-CSD",]$deltaG_ins,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "CSPB-CSD",]$scaled_fitness_ins,
         method = "pearson")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "CSPB-CSD",]$deltaG_dels,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "CSPB-CSD",]$scaled_fitness_dels,
         method = "pearson")

plot_validation("BL17-NTL9")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "BL17-NTL9",]$deltaG_ins,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "BL17-NTL9",]$scaled_fitness_ins,
         method = "pearson")
cor.test(cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "BL17-NTL9",]$deltaG_dels,
         cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == "BL17-NTL9",]$scaled_fitness_dels,
         method = "pearson")


############################################################################
######### distributions of 1aa indel effects and calculation of % deleterious

## plot distribution of  1aa indel effects

ggplot(scaled_variants_aPCA[scaled_variants_aPCA$type %in% c("singleINS", "singleDEL")], 
       aes(x = scaled_fitness, y = mut_type))+
  geom_density_ridges(aes(color=mut_type,
                          scale = 1.6,
                          fill=mut_type,
                          alpha="#FF0000A0"))+
  geom_vline(xintercept=1,
             linetype="dotted")+
  geom_vline(xintercept=0, 
             color="red", 
             linetype="dotted")+
  scale_y_discrete(expand = c(0, 0))+
  scale_color_manual(values=c("#7AD151FF","#2A788EFF"),
                     labels = c("1 aa dels", "1 aa ins"))+
  scale_fill_manual(values = c("#7AD151FF", "#2A788EFF"),
                    labels = c("1 aa dels", "1 aa ins"))+
  ylab("density")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        axis.text.x = element_text(size=22),
        legend.title = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        legend.text = element_text(size=22))+
  guides(alpha=FALSE,
         colour = guide_legend(reverse=T),
         fill = guide_legend(reverse=T))

## calculate the precentage deleterious: insertions
# find the mode and the interval of the mode: insertions
scaled_fitness <-scaled_variants_aPCA[scaled_variants_aPCA$type == c( "singleINS")]$scaled_fitness
# estimate the density using kernel density estimation
density_estimate <- density(scaled_fitness)
# find the modes of the density estimate
mode_values <- density_estimate$x[which(density_estimate$y == max(density_estimate$y))]
# calculate the weighted mean of the kernel density estimate
weighted_mean <- sum(density_estimate$x * density_estimate$y) / sum(density_estimate$y)
# calculate the variance of the kernel density estimate using weighted values
weighted_variance <- sum((density_estimate$x - weighted_mean)^2 * density_estimate$y) / sum(density_estimate$y)
# calculate the standard deviation as the square root of the weighted variance
std_dev_at_mode <- sqrt(weighted_variance)

# Define the interval around the mode(s)
interval_width <- 1 * std_dev_at_mode  # Adjust the multiplier for desired confidence

# Calculate the lower and upper bounds of the interval
lower_bound <- mode_values - interval_width
upper_bound <- mode_values + interval_width

# Count the observations within the interval
observations_in_interval <- sum(scaled_fitness >= lower_bound & scaled_fitness <= upper_bound)

# Calculate the total number of observations
total_observations <- length(scaled_fitness)

# Calculate the fraction of observations within the interval
fraction_in_interval <- observations_in_interval / total_observations

# Calculate the percentage by multiplying by 100
percentage_in_interval <- fraction_in_interval * 100

# Print the results
cat("Number of observations in the interval:", observations_in_interval, "\n")
cat("Total number of observations:", total_observations, "\n")
cat("Fraction of observations in the interval:", fraction_in_interval, "\n")
cat("Percentage of observations in the interval:", percentage_in_interval, "%\n")

## calculate the precentage deleterious: deletions
# find the mode and the interval of the mode: deletions
scaled_fitness <-scaled_variants_aPCA[scaled_variants_aPCA$type == c( "singleDEL")]$scaled_fitness
# Estimate the density using kernel density estimation
density_estimate <- density(scaled_fitness)
# Find the modes of the density estimate
mode_values <- density_estimate$x[which(density_estimate$y == max(density_estimate$y))]
# Calculate the weighted mean of the kernel density estimate
weighted_mean <- sum(density_estimate$x * density_estimate$y) / sum(density_estimate$y)
# Calculate the variance of the kernel density estimate using weighted values
weighted_variance <- sum((density_estimate$x - weighted_mean)^2 * density_estimate$y) / sum(density_estimate$y)
# Calculate the standard deviation as the square root of the weighted variance
std_dev_at_mode <- sqrt(weighted_variance)

# Define the interval around the mode(s)
interval_width <- 1 * std_dev_at_mode  # Adjust the multiplier for desired confidence

# Calculate the lower and upper bounds of the interval
lower_bound <- mode_values - interval_width
upper_bound <- mode_values + interval_width

# Count the observations within the interval
observations_in_interval <- sum(scaled_fitness >= lower_bound & scaled_fitness <= upper_bound)

# Calculate the total number of observations
total_observations <- length(scaled_fitness)

# Calculate the fraction of observations within the interval
fraction_in_interval <- observations_in_interval / total_observations

# Calculate the percentage by multiplying by 100
percentage_in_interval <- fraction_in_interval * 100

# Print the results
cat("Number of observations in the interval:", observations_in_interval, "\n")
cat("Total number of observations:", total_observations, "\n")
cat("Fraction of observations in the interval:", fraction_in_interval, "\n")
cat("Percentage of observations in the interval:", percentage_in_interval, "%\n")


############################################################################
######### correlations of 1-3 aa indel effects. 

###correlation between insertions 1-3aa
ins_single<-scaled_variants_aPCA[scaled_variants_aPCA$type == "singleINS",c(1:2,6:7)]
colnames(ins_single)[3:4]<-c("scaled_fitness_singleins", "scaled_sigma_singleins")

ins_double<-scaled_variants_aPCA[scaled_variants_aPCA$type == "doubleINS",c(1:2,6:7)]
colnames(ins_double)[3:4]<-c("scaled_fitness_doubleins", "scaled_sigma_doubleins")

ins_triple<-scaled_variants_aPCA[scaled_variants_aPCA$type == "tripleINS",c(1:2,6:7)]
colnames(ins_triple)[3:4]<-c("scaled_fitness_tripleins", "scaled_sigma_tripleins")

# make a df all insertions
all_doms_instypes<-merge(ins_single,
                         ins_double,
                         by=c("Pos","domain"))

all_doms_instypes<-merge(all_doms_instypes,
                         ins_triple,
                         by=c("Pos","domain"))

###correlation between deletions 1-3aa
dels_single<-scaled_variants_aPCA[scaled_variants_aPCA$type == "singleDEL",c(1:2,6:7)]
colnames(dels_single)[3:4]<-c("scaled_fitness_singledels", "scaled_sigma_singledels")

dels_double<-scaled_variants_aPCA[scaled_variants_aPCA$type == "doubleDEL",c(1:2,6:7)]
colnames(dels_double)[3:4]<-c("scaled_fitness_doubledels", "scaled_sigma_doubledels")

dels_triple<-scaled_variants_aPCA[scaled_variants_aPCA$type == "tripleDEL",c(1:2,6:7)]
colnames(dels_triple)[3:4]<-c("scaled_fitness_tripledels", "scaled_sigma_tripledels")

# make a df all deletions
all_doms_delstypes<-merge(dels_single,
                          dels_double,
                          by=c("Pos","domain"))

all_doms_delstypes<-merge(all_doms_delstypes,
                          dels_triple,
                          by=c("Pos","domain"))

#### merge deletions and insertions
all_doms_insVSdel<-merge(all_doms_instypes,
                         all_doms_delstypes,
                         by=c("Pos", "domain"))

## 1 aa indels correlation plot: scatter
ggplot(data=all_doms_insVSdel,
       aes(x=scaled_fitness_singleins,
           y=scaled_fitness_singledels))+
  geom_pointrange(aes(xmin=scaled_fitness_singleins-scaled_sigma_singleins,
                      xmax=scaled_fitness_singleins+scaled_sigma_singleins,
                      alpha=0.4))+
  geom_pointrange(aes(ymin=scaled_fitness_singledels-scaled_sigma_singledels,
                      ymax=scaled_fitness_singledels+scaled_sigma_singledels,
                      alpha=0.4))+
  geom_point()+
  scale_color_manual(values=c("gray90", "gray50"))+
  theme_classic()+
  #ggtitle("domain scan")+
  ylab("single deletion")+
  xlab("single insertion")+
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.title = element_blank(),
        legend.position = "none")

cor.test(all_doms_insVSdel$scaled_fitness_singledels,
         all_doms_insVSdel$scaled_fitness_singleins,
         method="pearson")

## 1 aa indels correlation plot: correlation of effects per domain
all_doms_insVSdel_1aa_cor<-c()
for (i in unique(all_doms_insVSdel$domain)){
  temp<-cor.test(all_doms_insVSdel[domain==i ,]$scaled_fitness_singleins,
                 all_doms_insVSdel[domain==i,]$scaled_fitness_singledels,method="pearson")
  temp2<-data.frame(domain=i,
                    cor=temp$estimate,
                    pvalue=temp$p.value)
  all_doms_insVSdel_1aa_cor<-all_doms_insVSdel_1aa_cor %>%
    rbind(temp2)
}

## check if cor is signifcant
all_doms_insVSdel_1aa_cor$sig<-""
for (i in 1:nrow(all_doms_insVSdel_1aa_cor)){
  if (all_doms_insVSdel_1aa_cor[i,]$pvalue>0.05){
    all_doms_insVSdel_1aa_cor[i,]$sig<-">0.05"
  }else if (all_doms_insVSdel_1aa_cor[i,]$pvalue<0.05){
    all_doms_insVSdel_1aa_cor[i,]$sig<-"<0.05"
  }
}


all_doms_insVSdel_1aa_cor$domain<- factor(all_doms_insVSdel_1aa_cor$domain, levels = c("CKS1",
                                                                                       "BL17-NTL9",
                                                                                       "CI2A-PIN1",
                                                                                       "GRB2-SH3",
                                                                                       "PSD95-PDZ3",
                                                                                       "CSPB-CSD",
                                                                                       "CSPA-CSD",
                                                                                       "VIL1-HP",
                                                                                       "FBP11-FF1"), ordered = TRUE)
# plot
ggplot(data=all_doms_insVSdel_1aa_cor,
       aes(x=cor,
           y=domain,
           fill=sig))+
  geom_vline(xintercept=0.5,
             linetype="dotted",
             color="darkred")+
  geom_vline(xintercept=-0.5,
             linetype="dotted",
             color="darkred")+
  geom_col()+
  scale_fill_manual(values = c("grey10", "grey70"))+
  xlab("Pearson's R")+
  xlim(-1,1)+
  theme_classic()+
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y =  element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  guides(fill=guide_legend(title="pvalue"))



## 2 aa indels correlation plot: scatter
ggplot(data=all_doms_insVSdel,
       aes(x=scaled_fitness_doubleins,
           y=scaled_fitness_doubledels))+
  geom_pointrange(aes(xmin=scaled_fitness_doubleins-scaled_sigma_doubleins,
                      xmax=scaled_fitness_doubleins+scaled_sigma_doubleins,
                      alpha=0.4))+
  geom_pointrange(aes(ymin=scaled_fitness_doubledels-scaled_sigma_doubledels,
                      ymax=scaled_fitness_doubledels+scaled_sigma_doubledels,
                      alpha=0.4))+
  geom_point()+
  scale_color_manual(values=c("gray90", "gray50"))+
  theme_classic()+
  #ggtitle("domain scan")+
  ylab("double deletion")+
  xlab("double insertion")+
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.title = element_blank(),
        legend.position = "none")

cor.test(all_doms_insVSdel$scaled_fitness_doubledels,
         all_doms_insVSdel$scaled_fitness_doubleins,
         method="pearson")

## 2 aa indels correlation plot: correlation per domain
all_doms_insVSdel_2aa_cor<-c()
for (i in unique(all_doms_insVSdel$domain)){
  temp<-cor.test(all_doms_insVSdel[domain==i ,]$scaled_fitness_doubleins,
                 all_doms_insVSdel[domain==i,]$scaled_fitness_doubledels,method="pearson")
  temp2<-data.frame(domain=i,
                    cor=temp$estimate,
                    pvalue=temp$p.value)
  all_doms_insVSdel_2aa_cor<-all_doms_insVSdel_2aa_cor %>%
    rbind(temp2)
}

## check if cor is signifcant
all_doms_insVSdel_2aa_cor$sig<-""
for (i in 1:nrow(all_doms_insVSdel_2aa_cor)){
  if (all_doms_insVSdel_2aa_cor[i,]$pvalue>0.05){
    all_doms_insVSdel_2aa_cor[i,]$sig<-">0.05"
  }else if (all_doms_insVSdel_2aa_cor[i,]$pvalue<0.05){
    all_doms_insVSdel_2aa_cor[i,]$sig<-"<0.05"
  }
}


all_doms_insVSdel_2aa_cor$domain<- factor(all_doms_insVSdel_2aa_cor$domain, levels = c("CKS1",
                                                                                       "BL17-NTL9",
                                                                                       "CI2A-PIN1",
                                                                                       "GRB2-SH3",
                                                                                       "PSD95-PDZ3",
                                                                                       "CSPB-CSD",
                                                                                       "CSPA-CSD",
                                                                                       "VIL1-HP",
                                                                                       "FBP11-FF1"), ordered = TRUE)
# plot
ggplot(data=all_doms_insVSdel_2aa_cor,
       aes(x=cor,
           y=domain,
           fill=sig))+
  geom_vline(xintercept=0.5,
             linetype="dotted",
             color="darkred")+
  geom_vline(xintercept=-0.5,
             linetype="dotted",
             color="darkred")+
  geom_col()+
  scale_fill_manual(values = c("grey10", "grey70"))+
  xlab("Pearson's R")+
  xlim(-1,1)+
  theme_classic()+
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y =  element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  guides(fill=guide_legend(title="pvalue"))

## 3 aa indels correlation plot: scatter
ggplot(data=all_doms_insVSdel,
       aes(x=scaled_fitness_tripleins,
           y=scaled_fitness_tripledels))+
  geom_pointrange(aes(xmin=scaled_fitness_tripleins-scaled_sigma_tripleins,
                      xmax=scaled_fitness_tripleins+scaled_sigma_tripleins,
                      alpha=0.4))+
  geom_pointrange(aes(ymin=scaled_fitness_tripledels-scaled_sigma_tripledels,
                      ymax=scaled_fitness_tripledels+scaled_sigma_tripledels,
                      alpha=0.4))+
  geom_point()+
  scale_color_manual(values=c("gray90", "gray50"))+
  theme_classic()+
  #ggtitle("domain scan")+
  ylab("triple deletion")+
  xlab("triple insertion")+
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.title = element_blank(),
        legend.position = "none")

cor.test(all_doms_insVSdel$scaled_fitness_tripledels,
         all_doms_insVSdel$scaled_fitness_tripleins,
         method="pearson")

## 3 aa indels correlation plot: correlation per domain
all_doms_insVSdel_3aa_cor<-c()
for (i in unique(all_doms_insVSdel$domain)){
  temp<-cor.test(all_doms_insVSdel[domain==i ,]$scaled_fitness_tripleins,
                 all_doms_insVSdel[domain==i,]$scaled_fitness_tripledels,method="pearson")
  temp2<-data.frame(domain=i,
                    cor=temp$estimate,
                    pvalue=temp$p.value)
  all_doms_insVSdel_3aa_cor<-all_doms_insVSdel_3aa_cor %>%
    rbind(temp2)
}

## check if cor is signifcant
all_doms_insVSdel_3aa_cor$sig<-""
for (i in 1:nrow(all_doms_insVSdel_3aa_cor)){
  if (all_doms_insVSdel_3aa_cor[i,]$pvalue>0.05){
    all_doms_insVSdel_3aa_cor[i,]$sig<-">0.05"
  }else if (all_doms_insVSdel_3aa_cor[i,]$pvalue<0.05){
    all_doms_insVSdel_3aa_cor[i,]$sig<-"<0.05"
  }
}


all_doms_insVSdel_3aa_cor$domain<- factor(all_doms_insVSdel_3aa_cor$domain, levels = c("CKS1",
                                                                                       "BL17-NTL9",
                                                                                       "CI2A-PIN1",
                                                                                       "GRB2-SH3",
                                                                                       "PSD95-PDZ3",
                                                                                       "CSPB-CSD",
                                                                                       "CSPA-CSD",
                                                                                       "VIL1-HP",
                                                                                       "FBP11-FF1"), ordered = TRUE)
# plot
ggplot(data=all_doms_insVSdel_3aa_cor,
       aes(x=cor,
           y=domain,
           fill=sig))+
  geom_vline(xintercept=0.5,
             linetype="dotted",
             color="darkred")+
  geom_vline(xintercept=-0.5,
             linetype="dotted",
             color="darkred")+
  geom_col()+
  scale_fill_manual(values = c("grey10", "grey70"))+
  xlab("Pearson's R")+
  xlim(-1,1)+
  theme_classic()+
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y =  element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  guides(fill=guide_legend(title="pvalue"))


############################################################################
######### correlations of 1-3 aa indel effects: insertions only

## 1 vs 2 aa insertions correlation plot: scatter
ggplot(data=all_doms_insVSdel,
       aes(x=scaled_fitness_singleins,
           y=scaled_fitness_doubleins))+
  geom_pointrange(aes(xmin=scaled_fitness_singleins-scaled_sigma_singleins,
                      xmax=scaled_fitness_singleins+scaled_sigma_singleins,
                      alpha=0.4))+
  geom_pointrange(aes(ymin=scaled_fitness_doubleins-scaled_sigma_doubleins,
                      ymax=scaled_fitness_doubleins+scaled_sigma_doubleins,
                      alpha=0.4))+
  geom_point()+
  scale_color_manual(values=c("gray90", "gray50"))+
  theme_classic()+
  #ggtitle("domain scan")+
  ylab("single insertion")+
  xlab("double insertion")+
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.title = element_blank(),
        legend.position = "none")

cor.test(all_doms_insVSdel$scaled_fitness_singleins,
         all_doms_insVSdel$scaled_fitness_doubleins,
         method="pearson")

## 1 vs 2 aa insertions correlation plot: correlation per domain
all_doms_insVSdel_1vs2aa_cor<-c()
for (i in unique(all_doms_insVSdel$domain)){
  temp<-cor.test(all_doms_insVSdel[domain==i ,]$scaled_fitness_singleins,
                 all_doms_insVSdel[domain==i,]$scaled_fitness_doubleins,method="pearson")
  temp2<-data.frame(domain=i,
                    cor=temp$estimate,
                    pvalue=temp$p.value)
  all_doms_insVSdel_1vs2aa_cor<-all_doms_insVSdel_1vs2aa_cor %>%
    rbind(temp2)
}

## check if cor is signifcant
all_doms_insVSdel_1vs2aa_cor$sig<-""
for (i in 1:nrow(all_doms_insVSdel_1vs2aa_cor)){
  if (all_doms_insVSdel_1vs2aa_cor[i,]$pvalue>0.05){
    all_doms_insVSdel_1vs2aa_cor[i,]$sig<-">0.05"
  }else if (all_doms_insVSdel_1vs2aa_cor[i,]$pvalue<0.05){
    all_doms_insVSdel_1vs2aa_cor[i,]$sig<-"<0.05"
  }
}


all_doms_insVSdel_1vs2aa_cor$domain<- factor(all_doms_insVSdel_1vs2aa_cor$domain, levels = c("CKS1",
                                                                                       "BL17-NTL9",
                                                                                       "CI2A-PIN1",
                                                                                       "GRB2-SH3",
                                                                                       "PSD95-PDZ3",
                                                                                       "CSPB-CSD",
                                                                                       "CSPA-CSD",
                                                                                       "VIL1-HP",
                                                                                       "FBP11-FF1"), ordered = TRUE)
# plot
ggplot(data=all_doms_insVSdel_1vs2aa_cor,
       aes(x=cor,
           y=domain,
           fill=sig))+
  geom_vline(xintercept=0.5,
             linetype="dotted",
             color="darkred")+
  geom_vline(xintercept=-0.5,
             linetype="dotted",
             color="darkred")+
  geom_col()+
  scale_fill_manual(values = c("grey10", "grey70"))+
  xlab("Pearson's R")+
  xlim(-1,1)+
  theme_classic()+
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y =  element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  guides(fill=guide_legend(title="pvalue"))

## 1 vs 3 aa insertions correlation plot: scatter
ggplot(data=all_doms_insVSdel,
       aes(x=scaled_fitness_singleins,
           y=scaled_fitness_tripleins))+
  geom_pointrange(aes(xmin=scaled_fitness_singleins-scaled_sigma_singleins,
                      xmax=scaled_fitness_singleins+scaled_sigma_singleins,
                      alpha=0.4))+
  geom_pointrange(aes(ymin=scaled_fitness_tripleins-scaled_sigma_tripleins,
                      ymax=scaled_fitness_tripleins+scaled_sigma_tripleins,
                      alpha=0.4))+
  geom_point()+
  scale_color_manual(values=c("gray90", "gray50"))+
  theme_classic()+
  #ggtitle("domain scan")+
  ylab("single insertion")+
  xlab("triple insertion")+
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.title = element_blank(),
        legend.position = "none")

cor.test(all_doms_insVSdel$scaled_fitness_singleins,
         all_doms_insVSdel$scaled_fitness_tripleins,
         method="pearson")

## 1 vs 3 aa insertions correlation plot: correlation per domain
all_doms_insVSdel_1vs3aa_cor<-c()
for (i in unique(all_doms_insVSdel$domain)){
  temp<-cor.test(all_doms_insVSdel[domain==i ,]$scaled_fitness_singleins,
                 all_doms_insVSdel[domain==i,]$scaled_fitness_tripleins,method="pearson")
  temp2<-data.frame(domain=i,
                    cor=temp$estimate,
                    pvalue=temp$p.value)
  all_doms_insVSdel_1vs3aa_cor<-all_doms_insVSdel_1vs3aa_cor %>%
    rbind(temp2)
}

## check if cor is signifcant
all_doms_insVSdel_1vs3aa_cor$sig<-""
for (i in 1:nrow(all_doms_insVSdel_1vs3aa_cor)){
  if (all_doms_insVSdel_1vs3aa_cor[i,]$pvalue>0.05){
    all_doms_insVSdel_1vs3aa_cor[i,]$sig<-">0.05"
  }else if (all_doms_insVSdel_1vs3aa_cor[i,]$pvalue<0.05){
    all_doms_insVSdel_1vs3aa_cor[i,]$sig<-"<0.05"
  }
}


all_doms_insVSdel_1vs3aa_cor$domain<- factor(all_doms_insVSdel_1vs3aa_cor$domain, levels = c("CKS1",
                                                                                             "BL17-NTL9",
                                                                                             "CI2A-PIN1",
                                                                                             "GRB2-SH3",
                                                                                             "PSD95-PDZ3",
                                                                                             "CSPB-CSD",
                                                                                             "CSPA-CSD",
                                                                                             "VIL1-HP",
                                                                                             "FBP11-FF1"), ordered = TRUE)
# plot
ggplot(data=all_doms_insVSdel_1vs3aa_cor,
       aes(x=cor,
           y=domain,
           fill=sig))+
  geom_vline(xintercept=0.5,
             linetype="dotted",
             color="darkred")+
  geom_vline(xintercept=-0.5,
             linetype="dotted",
             color="darkred")+
  geom_col()+
  scale_fill_manual(values = c("grey10", "grey70"))+
  xlab("Pearson's R")+
  xlim(-1,1)+
  theme_classic()+
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y =  element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  guides(fill=guide_legend(title="pvalue"))


############################################################################
######### correlations of 1-3 aa indel effects: deletions only

## 1 vs 2 aa deletions correlation plot: scatter
ggplot(data=all_doms_insVSdel,
       aes(x=scaled_fitness_singledels,
           y=scaled_fitness_doubledels))+
  geom_pointrange(aes(xmin=scaled_fitness_singledels-scaled_sigma_singledels,
                      xmax=scaled_fitness_singledels+scaled_sigma_singledels,
                      alpha=0.4))+
  geom_pointrange(aes(ymin=scaled_fitness_doubledels-scaled_sigma_doubledels,
                      ymax=scaled_fitness_doubledels+scaled_sigma_doubledels,
                      alpha=0.4))+
  geom_point()+
  scale_color_manual(values=c("gray90", "gray50"))+
  theme_classic()+
  #ggtitle("domain scan")+
  ylab("single deletion")+
  xlab("double deletion")+
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.title = element_blank(),
        legend.position = "none")

cor.test(all_doms_insVSdel$scaled_fitness_singledels,
         all_doms_insVSdel$scaled_fitness_doubledels,
         method="pearson")

## 1 vs 2 aa deletions correlation plot: correlation per domain
all_doms_insVSdel_1vs2aa_cor<-c()
for (i in unique(all_doms_insVSdel$domain)){
  temp<-cor.test(all_doms_insVSdel[domain==i ,]$scaled_fitness_singledels,
                 all_doms_insVSdel[domain==i,]$scaled_fitness_doubledels,method="pearson")
  temp2<-data.frame(domain=i,
                    cor=temp$estimate,
                    pvalue=temp$p.value)
  all_doms_insVSdel_1vs2aa_cor<-all_doms_insVSdel_1vs2aa_cor %>%
    rbind(temp2)
}

## check if cor is signifcant
all_doms_insVSdel_1vs2aa_cor$sig<-""
for (i in 1:nrow(all_doms_insVSdel_1vs2aa_cor)){
  if (all_doms_insVSdel_1vs2aa_cor[i,]$pvalue>0.05){
    all_doms_insVSdel_1vs2aa_cor[i,]$sig<-">0.05"
  }else if (all_doms_insVSdel_1vs2aa_cor[i,]$pvalue<0.05){
    all_doms_insVSdel_1vs2aa_cor[i,]$sig<-"<0.05"
  }
}


all_doms_insVSdel_1vs2aa_cor$domain<- factor(all_doms_insVSdel_1vs2aa_cor$domain, levels = c("CKS1",
                                                                                             "BL17-NTL9",
                                                                                             "CI2A-PIN1",
                                                                                             "GRB2-SH3",
                                                                                             "PSD95-PDZ3",
                                                                                             "CSPB-CSD",
                                                                                             "CSPA-CSD",
                                                                                             "VIL1-HP",
                                                                                             "FBP11-FF1"), ordered = TRUE)
# plot
ggplot(data=all_doms_insVSdel_1vs2aa_cor,
       aes(x=cor,
           y=domain,
           fill=sig))+
  geom_vline(xintercept=0.5,
             linetype="dotted",
             color="darkred")+
  geom_vline(xintercept=-0.5,
             linetype="dotted",
             color="darkred")+
  geom_col()+
  scale_fill_manual(values = c("grey10", "grey70"))+
  xlab("Pearson's R")+
  xlim(-1,1)+
  theme_classic()+
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y =  element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  guides(fill=guide_legend(title="pvalue"))

## 1 vs 3 aa deletions correlation plot: scatter
ggplot(data=all_doms_insVSdel,
       aes(x=scaled_fitness_singledels,
           y=scaled_fitness_tripledels))+
  geom_pointrange(aes(xmin=scaled_fitness_singledels-scaled_sigma_singledels,
                      xmax=scaled_fitness_singledels+scaled_sigma_singledels,
                      alpha=0.4))+
  geom_pointrange(aes(ymin=scaled_fitness_tripledels-scaled_sigma_tripledels,
                      ymax=scaled_fitness_tripledels+scaled_sigma_tripledels,
                      alpha=0.4))+
  geom_point()+
  scale_color_manual(values=c("gray90", "gray50"))+
  theme_classic()+
  #ggtitle("domain scan")+
  ylab("single deletion")+
  xlab("triple deletion")+
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.title = element_blank(),
        legend.position = "none")

cor.test(all_doms_insVSdel$scaled_fitness_singledels,
         all_doms_insVSdel$scaled_fitness_tripledels,
         method="pearson")

## 1 vs 3 aa deletions correlation plot: correlation per domain
all_doms_insVSdel_1vs3aa_cor<-c()
for (i in unique(all_doms_insVSdel$domain)){
  temp<-cor.test(all_doms_insVSdel[domain==i ,]$scaled_fitness_singledels,
                 all_doms_insVSdel[domain==i,]$scaled_fitness_tripledels,method="pearson")
  temp2<-data.frame(domain=i,
                    cor=temp$estimate,
                    pvalue=temp$p.value)
  all_doms_insVSdel_1vs3aa_cor<-all_doms_insVSdel_1vs3aa_cor %>%
    rbind(temp2)
}

## check if cor is signifcant
all_doms_insVSdel_1vs3aa_cor$sig<-""
for (i in 1:nrow(all_doms_insVSdel_1vs3aa_cor)){
  if (all_doms_insVSdel_1vs3aa_cor[i,]$pvalue>0.05){
    all_doms_insVSdel_1vs3aa_cor[i,]$sig<-">0.05"
  }else if (all_doms_insVSdel_1vs3aa_cor[i,]$pvalue<0.05){
    all_doms_insVSdel_1vs3aa_cor[i,]$sig<-"<0.05"
  }
}


all_doms_insVSdel_1vs3aa_cor$domain<- factor(all_doms_insVSdel_1vs3aa_cor$domain, levels = c("CKS1",
                                                                                             "BL17-NTL9",
                                                                                             "CI2A-PIN1",
                                                                                             "GRB2-SH3",
                                                                                             "PSD95-PDZ3",
                                                                                             "CSPB-CSD",
                                                                                             "CSPA-CSD",
                                                                                             "VIL1-HP",
                                                                                             "FBP11-FF1"), ordered = TRUE)
# plot
ggplot(data=all_doms_insVSdel_1vs3aa_cor,
       aes(x=cor,
           y=domain,
           fill=sig))+
  geom_vline(xintercept=0.5,
             linetype="dotted",
             color="darkred")+
  geom_vline(xintercept=-0.5,
             linetype="dotted",
             color="darkred")+
  geom_col()+
  scale_fill_manual(values = c("grey10", "grey70"))+
  xlab("Pearson's R")+
  xlim(-1,1)+
  theme_classic()+
  theme(axis.title.y  = element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y =  element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  guides(fill=guide_legend(title="pvalue"))



