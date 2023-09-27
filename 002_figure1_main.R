############################################################################
################ script to re-produce the Figure 1  plots #################
############################################################################

## load packages

library(ggplot2)
library(ggridges)

#########################################################################################################
### correlation aPCA and in vitro ddG Figure 1

## extract domain name
validated_domains_in_vitro_ddGs_and_growthrates$domain<-substr(validated_domains_in_vitro_ddGs_and_growthrates$variant_ID,1,14)

## change the domain name
validated_domains_in_vitro_ddGs_and_growthrates<- validated_domains_in_vitro_ddGs_and_growthrates %>% mutate(validated_domains_in_vitro_ddGs_and_growthrates, domain = ifelse(domain == "O75400_PF01846", "FBP11-FF1",
                                                                                                                                                                              ifelse(domain == "P01053_PF00280", "CI2A-PIN1",
                                                                                                                                                                                     ifelse(domain == "P02417_PF01281","BL17-NTL9" ,
                                                                                                                                                                                            ifelse(domain == "P02640_PF02209", "VIL1-HP",
                                                                                                                                                                                                   ifelse(domain == "P0A9X9_PF00313", "CSPA-CSD",
                                                                                                                                                                                                          ifelse(domain == "P32081_PF00313", "CSPB-CSD",
                                                                                                                                                                                                                        ifelse(domain == "P61024_PF01111", "CKS1",NA))))))))

## remove NA
doms_invitro_ddG_vs_gr<-validated_domains_in_vitro_ddGs_and_growthrates[!is.na(validated_domains_in_vitro_ddGs_and_growthrates$domain),]

## set color
my_palette<-"grey30"

## plot
plot_cor_aPCA_invitroddG(doms_invitro_ddG_vs_gr[doms_invitro_ddG_vs_gr$domain=="FBP11-FF1",],my_palette)
plot_cor_aPCA_invitroddG(doms_invitro_ddG_vs_gr[doms_invitro_ddG_vs_gr$domain=="VIL1-HP",],my_palette)
plot_cor_aPCA_invitroddG(doms_invitro_ddG_vs_gr[doms_invitro_ddG_vs_gr$domain=="CSPA-CSD",],my_palette)
plot_cor_aPCA_invitroddG(doms_invitro_ddG_vs_gr[doms_invitro_ddG_vs_gr$domain=="CSPB-CSD",],my_palette)
plot_cor_aPCA_invitroddG(doms_invitro_ddG_vs_gr[doms_invitro_ddG_vs_gr$domain=="CI2A-PIN1",],my_palette)
plot_cor_aPCA_invitroddG(doms_invitro_ddG_vs_gr[doms_invitro_ddG_vs_gr$domain=="BL17-NTL9",],my_palette)
plot_cor_aPCA_invitroddG(doms_invitro_ddG_vs_gr[doms_invitro_ddG_vs_gr$domain=="CKS1",],my_palette)


#########################################################################################################
#### density plots Figure 1

## split main df into indel and substitution+indel df
scaled_variants_aPCA_indels<-scaled_variants_aPCA[scaled_variants_aPCA$mut_type %in% c("insertions", "deletions"),]
scaled_variants_aPCA_substitutions_indels<-scaled_variants_aPCA[scaled_variants_aPCA$mut_type %in% c("insertions", "deletions","substitutions"),]

## FBP11-FF1
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "FBP11-FF1",])

## VIL1-HP
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "VIL1-HP",])

## CSPA-CSD
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "CSPA-CSD",])

## CSPB-CSD
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "CSPB-CSD",])

## CI2A-PIN1
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "CI2A-PIN1",])

## BL17-NTL9
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "BL17-NTL9",])

## CKS1
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "CKS1",])

## PSD95-PDZ3
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "PSD95-PDZ3",])+
  theme(legend.text = element_text(size=18))

## GRB2-SH3
plot_density_single_indels(scaled_variants_aPCA_substitutions_indels[scaled_variants_aPCA_substitutions_indels$domain == "GRB2-SH3",])+
  theme(legend.text = element_text(size=18))

#########################################################################################################
#### heatmaps Figure 1

############ GRB2-SH3: run function to plot the heatmap
result<-plot_heatmap_aPCA_grb2()
grb2_fold_hm_avg <- result
## grb2 wt sequence
WT<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
grb2_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:52,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 5.5,
           size=6)

############ PSD95-PDZ3: run function to plot the heatmap
result<-plot_heatmap_aPCA_pdz3()
pdz3_fold_hm_avg <- result
## pdz3 wt sequence
WT<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
pdz3_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:52,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 5.5,
           size=6)


############ FBP11-FF1: run function to plot the heatmap
result<-plot_heatmap_aPCA_fbp11ff1()
fbp11ff1_fold_hm_avg <- result
WT<-"EAKQAFKELLKEKRVPSNASWEQAMKMIINDPRYSALAKLSEKKQAFNAY"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
fbp11ff1_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:50,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 5.5,
           size=6)

############ CSPA-CSD: run function to plot the heatmap
result<-plot_heatmap_aPCA_cspacsd()
cspacsd_fold_hm_avg <- result
WT<-"MTGIVKWFNADKGFGFITPDDGSKDVFVHFSAIQNDGYKSLDEGQKVSFTIESGAKGPAAGNVTS"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
cspacsd_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:65,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 5.5,
           size=6)

############ CI2A-PIN1: run function to plot the heatmap
result<-plot_heatmap_aPCA_ci2apin1()
ci2apin1_fold_hm_avg <- result
WT<-"KTEWPELVGKSVEEAKKVILQDKPEAQIIVLPVGTIVTMEYRIDRVRLFVDKLDNIAQVPRVG"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
ci2apin1_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:63,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 4.5,
           size=6)

############ BL17-NTL9: run function to plot the heatmap
result<-plot_heatmap_aPCA_bl17ntl9()
bl17ntl9_fold_hm_avg <- result
WT<-"MKVIFLKDVKGKGKKGEIKNVADGYANNFLFKQGLAIEATPANLKAL"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
bl17ntl9_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:47,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 5.5,
           size=6)

############ VIL1-HP: run function to plot the heatmap
result<-plot_heatmap_aPCA_vil1hp()
vil1hp_fold_hm_avg <- result
WT<-"HLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
vil1hp_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:36,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 5.5,
           size=6)

############ CSPB-CSD: run function to plot the heatmap
result<-plot_heatmap_aPCA_cspbcsd()
cspbcsd_fold_hm_avg <- result
WT<-"EGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQAVSFEIVEGNRGPQAANVTK"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
cspbcsd_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:63,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 5.5,
           size=6)


############ CKS1: run function to plot the heatmap
result<-plot_heatmap_aPCA_cks1()
cks1_fold_hm_avg <- result
WT<-"IYYSDKYDDEEFEYRHVMLPKDIAKLVPKTHLMSESEWRNLGVQQSQGWVHYMIHEPEPHILLFRRP"
WT.Vectorised <- sapply(seq(from=1, to=nchar(WT), by=1), function(i) substr(WT, i, i))

## plot stacked heatmap.
cks1_fold_hm_avg+
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(clip = "off") +
  annotate(geom = "text",
           x = 1:67,
           y = "singleDEL",
           label = WT.Vectorised,
           vjust = 5.5,
           size=6)

#########################################################################################################
#### geom_line plots Figure 1

# make dfs with single effects only:
scaled_variants_aPCA_singleINS<-scaled_variants_aPCA %>% filter(type == "singleINS")
scaled_variants_aPCA_singleDEL<-scaled_variants_aPCA %>% filter(type == "singleDEL")

############ FBP11-FF1
## define the breaks based on secondary structure information in stride. 
# the breaks represent where helix/sheet starts and ends
data_breaks <- data.frame(start = c(1,21,34,44),  
                          end = c(13,30,39,46),
                          colors = factor(c(1,1,1,1)))

make_lineplots_single_effects(scaled_variants_aPCA_singleINS  %>% filter(domain == "FBP11-FF1"),
                              scaled_variants_aPCA_singleDEL  %>% filter(domain == "FBP11-FF1"),
                              scaled_variants_aPCA_singleINS  %>% filter(domain == "FBP11-FF1"))

############ VIL1-HP
data_breaks <- data.frame(start = c(4,15,24),  # Create data with breaks
                          end = c(11,20,33),
                          colors = factor(c(1,1,1)))

make_lineplots_single_effects(scaled_variants_aPCA_singleINS  %>% filter(domain == "VIL1-HP"),
                              scaled_variants_aPCA_singleDEL  %>% filter(domain == "VIL1-HP"),
                              scaled_variants_aPCA_singleINS  %>% filter(domain == "VIL1-HP"))

############ CSPA-CSD
data_breaks <- data.frame(start = c(1,14,26,30,45,59),  # Create data with breaks
                          end = c(9,19,29,32,52,65),
                          colors = factor(c(2,2,2,1,2,2)))

make_lineplots_single_effects(scaled_variants_aPCA_singleINS  %>% filter(domain == "CSPA-CSD"),
                              scaled_variants_aPCA_singleDEL  %>% filter(domain == "CSPA-CSD"),
                              scaled_variants_aPCA_singleINS  %>% filter(domain == "CSPA-CSD"))


############ CSPB-CSD
data_breaks <- data.frame(start = c(1,13,24,28,43,55),  # Create data with breaks
                          end = c(8,17,27,30,52,63),
                          colors = factor(c(2,2,2,1,2,2)))

make_lineplots_single_effects(scaled_variants_aPCA_singleINS  %>% filter(domain == "CSPB-CSD"),
                              scaled_variants_aPCA_singleDEL  %>% filter(domain == "CSPB-CSD"),
                              scaled_variants_aPCA_singleINS  %>% filter(domain == "CSPB-CSD"))

############ CI2A-PIN1
data_breaks <- data.frame(start = c(5,10,12,27,45,61),  # Create data with breaks
                          end = c(7,11,22,32,50,62),
                          colors = factor(c(1,2,1,2,2,2)))

make_lineplots_single_effects(scaled_variants_aPCA_singleINS  %>% filter(domain == "CI2A-PIN1"),
                              scaled_variants_aPCA_singleDEL  %>% filter(domain == "CI2A-PIN1"),
                              scaled_variants_aPCA_singleINS  %>% filter(domain == "CI2A-PIN1"))

############ BL17-NTL9
data_breaks <- data.frame(start = c(2,17,23,36,41),  # Create data with breaks
                          end = c(5,20,33,38,47),
                          colors = factor(c(2,2,1,2,1)))

make_lineplots_single_effects(scaled_variants_aPCA_singleINS  %>% filter(domain == "BL17-NTL9"),
                              scaled_variants_aPCA_singleDEL  %>% filter(domain == "BL17-NTL9"),
                              scaled_variants_aPCA_singleINS  %>% filter(domain == "BL17-NTL9"))

############ CKS1
data_breaks <- data.frame(start = c(2,7,12,21,35,50,61),  # Create data with breaks
                          end = c(3,8,18,25,41,53,67),
                          colors = factor(c(2,2,2,1,1,2,2)))

make_lineplots_single_effects(scaled_variants_aPCA_singleINS  %>% filter(domain == "CKS1"),
                              scaled_variants_aPCA_singleDEL  %>% filter(domain == "CKS1"),
                              scaled_variants_aPCA_singleINS  %>% filter(domain == "CKS1"))


############ for GRB2-SH3 and PSD95-PDZ3 we will be adding the substitutions as mean effect/position
# for grb2
grb2_meansub_pos <- calculate_meansub_pos("GRB2-SH3")

# for pdz3
pdz3_meansub_pos <- calculate_meansub_pos("PSD95-PDZ3")

############ GRB2-SH3
data_breaks <- data.frame(start = c(2,23,35,43,49,52),  # Create data with breaks
                          end = c(5,29,40,48,51,52),
                          colors = factor(c(2,2,2,2,1,2)))

make_lineplots_single_effects_subs(scaled_variants_aPCA_singleINS  %>% filter(domain == "GRB2-SH3"),
                                   scaled_variants_aPCA_singleDEL  %>% filter(domain == "GRB2-SH3"),
                                   grb2_meansub_pos,
                                   scaled_variants_aPCA_singleINS  %>% filter(domain == "GRB2-SH3"))

############ PSD95-PDZ3
data_breaks <- data.frame(start = c(2,15,26,36,47),  # Create data with breaks
                          end = c(7,19,31,40,52),
                          colors = factor(c(2,2,2,1,2)))

make_lineplots_single_effects_subs(scaled_variants_aPCA_singleINS  %>% filter(domain == "PSD95-PDZ3"),
                                   scaled_variants_aPCA_singleDEL  %>% filter(domain == "PSD95-PDZ3"),
                                   pdz3_meansub_pos,
                                   scaled_variants_aPCA_singleINS  %>% filter(domain == "PSD95-PDZ3"))

############################################################################
######### box plot of 1-3aa indel effects

#### isolate deletions
data <- scaled_variants_aPCA[scaled_variants_aPCA$mut_type == "deletions",]
# Define custom colors with fading
color1<-"#7AD151FF"
color2<-adjustcolor("#7AD151FF", alpha.f = 0.8)
color3<-adjustcolor("#7AD151FF", alpha.f = 0.4)

custom_colors <- c(color2, color1, color3)

# create the ggplot
ggplot(data, aes(x = factor(type, levels = c("singleDEL", "doubleDEL", "tripleDEL")), y = scaled_fitness, fill = type)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +  # Use the custom colors
  theme_classic() +  # Customize the theme as needed
  ylab("abundance: deletion")+
  scale_x_discrete(limits=c("singleDEL", "doubleDEL", "tripleDEL"),
                   labels=c("1x", "2x", "3x"))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20), legend.position = "none")

## test for signficant differences using wilcox test and bonferroni multiple testing correction
p_values <- pairwise.wilcox.test(data$scaled_fitness, data$type, p.adjust.method = "bonferroni")

# Print the Wilcoxon test results
cat("Wilcoxon Test P-values (Bonferroni-adjusted):\n")
print(p_values$p.value)

#### isolate insertions
data <- scaled_variants_aPCA[scaled_variants_aPCA$mut_type == "insertions",]
# Define custom colors with fading
color1<-"#2A788EFF"
color2<-adjustcolor("#2A788EFF", alpha.f = 0.8)
color3<-adjustcolor("#2A788EFF", alpha.f = 0.4)

custom_colors <- c(color2, color1, color3)

# Create the ggplot
ggplot(data, aes(x = factor(type, levels = c("singleINS", "doubleINS", "tripleINS")), y = scaled_fitness, fill = type)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +  # Use the custom colors
  theme_classic() +  # Customize the theme as needed
  ylab("abundance: CCC ins")+  
  scale_x_discrete(limits=c("singleINS", "doubleINS", "tripleINS"),
                   labels=c("1x", "2x", "3x"))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20), legend.position = "none")

## test for signficant differences using wilcox test and bonferroni multiple testing correction
p_values <- pairwise.wilcox.test(data$scaled_fitness, data$type, p.adjust.method = "bonferroni")

# Print the Wilcoxon test results
cat("Wilcoxon Test P-values (Bonferroni-adjusted):\n")
print(p_values$p.value)







