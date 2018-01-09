####################################################
### Beta diversity measures, graphs, and stats ####
####################################################

setwd("/Users/Becca/Dropbox/Microbiome_analysis/mal/")

library('phyloseq'); packageVersion('phyloseq')
library("ggplot2")
library('vegan'); packageVersion('vegan')
library("ggsignif")
library("plyr")
library("car")
library("ape")
library("cowplot")

# First run import_qiime_data_to_phyloseq.R to get the qd object in the R environment


##
# Rarefy the data to an even sampling depth, based on the sample with the smallest # of reads
min_lib <- min(sample_sums(qiimedata))
qd <- rarefy_even_depth(qiimedata,sample.size = min_lib, verbose = FALSE, replace = TRUE)

# Check if your tree is rooted
is.rooted(phy_tree(qd))
# Should say TRUE

## Data qd is now ready for use in phyloseq

#--------------------------------------------------------------------------

##
# Calculate all four distance matrices: weighted unifrac, unweighted unifrac, bray curtis, binary jaccard
qd_wu <- distance(qd, method = "wunifrac")
qd_un <- distance(qd, method = "unifrac")
qd_bc <- distance(qd, method = "bray")
qd_bj <- distance(qd, method = "jaccard")
head(qd_bc)

# Ordinations for PCoAs
ord_wu <- ordinate(qd, "PCoA", distance = qd_wu)
ord_un <- ordinate(qd, "PCoA", distance = qd_un)
ord_bc <- ordinate(qd, "PCoA", distance = qd_bc)
ord_bj <- ordinate(qd, "PCoA", distance = qd_bj)

# Example code for plotting an ordination
# plot_ordination(qd, ord_wu, color="interaction", title = "Bray Curtis")
# po1 <- plot_ordination(qd, ord_wu, color="temp", title = "Weighted Unifrac")
# po2 <- plot_ordination(qd, ord_un, color="temp", title = "Unweighted Unifrac")
# po3 <- plot_ordination(qd, ord_bc, color="temp", title = "Bray Curtis")
# po4 <- plot_ordination(qd, ord_bj, color="temp", title = "Binary Jaccard")
# multiplot(po1,po2,po3,po4, cols=2) 
# multiplot is a function you have to run separately before use
##


##
# PERMANOVA's with Adonis
# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(qd))

# I can do Adonis test for between group diversity, 
# The independent variables can be tested one at a time or as a formula
# example
# adonis(qd_bj ~ interaction, data = sampledf)
adonis(qd_bc ~ nutrient*corallivory*temp+colony+tank, data = sampledf)
adonis(qd_bj ~ nutrient*corallivory*temp+colony+tank, data = sampledf) 
adonis(qd_un ~ nutrient*corallivory*temp+colony+tank, data = sampledf) 
adonis(qd_wu ~ nutrient*corallivory*temp+colony+tank, data = sampledf) 
##


##
# Betadisper or the homogeneity of dispersion test
# change the sampledf$X to get each category
betabc <- betadisper(qd_bc, sampledf$interaction, bias.adjust=TRUE)
betabj <- betadisper(qd_bj, sampledf$interaction, bias.adjust=TRUE)
betaun <- betadisper(qd_un, sampledf$interaction, bias.adjust=TRUE)
betawu <- betadisper(qd_wu, sampledf$interaction, bias.adjust=TRUE)

beta5 <- betadisper(qd_un, sampledf$nutrient)
beta6 <- betadisper(qd_un, sampledf$corallivory)
beta7 <- betadisper(qd_bj, sampledf$temp)
beta8 <- betadisper(qd_un, sampledf$tank)
beta9 <- betadisper(qd_un, sampledf$colony)

# Add the pairwise=TRUE to get a pairwise comparision of all the results
permutest(betabc, pairwise=TRUE)
permutest(betabj, pairwise=TRUE)
permutest(betaun, pairwise=TRUE)
permutest(betauw, pairwise=TRUE)
# These pvalues are used to designate significance on the graphs of all twelve treatments (plots 1,2,3)

# Extra stuff
# Perform test
#anova(beta4)
#permutest.betadisper(qd_bc, sampledf$interaction)
# Tukey's Honest Significant Differences
#beta.HSD <- TukeyHSD(beta)
#plot(beta.HSD)
# Plot the groups and distances to centroids on the first two PCoA axes
#plot(beta)
# Draw a boxplot of the distances to centroid for each group
# boxplot(beta)
# List of the coordinates for each sample in coordinate space
# scores(beta, display = c("sites","centroid"))

# Gives the distances of each sample to its centroid by interaction
distbc <- betabc$distances
distbj <- betabj$distances
distun <- betaun$distances
distwu <- betawu$distances

# IMPORTANT STEP: Merge the distance to centroid for each distance measure by SampleID
# Put them into the map file for use in the graphs

# Makes a dataframe of the distance to centroid of one distance measure
# Just necessary to get the SampleID labels
df_bc <- data.frame(betabc$distances) 
# Extracts the SampleIDs from the dataframe
SampleID <- row.names(df_bc) 
# makes a new dataframe with the specified columnes
distances <- data.frame(SampleID, distbc, distbj, distun, distwu) 
# Makes metadata into a df to work with
s <- data.frame(sample_data(qd)) 
# merges metadata df and distances df
beta_data <- merge(distances, s, by = "SampleID") 
##


##
# Linear regression and anova with distance to centroid as the dependent variable
reg1 <- lm(distwu ~ nutrient*temp*corallivory + tank + colony, data = beta_data)
reg2 <- lm(distun ~ nutrient*temp*corallivory + tank + colony, data = beta_data)
reg3 <- lm(distbc ~ nutrient*temp*corallivory + tank + colony, data = beta_data)
reg4 <- lm(distbj ~ nutrient*temp*corallivory + tank + colony, data = beta_data)

summary(reg1)
anova(reg1)
summary(reg2)
anova(reg2)
summary(reg3)
anova(reg3)
summary(reg4)
anova(reg4)
# These pvalues are used to designate significance on the graphs by treatment (plot 4)
##


## Now we have a .csv file to make the boxplots from


##
# Plots to include:
# 1 # Four plots of all twelve treatments, one for each of the four distance matrices:
# Bray Curtis (bc), Binary Jaccard (bj), Weighted Unifrac (wu), Unweighted Unifrac (un)

# 2 # Single plot of twelve treatments, based on Weighted Unifrac

# 3 # Plot of single stressors: control, high temp, NO3, NH4, scarred

# 4 # Three plots by nutrients, temperature, and scarring
# End of plot discriptions, begin plot code

# 1 # Plot of four distance matrices with distance to centroid.
# Beta diversity distance to centroid boxplot of all the interactions
breakss <- c("C.26.I","A.26.I","N.26.S","A.26.S","A.30.S","N.30.S","C.30.S","A.30.I",
            "N.26.I","N.30.I","C.26.S","C.30.I")
labelss <- c("control","NH3+","NO3-, scarred", "NH3+, scarred","NH3+, scarred, 29","NO3-, scarred, 29","scarred, 29","NH3+, 29",
            "NO3-","NO3-, 29","scarred","29")

# Significance codes determined by betadisper
sigwu <- c("A", "A,B,C", "A,B,C", "B", "B,C", "B,C", "B,C", "B,C", "B,C", "C", "B,C", "B,C")
sigbc <- c("A", "B,C,D", "B", "B,C,D", "B,C,D", "B,C", "A,B,C,D", "B,C,D", "B", "C,D", "D", "D")
sigbj <- c("A", "B,C", "B", "B,C", "B", "B", "B,C", "A,B,C", "B,C", "C", "C", "B,C")
# All four distances are included in this data sheet
beta_data <- read.csv(file = "/Users/Becca/Dropbox/Microbiome_analysis/mal/beta/beta_dist_map.csv")

# Code for a single graph ordered from lowest to highest mean
p1 <- ggplot(beta_data, aes(x=reorder(interaction, distwu,FUN = median), y=distwu, fill=interaction)) +
  geom_boxplot() + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x = element_blank()) + scale_x_discrete(breaks=breakss, labels=labelss) + ylim(0,1.0) + 
  stat_summary(geom = 'text', label = sig1, fun.y = max, vjust = -1, size = 3) +
  labs(title="Weighted Unifrac", y = "Distance to centroid") 

p2 <- ggplot(beta_data, aes(x=reorder(interaction, distbc, FUN = median), y=distbc, fill=interaction)) +
  geom_boxplot()+ theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(breaks=breakss, labels=labelss) + 
  stat_summary(geom = 'text', label = sigbc, fun.y = max, vjust = -1, size = 3) +
  labs(title="Bray Curtis") + ylim(0,1.0)

p3 <- ggplot(beta_data, aes(x=reorder(interaction, distun,FUN = median), y=distun, fill=interaction)) +
  geom_boxplot()+ theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(breaks=breakss, labels=labelss) + ylim(0,1.0) + 
  labs(title="Unweighted Unifrac", y = "Distance to centroid", x="Treatment interaction")

p4 <- ggplot(beta_data, aes(x=reorder(interaction, distbj,FUN = median), y=distbj, fill=interaction)) +
  geom_boxplot()+ theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.y = element_blank()) + scale_x_discrete(breaks=breakss, labels=labelss) + ylim(0,1.0) + 
  stat_summary(geom = 'text', label = sig4, fun.y = max, vjust = -1, size = 3) +
  labs(title="Binary Jaccard", x="Treatment interaction") 

theme_set(theme_cowplot(font_size=10)) # reduce default font size
plot_grid(p1,p2,p3,p4, labels = "auto", nrow = 2, align = 'v')

# 2 #
ggplot(beta_data, aes(x=reorder(interaction, distwu,FUN = median), y=distwu, fill=interaction)) +
  geom_boxplot()+ theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(breaks=breakss, labels=labelss) + ylim(0,1.0) + 
  stat_summary(geom = 'text', label = sig1, fun.y = max, vjust = -1, size = 3) +
  labs(title="Weighted Unifrac",x="Treatment interaction", y = "Distance to centroid") 

# 3 #
# Single stressors versus the control plots
# Subset data to include only the single stressors
beta_data_single <- subset(beta_data, interaction=="C.26.I" | interaction=="C.30.I" | interaction=="C.26.S" |  interaction=="A.26.I" | interaction=="N.26.I")

ggplot(beta_data_single, aes(x=interaction, distwu, fill=interaction)) + 
  geom_boxplot()+ theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values=c("grey","tomato2","orange","gold3","olivedrab4")) +
  scale_x_discrete(breaks=breakss, labels=labelss) + labs(x="Temperature", y = "Distance to centroid") + 
  ylim(0,0.6) + geom_signif(comparisons = list(c("control","30")), annotation="*", tip_length = 0, map_signif_level=TRUE)

# 4 #
nut <- ggplot(beta_data, mapping = aes(x = nutrient, y = distwu, fill = nutrient)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.x = element_blank()) +
  scale_fill_manual(values=c("gold3","grey","olivedrab4")) + labs(x = "Nutrient Regime", y = "Chao1") 

scar <- ggplot(beta_data, mapping = aes(x = corallivory, y = distwu, fill = corallivory)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  scale_fill_manual(values=c("grey", "orange")) + labs(x = "Scarring")  +
  scale_x_discrete(labels = c("Control","Scarred"))

temp <- ggplot(beta_data, mapping = aes(x = temp, y = distwu, fill = temp)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.y = element_blank(), axis.title.x = element_blank()) + 
  scale_fill_manual(values=c("grey","tomato2")) + labs(x = "Temperature")  +
  scale_x_discrete(labels = c("Control","29"))

theme_set(theme_cowplot(font_size=10)) # reduce default font size
faith <- plot_grid(nut, scar, temp, nrow = 1, align = 'h', rel_widths = c(1.5,1,1))
##
