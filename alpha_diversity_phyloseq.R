####################################################
### Alpha diversity measures, graphs, and stats ####
####################################################

# HELLO WORLD

setwd("/Users/Becca/Dropbox/Microbiome_analysis/mal/")

library('phyloseq'); packageVersion('phyloseq')
library("ggplot2")
library('vegan'); packageVersion('vegan')
library("ggsignif")
library("plyr")
library("car")
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(picante)
library(cowplot)

# First run import_qiime_data_to_phyloseq.R to get the qd object in the R environment

# Modified from: (http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#alpha_diversity)


##
# Minimum number of reads in a sample: establishing the depth to rarefy each sample to for a standard sampling effort
# Subsample to the minimum number of reads (min_lib). We will repeat this 100 times and average the diversity 
# estimates from each trial
min_lib <- min(sample_sums(qd))
##


##
# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(qd)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(qd)

evenness <- matrix(nrow =nsamp, ncol = trials)
row.names(evenness) <- sample_names(qd)

faithPD <- matrix(nrow = nsamp, ncol = trials)
row.names(faithPD) <- sample_names(qd)
##


##
# It is always important to set a seed when you subsample so your result is replicable
set.seed(3)

# Rarefy the data to an even sampling depth, and calculate different measures of alpha diversity
# Options for measures = ("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(qd,sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(subset(estimate_richness(r, measures = "Chao1"), select = c(1))))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "Simpson")))
  evenness[ ,i] <- even
  
  # Calculate Faith's PD
  faith <- as.numeric(as.matrix(subset(estimate_pd(r), select = c(1))))  # estimate_pd is a function assigned at the end of the script
  faithPD[ ,i] <- faith
}
# Included the subset in "rich" because the Chao1 measurement outputs two measures per sample (Chao1 and se.chao1)
# and we only want Chao1, so we select for the first column
##

# test


##
# Calculate the mean and standard deviation per sample for observed richness and inverse simpson's index
# and store those values in a datafram.

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Evenness", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of Faith's PD estimates
SampleID <- row.names(faithPD)
mean <- apply(faithPD, 1, mean)
sd <- apply(faithPD, 1, sd)
measure <- rep("FaithPD", nsamp)
faith_stats <- data.frame(SampleID, mean, sd, measure)

# Combine our estimates for richness and evenness into one dataframe
alpha <- rbind(rich_stats, even_stats, faith_stats)

# Add the sample metadata into this dataframe using the merge() command
# I had to change my column name from X.SampleID to SampleID to match alpha
colnames(sample_data(qd))[1] <- "SampleID" 
s <- data.frame(sample_data(qd))
alphadiv <- merge(alpha, s, by = "SampleID")
##


##
# Plotting
# Subset richess and evenness
alpha_rich <- subset(alphadiv, measure == "Richness")
alpha_even <- subset(alphadiv, measure == "Evenness")
alpha_faith <- subset(alphadiv, measure == "FaithPD")

breakss <- c("C.26.I","A.26.I","N.26.S","A.26.S","A.30.S","N.30.S","C.30.S","A.30.I",
             "N.26.I","N.30.I","C.26.S","C.30.I")
labelss <- c("control","NH3+","NO3-, scarred", "NH3+, scarred","NH3+, scarred, 29","NO3-, scarred, 29","scarred, 29","NH3+, 29",
             "NO3-","NO3-, 29","scarred","29")

# Three stacked plots of chao1, Simpson, and Faith's PD by treatment and increasing means
chao <- ggplot(alpha_rich, aes(x = reorder(interaction, mean, FUN = median), y = mean, fill = interaction)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x = element_blank()) + scale_x_discrete(breaks=breakss, labels=labelss) + labs(y = "Chao1") 

simp <- ggplot(alpha_even, aes(x = reorder(interaction, mean, FUN = median), y = mean, fill = interaction)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x = element_blank()) + scale_x_discrete(breaks=breakss, labels=labelss) + labs(y = "Simpson")  

faith <- ggplot(alpha_faith, aes(x = reorder(interaction, mean, FUN = median), y = mean, fill = interaction)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), 
  axis.title.x = element_blank()) + scale_x_discrete(breaks=breakss, labels=labelss) + labs(y = "Faith's PD") 

theme_set(theme_cowplot(font_size=10)) # reduce default font size
plot_grid(chao,simp,faith, labels = "auto", ncol = 1, align = 'v')

# Three stacked plots of chao1, simpson, and Faith's PD by nutrient, scarring, and temperature
# Faith's PD
faith1 <- ggplot(alpha_faith, mapping = aes(x = nutrient, y = mean, fill = nutrient)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("gold3","grey","olivedrab4")) + labs(x = "Nutrient Regime", y = "Faith's PD") 

faith2 <- ggplot(alpha_faith, mapping = aes(x = corallivory, y = mean, fill = corallivory)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.y = element_blank()) + scale_fill_manual(values=c("grey","orange")) + labs(x = "Scarring") +
  scale_x_discrete(labels = c("Control","Scarred"))

faith3 <- ggplot(alpha_faith, mapping = aes(x = temp, y = mean, fill = temp)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.y = element_blank()) + scale_fill_manual(values=c("grey","tomato2")) + labs(x = "Temperature") +
  scale_x_discrete(labels = c("Control","29"))

# Group three graphs
theme_set(theme_cowplot(font_size=10)) # reduce default font size
faith <- plot_grid(faith1,faith2,faith3, nrow = 1, align = 'h', rel_widths = c(1.5,1,1))

# Evenness
even1 <- ggplot(alpha_even, mapping = aes(x = nutrient, y = mean, fill = nutrient)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.x = element_blank()) +
  scale_fill_manual(values=c("gold3","grey","olivedrab4")) + labs(x = "Nutrient Regime", y = "Simpson") 

even2 <- ggplot(alpha_even, mapping = aes(x = corallivory, y = mean, fill = corallivory)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  scale_fill_manual(values=c("grey", "orange")) + labs(x = "Scarring")  +
  scale_x_discrete(labels = c("Control","Scarred"))

even3 <- ggplot(alpha_even, mapping = aes(x = temp, y = mean, fill = temp)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  scale_fill_manual(values=c("grey","tomato2")) + labs(x = "Temperature")  +
  scale_x_discrete(labels = c("Control","29"))

# Group three evenness graphs
even <- plot_grid(even1,even2,even3, nrow = 1, align = 'h', rel_widths = c(1.5,1,1))

# Richness
rich1 <- ggplot(alpha_rich, mapping = aes(x = nutrient, y = mean, fill = nutrient)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.x = element_blank()) +
  scale_fill_manual(values=c("gold3","grey","olivedrab4")) + labs(x = "Nutrient Regime", y = "Chao1") 

rich2 <- ggplot(alpha_rich, mapping = aes(x = corallivory, y = mean, fill = corallivory)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  scale_fill_manual(values=c("grey", "orange")) + labs(x = "Scarring")  +
  scale_x_discrete(labels = c("Control","Scarred"))

rich3 <- ggplot(alpha_rich, mapping = aes(x = temp, y = mean, fill = temp)) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.title.y = element_blank(), axis.title.x = element_blank()) + 
  scale_fill_manual(values=c("grey","tomato2")) + labs(x = "Temperature")  +
  scale_x_discrete(labels = c("Control","29"))

# Group three richness graphs
rich <- plot_grid(rich1,rich2,rich3, nrow = 1, align = 'h', rel_widths = c(1.5,1,1))

# Group Faith's PD, evenness, and richness graphs
plot_grid(rich, even, faith, ncol = 1, align = 'w', labels = "auto")
##


##
# Statistical tests for significance 

# Pairwise t-tests
pairwise.t.test(alpha_even$mean, alpha_even$interaction, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_even$mean, alpha_even$nutrient, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_even$mean, alpha_even$temp, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_even$mean, alpha_even$corallivory, p.adjust.method = "bonferroni")

pairwise.t.test(alpha_rich$mean, alpha_rich$interaction, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_rich$mean, alpha_rich$nutrient, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_rich$mean, alpha_rich$temp, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_rich$mean, alpha_rich$corallivory, p.adjust.method = "bonferroni")

pairwise.t.test(alpha_rich$mean, alpha_faith$interaction, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_rich$mean, alpha_faith$nutrient, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_rich$mean, alpha_faith$temp, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_rich$mean, alpha_faith$corallivory, p.adjust.method = "bonferroni")

# Linear regression and anova
reg1 <- lm(mean ~ nutrient*temp*corallivory + tank + colony, data = alpha_even)
reg2 <- lm(mean ~ nutrient*temp*corallivory + tank + colony, data = alpha_rich)
reg3 <- lm(mean ~ nutrient*temp*corallivory + tank + colony, data = alpha_faith)

summary(reg1)
anova(reg1)

summary(reg2)
anova(reg2)

summary(reg3)
anova(reg3)
##


##
# Faith's Phylogenetic Distance function

# Use this function to calculate Faith's Phylogenetic Distance from a Phyloseq object
# After running this function, I input the function into the code above.

# @ Thomas W. Battaglia

#' Estimate phylogenetic diversity from phyloseq object
#'
#' Estimate the Faiths phylogenetic diverstiy from an OTU table and phylogenetic tree
#'
#' @param phylo A phyloseq object with an OTU table and phylogenetic tree slot.
#' @return A data.frame of phylogenetic diversity metrics.
#' @export
# Estimate PD-whole tree
estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}
# estimate_pd(r)
# Returns a data frame of PD (phylogenetic distance) and SR (species richness)
##
