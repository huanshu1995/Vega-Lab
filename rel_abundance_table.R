########################################################
###      Change in Relative Abundance Table          ###
########################################################

rm(list=ls())

library(plotrix)


# Load qd data before rarefying and conglomerate the data into families
family_table <- tax_glom(biom, taxrank ="Family")


otus-family <- make_biom()

kwres <- read.csv(file = "/Users/Becca/Dropbox/Microbiome_analysis/mal/Analyses/kruskalwallis/kw_ocs_og.csv")
biom <- read.csv(file = "/Users/Becca/Dropbox/Microbiome_analysis/mal/otu_table/otu_table_mc2_w_tax_no_pynast_failures_newnames_clean_o100_s1000_noneg_filtagain.csv",
                        skip = 1) # skip option removes the first row "# Constructed from a biom file"


# Sum across all rows and add to the biom dataframe
biom_no_ids <- biom[,-1]
biom_sum <- as.vector(rowSums(biom1, dims = 1))
biom_ids_only <- biom[,1, drop = FALSE]
colnames(biom_ids_only) [1] <- "OTU"
biom_ids_only$sums <- biom_sum
colnames(biom_ids_only) [2] <- "sums"

# Merge to the Kruskal Wallis df
kw_res <- merge(kwres, biom_ids_only, by = "OTU")
colnames(kw_res)
kw_res <- kw_res[kw_res$"P"<0.05,]

# Calculate relative abundance
kw_res$rel_abundance_control <- kw_res$"Control_mean" / kw_res$"sums"
kw_res$rel_abundance_NH4 <- kw_res$"NH4._mean" / kw_res$"sums"
kw_res$rel_abundance_NO3 <- kw_res$"NO3._mean" / kw_res$"sums"
# Calculate % change in relative abundance
kw_res$change_NH4 <- kw_res$"rel_abundance_NH4" - kw_res$"rel_abundance_control"
kw_res$per_change_NH4 <- kw_res$"change_NH4" * 100
kw_res$change_NO3 <- kw_res$"rel_abundance_NO3" - kw_res$"rel_abundance_control"
kw_res$per_change_NO3 <- kw_res$"change_NO3" * 100
kw <- kw_res[,c("per_change_NO3", "per_change_NH4")]
row.names(kw) <- kw_res$OTU

mm <- as.matrix(kw, ncol= 2, byrow = FALSE)
cols <- color.scale(mm, extremes = c("orange", "purple"))

par(mar = c(0.5, 1, 2, 0.5))
# create empty plot
plot(1:100, axes = FALSE, xlab = "", ylab = "", type = "n")
# add table
addtable2plot(x = 1, y = 1, table = kw,
              bty = "o", display.rownames = TRUE,
              hlines = TRUE, vlines = TRUE,
              bg = cols,
              xjust = 2, yjust = 1, cex = 3)

par(mar = c(0.5, 8, 3.5, 0.5))
color2D.matplot(kw, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 1,
                vcol = "black",
                extremes = c("blue", "purple"))
axis(3, at = seq_len(ncol(kw)) - 0.5,
     labels = c("NO3-", "NH4+"), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(kw)) -0.5,
     labels = rev(rownames(kw)), tick = FALSE, las = 1, cex.axis = 1)
