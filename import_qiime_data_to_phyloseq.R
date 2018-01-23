########################################################
### Importing microbiome data from Qiime to Phyloseq ###
########################################################

rm(list=ls())


library(phyloseq)
library(ape)


##
# Importing data and transformations ##
# Import qiime mapping file, biom otu table, and tree
mapfile = "/Users/Becca/Dropbox/Microbiome_analysis/mal/map_int_noneg.txt"
map = import_qiime_sample_data(mapfile)
class(map)

# Import tree
tree = read_tree("/Users/Becca/Dropbox/Microbiome_analysis/mal/rep_set.tre")

# This biom file is a result of the pick_open_reference_otu command in qiime1,
# In qiime1, I already removed mitochondria and chloroplasts, removed otus with < 100 occurrences
# and samples with < 1000 counts.
biomfile = "/Users/Becca/Dropbox/Microbiome_analysis/mal/otu_table_mc2_w_tax_no_pynast_failures_newnames_clean_o100_s1000_noneg_filtagain.biom"
biom = import_biom(biomfile, parseFunction = parse_taxonomy_default)

qd = merge_phyloseq(map,tree,biom) 
##


##
# Various necessary adjustments to the data (Specific to your data)
# Change rank names from Rank1, Rank2, ... to Kingdom, Phylum, etc.
colnames(tax_table(qd)) = c(k="Kingdom", p="Phylum", c="Class", o="Order",f="Family", g="Genus", s="Species")
# check if it worked 
rank_names(qd)
# Make a numeric factor into a categorical
# In this case, i had temperature as 26 or 29, R considers this numerical but I want to consider it categorical
sample_data(qd)$temp=factor(get_variable(qd,"temp"))
##

## The phyloseq object qd is now ready to use in diversity analyses
