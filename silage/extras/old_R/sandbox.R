# assumes already have clustered/picked OTUs
# (eg from QIIME)
# in XYZ format

# corn plant matter/grain is innoculated with microorganisms
# potentially treated with chemicals
# and fermented in a lab silo
# some proceesing takes place, which should result in a faithful extraction of microbial DNA
# which is then amplified to get V4 sequences

# how to account for sequences that a removed durring filtering, flashing, aligning, mapping, otu picking...

# phyloseq uses (vegan, ade4, ape, picante), (ggplot2)

# running uclust in QIIME
# what about rdp... krona plots (or blast?)
# mark's hypothesis:  shows all possible amtching taxa... not best single hypothesis???

# what is biom?

# distance
# what is ordination (related to dimensionaliyt reduction?)

library("phyloseq")
packageVersion("phyloseq")
## [1] '1.20.0'
library("ggplot2")
packageVersion("ggplot2")

# Caporaso, J. G., et al. (2011). Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample.
#   PNAS, 108, 4516-4522. PMCID: PMC3063599

theme_set(theme_bw())


data(GlobalPatterns)
# data(esophagus)
# data(enterotype)
# data(soilrep)

tree = phy_tree(GlobalPatterns)
tax  = tax_table(GlobalPatterns)
otu  = otu_table(GlobalPatterns)
sam  = sample_data(GlobalPatterns)
otutax = phyloseq(otu, tax)
otutax

GP2 = merge_phyloseq(otutax, sam, tree)
identical(GP2, GlobalPatterns)
