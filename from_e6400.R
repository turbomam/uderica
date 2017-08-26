library("phyloseq")
packageVersion("phyloseq")
## [1] '1.16.2'
data(GlobalPatterns)
library("ggplot2")
packageVersion("ggplot2")
## [1] '2.1.0'
library("plyr")
packageVersion("plyr")
## [1] '1.8.4'

theme_set(theme_bw())

GP = GlobalPatterns
wh0 = genefilter_sample(GP, filterfun_sample(function(x)
  x > 5), A = 0.5 * nsamples(GP))
GP1 = prune_taxa(wh0, GP)

GP1 = transform_sample_counts(GP1, function(x)
  1E6 * x / sum(x))

phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm =
                      TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

human = get_variable(GP1, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
sample_data(GP1)$human <- factor(human)


GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1,
                     GP.ord,
                     type = "taxa",
                     color = "Phylum",
                     title = "taxa")
print(p1)

R.version
# _
# platform       x86_64-pc-linux-gnu
# arch           x86_64
# os             linux-gnu
# system         x86_64, linux-gnu
# status
# major          3
# minor          4.1
# year           2017
# month          06
# day            30
# svn rev        72865
# language       R
# version.string R version 3.4.1 (2017-06-30)
# nickname       Single Candle

temp <- installed.packages()
temp[rownames(temp) %in% c("phyloseq", "vegan") , c("Package", "Version")]
