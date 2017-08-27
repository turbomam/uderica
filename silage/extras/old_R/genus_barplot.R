library("phyloseq")
library("ggplot2")
library("ape")
library("RColorBrewer")
library("tidyverse")
library("highcharter")

library("plyr")
packageVersion("plyr")

# how to get filenames and other parameters from calling script?
# command line args?
# environment variables?
# config file?

# run as script, dumping output to files?
# or run as markdown knitr, constructing monolithic document

# script provided in word file "Microbiome analysis.docx" filters > 0
# but hasn't filtering already been performed at > 2
# in filter_otus_from_otu_table.py ?
min.prune.count <- 3
# in word file "Microbiome analysis.docx", was saying “top ten”, was setting filter at 20
top.n.choice <-  30



otutable <-
  import_biom(
    BIOMfilename =
      '/home/mark/silage/16s/pairedReads/flashed/extended/flash_trim_cat_pick/otu_table_json.biom',
    treefilename =
      '/home/mark/silage/16s/pairedReads/flashed/extended/flash_trim_cat_pick/rep_set.tre',
    parseFunction = parse_taxonomy_greengenes
  )

warnings()

mapping <-
  import_qiime_sample_data(mapfilename =
                             "/home/mark/silage/16s/map16S_MAM_some_factors.txt")
phylo <- merge_phyloseq(otutable, mapping)
phylo <- prune_taxa(taxa_sums(phylo) > min.prune.count, phylo)
phylo
top.taxa.ids <- names(sort(taxa_sums(phylo), TRUE)[1:top.n.choice])
top.counts <- prune_taxa(top.taxa.ids, phylo)
# p = plot_bar(ent10, "treatment", fill = "Family", facet_grid = ~ day)
# p + geom_bar(aes(color = Family, fill = Family),
#              stat = "identity",
#              positon = "stack")

# removing day from factors...
# day 0 and tx control are unique, singular and identicial
# why color and fill?
# Warning: Ignoring unknown parameters: positon
# spelling error... also using old ggplot syntax?
# not even necessary in this example?
# NA taxa from open reference alignemnt?...  good, doesn't hide unmapped reads
# altough data may have been lost or masked at other steps

# what rank to report? was usign family.  so far, genus seems good to me.  lots of NAs in species.
# who and how do thsoe krona plots include species?  blast/rdp?
# show absolute counts or relative abundance?
# log scale ?
p <- plot_bar(top.counts, "treatment", fill = "Genus")
p + geom_bar(aes(color = Genus), stat = "identity")

# why b/w, not color?
p <-
  plot_richness(
    phylo,
    x = "treatment",
    measures = c("Observed", "Chao1", "Shanon", "Simpson", "InvSimpson")
  )
p + geom_point(size = 5, alpha = 0.7)

alpha.diversity <-
  estimate_richness(phylo,
                    measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))

# crude assumption:
#   one character (X) has been prepended onto each sampkle numeber
alpha.diversity$sample.id <- substring(rownames(alpha.diversity), 2)
rownames(alpha.diversity) <- NULL
alpha.diversity <-
  alpha.diversity[order(alpha.diversity$sample.id),]
# use a prettier printing method/function
alpha.diversity

###

theme_set(theme_bw())
wh0 <- genefilter_sample(phylo, filterfun_sample(function(x)
  x > 5), A = 0.5 * nsamples(phylo))
GP1 <- prune_taxa(wh0, phylo)
GP1 <- transform_sample_counts(GP1, function(x)
  1E6 * x / sum(x))

dist <- "bray"
ord_meths <- c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist <- llply(as.list(ord_meths), function(i, phylo, dist) {
  ordi = ordinate(phylo, method = i, distance = dist)
  plot_ordination(phylo, ordi, "sampletype", color = "sample")
}, GP1, dist)
names(plist) <- ord_meths
pdataframe <- ldply(plist, function(x) {
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] <- "method"
p <- ggplot(pdataframe,
            aes(Axis_1, Axis_2, color = treatment, fill = treatment))
p <- p + geom_point(size = 4) + geom_polygon()
p <- p + facet_wrap( ~ method, scales = "free")
p <- p + scale_fill_brewer(type = "qual", palette = "Set1")
p <- p + scale_colour_brewer(type = "qual", palette = "Set1")
p
p <- plist[[6]]
