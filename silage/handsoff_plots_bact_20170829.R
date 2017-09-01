#!/usr/bin/env Rscript

library("optparse")
library("phyloseq")
library("ggplot2")
library("plyr")
library("DESeq2")

packageVersion("optparse")
packageVersion("phyloseq")
packageVersion("ggplot2")
packageVersion("plyr")
packageVersion("DESeq2")

the.N <- 20
selected.rank <- "Family"
exp.factor <- "treatment"

re.min.filt.val <- 3

ordination.dist <- "bray"
ord_meths  <-
  c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")

# useful to report all of these diversity emasures?  how do you knwo which is most meaningful?
# cramped plot!
alpha.div.measures <-
  c("Observed", "Chao1", "Shanon", "Simpson", "InvSimpson")

today <- Sys.time()
time.stamp <- format(today, format = "%Y%m%d%H%M")

# aesthetics
theme_set(theme_bw())
alpha.div.plot.alab.fontsize <- 7

# not currenttly asserting an output device (or asing the user for an output filename or path)
# so will go to  in the current working directory

# account for discarded and unmapped seqs
# how may warnings to show, of what type?

# good news... easy to reanalyze with different parameters
#  parameter decisions, desired conclusions/hypothesis/output...
# include time course?  beta diversity (why not included?)
# min count for filtering
# distnace and ordination method (how to quanititate usefulness?)


# desired output... think about automation, scientific conclusion, paramters to vary
# plot_bar (abundance of taxa by treatment barplot)
# plot_richness (alpha diversity by treatment dot/whisker)
# alpha diversity table export
# ordination/PCoA trellis and featured plot
# differential abundance DEseq2

# in general, discovery/visual exploration or hypothesis testing?
# beta diversity

# variables to manipulate:
# filtering
# dist, ordination method

# what quantititavie output of ordination to measure?

# user experience?
# knit/markdown -> document (automated forma command line?) VS
# rscript?  littler? -> files
# arguments from envvars? command line arguments?


# assumptions
# otufile already converted from HDF5(?) to json
# already filtered for counts > 2 ???

if (interactive()) {
  study.kingdom <- "V-HMC/Bacteria"
  otufile <-
    paste0(
      '/home/mark/gitstage/uderica/silage/',
      study.kingdom,
      '/raw_data_bacteria/flashed/extended/flash_trim_cat_pick/otu_table_json.biom'
    )
  trefile <-
    paste0(
      '/home/mark/gitstage/uderica/silage/',
      study.kingdom,
      '/raw_data_bacteria/flashed/extended/flash_trim_cat_pick/rep_set.tre'
    )
  mapfile <-
    paste0('/home/mark/gitstage/uderica/silage/',
           study.kingdom,
           '/map.txt')
  proj.gene.prefix <- "interactive"
  
} else {
  option_list <- list(
    make_option(
      c("-o", "--otufile"),
      type = "character",
      default = NULL,
      help = "full path to otu file e.g. otu_table_json.biom",
      metavar = "character"
    ),
    
    make_option(
      c("-t", "--trefile"),
      type = "character",
      default = NULL,
      help = "full path to tree file e.g. rep_set.tre",
      metavar = "character"
    ),
    
    make_option(
      c("-m", "--mapfile"),
      type = "character",
      default = NULL,
      help = "full path to QIIME mapping file",
      metavar = "character"
    ),
    
    make_option(
      c("-p", "--prefix"),
      type = "character",
      default = NULL,
      help = "project- and gene-specific prefix for output",
      metavar = "character"
    )
  )
  
  opt_parser <- OptionParser(option_list = option_list)
  
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$otufile) ||
      is.null(opt$trefile) ||
      is.null(opt$mapfile) ||
      is.null(opt$prefix)) {
    print_help(opt_parser)
    stop(
      "example: Rscript handsoff_plots.R --otufile /path/to/otu_table_json.biom --trefile /path/to/rep_set.tre --mapfile /path/to/<QIIME mapping file>.txt --prefix octospring_16S\n",
      call. =
        FALSE
    )
  }
  
  otufile <- opt$otufile
  trefile <- opt$trefile
  mapfile <- opt$mapfile
  proj.gene.prefix  <- opt$prefix
  
}

print(otufile)
print(trefile)
print(mapfile)

otutable <-
  import_biom(
    BIOMfilename = otufile,
    treefilename = trefile,
    parseFunction = parse_taxonomy_greengenes
  )

mapping <-
  import_qiime_sample_data(mapfilename = mapfile)

# There were 50 or more warnings (use warnings() to see the first 50)
# how many warnings are there really? are they all the same?

x <- summary(warnings())


# In parseFunction(i$metadata$taxonomy) : No greengenes prefixes were found.
# Consider using parse_taxonomy_default() instead if true for all OTUs.
# Dummy ranks may be included among taxonomic ranks now.

phylo <- merge_phyloseq(otutable, mapping)

# get a summary (how does R know what default print function to apply when an object is named?)
# phylo
print(phylo)
# str not really helpful here
# str(qiimedata)
# class(phylo)
# [1] "phyloseq"
# attr(,"package")
# bare bones confrimation that phyloseq was able to process the suplied QIIME data
# not useful without some top-n filtering?
# plot_bar(phylo, x = exp.factor, fill = "Phylum")

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3688 taxa and 28 samples ]
# sample_data() Sample Data:       [ 28 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 3688 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3688 tips and 3686 internal nodes ]

# # scientifically uninteresting reality check:
# # is data avaialble for all samples?
# # does not necessarily mean that useful tax ranks have been assigned or that
##   sample have useful metadata factors
#
# View(as.data.frame(phylo@sam_data@.Data))
#
# View(as.data.frame(phylo@tax_table@.Data))
# # not interested in unidenfiitfable taxa and or with low counts?
# View(as.data.frame(phylo@otu_table@.Data))

# look thorugh object phylo
# only keep OTUs with more than 2 counts in at least 30% of the samples?
re.min.filtered <-
  genefilter_sample(phylo, filterfun_sample(function(x)
    x > 2), A = 0.3 * nsamples(phylo))
phylo.manip <- prune_taxa(re.min.filtered, phylo)
plot_bar(phylo.manip, fill = "Genus")
plot_bar(phylo.manip, exp.factor, fill = "Genus")


###   ###   ###

# monitor the rank with the most change?  lowest NA? number of levels at a particular rank?
# how to tell best combination of phyogentic rank and top N?

# 2017-08-29 I just dont get how this works
# couldn't these OTUs/taxa be of any rank/mixed ranks?
# not necessarily gettign teh top N at a particular rank?
sorted.obs.per.taxon <- sort(taxa_sums(phylo), TRUE)
topN <- names(sorted.obs.per.taxon)[1:the.N]

# is temp being used in amy way?  did I make it up?
temp <- sort(taxa_sums(phylo), TRUE)
temp <- cbind.data.frame(names(temp), as.numeric(temp))
names(temp) <- c('taxon', 'counts')

entN <- prune_taxa(topN, phylo)

message(" \n basic bar plot \n ")

p <- plot_bar(entN, exp.factor, fill = selected.rank)
print(p)

# don't even need any of the fololowing?
# p + geom_bar(aes(color = Family, fill = Family), stat = "identity")
# # removed: position = "stack" (warning) and facet on day (day 0 indistinguishable from tx control)

message(" \n richness plot \n ")

p <-
  plot_richness(phylo,
                x = exp.factor,
                color = exp.factor,
                measures = alpha.div.measures)

# parameterize
p + geom_point(size = 5, alpha = 0.7) + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 1,
  size = alpha.div.plot.alab.fontsize
))

# theme(axis.text.x = element_text(angle = 90, hjust = 1))

# print(p)

message(" \n create alpha diversity table \n ")

alpha.diversity <-
  estimate_richness(phylo, measures = alpha.div.measures)

# would be nice to put somethign meaningful in the fiilenames of dumped tables
write.table(alpha.diversity,
            paste0(proj.gene.prefix, "_alpha_div_",  time.stamp, ".txt"))


#
# PUT IN A PRETTY TABLE HEAD WITH SOEMTHING LIKE XTABLES (IF KNITTING)
#

message(" \n start ordination analysis \n ")


re.min.filtered <- genefilter_sample(phylo,
                                     filterfun_sample(function(x)
                                       x > re.min.filt.val),
                                     A = 0.5 * nsamples(phylo))
phylo.manip <- prune_taxa(re.min.filtered, phylo)
phylo.manip <- transform_sample_counts(phylo.manip, function(x)
  1E6 * x / sum(x))


## what does this add over the multi-method block below?
# expdata <- phylo.manip
# mini_ord_meths <- c("PCoA")
# plist <- llply(as.list(ord_meths), function(i) {
#   print(i)
#   ordi <- ordinate(physeq = expdata,
#                    method = i,
#                    distance = ordination.dist)
#   p <-
#     plot_ordination(expdata, ordi, "samples", color = exp.factor)
#   print(p)
# })

ord_meths <-
  c("CCA", "DCA", "DPCoA", "MDS", "NMDS", "PCoA", "RDA")

expdata <- phylo.manip
dist <- ordination.dist

ord_meths <- c("CCA", "DCA", "DPCoA", "MDS", "NMDS", "PCoA", "RDA")

# messages when including DPCoA
# 1) Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)
# 2) non-unique values when setting 'row.names':

# DPCoA = Double Principle Coordinate Analysis using a (corrected, if necessary) phylogenetic/patristic distance between species. The calculation is performed by DPCoA(), which ultimately uses dpcoa after making the appropriate accessions/corrections of the data.

ord_meths <- setdiff(ord_meths, "DPCoA")

plist <- llply(as.list(ord_meths), function(i, expdata, dist) {
  message(i)
  ordi <- ordinate(expdata, method = i, distance = dist)
  plot_ordination(expdata, ordi, "samples", color = exp.factor)
}, expdata, dist)

names(plist) <- ord_meths


pdataframe <- ldply(plist, function(x) {
  df <- x$data[, 1:2]
  colnames(df) <- c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})

names(pdataframe)[1] <- "method"


p <- ggplot(pdataframe,
            aes(Axis_1, Axis_2, color = treatment, fill = treatment))
p <- p + geom_point(size = 4) + geom_polygon()
p <- p + facet_wrap(~ method, scales = "free")
p <- p + scale_fill_brewer(type = "qual", palette = "Set1")
p <- p + scale_colour_brewer(type = "qual", palette = "Set1")
p

# user may desire an individual plot on a page by itself
# so print all on a page by themselves?
# how did Erica detmine that p = plist[[6]] (MDS) was the best ordination?


#
# these analyses are being run at different levels of filtering
# initila barplot and alpha deiversity at 2 (from Bash script)
# ordination refiltered to 3
# diff abundance back to 2
#

message(" \n start differential analysis \n ")


# DESeq2.data <- phylo
# DESeq2.data <- subset_samples(DESeq2.data, mapping != exp.factor)
# DESeq2.data
# head(sample_data(DESeq2.data)$treatment, 8)
# # http://joey711.github.io/phyloseq-extensions/DESeq2.html warns of "none" diagnoses
# # I don't think that kind of filetering is relevant here
# sample_data(DESeq2.data)$treatment

# DESeq2.data <- subset_samples(DESeq2.data, mapping != exp.factor)

# diagdds <- phyloseq_to_deseq2(DESeq2.data, ~ treatment)
#   V-HMC, without modificatiosn below:  
#      estimating size factors Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#      every gene contains at least one zero, cannot compute log geometric means
# DESeq2.rowSums <- rowSums(counts(diagdds))
# DESeq2.rowSums <- cbind.data.frame(names(DESeq2.rowSums),as.numeric(DESeq2.rowSums))
#
# DESeq2.colSums <- colSums(counts(diagdds))
# DESeq2.colSums <- cbind.data.frame(names(DESeq2.colSums),as.numeric(DESeq2.colSums))

# dezeroed  = diagdds[ rowSums(counts(diagdds)) > 1000 , ]
# dezeroed

# x <- estimateSizeFactors(diagdds, type="iterate")
# idx <- rowSums( counts(x, normalized=TRUE) >= 5 ) >= 3
# x <- x[idx,]
# x <- DESeq(x)

DESeq2.data <- phylo.manip
diagdds <- phyloseq_to_deseq2(DESeq2.data, ~ treatment)
sample_data(DESeq2.data)$treatment

cts <- counts(diagdds)

DESeq2.rowSums <- rowSums(cts)
DESeq2.rowSums <-
  cbind.data.frame(names(DESeq2.rowSums), as.numeric(DESeq2.rowSums))

DESeq2.colSums <- colSums(cts)
DESeq2.colSums <-
  cbind.data.frame(names(DESeq2.colSums), as.numeric(DESeq2.colSums))

geoMeans <-
  apply(cts, 1, function(row)
    if (all(row == 0))
      0
    else
      exp(mean(log(row[row != 0]))))
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)

diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")

res <- results(diagdds, cooksCutoff = FALSE)
alpha <- 0.01
sigtab <- res[which(res$padj < alpha),]
sigtab <-
  cbind(as(sigtab, "data.frame"), as(tax_table(DESeq2.data)[rownames(sigtab),], "matrix"))
# this shows a single significatnce and abundance fold change by taxon
# what condition is being compared to what? (there are more than two treatmenrts in LalStress)
head(sigtab)

dim(sigtab)

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x)
  max(x))
x <- sort(x, TRUE)
sigtab$Phylum <-
  factor(as.character(sigtab$Phylum), levels = names(x))
# Family order
x <- tapply(sigtab$log2FoldChange, sigtab$Family, function(x)
  max(x))
x <- sort(x, TRUE)
sigtab$Family <-
  factor(as.character(sigtab$Family), levels = names(x))
p <-
  ggplot(sigtab, aes(x = Family, y = log2FoldChange, color = Phylum)) + geom_point(size = 6) +
  theme(axis.text.x = element_text(
    angle = -90,
    hjust = 0,
    vjust = 0.5
  ))
p

# dev.off()
