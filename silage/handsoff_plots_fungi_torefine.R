library("optparse")
library("phyloseq")
library("ggplot2")
library("plyr")
library("DESeq2")
library("ggthemes")

packageVersion("optparse")
packageVersion("phyloseq")
packageVersion("ggplot2")
packageVersion("plyr")
packageVersion("DESeq2")
packageVersion("ggthemes")

theme_set(theme_bw())

# filtering method 1
# if the counts for a OTU are <= (max otu count)/max.frac, then discard
# see also the calulation for the minumum number of samples in which this condition must be satisifed
max.frac <- 1000

# fitlering method 2
# only keep the top N OTUs, which could have differnt depths
#  and could all come from a singe upper rank
#  ie, cumulative counts at a higher rank aren't considered here
the.N <- 20
# at what rank should the results of this filtering be plotted
top.N.plot.rank <- "Genus"

# rank selected for subsequent plotting
selected.rank <- "Family"

# (discrete?) factor for differential abundance comparison
# add continuous or complex modelling later
exp.factor <- "treatment"

alpha.div.measures <-
  c("Observed", "Chao1", "Shanon", "Simpson", "InvSimpson")
alpha.div.plot.alab.fontsize <- 6
alpha.point.size <- 5
alpha.point.alpha <- 0.7

ordination.dist <- "bray"
ord_meths  <-
  c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")

# see below!
# ord_meths <- setdiff(ord_meths, "DPCoA")

# max significance for reporting and plotting in DESeq2
DESeq2.alpha <- 0.01

# what two ranks should eb sued to group the OTUs
higher.rank.name <- "Phylum"
lower.rank.name <- "Family"

# what values from DESeq analysis should be plotted
y.val.name <- "log2FoldChange"

experiment <- "LalStress"
# experiment <- "V-HMC"


# THIS WILL DEFINEITLY \BE DIFFERENT FOR EACH AXPERIMETN
#   what leves of the exp.factor should be compared
#  may wnat to do multiple
# numerator... top portion of ratio
constrat.num <- "C"

# denominator... bottom portion of ratio
if (experiment == "LalStress") {
  constrat.den <- "LB500.LATE"
  
  otufile <-
    paste0(
      '/home/mark/gitstage/uderica/silage/LalStress/Bacteria/raw_data_bacteria/flashed/extended/flash_trim_cat_pick/otu_table_no_sing_doub.biom'
    )
  
} else if (experiment == "V-HMC") {
  constrat.den <- "LB"
  
  otufile <-
    paste0(
      '/home/mark/gitstage/uderica/silage/V-HMC/Bacteria/raw_data_bacteria/flashed/extended/picked_temp/otu_table_mc2_w_tax_no_pynast_failures.biom'
    )
  
}




# min # counts per sample relvant for qiime corediv, but not here (at least yet)
# no "sing/doub" fitlering
# also, no json conversion!
# but have I chosen the most appropriate biom file?





trefile <-
  paste0(
    '/home/mark/gitstage/uderica/silage/',
    experiment ,
    '/Bacteria',
    '/raw_data_bacteria/flashed/extended/flash_trim_cat_pick/rep_set.tre'
  )
mapfile <-
  paste0('/home/mark/gitstage/uderica/silage/',
         experiment ,
         '/Bacteria/',
         'map.txt')

# create a pdf devcie anmed after the experimetn witha timestamp

proj.gene.prefix <- "interactive_careful_cat"







otutable <-
  import_biom(
    BIOMfilename = otufile,
    treefilename = trefile,
    parseFunction = parse_taxonomy_greengenes
  )

otutab_warns <- summary(warnings())

mapping <-
  import_qiime_sample_data(mapfilename = mapfile)


phylo <- merge_phyloseq(otutable, mapping)

print(phylo)




otu.dat <- otutable@otu_table@.Data

unlisted.otu.counts <- unlist(otu.dat)

# histogram, of unaggreagted counts, log scale
# would it be OK to fitler out taxa with observatiosn that are less than 1/x of teh max
# where x is soemthing like 1000?


max.count <- max(unlisted.otu.counts)
min.proposal <- max.count / max.frac

print(max.count)
print(min.proposal)

old.mar <- par("mar")
par("mar" = c(5, 4, 3, 2))

hist(
  log10(unlisted.otu.counts),
  breaks = 100,
  main = "Histogram of reads per sample/OTU",
  xlab = "log10(read counts)"
)

par("mar" = old.mar)

head(otu.dat)

otu.rowsums <- rowSums(otu.dat)
otu.rowsums <-
  cbind.data.frame(names(otu.rowsums), as.numeric(otu.rowsums))
names(otu.rowsums) <- c("OTU", "count.sum")


otu.colsums <- colSums(otu.dat)
otu.colsums <-
  cbind.data.frame(names(otu.colsums), as.numeric(otu.colsums))
names(otu.colsums) <- c("sample", "count.sum")
otu.colsums[] <- lapply(otu.colsums[], as.character)
otu.colsums[] <- lapply(otu.colsums[], as.numeric)

# is new.cleanup.reference OTU a problem?
# might still be mapped to a lineage
# what about unmapped?
# what has been lost before this?
#   un-flashed (really un-extended)
#   discarded because of quality by qiime split?
#   so some sequences not even make it into the biom fiel at all for soem reason?



map.dat <- as.data.frame(mapping@.Data)
names(map.dat) <- mapping@names

exp.design.matrix <- map.dat[, c("X.SampleID", "treatment", "day")]



if (sum(grepl(pattern = "h", x = exp.design.matrix$day))) {
  # written for V-HMC, but it may not be correct or useful
  # the hours portion doesn't refer to additonal anaeorib time
  # rather to air exposure after fermentation
  
  names(exp.design.matrix) <-
    c("SampleID", "treatment", "days.string")
  
  temp <-
    strsplit(as.character(exp.design.matrix$days.string), "\\+|h")
  temp <- ldply(temp, rbind)
  temp <- as.data.frame(temp)
  temp[] <- lapply(temp[], as.character)
  temp[] <- lapply(temp[], as.numeric)
  names(temp) <- c("days", "hours")
  temp$hours[is.na(temp$hours)] <- 0
  exp.design.matrix <- cbind.data.frame(exp.design.matrix, temp)
  exp.design.matrix$days.hours <- temp$days + (temp$hours / 24)
  
  exp.design.tabulation <-
    table(exp.design.matrix$treatment, exp.design.matrix$days.hours)
  
} else {
  exp.design.tabulation <-
    table(exp.design.matrix$treatment, exp.design.matrix$day)
}


exp.design.tabulation <- as.data.frame.matrix(exp.design.tabulation)
names(exp.design.tabulation) <-
  round(as.numeric(as.character(names(
    exp.design.tabulation
  ))), digits = 1)
exp.design.tabulation[exp.design.tabulation == 0] <- NA
temp.values <- unlist(exp.design.tabulation)
temp.flag <- complete.cases(temp.values)

table(temp.values[temp.flag])

par("mar" = c(5, 4, 3, 2))

hist(
  temp.values[temp.flag],
  breaks = length(unique(temp.values[temp.flag])),
  labels = as.character(min(temp.values[temp.flag]):max(temp.values[temp.flag])),
  main = "Histogram of repeats per COMBINATION of factor levels",
  xlab = NULL,
  axes = FALSE,
  ylim = c(0, (max(temp.values[temp.flag]) + 15))
)

par("mar" = old.mar)

min.replication <- min(temp.values[temp.flag])

tax.dat <- as.data.frame(otutable@tax_table)
temp <- !is.na(tax.dat)
taxonomic.specificity <- rowSums(temp)
tax.dat$taxonomic.specificity <- taxonomic.specificity

max.specificity <- max(tax.dat$taxonomic.specificity)
temp.frame <-
  cbind.data.frame((1:max.specificity), names(tax.dat)[1:max.specificity])
names(temp.frame) <- c("rank", "label")

max.count <- hist(tax.dat$taxonomic.specificity, plot = FALSE)
max.count <- max.count$counts
max.count <-  max(max.count)

par("mar" = c(5, 4, 3, 2))

hist(
  tax.dat$taxonomic.specificity,
  labels = as.character(temp.frame$label),
  breaks = 0:max.specificity,
  main = "Histogram of most specific rank",
  xlab = "Taxonomic rank",
  ylim = c(0, (max.count + 20))
)

par("mar" =  old.mar)



re.min.filtered <-
  genefilter_sample(phylo,
                    filterfun_sample(function(x)
                      x > min.proposal),
                    A = min.replication)


# re.min.filtered <-
#   genefilter_sample(phylo,
#                     filterfun_sample(function(x)
#                       x > 2),
#                     A = 2)


phylo.manip <- prune_taxa(re.min.filtered, phylo)

print(phylo.manip)

# absolute or relative?
# plot bar is absolute reads
# how to get rid of black dividers
# not seeing acetobacter any more

# ungrouped... one way to make sure that all samples made it to this phase
# plot_bar(phylo.manip, fill = "Family")

# could use + scale_colour_brewer(palette = "Set1"), but that adds in orange seperators !?
# set3 is paler but ahs 4 more levels
# paired also ahs more levels but is "paried"!
# + scale_colour_manual(values = colorRamps::primary.colors(n = 10))

p <- plot_bar(phylo.manip, exp.factor, fill = "Family")

# this gets rid of the black otu seperators
# add faceting?
p + geom_bar(aes(color = Family, fill = Family),
             stat = "identity")


# erica starts with top n

taxa.named.vector <- sort(taxa_sums(phylo), TRUE)
taxa.ids <- names(taxa.named.vector)
top.taxa <- taxa.ids[1:the.N]

phylo.pruned <- prune_taxa(top.taxa, phylo)

p <- plot_bar(phylo.pruned, exp.factor, fill = top.N.plot.rank)
p + geom_bar(aes_(
  color = as.name(top.N.plot.rank),
  fill = as.name(top.N.plot.rank)
), stat = "identity")

# plot richness without fitlering?
p <-
  plot_richness(phylo,
                x = exp.factor,
                color = exp.factor,
                measures = alpha.div.measures)

# parameterize
p + geom_point(size = alpha.point.size, alpha = alpha.point.alpha) + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 1,
  size = alpha.div.plot.alab.fontsize
))


alpha.diversity <-
  estimate_richness(phylo,
                    measures = alpha.div.measures)
alpha.diversity$sample.id.string <-
  as.numeric(sub("^X", "", rownames(alpha.diversity)))
rownames(alpha.diversity) <- NULL
non.sample.cols <-
  sort(setdiff(names(alpha.diversity), "sample.id.string"))
alpha.diversity <-
  alpha.diversity[order(alpha.diversity$sample), c("sample.id.string", non.sample.cols)]
# write.table(alpha.diversity, "safeearly.txt")
# dump to xtable if knitting


###   ###   ###

# messages when including DPCoA
# 1) Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)
# 2) non-unique values when setting 'row.names':

# DPCoA = Double Principle Coordinate Analysis using a (corrected, if necessary) phylogenetic/patristic distance between species. The calculation is performed by DPCoA(), which ultimately uses dpcoa after making the appropriate accessions/corrections of the data.

ord_meths <- setdiff(ord_meths, "DPCoA")

plist <- llply(as.list(ord_meths), function(i, phylo.manip, dist) {
  message(i)
  ordi <- ordinate(phylo.manip, method = i, distance = dist)
  plot_ordination(phylo.manip, ordi, "samples", color = exp.factor)
}, phylo.manip, ordination.dist)

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
# these analyses WERE being run at different levels of filtering
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

# DESeq2.data <- phylo.manip
# make a DESEQ2 object


# review the treatments (see legacy fitlering above)
# also shows the order, whcih can be relevled
# sample_data(phylo.manip)$treatment <- relevel(sample_data(phylo.manip)$treatment, "LB-AS")
# the first is the reference
sample_data(phylo.manip)$treatment


###   ###   ###

diagdds <- phyloseq_to_deseq2(phylo.manip, ~ treatment)

###   ###   ###

class(diagdds)

# create a matrix with direct access to the counts (whicha are a slot of an S4 object)
cts <- counts(diagdds)
class(cts)

# do rowwise and colwise sums... was handy when figurting out DESeq2 calulation failure
# would also be applciible to intiial filtering
DESeq2.rowSums <- rowSums(cts)
DESeq2.rowSums <-
  cbind.data.frame(names(DESeq2.rowSums), as.numeric(DESeq2.rowSums))

DESeq2.colSums <- colSums(cts)
DESeq2.colSums <-
  cbind.data.frame(names(DESeq2.colSums), as.numeric(DESeq2.colSums))

# this was necessary for the DESeq calulation when non-zero values were extremely sparse
geoMeans <-
  apply(cts, 1, function(row)
    if (all(row == 0))
      0
    else
      exp(mean(log(row[row != 0]))))
class(geoMeans)

diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
class(diagdds)

###   ###   ###

# should I make the model explicit here
diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")
class(diagdds)

rnms <- resultsNames(diagdds)
rnms



res <-
  results(
    diagdds,
    contrast = c("treatment", constrat.num, constrat.den),
    cooksCutoff = FALSE
  )
class(res)

# contrast
# this argument specifies what comparison to extract from the object to build a results table. one of either:
#   a character vector with exactly three elements:
#       the name of a factor in the design formula, the name of the numerator level for the fold change,
#       and the name of the denominator level for the fold change (simplest case)
# a list of 2 character vectors:
#       the names of the fold changes for the numerator, and the names of the fold changes for the denominator.
#       these names should be elements of resultsNames(object).
#       if the list is length 1, a second element is added which is the empty character vector, character().
#       (more general case, can be to combine interaction terms and main effects)
# a numeric contrast vector with one element for each element in resultsNames(object) (most general case)
# If specified, the name argument is ignored.
#
# name
#     the name of the individual effect (coefficient) for building a results table.
#     Use this argument rather than contrast for continuous variables,
#     individual effects or for individual interaction terms.
#     The value provided to name must be an element of resultsNames(object).

###   ###   ###


# get table with adjsuted pvalues below a user-speciifed cutoff
sigtab <- res[which(res$padj < DESeq2.alpha),]
class(sigtab)

sigtab <-
  cbind(as(sigtab, "data.frame"), as(tax_table(phylo.manip)[rownames(sigtab), ], "matrix"))

# this shows a single significatnce and abundance fold change by taxon
# what condition is being compared to what? (there are more than two treatmenrts in LalStress)
# see the phyloseq to deseq2 converter which has a parameter for the model

head(sigtab)

dim(sigtab)


# volcano?

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# use color an x position to render two different ranks
# get the rank labels in x
# make sure they're ordered factors



higher.rank.pos <- which(names(sigtab) == higher.rank.name)
lower.rank.pos <- which(names(sigtab) == lower.rank.name)


# higher order
x <-
  tapply(sigtab$log2FoldChange, sigtab[, higher.rank.pos], function(x)
    max(x))
x <- sort(x, TRUE)
sigtab[, higher.rank.pos] <-
  factor(as.character(sigtab[, higher.rank.pos]), levels = names(x))


# lower order
x <-
  tapply(sigtab$log2FoldChange, sigtab[, lower.rank.pos], function(x)
    max(x))
x <- sort(x, TRUE)
sigtab[, lower.rank.pos] <-
  factor(as.character(sigtab[, lower.rank.pos]), levels = names(x))

# create plot
p <-
  ggplot(sigtab,
         aes_(
           x = as.name(lower.rank.name),
           y = as.name(y.val.name),
           color = as.name(higher.rank.name)
         )) +
  geom_point(size = 6) +
  theme(axis.text.x = element_text(
    angle = -90,
    hjust = 0,
    vjust = 0.5
  )) + labs(title = paste0(
    "DIfferential abundance ratio, ",
    constrat.num ,
    " / ",
    constrat.den
  ))
p

# dev.off()

# dump sigtab etoa fild
