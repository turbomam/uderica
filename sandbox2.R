# googd news... easy to reanalyze with different parameters
# accounting, decisions, desired conclusions/hypothesis/output... which alpha diversity measure... how do you know which to pick?
# time course?  beta diversity
# min count for filtering
# top N
# distnace and ordination method (how to quanititate usefulness?)
# analyze @ rank = ?
# aesthetics


# desired output... think about automation, scientific conclusion, paramters to vary
# plot_bar (abundance of taxa by treatment barplot)
# plot_richness (alpha diversity by treatment dot/whisker)
# alpha diversity table export
# ordination/PCoA trellis and featured plot
# differential abundance DEseq2

# in general, discovery/visual exploration or hypothesis testing?
# time course ?
# beta diversity (why not included?)



# variables to manipulate:
# filtering
# top N
# dist, ordination method
# rank to show



# what quantititavie output of ordination to measure?

# selection of alpha diversity measure?


# user experience?
# knit/markdown -> document (automated forma command line?) VS 
# rscript?  littler? -> files
# arguments from envvars? command line arguments?


# assumptions
# already converted from XXX to json
# already filtered for counts > 2 ???

otufile = "c:\\Users\\Mark Miller\\Desktop\\MAM_16S_POC.tar\\otu_table_json.biom"
trefile = "c:\\Users\\Mark Miller\\Desktop\\MAM_16S_POC.tar\\rep_set.tre"

otutable <-
  import_biom(
    BIOMfilename = otufile,
    treefilename = trefile,
    parseFunction = parse_taxonomy_greengenes
  )

# how may warnings to show of what tuype?


mapfile = "c:\\Users\\Mark Miller\\Desktop\\MAM_16S_POC.tar\\map16S_MAM_some_factors.txt"

mapping <-
  import_qiime_sample_data(mapfilename = mapfile)

phylo <- merge_phyloseq(otutable, mapping)

# get a summary (how does R know what function to apply when an object is named?)
# phylo
print(phylo)
# str not really helpful here
# str(qiimedata)
class(phylo)
# [1] "phyloseq"
# attr(,"package")
# not useful without some top-n filtering?
plot_bar(phylo, x = "treatment", fill = "Phylum")

# monitor level with the most change?  lowest NA? number of levels at a particular rank?
# how to tell best combination of phyogentic rank and top N?

the.N <- 20
topN <- names(sort(taxa_sums(phylo), TRUE)[1:the.N])

temp <- sort(taxa_sums(phylo), TRUE)
temp <- cbind.data.frame(names(temp), as.numeric(temp))
names(temp) <- c('taxon', 'counts')

entN <- prune_taxa(topN, phylo)

p = plot_bar(entN, "treatment", fill = "Family")
p + geom_bar(aes(color = Family, fill = Family), stat = "identity")
# removed: position = "stack" (warning) and facet on day (day 0 indistinguishable from tx control)

p <-
  plot_richness(
    phylo,
    x = "treatment",
    color = "treatment",
    measures = c("Observed", "Chao1", "Shanon", "Simpson", "InvSimpson")
  )
p + geom_point(size = 5, alpha = 0.7)

alpha.diversity <-
  estimate_richness(phylo,
                    measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))
# write.table(alpha.diversity, "safeearly.txt")

library("plyr")
packageVersion("plyr")
theme_set(theme_bw())
wh0 = genefilter_sample(phylo,
                        filterfun_sample(function(x)
                          x > 3),
                        A = 0.5 * nsamples(phylo))
GP1 = prune_taxa(wh0, phylo)
GP1 = transform_sample_counts(GP1, function(x)
  1E6 * x / sum(x))

expdata = phylo
ord_meths = c("PCoA")
dist = "bray"
plist = llply(as.list(ord_meths), function(i) {
  print(i)
  ordi = ordinate(physeq = expdata,
                  method = i,
                  distance = dist)
  p <-
    plot_ordination(expdata, ordi, "samples", color = "treatment")
  print(p)
})



ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
expdata = phylo
dist = "bray"
plist = llply(as.list(ord_meths), function(i, expdata, dist) {
  ordi = ordinate(expdata, method = i, distance = dist)
  plot_ordination(expdata, ordi, "samples", color = "treatment")
}, GP1, dist)
names(plist) <- ord_meths



pdataframe = ldply(plist, function(x) {
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})

names(pdataframe)[1] = "method"

p = ggplot(pdataframe,
           aes(Axis_1, Axis_2, color = treatment, fill = treatment))
p = p + geom_point(size = 4) + geom_polygon()
p = p + facet_wrap(~ method, scales = "free")
p = p + scale_fill_brewer(type = "qual", palette = "Set1")
p = p + scale_colour_brewer(type = "qual", palette = "Set1")
p

p = plist[[6]]
p = plist[[6]] + scale_colour_brewer(type = "qual", palette = "Set1")
p = p + scale_fill_brewer(type = "qual", palette = "Set1")
p = p + geom_point(size = 5) + geom_polygon(aes(fill = treatment))
p

# End of Erica's PCoA

# start of Differential Abundance, DESeq2

library("DESeq2")
kostic = phylo
kostic = subset_samples(kostic, mapping != "treatment")
kostic
head(sample_data(kostic)$treatment, 8)

kostic = subset_samples(kostic, mapping != "treatment")
diagdds = phyloseq_to_deseq2(kostic, ~ treatment)
diagdds = DESeq(diagdds, test = "Wald", fitType = "parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x)
  max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels = names(x))
# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x)
  max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels = names(x))
p = ggplot(sigtab, aes(x = Family, y = log2FoldChange, color = Phylum)) + geom_point(size = 6) +
  theme(axis.text.x = element_text(
    angle = -90,
    hjust = 0,
    vjust = 0.5
  ))
p

# display/export sigtab table