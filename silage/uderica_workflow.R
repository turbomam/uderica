# package versions
# warnings about phyloseq biom inport
# omitted snippets:
# could show snippets of otu counts, sample data, taxon data...

# my.exp <- "Aer_Comp"
# my.fact <- "Treatment"
# my.ref <- "C/0"
# my.other <- "S2-3H-LATE/56"

my.exp <- "LalStress"
my.fact <- "treatment"
my.ref <- "C"
my.other <- "LB500-EARLY"

# my.exp <- "V-HMC"
# my.fact <- "Treatment"
# my.ref <- "C/0"
# # fungi
# # my.other <- "LB-AS/90+50h"
# # bacteria
# my.other <- "LB/30"


# my.domain <- "Fungi"
my.domain <- "Bacteria"

my.timestamp <- format(Sys.time(), format = "%Y%m%d%H%M")

sinkfile <- paste0(my.exp,
                   "_",
                   my.domain,
                   "_",
                   my.timestamp,
                   ".txt")

con <- file(sinkfile)
sink(con, append = TRUE)
sink(con, append = TRUE, type = "message")

pdf(paste0(my.exp,
           "_",
           my.domain,
           "_",
           my.timestamp,
           ".pdf"))

setwd('/home/mark/gitstage/uderica/silage/')

source("uderica_fxns.R")





alpha.div.measures <-
  c("Observed", "Chao1", "Shanon", "Simpson", "InvSimpson")
alpha.div.plot.alab.fontsize <- 6
alpha.point.size <- 5
alpha.point.alpha <- 0.7

# first load data
phylo <- uderica.load(the.exp = my.exp, the.domain = my.domain)


###   ###   ###
# make me a function
otu.dat <- phylo@otu_table
dim(otu.dat)

# could show snippets of otu counts
# otu.tab.excerpt <- otu.dat[order(rownames(otu.dat)),]
# otu.tab.excerpt <- otu.tab.excerpt[,order(as.numeric(colnames(otu.dat)))]
# pander(otu.tab.excerpt[1:9,1:9])

unlisted.otu.counts <- unlist(as.matrix.data.frame(otu.dat))

old.mar <- par("mar")
par("mar" = c(5, 4, 3, 2))

hist(
  log10(unlisted.otu.counts),
  breaks = 100,
  main = "Histogram of reads per sample/OTU",
  xlab = "log10(read counts)"
)
par("mar" = old.mar)
###   ###   ###


###   ###   ###
# make me a function
map.dat <- as.data.frame(phylo@sam_data@.Data)
names(map.dat) <- phylo@sam_data@names
mdncol <- ncol(map.dat)
# assumes a particular map.txt structure (has been consistent so far)
exp.design.matrix <- map.dat[, 5:(mdncol - 1)]

print(pander(table(exp.design.matrix)))
###   ###   ###

###   ###   ###
# make me a function
tax.dat <- as.data.frame(phylo@tax_table@.Data)
temp <- as.matrix(tax.dat)
temp[grepl(pattern = "^uncultured", temp)] <- NA
temp[temp == "unidentified"] <- NA
temp <- !is.na(temp)
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
  main = "Histogram of INFERRED most specific rank",
  xlab = 'Number of non "uncultured" and "unidentified" phylogenetic levels',
  ylim = c(0, (max.count + 20))
)

par("mar" =  old.mar)
###   ###   ###



# optionally filter
phylo.top.N <- top.N.filter(the.phylo = phylo, the.N = 10)
phylo.k.of.A <- k.of.A.filter(the.phylo = phylo, k.factor = 1000)

# barplot
uderica.barplot(the.phylo = phylo.top.N,
                the.factor = my.fact,
                the.rank = "Family")
uderica.barplot(the.phylo = phylo.top.N,
                the.factor = my.fact,
                the.rank = "Genus")
uderica.barplot(the.phylo = phylo.k.of.A,
                the.factor = my.fact,
                the.rank = "Family")
uderica.barplot(the.phylo = phylo.k.of.A,
                the.factor = my.fact,
                the.rank = "Genus")

# alpha/richness
uderica.alpha(the.phylo = phylo,
              the.factor = my.fact,
              the.measures = alpha.div.measures)

# ordination/"CoA"
ord_meths  <- c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
uderica.ordinate(
  the.phylo = phylo,
  the.meths = ord_meths,
  the.dist = "bray",
  the.fact = my.fact,
  connect = TRUE
)

# other distances, some require addirtional arguments
# bray gower (daisy?) and jsd seem to work
# "unifrac"
# Original (unweighted) UniFrac distance, UniFrac
# "wunifrac"
# weighted-UniFrac distance, UniFrac
# "dpcoa"
# sample-wise distance used in Double Principle Coordinate Analysis, DPCoA
# "jsd"


# heatmap, most distant from control
constrast.helper(the.phylo = phylo,
                 the.fact = my.fact,
                 the.ref = my.ref)

# differentail analysis
uderica.diffy(
  the.phylo = phylo,
  the.factor = my.fact,
  the.contrast.num = my.ref,
  the.contrast.denom =  my.other,
  the.alpha = .05,
  the.high.rank = "Phylum",
  the.low.rank = "Class",
  the.y.axis = "log2FoldChange"
)


dev.off()

sink()
sink(type = "message")
