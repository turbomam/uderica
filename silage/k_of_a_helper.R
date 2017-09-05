# k of A "helper"

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

exp.factor <- "Treatment"

# max significance for reporting and plotting in DESeq2
DESeq2.alpha <- 0.05

# what two ranks should eb sued to group the OTUs
higher.rank.name <- "Family"
lower.rank.name <- "Genus"

# what values from DESeq analysis should be plotted
y.val.name <- "log2FoldChange"


experiment <- "Aer_Comp"
# experiment <- "LalStress"
# experiment <- "V-HMC"

# denominator... bottom portion of ratio
if (experiment == "LalStress") {
  constrat.num <- "C"
  constrat.den <- "LB500.LATE"
  
} else if (experiment == "V-HMC") {
  constrat.num <- "C"
  constrat.den <- "LB"
  
} else if (experiment == "Aer_Comp") {
  constrat.num <- "C"
  constrat.den <- "LB"
  
}

otufile <-
  paste0(
    '/home/mark/gitstage/uderica/silage/',
    experiment ,
    '/Fungi/raw_data_fungi/flashed/extended/flash_trim_cat_pick/otu_table_mc2_w_tax.biom'
  )

mapfile <-
  paste0('/home/mark/gitstage/uderica/silage/',
         experiment ,
         '/Fungi/',
         'map.txt')

otutable <-
  import_biom(BIOMfilename = otufile,
              parseFunction = parse_taxonomy_greengenes)

mapping <-
  import_qiime_sample_data(mapfilename = mapfile)

phylo <- merge_phyloseq(otutable, mapping)

print(phylo)

otu.dat <- otu_table(phylo)
otu.dat <- otu.dat@.Data

otu.melt <- reshape2::melt(otu.dat)
names(otu.melt) <-  c("OTU", "sample", "count")


otu.melt$log10count <- log10(otu.melt$count)
otu.melt <- otu.melt[otu.melt$log10count > 0 , ]

lastfloor <- floor(max(otu.melt$log10count))
floor.cutoffs <- seq(from = 0, to = lastfloor, by = 0.5)

success.list <-
  lapply(floor.cutoffs, function(cutoff) {
    print
    otu.melt$log10count > cutoff
  })
names(success.list) <- floor.cutoffs
success.list <- do.call(cbind.data.frame, success.list)

otu.bound <- cbind.data.frame(otu.melt, success.list)

keeper.cols <-
  setdiff(names(otu.bound), c("count", "log10count"))
otu.bound <- otu.bound[, keeper.cols]
otu.bound <- reshape2::melt(otu.bound, id = c("OTU", "sample"))
otu.bound <- otu.bound[otu.bound$value,]
otu.bound$value <- as.numeric(otu.bound$value)

names(otu.bound) <- c("OTU", "sample", "count.floor", "value")

otu.bound.agg <-
  aggregate(
    otu.bound$value,
    by = list(otu.bound$OTU, otu.bound$count.floor),
    FUN = sum
  )

names(otu.bound.agg) <- c("OTU", "count.floor", "num.samples")

otu.re.agg <-
  aggregate(
    otu.bound.agg$OTU,
    by = list(otu.bound.agg$num.samples, otu.bound.agg$count.floor),
    FUN = length
  )

names(otu.re.agg) <- c("num.samples", "count.floor", "num.OTUs")
otu.re.agg$log10.num.OTUs <- log10(otu.re.agg$num.OTUs)

p <-
  ggplot(otu.re.agg, aes(num.samples, count.floor)) +
  geom_tile(aes(fill = log10.num.OTUs), colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue")

p
