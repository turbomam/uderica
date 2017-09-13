# load data

library("optparse")
library("phyloseq")
library("ggplot2")
library("plyr")
library("DESeq2")
library("ggthemes")
library("pander")

packageVersion("optparse")
packageVersion("phyloseq")
packageVersion("ggplot2")
packageVersion("plyr")
packageVersion("DESeq2")
packageVersion("ggthemes")
packageVersion("pander")

theme_set(theme_bw())

# list experiments and avialbe domains?  check filesystem?


###   ###   ###

# just iterate?

# experiment <- "Aer_Comp"
# exp.fact <- "Treatment"
# reference.level <- "C/0"

# experiment <- "LalStress"
# exp.fact <- "treatment"
# reference.level <- "C"

# experiment <- "V-HMC"
# exp.fact <- "Treatment"
# reference.level <- "C/0"

###   ###   ###

# domain <- "Bacteria"

domain <- "Fungi"

###   ###   ###

uderica.load <- function(the.exp, the.domain) {
  raw_data_folder <- paste0("raw_data_", tolower(the.domain))
  
  if (the.domain == "Fungi") {
    biom.file.basename <- "otu_table_mc2_w_tax.biom"
  } else if (the.domain == "Bacteria") {
    biom.file.basename <- "otu_table_mc2_w_tax_no_pynast_failures.biom"
    trefile <-
      paste0(
        '/home/mark/gitstage/uderica/silage/',
        the.exp ,
        '/Bacteria',
        '/raw_data_bacteria/flashed/extended/flash_trim_cat_pick/rep_set.tre'
      )
  }
  
  otufile <-
    paste0(
      '/home/mark/gitstage/uderica/silage/',
      the.exp ,
      "/",
      the.domain,
      "/",
      raw_data_folder,
      '/flashed/extended/flash_trim_cat_pick/',
      biom.file.basename
    )
  
  mapfile <-
    paste0('/home/mark/gitstage/uderica/silage/',
           the.exp ,
           '/Fungi/',
           'map.txt')
  
  
  if (the.domain == "Fungi") {
    otutable <-
      import_biom(BIOMfilename = otufile,
                  parseFunction = parse_taxonomy_greengenes)
  } else if (the.domain == "Bacteria") {
    otutable <-
      import_biom(
        BIOMfilename = otufile,
        treefilename = trefile,
        parseFunction = parse_taxonomy_greengenes
      )
  }
  
  mapping <-
    import_qiime_sample_data(mapfilename = mapfile)
  
  map.dat <- as.data.frame(mapping@.Data)
  names(map.dat) <- mapping@names
  
  # confirm factors and levels
  last.col <- ncol(map.dat)
  names(map.dat)[5:(last.col - 1)]
  print(table(map.dat[, 5:(last.col - 1)]))
  
  phylo <- merge_phyloseq(otutable, mapping)
  
  print(phylo)
  
  phylo
  
}

top.N.filter <- function(the.phylo, the.N) {
  taxa.named.vector <- sort(taxa_sums(the.phylo), TRUE)
  taxa.ids <- names(taxa.named.vector)
  top.taxa <- taxa.ids[1:the.N]
  
  prune_taxa(top.taxa, the.phylo)
  
}

k.of.A.filter <- function(the.phylo, k.factor) {
  one.plus <- as.matrix(the.phylo@otu_table@.Data)
  one.plus[one.plus == 0] <- NA
  one.plus <- as.data.frame(one.plus)
  one.plus <- unlist(one.plus)
  
  # not even close to normal distribution
  # one.plus.mean <- mean(one.plus, na.rm = TRUE)
  one.plus.max <- max(one.plus, na.rm = TRUE)
  # one.plus.sd <- sd(one.plus, na.rm = TRUE)
  
  proposed.k <- one.plus.max / k.factor
  
  level.counter <- as.data.frame(the.phylo@sam_data@.Data)
  names(level.counter) <- the.phylo@sam_data@names
  lc.col.count <- ncol(level.counter)
  level.counter <- level.counter[, 5:(lc.col.count - 1)]
  level.counter <- table(level.counter)
  print(level.counter)
  lc.non.zero <- level.counter[level.counter > 0]
  if (length(unique(lc.non.zero)) > 1) {
    table(lc.non.zero)
  } else {
    print(paste0("all levels have ", unique(lc.non.zero), " repeats"))
  }
  proposed.A <- min(lc.non.zero)
  
  prune.func <-
    genefilter_sample(the.phylo,
                      filterfun_sample(function(x)
                        x > proposed.k),
                      A = proposed.A)
  
  
  prune_taxa(prune.func, the.phylo)
}

uderica.barplot <- function(the.phylo, the.factor, the.rank) {
  p <- plot_bar(the.phylo, the.factor, fill = the.rank)
  
  # p + geom_bar(aes(color = Family,
  #                  fill = Family),
  #              stat = "identity")
  
  
  p + geom_bar(aes_(color = as.name(the.rank),
                    fill = as.name(the.rank)),
               stat = "identity")
  
  
}

uderica.alpha <- function(the.phylo, the.factor, the.measures) {
  print(my.exp)
  print(my.domain)
  
  # plot richness without filtering?
  p <-
    plot_richness(the.phylo,
                  x = the.factor,
                  color = the.factor,
                  measures = the.measures)
  
  # parameterize
  p <-
    p + geom_point(size = alpha.point.size, alpha = alpha.point.alpha) + theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = alpha.div.plot.alab.fontsize
    ))
  
  print(p)
  
  alpha.diversity <-
    estimate_richness(phylo,
                      measures = the.measures)
  alpha.diversity$sample.id.string <-
    as.numeric(sub("^X", "", rownames(alpha.diversity)))
  rownames(alpha.diversity) <- NULL
  non.sample.cols <-
    sort(setdiff(names(alpha.diversity), "sample.id.string"))
  alpha.diversity <-
    alpha.diversity[order(alpha.diversity$sample), c("sample.id.string", non.sample.cols)]
  
  timestamp <- format(Sys.time(), format = "%Y%m%d%H%M")
  
  write.table(
    alpha.diversity,
    paste0(my.exp, "_", my.domain, "_", timestamp, "_alpha.txt"),
    row.names = FALSE
  )
  
  # dump to xtable if knitting
  
}

constrast.helper <- function(the.phylo, the.fact, the.ref) {
  # the.phylo <- phylo
  # the.fact <- "treatment"
  # the.ref <- "C"
  
  print(the.phylo)
  print(the.fact)
  print(the.ref)
  
  
  map.dat <- as.data.frame(phylo@sam_data@.Data)
  names(map.dat) <- phylo@sam_data@names
  
  group.dist <- phyloseq::distance(the.phylo, "bray")
  group.dist <- reshape2::melt(as.matrix(group.dist))
  group.dist <- merge(
    x = group.dist,
    y = map.dat[, c("X.SampleID", the.fact)],
    by.x = "Var1",
    by.y = "X.SampleID"
  )
  
  group.dist <- merge(
    x = group.dist,
    y = map.dat[, c("X.SampleID", the.fact)],
    by.x = "Var2",
    by.y = "X.SampleID",
    suffixes = "Var2"
  )
  
  group.dist.agg <-
    aggregate(
      group.dist[, c("value")],
      by = list(group.dist[, 4], group.dist[, 5]),
      FUN = mean,
      na.rm = TRUE
    )
  
  names(group.dist.agg) <- c("Group.1", "Group.2", "distance")
  
  to.report <-
    group.dist.agg[group.dist.agg$Group.1 == the.ref , ]
  to.report <-
    to.report[order(to.report$distance, decreasing = TRUE),]
  print(to.report)
  
  # add titel:  study and domain
  # unfiltered so otu y axis is very busy
  # try different ordiantions or distance matrics?
  # manually group by treatment?
  p <- plot_heatmap(the.phylo, sample.label = the.fact)
  p <- p + ggtitle(paste0(my.exp, ", ", my.domain))
  print(p)
  
}

# refactor... repeating slow DESeq2 step

uderica.diffy <-
  function(the.phylo,
           the.factor,
           the.contrast.num,
           the.contrast.denom,
           the.alpha,
           the.high.rank,
           the.low.rank,
           the.y.axis) {
    the.phylo = phylo
    the.factor = my.fact
    the.contrast.num = my.ref
    the.contrast.denom =  my.other
    the.alpha = 0.05
    the.high.rank = "Phylum"
    the.low.rank = "Class"
    the.y.axis = "log2FoldChange"
    
    
    print(the.phylo)
    
    ###   ###   ###
    
    sample_data(the.phylo)$treatment
    
    ###   ###   ###
    
    # set formula programatically
    
    the.form <- as.formula(paste0(" ~ ", the.factor))
    diagdds <- phyloseq_to_deseq2(the.phylo, the.form)
    
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
    
    diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")
    class(diagdds)
    
    rnms <- resultsNames(diagdds)
    rnms <- sub(pattern = the.factor, replacement = "", rnms)
    print(rnms)
    
    res <-
      results(
        diagdds,
        contrast = c(the.factor, the.contrast.num, the.contrast.denom),
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
    sigtab <- res[which(res$padj < the.alpha),]
    class(sigtab)
    
    sigtab <-
      cbind(as(sigtab, "data.frame"), as(tax_table(the.phylo)[rownames(sigtab), ], "matrix"))
    
    to.report <- sigtab[, c(
      "log2FoldChange",
      "pvalue",
      "padj",
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus"
    )]
    
    print(to.report)
    
    # volcano?
    
    # scale_fill_discrete <- function(palname = "Set1", ...) {
    #   scale_fill_brewer(palette = palname, ...)
    # }
    
    
    # use color an x position to render two different ranks
    # get the rank labels in x
    # make sure they're ordered factors
    
    # AUTOMATE SELECTION OF HIGH AND LOW RANKS?
    # GET AS SPECIFIC AS POSSIBLE WITHOUT INCLUDING A LOT OF UNIDIENFIEDS OR UNCULTUREDS
    
    higher.rank.pos <- which(names(sigtab) == the.high.rank)
    lower.rank.pos <- which(names(sigtab) == the.low.rank)
    
    
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
               x = as.name(the.low.rank),
               y = as.name(the.y.axis),
               color = as.name(the.high.rank)
             )) +
      geom_point(size = 6) +
      theme(axis.text.x = element_text(
        angle = -90,
        hjust = 0,
        vjust = 0.5
      )) + labs(
        title = paste0(
          "Differential abundance ratio, ",
          the.contrast.num ,
          " vs. ",
          the.contrast.denom
        )
      )
    print(p)
    
    timestamp <- format(Sys.time(), format = "%Y%m%d%H%M")
    
    sigtab <- cbind.data.frame(rownames(sigtab), sigtab)
    names(sigtab)[1] <- "OTU"
    
    write.table(
      sigtab,
      file = paste0(
        my.exp,
        "_",
        my.domain,
        "_",
        timestamp,
        "_",
        make.names(the.contrast.num),
        "_vs_",
        make.names(the.contrast.denom),
        ".txt"
      ),
      row.names = FALSE
    )
  }

uderica.ordinate <-
  function(the.phylo,
           the.meths,
           the.dist,
           the.fact,
           connect) {
    # messages when including DPCoA
    # 1) Species coordinates not found directly in ordination object.
    #    Attempting weighted average (`vegan::wascores`)
    # 2) non-unique values when setting 'row.names':
    
    # DPCoA = Double Principle Coordinate Analysis using a (corrected, if necessary)
    #   phylogenetic/patristic distance between species.
    # The calculation is performed by DPCoA(), which ultimately uses dpcoa
    #   after making the appropriate accessions/corrections of the data.
    
    the.meths <- setdiff(the.meths, "DPCoA")
    
    plist <-
      llply(as.list(the.meths), function(i, the.phylo, the.dist) {
        message(i)
        ordi <- ordinate(the.phylo, method = i, distance = the.dist)
        plot_ordination(the.phylo, ordi, "samples", color = the.fact)
      }, the.phylo, the.dist)
    
    names(plist) <- the.meths
    
    pdataframe <- ldply(plist, function(x) {
      df <- x$data[, 1:2]
      colnames(df) <- c("Axis_1", "Axis_2")
      return(cbind(df, x$data))
    })
    
    names(pdataframe)[1] <- "method"
    
    a1 <- "Axis_1"
    a2 <- "Axis_2"
    
    p <-
      ggplot(pdataframe,
             aes_(
               as.name(a1),
               as.name(a2),
               color = as.name(the.fact),
               fill = as.name(the.fact)
             ))
    p <- p + geom_point(size = 2)
    if (connect) {
      p <- p + geom_polygon()
    }
    
    p <- p + facet_wrap(~ method, scales = "free")
    p <- p + scale_fill_brewer(type = "qual", palette = "Set1")
    p <- p + scale_colour_brewer(type = "qual", palette = "Set1")
    print(p)
    
    # user may desire an individual plot on a page by itself
    # so print all on a page by themselves?
    # how did Erica determine that p = plist[[6]] (MDS) was the best ordination?
    
    # these analyses WERE being run at different levels of filtering
    # initial barplot and alpha diversity at 2 (from Bash script)
    # ordination refiltered to 3
    # diff abundance back to 2
    #
    
  }
