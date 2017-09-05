tab <- data.frame(logFC = res$log2FoldChange,
                 negLogPval = -log10(res$pvalue))
# head(tab)
par(mar = c(5, 4, 4, 4))
plot(
  tab,
  pch = 16,
  cex = 0.6,
  xlab = expression(log[2] ~ fold ~ change),
  ylab = expression(-log[10] ~ pvalue)
)

lfc <- 2
pval <- 0.01

signGenes <- (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes,],
       pch = 16,
       cex = 0.8,
       col = "red")
abline(h = -log10(pval),
       col = "green3",
       lty = 2)
abline(v = c(-lfc, lfc),
       col = "blue",
       lty = 2)
mtext(
  paste("pval =", pval),
  side = 4,
  at = -log10(pval),
  cex = 0.8,
  line = 0.5,
  las = 1
)
mtext(
  c(paste("-", lfc, "fold"), paste("+", lfc, "fold")),
  side = 3,
  at = c(-lfc, lfc),
  cex = 0.8,
  line = 0.5
)
