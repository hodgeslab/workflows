#!/usr/bin/env Rscript

gmtFile <- "/s1/share/msigdb-6.2/msigdb_v6.2_GMTs/custom/h_c2_c6.all.v6.2.symbols.gmt"
outDir <- "gsea"
padjThresh <- 0.05
topN <- 40

library(fgsea)
library(ggplot2)

dir.create(file.path(outDir), showWarnings = F)

# load ranked data
rnkFile <- "gseaData.rnk"
ranks <- read.table(rnkFile,
                    header=F, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$V2, ranks$V1)
str(ranks)

# load gmt file
pathways <- gmtPathways(gmtFile)
str(head(pathways))

# perform gsea
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1e6, nproc=4)
fn <- paste0(outDir,"/summary.txt")

fgseaResOut <- fgseaRes
fgseaResOut$leadingEdge <- lapply(fgseaResOut$leadingEdge,toString)
write.table(as.matrix(fgseaResOut),fn,quote=F,sep="\t",row.names=F)

head(fgseaRes[order(pval), ])

# count padj results
sum(fgseaRes[, padj < padjThresh])

# collapse pathways to remove redundant pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < padjThresh], 
                                      pathways, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                         order(-NES), pathway]

# topPathwaysUp <- fgseaRes[ES > 0][head(order(pval,-NES), n=topN), pathway]
topPathwaysUp <- fgseaRes[ES > 0][head(order(NES, decreasing=T), n=topN), pathway]
# topPathwaysDown <- fgseaRes[ES < 0][head(order(pval,NES), n=topN), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(NES, decreasing=F), n=topN), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
topPathways <- topPathways[topPathways %in% mainPathways]

fn <- paste0(outDir,"/summary.pdf")
pdf(fn, height=8, width=10 , useDingbats=F)
plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5)
dev.off()

for(sig in topPathways) {
  fn <- paste0(outDir,"/",sig,".pdf")
  pdf(fn,height=2.5, width=3.5,useDingbats=F)
  p1 <- plotEnrichment(pathways[[sig]], ranks, gseaParam = 1, ticksSize=0.1) +
    labs(title=sig) +
    theme(plot.title = element_text(size=6))
  print(p1)
  dev.off()
}
