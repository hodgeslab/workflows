#!/usr/bin/env Rscript

library(ggplot2)
library(hodgeslabR)

args <- commandArgs(TRUE)
inFn <- args[1]
outFn <- args[2]

df <- read.table(inFn, header=F, stringsAsFactors=F, sep="\t")
colnames(df) <- c("len","count")
df$count <- df$count/max(df$count)

p1 <- ggplot(df,aes(x=len,y=count)) +
  theme_hodgeslab_basic() +
  geom_step(size=0.5*linescale) +
  xlim(c(0,1000)) +
  xlab("Fragment size (bp)") +
  ylab("Frequency (norm.)")

pdf(outFn,height=1.5,width=1.75,useDingbats=F)
print(p1)
dev.off()
