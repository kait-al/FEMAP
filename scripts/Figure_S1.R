# aldex effect volcano plots from metabolic responder comparisons
library(ggplot2)
library(ggrepel)
library(gridExtra)
 
pre.sv<-read.table("aldex_output/metR_vsNR_preSV.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
#colours based on significance!
# add a column
pre.sv$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
pre.sv$diffexpressed[pre.sv$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
pre.sv$diffexpressed[pre.sv$effect < -0.5] <- "DOWN"
#labels based on significance!
pre.sv$delabel <- NA
pre.sv$delabel[pre.sv$diffexpressed != "NO"] <- rownames(pre.sv)[pre.sv$diffexpressed != "NO"]

pre.sv.v <- ggplot(data=pre.sv, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 30)
pre.sv.v

post.sv<-read.table("aldex_output/metR_vsNR_postSV.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
#colours based on significance!
# add a column
post.sv$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
post.sv$diffexpressed[post.sv$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
post.sv$diffexpressed[post.sv$effect < -0.5] <- "DOWN"
#labels based on significance!
post.sv$delabel <- NA
post.sv$delabel[post.sv$diffexpressed != "NO"] <- rownames(post.sv)[post.sv$diffexpressed != "NO"]

post.sv.v <- ggplot(data=post.sv, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 30)
post.sv.v

pre.ec<-read.table("aldex_output/metR_vsNR_preECs.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
#colours based on significance!
# add a column
pre.ec$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
pre.ec$diffexpressed[pre.ec$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
pre.ec$diffexpressed[pre.ec$effect < -0.5] <- "DOWN"
#labels based on significance!
pre.ec$delabel <- NA
pre.ec$delabel[pre.ec$diffexpressed != "NO"] <- rownames(pre.ec)[pre.ec$diffexpressed != "NO"]

pre.ec.v <- ggplot(data=pre.ec, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 30)
pre.ec.v

post.ec<-read.table("aldex_output/metR_vsNR_postEC.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
#colours based on significance!
# add a column
post.ec$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
post.ec$diffexpressed[post.ec$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
post.ec$diffexpressed[post.ec$effect < -0.5] <- "DOWN"
#labels based on significance!
post.ec$delabel <- NA
post.ec$delabel[post.ec$diffexpressed != "NO"] <- rownames(post.ec)[post.ec$diffexpressed != "NO"]

post.ec.v <- ggplot(data=post.ec, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 30)
post.ec.v

pre.p<-read.table("aldex_output/metR_vsNR_prePathways.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
#colours based on significance!
# add a column
pre.p$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
pre.p$diffexpressed[pre.p$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
pre.p$diffexpressed[pre.p$effect < -0.5] <- "DOWN"
#labels based on significance!
pre.p$delabel <- NA
pre.p$delabel[pre.p$diffexpressed != "NO"] <- rownames(pre.p)[pre.p$diffexpressed != "NO"]

pre.p.v <- ggplot(data=pre.p, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 30)
pre.p.v

post.p<-read.table("aldex_output/metR_vsNR_postPathways.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
#colours based on significance!
# add a column
post.p$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
post.p$diffexpressed[post.p$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
post.p$diffexpressed[post.p$effect < -0.5] <- "DOWN"
#labels based on significance!
post.p$delabel <- NA
post.p$delabel[post.p$diffexpressed != "NO"] <- rownames(post.p)[post.p$diffexpressed != "NO"]

post.p.v <- ggplot(data=post.p, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 30)
post.p.v

grid.arrange(pre.sv.v, post.sv.v, pre.ec.v, post.ec.v, pre.p.v, post.p.v, nrow=3)
