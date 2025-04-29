#aldex ttests

# load required library 
library(ALDEx2)

# taxonomy SV level (or aggregate to higher taxonomic levels)
d<-read.table("data/counts_01abun_100filt.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

# EC level
ec<-read.table("PICRUSt2_output_directory/EC_metagenome_out/pred_metagenome_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
ec<-round(ec, digits=0)

# pathway level
p<-read.table("PICRUSt2_output_directory/pathways_out/path_abun_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
p<-round(p, digits=0)


pre<-c("001_pre","002_pre","005_pre","006_pre","007_pre","008_pre","009_pre","010_pre","011_pre","012_pre")
post<-c("001_post","002_post","005_post","006_post","007_post","008_post","009_post","010_post","011_post","012_post")

metR_pre<-c("001_pre","002_pre","005_pre","007_pre","009_pre","011_pre")
metR_post<-c("001_post","002_post","005_post","007_post","009_post","011_post")

metNR_pre<-c("006_pre","008_pre","010_pre","012_pre")
metNR_post<-c("006_post","008_post","010_post","012_post")


#####################################################################
# pre vs post all samples
aldex.in<-d[,c(pre, post)]
conds<-c(rep("pre", length(pre)), rep("post", length(post)))
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)
x.tt <- aldex.ttest(x, paired.test=TRUE)
x.effect <- aldex.effect(x)

#merge the data
x.all <- data.frame(x.tt, x.effect)

taxa<-read.table("data/cutadapt_tax2.txt", sep="\t", quote="", header=TRUE, row.names=1, comment.char="")
x.all.tax.sv<-merge(x.all, taxa[,2:7], by='row.names', all.y=FALSE)
#test <- merge(x.all, taxa, by=0)

#write a .txt with the results
#write.table(x.all.tax, file="aldex_output/pre_vspost_SV.txt", sep="\t", quote=F, col.names=NA)

#####################################################################
# make a heatmap of these results

aldex_out<-read.table("aldex_output/pre_vspost_SV.txt", sep="\t", quote="", header=T, row.names=1)
filt_effect <- aldex_out[abs(aldex_out$effect) > 0.5, ]

# relative abundance
abun <- apply(aldex.in, 2, function(x) {x/sum(x)})
colSums(abun) #check all sum to 1

#transpose to samples as rows
abun.t<-t(abun)

#samples must be as rows
d.czm <- cmultRepl(abun.t,  label=0, method="CZM", z.warning = 0.95)

# The table needs to be transposed again (samples as COLUMNS)
d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))}))

#only get the differentially abundant SV
d.clr.diff<-d.clr[,colnames(d.clr) %in% filt_effect$Row.names]

library(reshape2)
clr_melt<-melt(d.clr.diff)
colnames(clr_melt)<-c("sample","SV","clr")

#add a participant and timepoint column
clr_melt$participant<-gsub("_.*","", clr_melt$sample)
clr_melt$timepoint<-gsub(".*_","",clr_melt$sample)
clr_melt$timepoint<-gsub("pre","1pre",clr_melt$timepoint)

filt_effect <- filt_effect %>%
  mutate(bugnames = ifelse(!is.na(Species), paste(Row.names, Genus, Species, sep = "_"), paste(Row.names, Genus, sep = "_")))

filt_effect<- filt_effect %>% arrange(effect)

clr_final<- clr_melt %>%
  left_join(filt_effect %>% select(Row.names, bugnames), by = c("SV" = "Row.names"))

clr_final <- clr_final %>%
  left_join(filt_effect %>% select(Row.names, effect), by =c("SV"="Row.names"))

clr_final$effect<-clr_final$effect * -1
clr_final <- clr_final %>%
  arrange(effect) %>%
  mutate(bugnames = factor(bugnames, levels = unique(bugnames)))  # Preserve order in factor

filt_effect$effect<-filt_effect$effect *-1

filt_effect<-filt_effect %>% mutate(bugnames = factor(bugnames, levels = unique(bugnames)))

heat <- ggplot(clr_final, aes(x = timepoint, y = bugnames, fill = clr)) +
  geom_tile() +
  facet_grid(~participant, scales = "free_x") +
  scale_fill_gradientn(
    colors = rev(c("#98033A","#D43747","#F76340","#FEA55C","#FEFFBA","white","#9CDA9C","#267AB1","#5a78b5","#554394")),
    values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), # Keeping the same values
    limits = range(clr_final$clr)) +
  scale_x_discrete(labels=c("1pre"="Pre", "post"="Post")) +
  theme_bw() +
  theme(axis.title = element_blank(), 
        legend.position = c(-0.2, 1.05), 
        legend.direction = "horizontal",
        axis.text.y=element_text(face="italic"),
        axis.text.x = element_text(size=8))
heat

bars <- ggplot(filt_effect, aes(x = effect, y = bugnames, fill = ifelse(effect > 0.5, "salmon", "#5a78b5"))) +
  geom_col() +
  scale_fill_identity() +  # Directly apply colors without expecting discrete categories
  geom_text(aes(label = round(effect, 2), x = ifelse(effect > 0.5, 0.1, -0.1)),  # Position text dynamically
            color = "black", hjust = ifelse(filt_effect$effect > 0.5, 0, 1), size = 3) +  # Adjust text alignment
  scale_y_discrete(limits = rev(unique(filt_effect$bugnames))) +  # Keep order consistent with heatmap
  theme_minimal() +
  xlab("Effect Size") +
  theme(
    axis.text.y = element_blank(),  # Hide y-axis text to avoid duplication
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank())

bars

final_plot<- heat + bars + plot_layout(widths = c(4, 1))
final_plot

#########################################
# pre R vs NR

aldex.in<-d[,c(metR_pre, metNR_pre)]
conds<-c(rep("metR_pre", length(metR_pre)), rep("metNR_pre", length(metNR_pre)))
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)
x.tt <- aldex.ttest(x, paired.test=FALSE)
x.effect <- aldex.effect(x)

#merge the data
x.all <- data.frame(x.tt, x.effect)

taxa<-read.table("data/cutadapt_tax2.txt", sep="\t", quote="", header=TRUE, row.names=1, comment.char="")
x.all.tax<-merge(x.all, taxa[,2:7], by='row.names', all.y=FALSE)

write.table(x.all, file="aldex_output/metR_vsNR_preSV.txt", sep="\t", quote=F, col.names=NA)
#repeat for EC and pathways at both pre and post

#########################################
# post R vs NR

aldex.in<-d[,c(metR_post, metNR_post)]
conds<-c(rep("metR_post", length(metR_post)), rep("metNR_post", length(metNR_post)))
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)
x.tt <- aldex.ttest(x, paired.test=FALSE)
x.effect <- aldex.effect(x)

#merge the data
x.all <- data.frame(x.tt, x.effect)

taxa<-read.table("data/cutadapt_tax2.txt", sep="\t", quote="", header=TRUE, row.names=1, comment.char="")
x.all.tax<-merge(x.all, taxa[,2:7], by='row.names', all.y=FALSE)

write.table(x.all.tax, file="aldex_output/metR_vsNR_postSV.txt", sep="\t", quote=F, col.names=NA)

#####################################################################
#####################################################################

### repeat with picrust outputs
#####################################################################
# EC level
ec<-read.table("PICRUSt2_output_directory/EC_metagenome_out/pred_metagenome_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
ec<-round(ec, digits=0)

# pre vs post all samples
aldex.in<-ec[,c(pre, post)]
conds<-c(rep("pre", length(pre)), rep("post", length(post)))
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)
x.tt <- aldex.ttest(x, paired.test=TRUE)
x.effect <- aldex.effect(x)

#merge the data
x.all <- data.frame(x.tt, x.effect)
write.table(x.all, file="aldex_output/prevspost_ECnumbers.txt", sep='\t', quote=F)

#####################################################################
# KO level
ko<-read.table("PICRUSt2_output_directory/KO_metagenome_out/pred_metagenome_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
ko<-round(ko, digits=0)

# pre vs post all samples
aldex.in<-ko[,c(pre, post)]
conds<-c(rep("pre", length(pre)), rep("post", length(post)))
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)
x.tt <- aldex.ttest(x, paired.test=TRUE)
x.effect <- aldex.effect(x)

#merge the data
x.all.ko <- data.frame(x.tt, x.effect)
write.table(x.all.ko, file="aldex_output/prevspost_KOnumbers.txt", sep='\t', quote=F)

#####################################################################
# pathways level
p<-read.table("PICRUSt2_output_directory/pathways_out/path_abun_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
p<-round(p, digits=0)

# pre vs post all samples
aldex.in<-p[,c(pre, post)]
conds<-c(rep("pre", length(pre)), rep("post", length(post)))
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)
x.tt <- aldex.ttest(x, paired.test=TRUE)
x.effect <- aldex.effect(x)

#merge the data
x.all <- data.frame(x.tt, x.effect)
write.table(x.all, file="aldex_output/prevspost_pathways.txt", sep='\t', quote=F)

################################################################################
################################################################################
################################################################################


# volcano plot of aldex results
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(zCompositions)

################################################################################
# ec numbers
ec_num<-read.table("aldex_output/prevspost_ECnumbers.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
ec_num$effect<-ec_num$effect * -1

#colours based on significance!
# add a column
ec_num$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
ec_num$diffexpressed[ec_num$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
ec_num$diffexpressed[ec_num$effect < -0.5] <- "DOWN"

#labels based on significance!
ec_num$delabel <- NA
#ec_num$delabel[ec_num$diffexpressed != "NO"] <- rownames(ec_num)[ec_num$diffexpressed != "NO"]
ec_num$delabel[ec_num$effect > 0.65 | ec_num$effect < -0.6] <- rownames(ec_num)[ec_num$effect > 0.65 | ec_num$effect < -0.65]

ec_volcano <- ggplot(data=ec_num, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 22)
ec_volcano	

################################################################################
# KO numbers
ko_num<-read.table("aldex_output/prevspost_KOnumbers.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
ko_num$effect<-ko_num$effect * -1 #multiply by -1 so those that are enriched post-acv are positive

#colours based on significance!
# add a column
ko_num$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
ko_num$diffexpressed[ko_num$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
ko_num$diffexpressed[ko_num$effect < -0.5] <- "DOWN"

#labels based on significance!
ko_num$delabel <- NA
#ko_num$delabel[ko_num$diffexpressed != "NO"] <- rownames(ko_num)[ko_num$diffexpressed != "NO"]
ko_num$delabel[ko_num$effect > 0.65 | ko_num$effect < -0.65] <- rownames(ko_num)[ko_num$effect > 0.65 | ko_num$effect < -0.65]

ko_volcano <- ggplot(data=ko_num, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 22)
ko_volcano	

################################################################################
# pathways 
d<-read.table("aldex_output/pre_vspost_pathways.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d$effect<-d$effect * -1

#colours based on significance!
# add a column
d$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
d$diffexpressed[d$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
d$diffexpressed[d$effect < -0.5] <- "DOWN"

#labels based on significance!
d$delabel <- NA
d$delabel[d$diffexpressed != "NO"] <- rownames(d)[d$diffexpressed != "NO"]

path_volcano <- ggplot(data=d, aes(x=effect, y=-log10(wi.ep), col=diffexpressed)) +
  geom_point(aes(alpha = ifelse(diffexpressed == "NO", 0.5, 1))) +
  theme_bw() +
  scale_color_manual(values=c("#5a78b5", "black", "salmon")) +
  theme(legend.position="none") +
  labs(x="Effect Size", y="Wilcoxon P (-log10)") +
  geom_text_repel(aes(label = delabel), size = 3, max.overlaps = 30)
path_volcano

##########################################
# Individual pathway plots
p<-read.table("PICRUSt2_output_directory/pathways_out/path_abun_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
p<-round(p, digits=0)

# relative abundance
p.f <- apply(p, 2, function(x) {x/sum(x)})
colSums(p.f) #check all sum to 1

to.remove<-c("003_post","004_pre","DNA_BLANK2","PCR_BLANK2")

p.f2<-p.f[,!colnames(p.f) %in% to.remove]

#transpose to samples as rows
p.f3<-t(p.f2)

#samples must be as rows
p.czm <- cmultRepl(p.f3,  label=0, method="CZM")
p.clr <- t(apply(p.czm, 1, function(x){log(x) - mean(log(x))}))

meta<-read.table("data/metadata.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

p.clr.m<-merge(p.clr, meta, by=0, all=FALSE)
p.clr.m$participant<-as.factor(p.clr.m$participant)

p1 <- ggplot(p.clr.m, aes(x = timepoint, y = p.clr.m[[163]], group = participant, colour = participant)) +
  geom_point(aes(fill = participant), size = 3, shape = 21, stroke = 1, alpha = 0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "PWY-5913", y = "CLR relative abundance") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  theme_bw() +
  annotate("text", x = 1.5, y = 1.25, label = "ES = -0.65", size = 3) +
  theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=10))

p2 <- ggplot(p.clr.m, aes(x = timepoint, y = p.clr.m[[250]], group = participant, colour = participant)) +
  geom_point(aes(fill = participant), size = 3, shape = 21, stroke = 1, alpha = 0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "PWY-7328", y = "CLR relative abundance") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  annotate("text", x = 1.5, y = 0.5, label = "ES = -0.61", size = 3) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=10))

p3 <- ggplot(p.clr.m, aes(x = timepoint, y = p.clr.m[[108]], group = participant, colour = participant)) +
  geom_point(aes(fill = participant), size = 3, shape = 21, stroke = 1, alpha = 0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "PWY-3781", y = "CLR relative abundance") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  annotate("text", x = 1.5, y = 0.5, label = "ES = -0.58", size = 3) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=10))

p4 <- ggplot(p.clr.m, aes(x = timepoint, y = p.clr.m[[80]], group = participant, colour = participant)) +
  geom_point(aes(fill = participant), size = 3, shape = 21, stroke = 1, alpha = 0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "P125-PWY", y = "CLR relative abundance") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  annotate("text", x = 1.5, y = 0.5, label = "ES = -0.57", size = 3) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=10))

p5 <- ggplot(p.clr.m, aes(x = timepoint, y = p.clr.m[[64]], group = participant, colour = participant)) +
  geom_point(aes(fill = participant), size = 3, shape = 21, stroke = 1, alpha = 0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "LACTOSECAT-PWY", y = "CLR relative abundance") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  annotate("text", x = 1.5, y = 0.5, label = "ES = -0.55", size = 3) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(size=10))

################################################################################
# combination panel
################################################################################

layout <- "
AAAACDE
BBBBFG#
"

ec_volcano + ko_volcano + p1 + p2 + p3 + p4 + p5 +
  plot_layout(design = layout)
