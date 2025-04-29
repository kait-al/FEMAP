# Script to make panel figure containing the barplot, pca, aitchison distance comparison, and alpha diversity

# Load required libraries
library(ggplot2)
library(RColorBrewer)
library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)
library(gridExtra)
library(stringr)
library(reshape2)
library(forcats)
library(plyr)
library(sjmisc)
library(rlist)
library(zCompositions)

#######################################################################################################
##### A, Barplot

# Read input data
d <- read.table("data/counts_01abun_100filt.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

# Extract and clean taxonomic information
taxa <- data.frame(str_split_fixed(d$tax.vector, ":", 6))
taxa$Genus <- paste(taxa$X2, taxa$X3, taxa$X4, taxa$X5, taxa$X6, sep="_")

# Merge cleaned taxa with count data
d1 <- cbind(d[,1:(ncol(d)-1)], Genus=taxa$Genus)

# Aggregate at the genus level
gen <- ddply(d1, "Genus", numcolwise(sum))
row.names(gen) <- gen$Genus
gen$Genus <- NULL

# Calculate relative abundance
gen.f <- apply(gen, 2, function(x) x / sum(x))

# Order by abundance
y1 <- gen.f[order(rowSums(gen.f), decreasing = TRUE),]

# Filter taxa above 1% abundance
abund <- 0.01
keep.taxa.index <- rownames(y1[rowMeans(y1) > abund,])

# Retain only the relevant taxa and calculate remainder
y3 <- as.data.frame(y1) %>% filter(rownames(.) %in% keep.taxa.index)
remainder <- colSums(y1[!rownames(y1) %in% keep.taxa.index, , drop = FALSE])
y3 <- rbind(y3, remainder)
rownames(y3)[nrow(y3)] <- "remainder"

# Remove unwanted columns
y3 <- y3 %>% select(-"003_post", -"004_pre", -"DNA_BLANK2", -"PCR_BLANK2")

# Convert to matrix and reshape for plotting
y3 <- as.matrix(y3)
melted <- melt(y3)
colnames(melted) <- c("Genus", "Sample", "value")

# Shorten genus names
gen_name <- data.frame(str_split_fixed(melted$Genus, "_", 5))
melted$Genus <- ifelse(gen_name$X5 == "", "Genera <1% abundance", gen_name$X5)

# Custom sample ordering
custom_order <- function(samples) {
  df <- data.frame(
    Sample = samples,
    Participant = gsub("_(pre|post)", "", samples),
    Visit = ifelse(grepl("_pre", samples), "A", "B")
  ) %>% arrange(Participant, Visit)
  return(df$Sample)
}
melted$Sample <- factor(melted$Sample, levels = custom_order(unique(melted$Sample)))

# Define color palettes
pal1 <- rep(c("lightgray", "#006F6B", "#00ADAB", "#ACDEE0", "#3E783A", "#76BB47", "#AAD486",
              "#CBC02D", "#FCF281", "#C36928", "#F58123", "#F8A96E", "#F9DDCA", "#851719",
              "#B01F24", "#ED2027", "#F3786D", "#4A2970", "#7E6BAD", "#B8ACD1", "#A71D47",
              "#ED1F6B"))
pal2 <- rlist::list.reverse(pal1)

# Plot with faceting
melted$Sample <- gsub("X", "", melted$Sample)
melted$Participant <- sub("_.*", "", melted$Sample)
melted$Timepoint <- sub(".*_", "", melted$Sample)
melted$Timepoint <- gsub("pre", "Pre", melted$Timepoint)
melted$Timepoint <- gsub("post", "Post", melted$Timepoint)
melted$Timepoint <- factor(melted$Timepoint, levels = c("Pre", "Post"))

bars <- ggplot(melted, aes(x=Timepoint, y=value, fill=fct_reorder(Genus, value))) +
  geom_bar(stat="identity", position="stack") + 
  scale_fill_manual(values=pal2) +
  facet_grid(~Participant, scales="free_x") +
  ylab("Relative Abundance") +
  guides(fill = guide_legend(ncol=2, reverse=T, title="Genus")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=8), legend.position="right", axis.title.x=element_blank(), legend.text=element_text(size=8, face="italic"), legend.title=element_text(size=10), axis.title=element_text(size=10))
bars

#######################################################################################################
##### B, PCA

tax <- d$tax.vector

# Split to the 6th taxonomic level -> genus (separated by :)
split6 <- sapply(strsplit(as.character(tax), ":"), "[", 6)
split6 <- as.data.frame(split6)
rownames(split6)<-rownames(d)
split6$sv<-rownames(d)
split6$sv_split6<-paste(split6$sv, split6$split6, sep='_')
rownames(split6)<-rownames(d)

#get only the count columns
dm <- d[,1:ncol(d)-1]

gen.f <- apply(dm, 2, function(x) {x/sum(x)})
colSums(gen.f) #check all sum to 1

# Apply compositional data transformation (CZM method) and CLR transformation
d.czm <- cmultRepl(t(gen.f),  label=0, method="CZM")
d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))}))

# Compute Aitchison distances
aitch_dists <- as.matrix(dist(d.clr))
#write.table(aitch_dists, file="data/aitchdist_filtered.txt", sep="\t", quote=F)

# Perform PCA
d.pcx <- prcomp(d.clr)

# Define parameters for PCA
sv_positions <- data.frame(d.pcx[["rotation"]])
# Merge on genus table
sv_pos <- merge(sv_positions, split6, by = 0)
# Calculate Euclidean distance
arrow_len <- function(x, y) {
  sqrt((x - 0)^2 + (y - 0)^2)
}
sv_pos_dist <- mapply(arrow_len, sv_pos$PC1, sv_pos$PC2)
sv_pos_dist <- as.data.frame(cbind(sv_pos_dist, sv_pos$Row.names, sv_pos$sv_split6, sv_pos$PC1, sv_pos$PC2))
colnames(sv_pos_dist) <- c("Distance", "SV", "Genus", "PC1", "PC2")
# Filter arrow based on distance
filter <- sv_pos_dist %>% filter(Distance >= 0.15)

d.mvar <- sum(d.pcx$sdev^2)
# Calculate the PC1 and PC2 variance
PC1 <- paste("PC1: ", round(sum(d.pcx$sdev[1]^2)/d.mvar, 3))
PC2 <- paste("PC2: ", round(sum(d.pcx$sdev[2]^2)/d.mvar, 3))

meta<-read.table("data/metadata.txt", header=T, row.names = 1, sep='\t', comment.char = "")
metadata<-tibble::rownames_to_column(meta, "SampleID")

loadings<- data.frame(Variables=rownames(d.pcx$rotation), d.pcx$rotation)
values<-merge(d.pcx$x[,c(1,2)], metadata[,c("SampleID","participant","timepoint","metabolic_responder")],
              by.x="row.names", by.y="SampleID", all=F)

# Remove unwanted rows
values2 <- values[!(values$Row.names %in% c("003_post", "004_pre","DNA_BLANK2","PCR_BLANK2")), ]
values2$time<-gsub("._","",values2$timepoint)

# Convert participant to a factor so it's treated as discrete
values2$participant <- as.factor(values2$participant)


# Create the plot
p <- ggplot(values2, aes(x = PC1, y = PC2)) +
  geom_segment(data = sv_pos, aes(x = 0, y = 0, xend = 50 * PC1, yend = 50 * PC2),
               arrow = arrow(length = unit(1/2, 'picas')),
               color = "grey69", alpha = 0.8, size = 0.15) + # Plot features
  geom_text(data = filter, aes(x = 50 * as.numeric(PC1), y = 50 * as.numeric(PC2)), 
            nudge_x = 0.25, nudge_y = 0.25, label = filter$Genus, 
            color = "grey69", size = 3, fontface = "italic") +
  geom_point(aes(colour = participant), size = 6) +  # Points colored by participant (Viridis)
  geom_text(aes(label = time), size = 2.5) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ", round(100 * (d.pcx$sdev[1]^2 / sum(d.pcx$sdev^2)), 1), "%")) +
  ylab(paste0("PC2: ", round(100 * (d.pcx$sdev[2]^2 / sum(d.pcx$sdev^2)), 1), "%")) +
  labs(colour="Participant") + 
  theme_bw() +
  theme(legend.position = c(0.1, 0.67), legend.background = element_blank(), legend.key = element_blank(), axis.title=element_text(size=10), legend.title = element_text(size=10))
p

#######################################################################################################
##### C, Pre vs post aitchison distance

dm<-read.table("data/aitchdist_filtered.txt", sep='\t', header=TRUE, row.names=1, quote="")
colnames(dm)<-gsub("X","", colnames(dm))

time_dist_values <- list()

# Loop through the dm matrix to capture the filtered values
for(i in 1:nrow(dm)) {
  for(j in 1:ncol(dm)) {
    currentrow <- unlist(strsplit(rownames(dm)[i], "_"))
    currentcol <- unlist(strsplit(colnames(dm)[j], "_"))
    
    # Apply the condition: same identifier (currentrow[1] == currentcol[1]), pre in row, not equal in the second part
    if(currentrow[1] == currentcol[1] 
       && currentrow[2] != currentcol[2] 
       && currentrow[2] == "pre") {
      # Store the value into the list
      time_dist_values[[currentrow[1]]] <- c(time_dist_values[[currentrow[1]]], dm[i, j])
    }
  }
}

# Convert the list into a data frame and set appropriate column names
time_dist_df <- as.data.frame(time_dist_values)
time_test<-t(time_dist_df)

colnames(time_test)<-"aitch"
rownames(time_test)<-gsub("X","", rownames(time_test))
rownames(time_test)<-paste0(rownames(time_test), "_post")

dist_merge<-merge(meta, time_test, by=0, all=FALSE)
dist_merge$participant<-as.factor(dist_merge$participant)


aitch_plot <- ggplot(dist_merge, aes(x = participant, y = aitch, fill = factor(participant))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Paired") +
  labs(y = "Pre vs. Post Aitchison Distance", x = "Participant") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 10), axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +  # Remove x-axis labels
  geom_text(aes(x = participant, y = 2, label = participant), 
            color = "white", size = 3, inherit.aes = FALSE)  # Adjust size as needed

aitch_plot

interindiv_pre_values <- list()
# Loop through the d1 matrix to capture the filtered values
for(i in 1:nrow(d1)) {
  for(j in 1:ncol(d1)) {
    currentrow <- unlist(strsplit(rownames(d1)[i], "_"))
    currentcol <- unlist(strsplit(colnames(d1)[j], "_"))
    
    # Apply the condition: same identifier (currentrow[1] == currentcol[1]), pre in row, not equal in the second part
    if(currentrow[1] != currentcol[1] 
       && currentrow[2] == currentcol[2] 
       && currentrow[2] == "pre") {
      # Store the value into the list
      interindiv_pre_values[[currentrow[1]]] <- c(interindiv_pre_values[[currentrow[1]]], d1[i, j])
    }
  }
}

interindiv_post_values <- list()
# Loop through the d1 matrix to capture the filtered values
for(i in 1:nrow(d1)) {
  for(j in 1:ncol(d1)) {
    currentrow <- unlist(strsplit(rownames(d1)[i], "_"))
    currentcol <- unlist(strsplit(colnames(d1)[j], "_"))
    
    # Apply the condition: same identifier (currentrow[1] == currentcol[1]), post in row, not equal in the second part
    if(currentrow[1] != currentcol[1] 
       && currentrow[2] == currentcol[2] 
       && currentrow[2] == "post") {
      # Store the value into the list
      interindiv_post_values[[currentrow[1]]] <- c(interindiv_post_values[[currentrow[1]]], d1[i, j])
    }
  }
}

# Convert lists to data frames
intra_df <- data.frame(Distance = unlist(intra_values), Category = "Intra")
inter_pre_df <- data.frame(Distance = unlist(interindiv_pre_values), Category = "Inter-Pre")
inter_post_df <- data.frame(Distance = unlist(interindiv_post_values), Category = "Inter-Post")

distance_data <- bind_rows(intra_df, inter_pre_df, inter_post_df) # Combine data frames

distance_data$Category <- factor(distance_data$Category, levels = c("Intra", "Inter-Pre", "Inter-Post"))

dist_comp<-ggplot(distance_data, aes(x = Category, y = Distance)) +
  geom_violin(trim = FALSE, alpha = 0.7) + 
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +  # Add boxplot inside violin
  theme_bw() +
  #scale_fill_manual(values = c("Intra" = "#1f77b4", "Inter-Pre" = "#ff7f0e", "Inter-Post" = "#2ca02c")) +  # Custom colors
  labs(y = "Aitchison Distance") +
  theme(legend.position = "none", text = element_text(size = 10), axis.title.x = element_blank())


#######################################################################################################
##### D and E, alpha diversity

div<-read.table("data/alpha_diversity.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
div$participant<-as.factor(div$participant)
div<-div[!((div$Row.names) %in% c("003_post", "004_pre")),]

shan<-ggplot(div, aes(x=timepoint, y=Shannon, group = participant, colour = participant)) +
  geom_point(aes(fill=participant), size=4, shape=21, stroke=1, alpha=0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette="Paired") +
  scale_fill_brewer(palette="Paired") +
  labs(y="Shannon's Index") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position ="none", axis.title=element_text(size=10))
shan

core<-ggplot(div, aes(x=timepoint, y=core_abundance, group = participant, colour = participant)) +
  geom_point(aes(fill=participant), size=4, shape=21, stroke=1, alpha=0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette="Paired") +
  scale_fill_brewer(palette="Paired") +
  labs(y="Core Abundance") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position="none", axis.title=element_text(size=10))
core

#######################################################################################################
#######################################################################################################
#### final panel layout

lay<-rbind(c(1,1,1,1),c(1,1,1,1),c(2,2,3,4),c(2,2,5,6))
grid.arrange(bars, p, aitch_plot, dist_comp, shan, core, layout_matrix=lay)
