#### Load data ####

counts <- read.table("data/cutadapt_counts2.txt", 
                     header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, 
                     quote = "", stringsAsFactors = FALSE)
tax <- read.table("data/cutadapt_tax2.txt", 
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, 
                  quote = "", stringsAsFactors = FALSE)

track <- read.table("data/cutadapt_track2.txt", header = TRUE, row.names = 1, sep = "\t", quote = "")                  

# Define the samples to keep
rows_to_keep <- c(
  "001_post", "001_pre", "002_post", "002_pre", "003_post", "004_pre", 
  "005_post", "005_pre", "006_post", "006_pre", "007_post", "007_pre", 
  "008_post", "008_pre", "009_post", "009_pre", "010_post", "010_pre", 
  "011_post", "011_pre", "012_post", "012_pre", "DNA_BLANK2", "POS", "PCR_BLANK2"
)

# Subset 'track' to only contain the samples of interest
track1<- track[rownames(track) %in% rows_to_keep, ]

print(track1)

#Looking a the FEMAP samples only (not including the technical controls)
femap_reads <- track1$Tax.Filtered[1:22]

mean(femap_reads)
median(femap_reads)
min(femap_reads)
max(femap_reads)

# Take a look at the other two files
counts[1:5,1:10]

tax[1:5,]

# what are the dimensions of the original counts file?
dim(counts)

# Subset the counts data frame based on row names
counts_filt <- counts[rownames(counts) %in% rows_to_keep, ]

# What are the filtered file dimensions?
dim(counts_filt)

# Combine the counts and tax tables into one tidier working file
t.counts<-t(counts_filt)
ct<-merge(t.counts, tax, by="row.names")

rownames(ct)<-ct$Row.names
ct$Row.names <- NULL
ct$Sequence <- NULL

library(tidyr)
ct<-unite(ct,"tax.vector", Kingdom:Species, sep=':', remove=TRUE)

dim(ct)
#[1] 2418   25
#write.table(ct, file='data/counts_unfilt.txt', sep='\t', quote=F)

#repeat with taxonomy Kingdom:Genus
gen<- merge(t.counts, tax, by="row.names")
gen<-unite(gen, "tax.vector", Kingdom:Genus, sep=':', remove=TRUE)
gen[,27:28]<-list(NULL)
rownames(gen)<-gen$Row.names
gen$Row.names<-NULL
#write.table(gen, file="data/counts_001_gen.txt", sep='\t', quote=F)

### Filtering and pruning table

#d <- read.table("data/counts_unfilt.txt", sep="\t", quote="", header=T, row.names=1)

#dim(d)
# 2418   28

tax <- ct$tax.vector
#get only the count columns
dm <- ct[,1:ncol(ct)-1]
dim(dm)
#2418	24

i <- (colSums(dm) <=1000)
d.s <- dm[, !i]
dim(d.s)
#2418	24
ncol(dm)-ncol(d.s)
#0
#all samples (columns) have colSums > 1000

d.freq <- apply(d.s, 2, function(x) {x/sum(x)})

#keep SVs > 0.1% in any sample
d.0 <- d.s[apply(d.freq, 1, max)>0.001,]
dim(d.0)
# 784  27

# make relative abundance table
d.freq0 <- data.frame(apply(d.0, 2, function(x){x/sum(x)}))

#add taxonomy back in and save the filtered counts file
d.0$tax.vector = d$tax.vector[match(rownames(d.0), rownames(d))]
#write.table(d.0, file="data/dada2_tax_counts_001filt.txt", sep='\t', quote=F)

# make relative abundance table
d.freq0 <- data.frame(apply(d.0, 2, function(x){x/sum(x)}))

#add taxonomy back in and save the filtered counts file
d.freq0$tax.vector = d$tax.vector[match(rownames(d.freq0), rownames(d))]
#write.table(d.freq0, file="data/dada2_tax_freq_001filt.txt", sep='\t', quote=F)

#filter SVs based on a read count cutoff
count = 100
d.2 <- data.frame(d.s[which(apply(d.s, 1, function(x){sum(x)}) > count),])
dim(d.2)
# 532  24

#combine d.freq0 and d.2 so that we have a filtered table that contains only SVs present a 0.1% abundance AND >100 reads across all samples
row_list<-intersect(rownames(d.0), rownames(d.2))
length(row_list)
#526

d_filt<-d.s[rownames(d.s) %in% row_list,]
#add taxonomy back in and save the filtered counts file
d_filt$tax.vector = gen$tax.vector[match(rownames(d_filt), rownames(gen))]
#write the file
write.table(d_filt, file="data/counts_01abun_100filt_2.txt", sep='\t', quote=F)

d_filt_freq<-d.s[rownames(d.s) %in% row_list,]
d_filt_freq<-data.frame(apply(d_filt_freq, 2, function(x){x/sum(x)}))
d_filt_freq$tax.vector = gen$tax.vector[match(rownames(d_filt_freq), rownames(gen))]
write.table(d_filt_freq, file="data/abundance_01abun_100filt.txt", sep='\t', quote=F)
