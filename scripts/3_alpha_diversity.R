########################################################################################
# alpha diversity

library(phyloseq)
library(microbiome)

raw_counts<-read.table("data/cutadapt_counts2.txt", sep='\t', header=TRUE, row.names=1, quote="")
tax<-read.table("data/cutadapt_tax2.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
tax<-as.matrix(tax) 

OTU = otu_table(raw_counts, taxa_are_rows = FALSE)
TAX = tax_table(tax)

physeq = phyloseq(OTU, TAX)
physeq

div<- estimate_richness(physeq, split = TRUE, measures = NULL)

#from microbiome R package, calculate the Berger-Parker dominance index
dom<- dominance(physeq)

#merge the data
div_all <- data.frame(div, dom)
rownames(div_all)<-gsub("X","", rownames(div_all))

div_merge<-merge(div_all, meta, by=0, all=FALSE)
div_merge2<-subset(div_merge, !is.na(timepoint))
div_merge2$participant<-as.factor(div_merge2$participant)
div_merge2<-div_merge2[!((div_merge2$Row.names) %in% c("003_post", "004_pre")),]

#write.table(div_merge2, file="data/alpha_diversity.txt", sep='\t', quote=F)
div_merge2<-read.table("data/alpha_diversity.txt", sep='\t', header=TRUE, quote="")

shan<-ggplot(div_merge2, aes(x=timepoint, y=Shannon, group = participant, colour = participant)) +
  geom_point(aes(fill=participant), size=4, shape=21, stroke=1, alpha=0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette="Paired") +
  scale_fill_brewer(palette="Paired") +
  labs(y="Shannon's Index") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position ="none")
shan


core<-ggplot(div_merge2, aes(x=timepoint, y=core_abundance, group = participant, colour = participant)) +
  geom_point(aes(fill=participant), size=4, shape=21, stroke=1, alpha=0.5) +
  geom_line(aes(group = participant), size = 1) +
  scale_color_brewer(palette="Paired") +
  scale_fill_brewer(palette="Paired") +
  labs(y="Core Abundance") +
  scale_x_discrete(labels = c("Pre", "Post")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position="none")
core

# calculate the difference post vs. pre

# Load necessary library
library(dplyr)

# Ensure the data is sorted by participant and timepoint
div_merge2 <- div_merge2 %>% 
  arrange(participant, timepoint)

# Extract numeric columns for difference calculation
numeric_cols <- sapply(div_merge2, is.numeric)
numeric_data <- div_merge2[, numeric_cols]

# Create a new dataframe with differences (post - pre)
diff_df <- div_merge2 %>% 
  mutate(participant_id = sub("_.*", "", Row.names)) %>% # Extract participant ID
  group_by(participant_id) %>% 
  summarise(across(where(is.numeric), 
                   ~ .[grepl("post", Row.names)] - .[grepl("pre", Row.names)], 
                   .names = "diff_{col}"), 
            metabolic_responder = first(metabolic_responder))

write.table(diff_df, file="data/pre_post_differences.txt", sep='\t', quote=F)

# perform a 2-way anova
# Create the 'participant_id' column to group by
div_merge2$participant_id <- sub("_.*", "", div_merge2$Row.names)
div_merge2$participant_id<-as.factor(div_merge2$participant_id)
aov2<-aov(Shannon ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov2)
# 
# Error: participant_id
# Df Sum Sq Mean Sq F value Pr(>F)
# metabolic_responder  1 0.0015 0.00149   0.005  0.945
# Residuals            8 2.3191 0.28989               
# 
# Error: participant_id:timepoint
# Df  Sum Sq Mean Sq F value Pr(>F)
# timepoint                      1 0.02981 0.02981   0.775  0.404
# timepoint:metabolic_responder  1 0.07666 0.07666   1.993  0.196
# Residuals                      8 0.30774 0.03847 

aov_obs<-aov(Observed ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_obs)
aov_chao<-aov(Chao1 ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_chao)
aov_ace<-aov(ACE ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_ace)
aov_simp<-aov(Simpson ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_simp)
aov_invsimp<-aov(InvSimpson ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_invsimp)
aov_fisher<-aov(Fisher ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_fisher)
aov_dbp<-aov(dbp ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_dbp)
aov_dmn<-aov(dmn ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_dmn)
aov_abs<-aov(absolute ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_abs)
aov_rel<-aov(relative ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_rel)
aov_simp2<-aov(simpson ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_simp2)
aov_core<-aov(core_abundance ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_core)
aov_gini<-aov(gini ~ timepoint * metabolic_responder + Error(participant_id/(timepoint * metabolic_responder)), data=div_merge2)
summary(aov_gini)

# pairwise t-tests

# Load necessary library
library(dplyr)

# Ensure data is sorted properly
div_merge2 <- div_merge2 %>% arrange(participant, timepoint)

# Select columns to test
columns_to_test <- colnames(div_merge2)[2:17]

# Initialize results dataframe
results <- data.frame(Variable = character(), t_statistic = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Perform pairwise t-tests
for (col in columns_to_test) {
  test_result <- t.test(div_merge2[[col]] ~ div_merge2$timepoint, paired = TRUE)
  results <- rbind(results, data.frame(Variable = col, t_statistic = test_result$statistic, p_value = test_result$p.value))
}

# Print results
print(results)

# Variable t_statistic   p_value p_adjusted
# t         Observed  0.80387016 0.4421860   0.986876
# t1           Chao1  0.79670624 0.4461286   0.986876
# t2        se.chao1 -0.38086323 0.7121391   0.986876
# t3             ACE  0.78996671 0.4498589   0.986876
# t4          se.ACE  0.61018987 0.5568275   0.986876
# t5         Shannon  0.83545689 0.4250832   0.986876
# t6         Simpson  0.44420026 0.6673864   0.986876
# t7      InvSimpson  0.01691174 0.9868760   0.986876
# t8          Fisher  0.93800121 0.3727319   0.986876
# t9             dbp -0.04685039 0.9636556   0.986876
# t10            dmn -0.22667915 0.8257381   0.986876
# t11       absolute  0.04546791 0.9647272   0.986876
# t12       relative -0.04685039 0.9636556   0.986876
# t13        simpson -0.44420026 0.6673864   0.986876
# t14 core_abundance -1.43804346 0.1842696   0.986876
# t15           gini -1.14823941 0.2804727   0.986876