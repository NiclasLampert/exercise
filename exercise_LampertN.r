# exercise_LampertN.r
# author: Niclas Lampert (niclas.lampert@gmail.com)
#######################################################################################
# The objective of this script is to identify individuals for which the levels of the 
# malignant isoform (m) significantly changed from August to December independently of 
# the changes in the gene(t).
#######################################################################################

## thresholds for p value and fold change. CHANGE IF NECESSARY ####
p.threshold <- 0.05 # p value below which the null hypothesis will be rejected
fc.threshold <- 1.5 # common fold change threshold considered biologically meaningful

## load libraries ####
library(dplyr)
library(tidyr)
library(ggpubr)
library(ibb)

## import raw data ####
counts <- as_tibble(read.csv("raw_data.tsv", sep = "\t"))
## The data is from a longitudinal study, with counts of malignant (m) and total (t) isoforms,
## four replicates each, and two time points (August, December).
glimpse(counts)
summary(counts)
## identify individuals
counts <- counts %>% mutate(Individual = 1:nrow(counts)) %>% select(Individual, August_1m:December_4t)
## Negative values will be disregarded
counts.pos <- counts %>% mutate_all(funs(replace(., .<0, NA)))
## infinite values will be disregarded
counts.good <- do.call(tibble, lapply(counts.pos, function(x) replace(x, is.infinite(x), NA)))
summary(counts.good)

## filter dataset for complete cases only
counts.compl <- counts[complete.cases(counts.good),]
summary(counts.compl)
glimpse(counts.compl)

## create matrices for malignant and total isoforms, respectively
mal.compl.mtr <- as.matrix(counts.compl %>% select(August_1m:December_4m))
tot.compl.mtr <- as.matrix(counts.compl %>% select(August_1t:December_4t))

## use the ibb package for beta binomial test:
group <- rep(c("August", "December"), each = 4) # define the time points to be compared
## Perform the inverted beta-binomial test for paired count data,
## since the counts 1-4 are paired replicates of total (t) and malignant (m) gene counts
compl.bb <- ibb.test(mal.compl.mtr, tot.compl.mtr, group, n.threads = 0)

## append the p-values to the counts table
counts.compl <- counts.compl %>% 
  mutate(p.values = as.vector(unlist(compl.bb$p.value)), fc = as.vector(unlist(compl.bb$fc)))
## select individuals where malignant isoforms underwent significant changes, 
## with p value and fold change thresholds as defined above.
significant.changes <- counts.compl %>% 
  filter(p.values < p.threshold & (fc >= fc.threshold | fc <= 1/fc.threshold))
write.table(significant.changes, "significant.changes.csv", sep = "\t", dec = ".", row.names = FALSE)
