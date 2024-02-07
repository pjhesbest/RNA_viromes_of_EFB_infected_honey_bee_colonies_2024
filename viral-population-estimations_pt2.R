#!/usr/bin/env Rscript

library(metagenomeSeq)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(data.table)
library(parallel)

theme_set(theme_minimal(base_size = 10, base_family = 'Source Sans Pro'))

set.seed(1234)

#########################
# Import text files
#########################
# Read idxstats file as TSV
idxstats.df <- read.delim("mapping.idxstats.tsv", 
                          header = F, sep = "\t",
                          col.names=c("sequencingBC", "datatype", "chr", "length", "counts", "totalReads"),) %>% subset(chr != "*")
# Read depth file as TSV - this file tends to be too big
depth.df <- read.delim("mapping.depth.tsv", 
                       header = F, sep = "\t", 
                       col.names=c("sequencingBC", "datatype", "chr", "locus", "counts", "totalReads"))

# Read taxonomy file as TSV
taxonomy.df <- read.delim("population.pool.kaiju.merge", 
                          header = F, sep = "\t", 
                          col.names=c("classification", "chr", "taxID",
						"lineage","V5","V6","V7","V8","V9","V10",
						"V11", "V12", "V13"))
drop <- c("V5","V6","V7","V8","V9","V10","V11", "V12", "V13")
taxonomy.df = taxonomy.df[,!(names(taxonomy.df) %in% drop)]
taxonomy.df <- separate_wider_delim(taxonomy.df, cols = lineage, delim = ";", 
                     names = c("superkingdom","phylum","order","class","family","genus","species","V8"),
                     too_many = "merge", too_few = "align_start")
taxonomy.df <- subset(taxonomy.df, select = -c(V8) )

# Ensure values are numeric:
idxstats.df$length <- as.numeric(idxstats.df$length)
idxstats.df$counts <- as.numeric(idxstats.df$counts)
idxstats.df$totalReads <- as.numeric(idxstats.df$totalReads)
depth.df$locus <- as.numeric(depth.df$locus)
depth.df$counts <- as.numeric(depth.df$counts)
depth.df$totalReads <- as.numeric(depth.df$totalReads)

#########################
# Calculate different abundance metrics
idxstats.df$normalised <- with(idxstats.df, ((counts)/(totalReads))/length)
idxstats.df$RPKM <- with(idxstats.df, (counts)/((length/1e3) * (totalReads/1e6)))
idxstats.df$Cov <- with(idxstats.df, ((counts * 200)/length))

counts.mat <- idxstats.df %>% select(chr, counts, sequencingBC) %>% subset(chr != "*") %>%
            pivot_wider(names_from = sequencingBC, values_from = counts) %>%
            ungroup() %>% remove_rownames %>% column_to_rownames(var="chr") %>% as.matrix()

normalised.mat <- idxstats.df %>% select(chr, normalised, sequencingBC) %>% 
            select(chr, normalised, sequencingBC) %>% subset(chr != "*") %>%
            pivot_wider(names_from = sequencingBC, values_from = normalised) %>%
            ungroup() %>% remove_rownames %>% column_to_rownames(var="chr") %>% as.matrix()

RPKM.mat <- idxstats.df %>% select(chr, RPKM, sequencingBC) %>% subset(chr != "*") %>%
            pivot_wider(names_from = sequencingBC, values_from = RPKM) %>%
            ungroup() %>% remove_rownames %>% column_to_rownames(var="chr") %>% as.matrix()
# Calculate average coverage per bp
depth.df$length <- idxstats.df$length[match (depth.df$chr, idxstats.df$chr)]
depth.df$Cov <- idxstats.df$Cov[match (depth.df$chr, idxstats.df$chr)]
depth.df$Cov_per_bp <- with(depth.df, (counts * 200)/(length))
head(depth.df)

depth.min.df <- depth.df %>%
  group_by(sequencingBC, chr, datatype) %>%
  summarise(meanCov_per_bp = mean(Cov_per_bp),
            minCov_per_bp = min(Cov_per_bp), 
            maxCov_per_bp = max(Cov_per_bp),
            Cov = mean(Cov),)

#MetagenomeSeq cumalitive normalisation
Coverage.MR <- newMRexperiment(counts.mat)
p <- cumNormStat(Coverage.MR, pFlag=TRUE)
Coverage.MR <- cumNorm(Coverage.MR, p=p)
normFactors(Coverage.MR)
Coverage.norm <- MRcounts(Coverage.MR, norm=T, log=F)
library(reshape2)
cumNorm.df <- as.data.frame(as.table(Coverage.norm))
colnames(cumNorm.df) <- c("chr", "sequencingBC", "cumNorm")

# Add taxonomy to mapping tables and filter out idxstats to remove annoying * (i.e no. of unmapped reads)
idxstats.df <- idxstats.df %>% subset(chr != "*")
idxstats.df <- merge(idxstats.df, taxonomy.df, by= c("chr", "chr"))
depth.df <- merge(depth.df, taxonomy.df, by= c("chr", "chr"))
cumNorm.df <- merge(cumNorm.df, taxonomy.df, by= c("chr", "chr"))

#export normalized tables
write.csv(Coverage.norm, paste("tables/MGSeq_Read_counts.matrix.csv",sep=""),quote=F)
write.csv(counts.mat, paste("tables/Read_counts.matrix.csv",sep=""),quote=F)
write.csv(normalised.mat, paste("tables/Normalised_counts.matrix.csv",sep=""),quote=F)
write.csv(RPKM.mat, paste("tables/RPKM_counts.matrix.csv",sep=""),quote=F)
write.csv(cumNorm.df, paste("tables/mapping.cumNorm.MGSeq.csv", sep=""),quote=F)
write.csv(depth.min.df, paste("tables/mapping.depth.summary.csv", sep=""),quote=F)

write.csv(depth.df, paste("tables/mapping.depth.merged.csv", sep=""), quote=F)
write.csv(idxstats.df, paste("tables/mapping.idxstats.merged.csv", sep=""), quote=F)

#######################################
con <- read.csv("tmp.ALLcontig", header=TRUE, sep="\t")
sam_con <- read.csv("tmp.ALLstats.list.count0", header=TRUE, sep="\t")
head(bas_1)
bas_1 <- read.csv("tmp.ALLstats.list.count1", header=TRUE, sep="\t")
bas_2 <- read.csv("tmp.ALLstats.list.count2", header=TRUE, sep="\t")
bas_3 <- read.csv("tmp.ALLstats.list.count3", header=TRUE, sep="\t")

df_merged <- merge(con, sam_con, by.x="contig", by.y="contig")
df_merged <- merge(bas_1, df_merged, by.x=c("reads","contig"), by.y=c("reads","contig"), all.y=TRUE)
df_merged <- merge(bas_2, df_merged, by.x=c("reads","contig"), by.y=c("reads","contig"), all.y=TRUE)
df_merged <- merge(bas_3, df_merged, by.x=c("reads","contig"), by.y=c("reads","contig"), all.y=TRUE)

## add average coverage depth per base for each contig (total coveage per base for contig divided by contig length), then times 100 to get the percent of the contig length in bases that passed the coverage threshold in each category - e.g. 50 % of contig had at least 5X coverage
df_merged <- transform(df_merged, percent_of_contig_len_with_over_5X_coverage = total_bases_over_5_coverage / contig_length * 100)
df_merged <- transform(df_merged, percent_of_contig_len_with_over_10X_coverage = total_bases_over_10_coverage / contig_length * 100)
df_merged <- transform(df_merged, percent_of_contig_len_with_over_100X_coverage = total_bases_over_100_coverage / contig_length * 100)

write.csv(df_merged, paste("tables/mapping.depth.coverages.csv", sep=""), quote=F)
