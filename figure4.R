#!/usr/bin/env Rscript

library(ggplot2)
library(tidyverse)
library(googlesheets4)
library(ggpmisc)
library(smplot2)
library(plotly)
library(patchwork)
library(metagenomeSeq)
library(ggh4x)

theme_set(theme_minimal(base_size = 12, base_family = 'Source Sans Pro'))

set.seed(1234)

setwd("/home/dcschroe/heske011/projects/BLUEBERRY/POPULATION_ABUNDANCE")

############################################################################################################################
############################################################################################################################
############################################################################################################################

metadata.df <- read.csv("data/")
idxstats.df <- read.csv("data/", 
                            header=FALSE, sep="\t", col.names=c("sequencingBC","datatype", "chr", "locus", "numReads", "totalReads"))
merged.df <- merge(idxstats.df, metadata.df, by = "sequencingBC")
taxonomy.df <- read.delim("data/", 
                       header = F, sep = "\t", 
                       col.names=c("classification", "chr", "taxID",
                                   "lineage","V5","V6","V7","V8","V9","V10",
                                   "V11", "V12", "V13"))

drop <- c("V5","V6","V7","V8","V9","V10","V11", "V12", "V13", "classification")
taxonomy.df = taxonomy.df[,!(names(taxonomy.df) %in% drop)]
taxonomy.df <- separate_wider_delim(taxonomy.df, cols = lineage, delim = ";", 
                                    names = c("superkingdom","phylum","order","class","family","genus","species","V8"),
                                    too_many = "merge", too_few = "align_start")
taxonomy.df <- subset(taxonomy.df, select = -c(V8) )
taxonomy.df$species <- trimws(taxonomy.df$species, which = "left")
taxonomy.df$genus <- trimws(taxonomy.df$genus, which = "left")
taxonomy.df$family <- trimws(taxonomy.df$family, which = "left")

merged.df <- merge(merged.df, taxonomy.df, by = "chr")


head(merged.df)
colnames(merged.df)

contig.meta.df <- read_sheet("https://docs.google.com/spreadsheets/d/15VMY_-fpY6yf4PCuAu02WDB94lU6K6stdS30o4IvMjI/edit?usp=sharing")
head(contig.meta.df)

merged.df <- merge(merged.df, contig.meta.df, by = "chr")

############################################################################################################################
############################################################################################################################
############################################################################################################################

# Plot cumSum normalised data heatmap

counts.mat <- merged.df %>% select(chr, numReads, sequencingBC) %>% subset(chr != "*") %>%
  pivot_wider(names_from = sequencingBC, values_from = numReads) %>%
  ungroup() %>% remove_rownames %>% column_to_rownames(var="chr") %>% as.matrix()
#MetagenomeSeq cumalitive normalisation
Coverage.MR <- newMRexperiment(counts.mat)
p <- cumNormStat(Coverage.MR, pFlag=TRUE)
Coverage.MR <- cumNorm(Coverage.MR, p=p)
normFactors(Coverage.MR)
Coverage.norm <- MRcounts(Coverage.MR, norm=T, log=F)
library(reshape2)
cumNorm.df <- as.data.frame(as.table(Coverage.norm))
colnames(cumNorm.df) <- c("chr", "sequencingBC", "cumNorm")

merged.norm.df <- merge(merged.df, cumNorm.df,
                                       by.x = c("chr", "sequencingBC"), 
                                       by.y = c("chr", "sequencingBC"), all.x = TRUE)

write.csv(merged.norm.df, "merged.norm.df.csv")
colour.man <- list(
                   Organism = c("Deformed Wing Virus A","Deformed Wing Virus B","Deformed Wing Virus",
                                "Lake Sinai Virus","Israel Acute Paralysis Virus","Black Queen Cell Virus",
                                "Sacbrood Virus"),
                   value = c("#d53e4f", "#3288bd", "#6a3d9a", "#66c2a5", "#abdda4", "#B2BEB5", "#fdae61"))
unique(merged.norm.df$Organism)


read.recruitment.summary <- merged.norm.df %>% group_by(sequencingBC, Variant) %>% summarise(meanCumNorm = mean(cumNorm))

############################################################################################################
############################################################################################################

p1 <- merged.norm.df %>% 
  subset(superkingdom != "Bacteria" & superkingdom != "" & superkingdom != "Eukaryota") %>%
  ggplot() + geom_bar(aes(x = sequencingBC, y = cumNorm, fill= Variant),
                      position="fill", stat="identity") +
  scale_fill_manual(values = c("#161565","#6a3d9a","#d53e4f", "#3288bd","#66c2a5","#fdae61",
                               "#fdae61","#fdae61","#fdae61","#fdae61","#f46d43")) +
  facet_nested(.~ Foraging + EFB_status, space = "free", scale = "free", switch = "y") +
  theme(strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"),
        strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside",
        axis.text.x=element_blank()) +
        labs(x = "", y = "")
p1

############################################################################################################
############################################################################################################

aves.merged <- merged.df %>% group_by(Foraging, EFB_status, Organism, Variant) %>%
              summarise(meanCumNorm = mean(cumNorm), 
                        sdCumNorm = sd(cumNorm)) %>% ungroup()

p2 <- ggplot(aves.merged, aes( x = Variant, y = meanCumNorm, fill = Variant )) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#161565","#6a3d9a","#d53e4f", "#3288bd","#66c2a5",
                               "#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#f46d43")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none") +
  scale_y_continuous(name="Cumilatively normalised reads", trans = "log10") +
  facet_nested(Foraging~ EFB_status, switch = "y")

p2


aves.merged2 <- merged.df %>% group_by(Foraging, EFB_status, Variant) %>%
  summarise(meanCumNorm = mean(cumNorm),
            medianCumNorm = median(cumNorm),
            sdCumNorm = sd(cumNorm)) %>% ungroup()

rank_abundance <- function(foraging, efb){
  ggplot(subset(aves.merged2,  
                Foraging == foraging &  
                  EFB_status == efb),
       aes( x = fct_reorder(Variant, meanCumNorm),
            y = meanCumNorm, color = Variant )) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#161565","#6a3d9a","#d53e4f", "#3288bd","#66c2a5",
                               "#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#f46d43")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          ) +
  scale_y_continuous(name="", 
                     trans = "log10", limits = c(1,1.2e5))

}

p2a <- rank_abundance(foraging = "Blueberry", efb = "Negative")
p2b <- rank_abundance(foraging = "Blueberry", efb = "Positive")
p2c <- rank_abundance(foraging = "Holding", efb = "Negative")
p2d <- rank_abundance(foraging = "Holding", efb = "Positive")

p2 <- (p2a + p2b) / (p2c + p2d)
p2

ggplot(subset(aves.merged2,  
                Foraging == "Blueberry" &  
                  EFB_status == "Negative"),
       aes( x = fct_reorder(short_name, meanCumNorm),
            y = meanCumNorm, color = Organism )) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#161565","#6a3d9a","#d53e4f", "#3288bd","#66c2a5",
                                "#fdae61","#f46d43")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(name="Mean mapped reads\n(normalised)", trans = "log10",
                     limits = c(1,1e6))


p2

############################################################################################################
############################################################################################################

colnames(merged.df)
merged.mat <- merged.df %>% 
  select(sequencingBC, chr, cumNorm) %>% 
  pivot_wider(names_from = chr, values_from = cumNorm) %>%
  remove_rownames %>% column_to_rownames(var="sequencingBC") %>%
  ungroup() %>% as.matrix()

community.mds <- vegan::metaMDS(comm = merged.mat, 
                                distance = "bray", 
                                trace = FALSE, autotransform = FALSE)
plot(community.mds$points)
MDS_xy <- data.frame(community.mds$points)
MDS_xy.df <- merge(MDS_xy, metadata.df, by.x = 0, by.y = "sequencingBC", all.x = TRUE)  # Merge by row names

p3 <- ggplot(MDS_xy.df, aes(x = MDS1, y = MDS2, color = EFB_status, shape = Foraging)) +
  geom_point(aes(size = bctrimmedreads)) + 
  theme_bw() + #theme(legend.position = "none") +
  scale_color_manual(values = c("#719400ff","#800000ff")) # Customize the plot appearance as desired

community.mds$stress

metadata.df2 <- metadata.df %>% remove_rownames %>% column_to_rownames(var="sequencingBC")

adonis(merged.mat ~ EFB_status*Foraging, data=metadata.df2, permutations=99)


p3

############################################################################################################
############################################################################################################

bottom <- (p2 + p3)
p1 / bottom
