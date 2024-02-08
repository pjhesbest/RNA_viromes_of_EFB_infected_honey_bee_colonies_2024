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

############################################################################################################################
############################################################################################################################

metadata.df <- read.csv("data/ApisM_MM_metadata.csv")
idxstats.df <- read.csv("data/ApisM_MM_viral-genome-read-recruitment.csv"))

############################################################################################################################
############################################################################################################################

p1 <- metadata.df %>% 
  ggplot() + geom_bar(aes(x = libraryID, y = cumulatively.normalise.reads, fill= Variant),
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

rank_abundance_plot <- function(foraging, efb){
              ggplot(subset(aves.merged2,  
                      Foraging == foraging &  
                      EFB_status == efb),
             aes( x = fct_reorder(Variant, meanCumNorm),
                  y = meanCumNorm, color = Variant )) +
            geom_point(size = 2) +
            scale_color_manual(values = c("#161565","#6a3d9a","#d53e4f", "#3288bd",
                                          "#66c2a5","#fdae61","#fdae61","#fdae61",
                                          "#fdae61","#fdae61","#f46d43")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            theme(legend.position = "none",
            axis.title.x=element_blank(),
            ) +
            scale_y_continuous(name="", 
                     trans = "log10", limits = c(1,1.2e5))
}

p2a <- rank_abundance_plot(foraging = "Blueberry", efb = "Negative")
p2b <- rank_abundance_plot(foraging = "Blueberry", efb = "Positive")
p2c <- rank_abundance_plot(foraging = "Holding", efb = "Negative")
p2d <- rank_abundance_plot(foraging = "Holding", efb = "Positive")

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

cumulatively.normalise.reads.mat <- idxstats.df %>% select(genomeID, cumulatively.normalise.reads, libraryID) %>%
                                    pivot_wider(names_from = libraryID, values_from = cumulatively.normalise.reads) %>%
                                    ungroup() %>% 
                                    remove_rownames %>% column_to_rownames(var="genomeID") %>% 
                                    as.matrix()

community.mds <- vegan::metaMDS(comm = cumulatively.normalise.reads.mat, 
                                distance = "bray", 
                                trace = FALSE, autotransform = FALSE)
plot(community.mds$points)
MDS_xy <- data.frame(community.mds$points)
MDS_xy.df <- merge(MDS_xy, metadata.df, by.x = 0, by.y = "sequencingBC", all.x = TRUE)  # Merge by row names

p3 <- ggplot(MDS_xy.df, 
             aes(x = MDS1, 
                 y = MDS2, 
                 color = EFB_status, 
                 shape = Foraging)) +
              geom_point(aes(size = bctrimmedreads)) + 
              theme_bw() + #theme(legend.position = "none")
              scale_color_manual(values = c("#719400ff","#800000ff")) # Customize the plot appearance as desired

community.mds$stress #stress value = 

metadata.df2 <- metadata.df %>% remove_rownames %>% column_to_rownames(var="sequencingBC")

adonis(merged.mat ~ EFB_status*Foraging, data=metadata.df2, permutations=99)


p3

############################################################################################################
############################################################################################################

bottom <- (p2 + p3)
p1 / bottom
