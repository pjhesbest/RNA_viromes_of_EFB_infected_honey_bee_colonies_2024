#!/usr/bin/env Rscript

library(ggplot2)
library(tidyverse)
library(googlesheets4)
library(ggpmisc)
library(smplot2)
library(patchwork)
library(vegan)
library(ggh4x)

theme_set(theme_minimal())

set.seed(1234)

############################################################################################################################
############################################################################################################################

metadata.df <- read.csv("data/ApisM_MM_metadata.csv")
idxstats.df <- read.csv("data/ApisM_MM_viral-genome-read-recruitment.csv")

############################################################################################################################
############################################################################################################################

p1 <- idxstats.df %>% 
  ggplot() + geom_bar(aes(x = libraryID, 
                          y = cumulatively.normalise.reads, 
                          fill= Variant),
                      position="fill", stat="identity") +
  scale_fill_manual(values = c("#161565","#6a3d9a","#d53e4f", "#3288bd","#66c2a5","#fdae61",
                               "#fdae61","#fdae61","#fdae61","#fdae61","#f46d43")) +
  facet_nested(.~ Yard + EFB, space = "free", scale = "free", switch = "y") +
  theme(strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"),
        strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside",
        axis.text.x=element_blank()) +
        labs(x = "", y = "")
p1

############################################################################################################
############################################################################################################

aves.merged <- idxstats.df %>% group_by(Yard, EFB, CommonName, Variant) %>%
  summarise(meanCumNorm = mean(cumulatively.normalise.reads),
            medianCumNorm = median(cumulatively.normalise.reads),
            sdCumNorm = sd(cumulatively.normalise.reads)) %>% 
  ungroup()

rank_abundance_plot <- function(foraging, efb){
              ggplot(subset(aves.merged,  
                      Yard == foraging &  
                      EFB == efb),
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

p2a <- rank_abundance_plot(foraging = "Yard A", efb = "Negative")
p2b <- rank_abundance_plot(foraging = "Yard A", efb = "Positive")
p2c <- rank_abundance_plot(foraging = "Yard B", efb = "Negative")
p2d <- rank_abundance_plot(foraging = "Yard B", efb = "Positive")

p2 <- (p2a + p2b) / (p2c + p2d)
p2

############################################################################################################
############################################################################################################

cumulatively.normalise.reads.mat <- idxstats.df %>% select(genomeID, 
                                                           cumulatively.normalise.reads, 
                                                           libraryID) %>%
                                    pivot_wider(names_from = genomeID, values_from = cumulatively.normalise.reads) %>%
                                    ungroup() %>% 
                                    remove_rownames %>% column_to_rownames(var="libraryID") %>% 
                                    as.matrix()

community.mds <- vegan::metaMDS(comm = cumulatively.normalise.reads.mat, 
                                distance = "bray", 
                                trace = FALSE, autotransform = FALSE)
plot(community.mds$points)
MDS_xy <- data.frame(community.mds$points)
MDS_xy.df <- merge(MDS_xy, metadata.df, by.x = 0, by.y = "libraryID", all.x = TRUE)  # Merge by row names

p3 <- ggplot(MDS_xy.df, 
             aes(x = MDS1, 
                 y = MDS2, 
                 color = EFB, 
                 shape = Yard)) +
              geom_point(aes(size = LibrarySize)) + 
              scale_color_manual(values = c("#719400ff","#800000ff")) # Customize the plot appearance as desired
p3

community.mds$stress #stress value = 0.140

# ADONIS

adanois.df <- idxstats.df %>% select(Yard, EFB, CommonName, cumulatively.normalise.reads) %>%
  subset(cumulatively.normalise.reads > 0)

adonis.res <- adonis(cumulatively.normalise.reads ~ EFB * Yard * CommonName, 
                     data = adanois.df, 
                     permutations = 99, 
                     method = "bray")
adonis.res

############################################################################################################
############################################################################################################

bottom <- (p2 + p3)
p1 / bottom
