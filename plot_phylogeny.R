library(ape)
library(ggtree)
library(ggrepel)
library(googlesheets4)
library(dplyr)
library(viridis)
library(TDbook)
library(Biostrings)
library(patchwork)

#####################################################################################################
#####################################################################################################

set.seed(1234)

  # Import data
ctg.metadata <- read.csv("data/ApisM_MM_contigs-metadata.csv", header=TRUE)

  # Custom color palette
color.man <- list(Region = c("Asia","Europe","Africa","North America","NA",
                               "Middle East","South America","Oceania",   
                               "EFB+", "EFB-"),
              cols = c("#f46d43","#fdae61","#fee08b","#e6f598","#36454F",
                       "#abdda4","#66c2a5","#3288bd",   
                       "#800000ff", "#719400ff"))

#####################################################################################################
#####################################################################################################

  # DWV phylogeny
dwv.metadata <- read.table("data/phylogeny/ApisM_MM_DWV-contigs.list", quote="\"", comment.char="", col.names = c("chr"))
dwv.metadata <- merge(dwv.metadata, ctg.metadata, by = "chr", all.x = TRUE)

  # Read in tree
dwv.tree <- read.tree("Bin_4_DWV-contigs.phylogeny.contree")

  # Root tree
dwv.tree <- root(dwv.tree, "AF092924.1_Sacbrood_virus", resolve.root=TRUE, edgelabel=TRUE)

  # Plot tree
dwv.plot <-
  ggtree(dwv.tree, linewidth=0.9, 
         aes(color = Region)) %<+% dwv.metadata +
        scale_color_manual(name = 'Region', 
                     values = c("#f46d43","#800000ff", "#719400ff",
                                "#fdae61","#fee08b","#e6f598",
                                "#abdda4","#66c2a5","#3288bd","#36454F")
        ) +
        geom_hilight(node = 69, alpha = 0.2, fill = "#d53e4f", align="right") + # DWV-A
        geom_hilight(node = 54, alpha = 0.2, fill = "#3288bd", align="right") + # DWV-B
        geom_hilight(node = c(66,67,65), alpha = 0.2, fill = "#6a3d9a", align="right") + # RECOMB
        geom_tiplab(aes(label=Tip_short), size = 3, color = "black", offset=0.01) + 
        geom_treescale(x=0, y=-5, width=0.1, color='black') +
        geom_point(aes(color=EFB, shape=Foraging), na.rm=TRUE, size = 3) +
        geom_text2(aes(label=label,  # add the bootstrap values
                 subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
                 nudge_x = -0.01, nudge_y = 0.04, check_overlap = TRUE, size = 3)

dwv.plot

ggsave(dwv.plot, filename = "dwv.plot.pdf", height = 10, width = 8, dpi = 1200)

#####################################################################################################
#####################################################################################################

  # SBV phylogeny
sbv.metadata <- read.table("Bin_5_SBV-contigs.phylogeny.aln.trim.list", quote="\"", comment.char="",
                           col.names = c("chr"))
sbv.metadata <- merge(sbv.metadata, ctg.metadata, by = "chr", all.x = TRUE)
head(sbv.metadata)
  #Read in tree
sbv.tree <- read.tree("Bin_5_SBV-contigs.phylogeny.contree")
  #Root tree
sbv.tree <- root(sbv.tree, "AY292384.1_DWV_isolate_PA",
                 resolve.root=TRUE, edgelabel=TRUE)

sbv.plot <-
  ggtree(sbv.tree, linewidth=0.9, 
         aes(color = Region), ) %<+% sbv.metadata +
  scale_color_manual(name = 'Region', 
                     values = c("#f46d43","#800000ff", "#719400ff",
                                "#fdae61","#e6f598","#abdda4")
  ) +
  geom_tiplab(aes(label=Tip_short), size = 3, color = "black", offset=0.01) + 
  geom_treescale(x=0, y=-5, width=0.1, color='black') +
  geom_point(aes(color=EFB, shape=Foraging), na.rm=TRUE, size = 3) +
  geom_text2(aes(label=label,  # add the bootstrap values
                 subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
             nudge_x = -0.01, nudge_y = 0.04, check_overlap = TRUE, size = 3)

sbv.plot

ggsave(sbv.plot, filename = "sbv.plot.unedited.pdf", height = 10, width = 8, dpi = 1200)

#####################################################################################################
#####################################################################################################

  # BQCV phylogeny
bqcv.metadata <- read.table("Bin_7_BQCV-contigs.phylogeny.aln.trim.list", quote="\"", comment.char="",
                           col.names = c("chr"))
bqcv.metadata <- merge(bqcv.metadata, ctg.metadata, by = "chr", all.x = TRUE)
head(bqcv.metadata)
#Read in tree
bqcv.tree <- read.tree("Bin_7_BQCV-contigs.phylogeny.contree")

bqcv.plot <-
  ggtree(bqcv.tree, linewidth=0.9,
         aes(color = Region), ) %<+% bqcv.metadata +
  scale_color_manual(name = 'Region', 
                     values = c("#d53e4f","#f46d43","#800000ff", "#719400ff",
                                "#fdae61","#fee08b","#e6f598",
                                "#abdda4","#66c2a5","#36454F")
  ) +
  geom_tiplab(aes(label=Tip_short), size = 3, color = "black", offset=0.005) + 
  geom_treescale(x=0, y=-5, width=0.1, color='black') +
  geom_point(aes(color=EFB, shape=Foraging), na.rm=TRUE, size = 3) +
  geom_text2(aes(label=label,  # add the bootstrap values
                 subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
             nudge_x = -0.01, nudge_y = 0.04, check_overlap = TRUE, size = 3)
bqcv.plot

msaplot(bqcv.plot, "Bin_7_BQCV-contigs.phylogeny.aln.trim.fasta", offset=0.5, width =  0.5)

ggsave(bqcv.plot, filename = "bqcv.plot.unedited.pdf", height = 10, width = 8, dpi = 1200)


#####################################################################################################
#####################################################################################################

  # Lake Siani Virus
lsv.metadata <- read.table("Bin_2_LSV-contigs.phylogeny.aln.trim.list", quote="\"", comment.char="",
                            col.names = c("chr"))
lsv.metadata <- merge(lsv.metadata, ctg.metadata, by = "chr", all.x = TRUE)
head(lsv.metadata)
#Read in tree
lsv.tree <- read.tree("Bin_2_LSV-contigs.phylogeny.contree")

lsv.plot <- ggtree(lsv.tree, linewidth=0.9,
                   aes(color = Region),
                   ) %<+% lsv.metadata + 
  scale_color_manual(name = 'Region', 
                     values = c("#d53e4f","#f46d43","#800000ff", "#719400ff",
                                "#fdae61","#fee08b","#e6f598",
                                "#abdda4","#66c2a5","#36454F")
  ) +
  geom_hilight(node = 104, alpha = 0.8, fill = "darkgrey", align="right") + # LSV8
  geom_hilight(node = 121, alpha = 0.8, fill = "lightgrey", align="right") + # LSV2
  geom_hilight(node = 93, alpha = 0.8, fill = "darkgrey", align="right") + # LSV4
  geom_hilight(node = 132, alpha = 0.8, fill = "lightgrey", align="right") + # LSV3
  geom_hilight(node = 76, alpha = 0.8, fill = "darkgrey", align="right") + # LSV3

  geom_tiplab(aes(label=Tip_short), 
    size = 3, color = "black", offset=0.005) + 
  geom_treescale(x=0, y=-5, width=0.1, color='black') +
  geom_point(aes(color=EFB, shape=Foraging), na.rm=TRUE, size = 3) +
  geom_text2(aes(label=label,  # add the bootstrap values
                 subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
             nudge_x = -0.01, nudge_y = 0.04, check_overlap = TRUE, size = 3) #+
  #geom_text(aes(label=node), hjust=-.3) # add the node numbers of find highl + ggtree(data aes(color = Region)) %<+% lsv.metadata
lsv.plot

msaplot(lsv.plot, "Bin_2_LSV-contigs.phylogeny.aln.trim.fasta", offset=0.5, width =  0.5)

ggsave(lsv.plot, filename = "lsv.plot.unedited.pdf", height = 10, width = 8, dpi = 1200)

#####################################################################################################
#####################################################################################################

# AKI complex viruses

iapv.metadata <- read.table("Bin_6_IAPV-contigs.phylogeny.aln.trim.list", quote="\"", comment.char="",
                            col.names = c("chr"))
iapv.metadata <- merge(iapv.metadata, ctg.metadata, by = "chr", all.x = TRUE)
head(iapv.metadata)
#Read in tree
iapv.tree <- read.tree("Bin_6_IAPV-contigs.phylogeny.contree")

iapv.plot <-
  ggtree(iapv.tree, linewidth=0.9,
         aes(color = Region), ) %<+% iapv.metadata +
  scale_color_manual(name = 'Region', 
                     values = c("#d53e4f","#f46d43","#800000ff", "#719400ff",
                                "#fdae61","#fee08b","#e6f598",
                                "#abdda4","#66c2a5","#36454F")
  ) +
  geom_tiplab(aes(label=Tip_short), size = 3, color = "black", offset=0.005) + 
  geom_treescale(x=0, y=-5, width=0.1, color='black') +
  geom_point(aes(color=EFB, shape=Foraging), na.rm=TRUE, size = 3) +
  geom_text2(aes(label=label,  # add the bootstrap values
                 subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
             nudge_x = -0.01, nudge_y = 0.04, check_overlap = TRUE, size = 3)
iapv.plot


ggsave(bqcv.plot, filename = "bqcv.plot.unedited.pdf", height = 10, width = 8, dpi = 1200)
