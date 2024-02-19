# EFB_and_Foraging_honeybee_colonies-RNA_viromes_2024
Script and data files relating to the publication:
### Sacbrood viruses and select Lake Sinai virus variants dominated <i>Apis mellifera</i> colonies symptomatic for European foulbrood (<i>in review</i> J. Vir)
#### Poppy J. Hesketh-Best<sup>1</sup>, Peter D. Fowler<sup>2</sup>, Nkechi M. Odogwu<sup>1</sup>, Meghan O.G. Milbrath<sup>2</sup>, Declan C. Schroeder<sup>1</sup>
##### <sup>1</sup>Department of Veterinary Population Medicine, University of Minnesota, St. Paul, MN 55108, USA
##### <sup>2</sup>Department of Entomology, Michigan State University, Pollinator Performance Center, 4090 N. College Road, RM 100, Lansing, MI 48910, USA

---------------

Within this repository, you will find scripts to recreate the analyses and figures described in the manuscript. Additionally are the alignment files and tree files, as well as whole-genomes of the honey bee RNA viruses.

I have attempted to make the bash files general purpose (as long as the dependables are available on your system), however, <code>viral-population-estimations_p1.sh</code> requires specific Kaiju databases to be available and the paths modified in the script. Please refer to this repository for instructions on how to build the required databases and conda environments to run these scripts: https://github.com/dmckeow/bioinf

### Data provided:
- <code>data/ApisM_MM_contigs.metadata.csv</code> : Metadata relating to viral metagenome-assembled genomes (vMAG), taxonomy, genome quality (checkV), genome origin, ect.
- <code>data/ApisM_MM_metadata.csv</code> : Study metadata and Sequence Read Archives (SRA) accession numbers for accessing metagenomic reads.
- <code>data/ApisM_MM_phylogeny.metadata.csv</code> : Metadata relating to genome (publically available and vMAGS generated in this study) used in phylogeny.
- <code>data/ApisM_MM_read_kaiju_counts.csv</code> : Count tables of taxonomic classification of metagenomic reads by Kaiju (NCBI-nr and RVDB)
- <code>data/ApisM_MM_viral-genome-read-recruitment.csv</code> : Count tables of read recruitment to study vMAGs.
- <code>data/ApisM_MM_viral-genome.fasta</code> : Whole-genomes of vMAGs generated in this study, consult <code>data/ApisM_MM_contigs.metadata.csv</code> for relevant metadata.
- <code>data/ApisM_MM_viral-genomes-NCBI.fa</code> : NCBI format of whole-genomes of vMAGs generated in this study, consult <code>data/ApisM_MM_contigs.metadata.csv</code> for relevant metadata.
