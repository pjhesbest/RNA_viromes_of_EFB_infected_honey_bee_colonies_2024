#!/bin/bash

echo -e "\e[34mCOMMAND LINE JOB SUBMISSSION:\n\tbash mapping_and_counts.sh$
echo -e "\e[32m\nThis script performs the following:\n\n\t1. Filters out reads with a minimum (-m) and m$

Help()
{
echo -e "This script will perform the following:\n\t\t\t(1) Filter reads (ONT for the original purpose of this script) base on size (-m, --min/-M , --max)\n\t\t\t(2) Perform alignment of reads to reference using minimap2\n\t\t\t(3) Use samtools to convert sam to bam and output tsv files of the number of reads mapped per contig and the depths."   
echo -e "-i, --input\tfasta file of reads for mapping \e[31m[Required]\e[0m"
echo -e "-p, --prefix\tprefix for read file \e[31m[Optional]\e[0m"
echo -e "-r, --reference\treference file for mapping \e[31m[Required]\e[0m"
echo -e "-m, --min\tMinimum read length for mapping [default: -m 200]"
echo -e "-M, --max\tMaximum read length for mapping \e[31m[Optional]\e[0m\n"
}

###########################################################################################################
###########################################################################################################

while getopts i:p:r:m:M:h option
do
    case "${option}" in
        i)input=${OPTARG};;
        p)prefix=${OPTARG};;
        r)reference=${OPTARG};;
        m)min=${OPTARG};;
        M)max=${OPTARG};;
    h)Help; exit;;
    esac
done

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${prefix}" ]]; then
    prefix=$(basename ${input} .gz | sed 's/\.fasta//g; s/\.fastq//g; s/\.fq//g; s/\.fa//g; s/\.fq//g');$
if [[ -z "${reference}" ]]; then echo "-r, --reference REQUIRED"; Help; exit; fi
if [[ -z "${min}" ]]; then min=200; fi

###########################################################################################################
###########################################################################################################

# Perform the mapping

if [[ -z "${max}" ]]; then
    seqkit seq --min-len ${min} ${input} | minimap2 -ax splice -uf -t $SLURM_CPUS_PER_TASK --secondary=n$
    samtools index -@ $SLURM_CPUS_PER_TASK -b ${prefix}.bam
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
        samtools depth ${prefix}.bam > ${prefix}.depth

    totalReads=$(seqkit seq --min-len ${min} ${input} | seqkit stats -T - | awk '{print $4}' | sed '1d')

elif [[ ! -z "${max}" ]]; then
    seqkit seq --min-len $min --max-len $max ${input} | minimap2 -ax splice -uf -t $SLURM_CPUS_PER_TASK $
    samtools index -@ $SLURM_CPUS_PER_TASK -b ${prefix}.bam
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
        samtools depth ${prefix}.bam > ${prefix}.depth

    totalReads=$(seqkit seq --min-len ${min} --max-len ${max} ${input} | seqkit stats -T - | awk '{print$

fi

###########################################################################################################
###########################################################################################################

# Modify the idxstats file to contains the prefix and remove suffix

awk -v prefix="$prefix" '{print prefix "\t" $0}' ${prefix}.idxstats | sed 's/\.gz//g; s/\.fasta//g; s/\.$
awk '{print $1,$2,$3,$4}' ${prefix}.idxstats.tmp1 > ${prefix}.idxstats.tmp2

        # Filter out mapping with 0 reads mapped to contigs
                awk '{if ($4 > 0) print}' ${prefix}.idxstats.tmp3 > ${prefix}.idxstats.tmp4
                mv ${prefix}.idxstats.tmp4 ${prefix}.idxstats
awk -v totalReads="$totalReads" '{print $0 "\t" totalReads}' ${prefix}.idxstats.tmp2 > ${prefix}.idxstat$

# Modify the depth file to contains the prefix and remove suffix
awk -v prefix="$prefix" '{print prefix "\t" $0}' ${prefix}.depth | sed 's/\.gz//g; s/\.fasta//g; s/\.fas$
mv ${prefix}.depth.tmp1 ${prefix}.depth

awk '{if ($4 > 0) print}' ${prefix}.idxstats.tmp3 > ${prefix}.idxstats.tmp4
mv ${prefix}.idxstats.tmp4 ${prefix}.idxstats

###########################################################################################################
###########################################################################################################

# Clean up intermediate files
rm ${prefix}.idxstats.tmp*
rm ${prefix}.depth.tmp*

sed -i 's/ /\t/g' ${prefix}.idxstats
sed -i 's/ /\t/g' ${prefix}.depth