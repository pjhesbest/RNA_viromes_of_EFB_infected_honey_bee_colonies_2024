#!/bin/bash

echo -e "\e[34mCOMMAND LINE JOB SUBMISSSION to SLURM:\n\tbash viral-population-estimation_p1.sh$

echo -e "\e[32m\nThis script performs the following:\n\n\t1. Filters out reads with a minimum (-m; defau"
echo -e "Arguments:"

Help()
{
echo -e "-i, --input\tFull path to directory containing assembly output (fasta/fastq/.gz) file of reads"
echo -e "-r, --reference\treference file of viral genomes \e[31m[Required]\e[0m"
echo -e "-m, --min\tMinimum read length for mapping [default: -m 200]"
echo -e "-M, --max\tMaximum read length for mapping \e[31m[Optional]\e[0m\n"
}

while getopts i:p:r:m:M:h option
do
    case "${option}" in
                i)input=${OPTARG};;
        r)reference=${OPTARG};;
        m)min=${OPTARG};;
        M)max=${OPTARG};;
    h)Help; exit;;
    esac
done

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${reference}" ]]; then echo "-r, --reference REQUIRED"; Help; exit; fi
if [[ -z "${min}" ]]; then min=200; fi

######################
