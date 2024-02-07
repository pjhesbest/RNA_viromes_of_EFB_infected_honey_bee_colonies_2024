#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=150GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

eval "$(conda shell.bash hook)"

conda activate phylogeny


echo -e "\e[34mCOMMAND LINE JOB SUBMISSSION:\n\tbash iq-tree_phylogeny.sh$

Help()  
{
echo -e "-i, --input\t/full/path/to/file.fasta for MSA and phylogeny \e[31m[Required]\e[0m"
echo -e "-b, --bootstraps\tNumber of boottraps for tree \e[31m[Default: -b 200 ]\e[0m"
echo -e "-R, --Retree\tRerun only the tree using existing alignment \e[31m[Optional]\e[0m"
}


while getopts i:b:R:h option
do
    case "${option}" in
        i)input=${OPTARG};;
                b)bootstraps=${OPTARG};;
                R)retree=true;;
    h)Help; exit;;
    esac

done

    	# check flags are used or set defaults
if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${bootstraps}" ]]; then bootstraps="200"; fi

        # create prefix variable for naming outputs
prefix=$(basename ${input} .gz | sed 's/\.fasta//g; s/\.fa//g; s/\.fq//g; s/\.fastq//g; s/\.fna//g')

# Check if the input file has a .gz extension and unzip it if necessary
if [[ "${input}" == *.gz ]]; then
  gunzip -c "${input}" >"${prefix}.fasta"
  input="${prefix}.fasta"
fi

# Check if the -R flag is used to determine whether to run mafft and iqtree
if [[ -n "${retree}" ]]; then

  sed -i 's/>_R_/>/g' ${prefix}.aln.trim.fasta
  iqtree -s "${prefix}.aln.fasta" -m GTR -T AUTO -ntmax "${SLURM_CPUS_PER_TASK}" -b "${bootstraps}" --pr$

else

  mafft --thread "${SLURM_CPUS_PER_TASK}" --adjustdirectionaccurately --auto "${input}" > "${prefix}.aln$
  trimal -in "${prefix}.aln.fasta" -out "${prefix}.aln.trim.fasta" -gt 0.9 -cons 60
  sed -i 's/>_R_/>/g' ${prefix}.aln.trim.fasta
  iqtree -s "${prefix}.aln.trim.fasta" -m GTR -T AUTO -ntmax "${SLURM_CPUS_PER_TASK}" -b "${bootstraps}"$
fi
