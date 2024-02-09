#!/bin/bash

# Requirements:
#	mafft https://mafft.cbrc.jp/alignment/software/
#	cd-hit http://cd-hit.org/
#	kaiju https://bioinformatics-centre.github.io/kaiju/
#	seqkit https://bioinf.shenwei.me/seqkit/
#	minimap2 https://github.com/lh3/minimap2
#	samtools https://www.htslib.org/

# To succesfully run this script you will be required to install the necessary databases. 
# Follow instructions from the following repository to do so: https://github.com/dmckeow/bioinf

##################################################################################################
##################################################################################################

echo -e "\e[34mCOMMAND LINE JOB SUBMISSSION to SLURM:\n\tsbatch /home/dcschroe/heske011/codebase/population_estimations.sh $@\e[0m"
echo -e "\e[32m\nThis script performs the following:\n\n\t1. Filters out reads with a minimum (-m; default 200 bp) and maximum (-M) length\n\t2. Reference genome are clustered at ≥95% ANI across ≥80% of their lengths - viral genome 'population' pool.\n\t3. Filtered reads are mapped to viral genome 'population' pool.\n\t4. Secondary mapped reads are discarded, retaining only primary.\n\t4. Read counts are extracted from mapping file (BAM) three files generated - number of reads mapped ({prefix}.idxstats), sequencing depth ({prefix}.depth) and coverage ({prefix}.coverage)\n\e[0m\n"
echo -e "Arguments:"

Help()
{
echo -e "-i, --input\tFull path to directory containing assembly output (fasta/fastq/.gz) file of reads for mapping \e[31m[Required]\e[0m"
echo -e "-r, --reference\treference file of viral genomes \e[31m[Required]\e[0m"
echo -e "-m, --min\tMinimum read length for mapping [default: -m 200]"
echo -e "-M, --max\tMaximum read length for mapping \e[31m[Optional]\e[0m\n"
echo -e "-t, --threads\tNumber of threads \e[31m[Default: -t 4]\e[0m\n"
}

##################################################################################
##################################################################################

while getopts i:p:r:m:M:t:h option
do 
    case "${option}" in
		i)input=${OPTARG};;
        r)reference=${OPTARG};;
        m)min=${OPTARG};;
        M)max=${OPTARG};;
	t)threads=${OPTARG};;
    h)Help; exit;;
    esac
done

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${reference}" ]]; then echo "-r, --reference REQUIRED"; Help; exit; fi
if [[ -z "${min}" ]]; then min=200; fi
if [[ -z "${threads}" ]]; then threads=4; fi

##################################################################################
##################################################################################

# Generate viral population pool clustering at 95% identity and >80% length
cd-hit-est -i ${reference} -c 0.95 -s 0.8 -n 4 -o population.pool.fasta

##################################################################################
##################################################################################

for reads in $input/*/*.{fastq.gz,fasta.gz,fq.gz,fa.gz,fastq,fasta,fq,fa}; do
	
prefix=$(basename ${reads} .gz | sed 's/\.fasta\?//g; s/\.fastq\?//g; s/\.fa\?//g; s/\.fq\?//g')

# Perform the mapping
if [[ -z "${max}" ]]; then
    seqkit seq --min-len ${min} ${reads} | minimap2 -ax map-ont -t ${threads} --secondary=no population.pool.fasta - | samtools sort -O BAM - > ${prefix}.bam
    samtools index -@ ${threads} -b ${prefix}.bam
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
	samtools depth ${prefix}.bam > ${prefix}.depth
	# store a variable of the total number of reads after filtering
    totalReads=$(seqkit seq --min-len ${min} ${reads} | seqkit stats -T - | awk '{print $4}' | sed '1d')

elif [[ ! -z "${max}" ]]; then
    seqkit seq --min-len $min --max-len $max ${reads} | minimap2 -ax map-ont -t ${threads} --secondary=no population.pool.fasta - | samtools sort -O BAM - > ${prefix}.bam
    samtools index -@ ${threads} -b ${prefix}.bam
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
	samtools depth ${prefix}.bam > ${prefix}.depth
	# store a variable of the total number of reads after filtering
    totalReads=$(seqkit seq --min-len ${min} --max-len ${max} ${reads} | seqkit stats -T - | awk '{print $4}' | sed '1d')
	
fi

    ## Modify the idxstats file to contains the prefix and remove suffix
awk -v prefix="$prefix" '{print prefix "\t" $0}' ${prefix}.idxstats > ${prefix}.idxstats.tmp1
awk '{print $1,$2,$3,$4}' ${prefix}.idxstats.tmp1 > ${prefix}.idxstats.tmp2
	# Filter out mapping with 0 reads mapped to contigs
		awk '{if ($3 != "*") print}' ${prefix}.idxstats.tmp2 > ${prefix}.idxstats.tmp3
		mv ${prefix}.idxstats.tmp3 ${prefix}.idxstats
awk -v totalReads="$totalReads" '{print $0 "\t" totalReads}' ${prefix}.idxstats.tmp2 > ${prefix}.idxstats

    ## Modify the depth file to contains the prefix and remove suffix
awk -v prefix="$prefix" '{print prefix "\t" $0}' ${prefix}.depth > ${prefix}.depth.tmp1
awk -v totalReads="$totalReads" '{print $0 "\t" totalReads}' ${prefix}.depth.tmp1 > ${prefix}.depth

    # Clean up intermediate files
rm ${prefix}.idxstats.tmp*
rm ${prefix}.depth.tmp*
    # make all the files tab deliminated
sed -i 's/ /\t/g' ${prefix}.idxstats
sed -i 's/ /\t/g' ${prefix}.depth

#Schroeder lab specific - make the read type a seperate collumn
sed -i 's/\.bctrimmedreads/\tbctrimmedreads/g; s/\.correctedReads/\tcorrectedReads/g; s/\.contigs/\tcontigs/g' ${prefix}.idxstats
sed -i 's/\.bctrimmedreads/\tbctrimmedreads/g; s/\.correctedReads/\tcorrectedReads/g; s/\.contigs/\tcontigs/g' ${prefix}.depth

done


##################################################################################
##################################################################################

    # Tidy the outputs

cat *.idxstats > mapping.idxstats.tsv
awk '{if ($3 != "*") print}' mapping.idxstats.tsv > mapping.idxstats.tsv.tmp
mv mapping.idxstats.tsv.tmp mapping.idxstats.tsv
cat *.depth > mapping.depth.tsv

mkdir mapping # clean up all the mapping files
mv *bam mapping/
rm *.depth
rm *.idxstats

##################################################################################
##################################################################################

for f in mapping/*.bam; do
N=$(basename $f)
samtools index $f
## get contig lengths
samtools idxstats $f | cut -f 1,2 | sed '/^\*/d' > ${f}.idxstats
## total coverage per base for each contig
samtools depth $f | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.samtoolsdepth
done

## contig length fix
rm -f tmp.ALLcontig; touch tmp.ALLcontig
cat mapping/*.idxstats >> tmp.ALLcontig
sort -Vu tmp.ALLcontig -o tmp.ALLcontig
## get total number of bases covered at MIN_COVERAGE_DEPTH or higher PER contig
rm -f tmp.ALLstats.list.count0 tmp.ALLstats.list.count1 tmp.ALLstats.list.count2 tmp.ALLstats.list.count3; touch tmp.ALLstats.list.count0 tmp.ALLstats.list.count1 tmp.ALLstats.list.count2 tmp.ALLstats.list.count3
for f in mapping/*.samtoolsdepth; do
    cut -f 1,2 $f | sort -Vu >> tmp.ALLstats.list.count0
    awk '$4 >= 5' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count1
    awk '$4 >= 10' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count2
    awk '$4 >= 100' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count3
done
### 1X coverage
rm -f tmp.ALLstats.list.count0.5; touch tmp.ALLstats.list.count0.5
for f in mapping/*.samtoolsdepth; do
    cut -f 1,2 $f | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count0.5
done

## add headers
## add headers
sed -i 's/\.bam\t/\t/g' tmp.ALLcontig
sed -i 's/\.bam\t/\t/g' tmp.ALLstats.list.count*
sed -i -z 's/^/contig\tcontig_length\n/1' tmp.ALLcontig
sed -i -z 's/^/reads\tcontig\n/1' tmp.ALLstats.list.count0
sed -i -z 's/^/reads\tcontig\ttotal_bases_over_1_coverage\n/1' tmp.ALLstats.list.count0.5
sed -i -z 's/^/reads\tcontig\ttotal_bases_over_5_coverage\n/1' tmp.ALLstats.list.count1
sed -i -z 's/^/reads\tcontig\ttotal_bases_over_10_coverage\n/1' tmp.ALLstats.list.count2
sed -i -z 's/^/reads\tcontig\ttotal_bases_over_100_coverage\n/1' tmp.ALLstats.list.count3

##################################################################################
##################################################################################

    # Taxonomy of contig population using kaiju
kaiju -t rvdb/nodes.dmp \
	-f rvdb/rvdb.fmi \
		-i population.pool.fasta -z ${threads} -o population.pool.kaiju.rvdb.out

sort -t $'\t' -V -k 2,2 population.pool.kaiju.rvdb.out -o population.pool.kaiju.rvdb.out.sorted
kaiju-addTaxonNames \
	-t rvdb/nodes.dmp \
		-n rvdb/names.dmp \
			-i population.pool.kaiju.rvdb.out.sorted \
				-o population.pool.kaiju.rvdb.names \
					-r superkingdom,phylum,order,class,family,genus,species

kaiju -t nr_euk/nodes.dmp \
	-f /panfs/jay/groups/27/dcschroe/shared/bioinfdb/KAIJU/nr_euk/nr_euk.fmi \
		-i population.pool.fasta -z ${threads} -o population.pool.kaiju.nr_euk.out

sort -t $'\t' -V -k 2,2 population.pool.kaiju.nr_euk.out -o population.pool.kaiju.nr_euk.sorted
kaiju-addTaxonNames \
	-t rvdb/nodes.dmp \
		-n rvdb/names.dmp \
			-i population.pool.kaiju.nr_euk.sorted \
				-o population.pool.kaiju.nr_euk.names \
					-r superkingdom,phylum,order,class,family,genus,species

awk '{print $0"\t0\tNA\tNA\tNA\tNA"}' population.pool.kaiju.*.names | cut -f 1-8 | sort -t $'\t' -k 2,2V -k 4,4nr | awk -F "\t" '!a[$2]++' | sed 's/0\tNA\tNA\tNA\tNA//g' > population.pool.kaiju.merge

##################################################################################
##################################################################################

    # Run the R script that outputs all the final tables

mkdir tables

Rscript viral-population-estimations_pt2.R
