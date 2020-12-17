#! /bin/bash

# ./Gatk.sh

FILE1="Reads/test_new_dup_dna_1.fq"
FILE2="Reads/test_new_dup_dna_2.fq"
PUT="Rezultat"
INDEX="Ref/Index"
FASTA_FILE="Ref/test.fa"
VARIANTS="Reads/test.dbsnp.vcf.gz"
# BAM_FILE= "Rezultat/Bowtie/test_new_dup_dna.bam"
# SORT_BAM_FILE= "Rezultat/Bowtie/test_new_dup_dna_sort.bam"


if [ ${FILE1%_*} == ${FILE2%_*} ]
then
   FILE="$(basename ${FILE1%_*})"
   echo $FILE
else
    echo -e "\n Imena fajlova nisu konzistentna.\n"
    exit
fi

sleep 1

OUT1="$PUT/Fastp/${FILE}_trim1.fq"
OUT2="$PUT/Fastp/${FILE}_trim2.fq"

mkdir $PUT/Fastp
./fastp  --html "$PUT/Fastp/fastp.html" --json "$PUT/Fastp/fastp.json" -i $FILE1 -I $FILE2 -o $OUT1 -O $OUT2

sleep 2

if ! [ -e "$INDEX" ]
then
    echo -e "\n Sada cemo napraviti index fajl.\n"
    sleep 2
    mkdir $INDEX
    bowtie2-build Ref/test.fa $INDEX/moj_index

else
    echo -e "\n Index fajl vec postoji.\n"
fi

sleep 2

if ! [ -e "$PUT/Bowtie/$FILE.bam" ]
then
    echo -e "\n Sada cemo napraviti BAM fajl.\n"
    sleep 2
    mkdir $PUT/Bowtie
    (bowtie2 -x $INDEX/moj_index --rg-id test --rg SM:test --rg LB:GRC --rg PL:ILLUMINA --rg DS:HiSeq2000 -1 $OUT1 -2 $OUT2 | samtools view -bS - > "Rezultat/Bowtie/$FILE.bam") 2>"Rezultat/Bowtie/file.log"


else
    echo -e "\n $FILE.bam fajl vec postoji\n"
fi

sleep 2


multiqc Rezultat -o Rezultat

#
# echo index fasta file
# samtools faidx $FASTA_FILE
# echo prepare FASTA genome sequence dictionary with Picard
# java -jar picard.jar CreateSequenceDictionary \
#       -R $FASTA_FILE
#
echo sorting bam file
samtools sort "Rezultat/Bowtie/$FILE.bam" > "Rezultat/Bowtie/$FILE.sort.bam"
echo marking of duplictes
mkdir $PUT/MarkDuplicates
java -jar picard.jar MarkDuplicates \
       -I Rezultat/Bowtie/$FILE.sort.bam \
       -O Rezultat/MarkDuplicates/$FILE_marked_duplicates.bam \
       -M Rezultat/MarkDuplicates/$FILE_marked_dup_metrics.txt

#
mkdir $PUT/BaseRecalibrator
echo base quality score recalibration
java -jar gatk.jar BaseRecalibrator \
      -I Rezultat/Bowtie/$FILE.sort.bam  \
      -R $FASTA_FILE \
      --known-sites $VARIANTS \
      -O Rezultat/BaseRecalibrator/recal_data.table
# #
# # java -jar gatk.jar ApplyBQSR \
# #       -R $FASTA_FILE \
# #       -I Rezultat/Bowtie/$FILE.bam \
# #       --bqsr-recal-file Rezultat/BaseRecalibrator/recal_data.table \
# #       -O Rezultat/BaseRecalibrator/output.bam
