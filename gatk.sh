#! /bin/bash

# ./Gatk.sh

# files gatk.jar and picard.jar must be places in a folder where script is

FASTQ_1="Reads/test_new_dup_dna_1.fq"
FASTQ_2="Reads/test_new_dup_dna_2.fq"
RES_DIR="results"
INDEX="Ref/Index"
FASTA_FILE="Ref/test.fa"
VARIANTS="Reads/test.dbsnp.vcf.gz"


mkdir $RES_DIR
if [ ${FASTQ_1%_*} == ${FASTQ_2%_*} ]
then
   FILE="$(basename ${FASTQ_1%_*})"
   echo "names of files are the same, continue"
   echo $FILE
else
    echo -e "\n File names are not good, exiting \n"
    exit
fi

sleep 1

TRIMM_FASTA_1="$RES_DIR/fastp/${FILE}_trim1.fq"
TRIMM_FASTA_2="$RES_DIR/fastp/${FILE}_trim2.fq"

mkdir $RES_DIR/fastp
./fastp  --html "$RES_DIR/fastp/fastp.html" --json "$RES_DIR/fastp/fastp.json" -i $FASTQ_1 -I $FASTQ_2 -o $TRIMM_FASTA_1 -O $TRIMM_FASTA_2

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

if ! [ -e "$RES_DIR/bowtie2/$FILE.bam" ]
then
    echo -e "\n Sada cemo napraviti BAM fajl.\n"
    sleep 2
    mkdir $RES_DIR/bowtie2
    (bowtie2 -x $INDEX/moj_index --rg-id test \
                                 --rg SM:test \
                                 --rg LB:GRC \
                                 --rg PL:ILLUMINA \
                                 --rg DS:HiSeq2000 \
                                 -1 $TRIMM_FASTA_1 -2 $TRIMM_FASTA_2 | samtools view -bS - > "results/bowtie2/$FILE.bam") 2> "results/bowtie2/file.log"


else
    echo -e "\n $FILE.bam file already exist\n"
fi

sleep 2


# multiqc results -o results
#
#
echo index fasta file
samtools faidx $FASTA_FILE
echo prepare FASTA genome sequence dictionary with Picard
java -jar picard.jar CreateSequenceDictionary \
      -R $FASTA_FILE


echo sorting bam file
samtools sort "results/bowtie2/$FILE.bam" > "results/bowtie2/$FILE.sort.bam"
#
#
echo marking of duplictes
mkdir $RES_DIR/MarkDuplicates
java -jar picard.jar MarkDuplicates \
       -I results/bowtie2/$FILE.sort.bam \
       -O results/MarkDuplicates/marked_duplicates.bam \
       -M results/MarkDuplicates/marked_dup_metrics.txt

# #
mkdir $RES_DIR/BaseRecalibrator
echo base quality score recalibration
java -jar gatk.jar BaseRecalibrator \
      -I results/bowtie2/$FILE.sort.bam  \
      -R $FASTA_FILE \
      --known-sites $VARIANTS \
      -O results/BaseRecalibrator/$FILE_recal_data.table

java -jar gatk.jar ApplyBQSR \
      -R $FASTA_FILE \
      -I results/bowtie2/$FILE.sort.bam \
      --bqsr-recal-file results/BaseRecalibrator/$FILE_recal_data.table \
      -O results/BaseRecalibrator/output.bam


echo variant calling using GATK
mkdir $RES_DIR/HaplotypeCaller
java -jar gatk.jar  HaplotypeCaller \
     -R $FASTA_FILE -I results/bowtie2/$FILE.sort.bam \
     -O results/HaplotypeCaller/output.gatk.vcf.gz
