#! /usr/bin/bash

FILE1="Projects/GATK/Reads/test_new_dup_dna_1.fq"
FILE2="Projects/GATK/Reads/test_new_dup_dna_2.fq"
PUT="Projects/GATK/Rezultat"
INDEX="$PUT/Bowtie/Index"

if [ ${FILE1%_*} == ${FILE2%_*} ]
then 
   FILE="$(basename ${FILE1%_*})" 
   echo $FILE
else
    echo -e "Imena fajlova nisu konzistentna.\n"
    exit
fi

sleep 1

OUT1="$PUT/Fastp/${FILE}_trim1.fq"
OUT2="$PUT/Fastp/${FILE}_trim2.fq"
echo $OUT1

fastp  --html "$PUT/Fastp/fastp.html" --json "$PUT/Fastp/fastp.json" -i $FILE1 -I $FILE2 -o $OUT1 -O $OUT2    


if ! [ -e "$INDEX" ]
then 
    echo -e "Sada cemo napraviti index fajl.\n"
    sleep 2
    mkdir $INDEX
    bowtie2-build Projects/GATK/Ref/test.fa $INDEX/moj_index

else
    echo -e "Index fajl vec postoji.\n"
fi

sleep 2 

if ! [ -e "$PUT/Bowtie/$FILE.bam" ]
then 
    echo -e "Sada cemo napraviti BAM fajl.\n"
    sleep 2 
    (bowtie2 -x $INDEX/moj_index -1 $OUT1 -2 $OUT2 |samtools view -bS - > "Projects/GATK/Rezultat/Bowtie/$FILE.bam") 2>"Projects/GATK/Rezultat/Bowtie/file.log"
    
else
    echo -e "$FILE.bam fajl vec postoji\n"
fi

sleep 2 


cd Projects/GATK/Rezultat
multiqc .

