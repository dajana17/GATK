#! /bin/bash

# ./Gatk/Gatk.sh

FILE1="Gatk/Reads/test_new_dup_dna_1.fq"
FILE2="Gatk/Reads/test_new_dup_dna_2.fq"
PUT="Gatk/Rezultat"
INDEX="Gatk/Ref/Index"

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


fastp  --html "$PUT/Fastp/fastp.html" --json "$PUT/Fastp/fastp.json" -i $FILE1 -I $FILE2 -o $OUT1 -O $OUT2    

sleep 2

if ! [ -e "$INDEX" ]
then 
    echo -e "\n Sada cemo napraviti index fajl.\n"
    sleep 2
    mkdir $INDEX
    bowtie2-build Gatk/Ref/test.fa $INDEX/moj_index

else
    echo -e "\n Index fajl vec postoji.\n"
fi

sleep 2 

if ! [ -e "$PUT/Bowtie/$FILE.bam" ]
then 
    echo -e "\n Sada cemo napraviti BAM fajl.\n"
    sleep 2 
    (bowtie2 -x $INDEX/moj_index -1 $OUT1 -2 $OUT2 |samtools view -bS - > "Gatk/Rezultat/Bowtie/$FILE.bam") 2>"Gatk/Rezultat/Bowtie/file.log"
    
else
    echo -e "\n $FILE.bam fajl vec postoji\n"
fi

sleep 2 


multiqc Gatk/Rezultat -o Gatk/Rezultat


