#!/bin/zsh

#split up arab proteome into different fasta files containing 100 proteins each

#cat cache/arab_file_spam/arab_proteome_standardaa.fasta | bioawk -c 'fastx' 'BEGIN{JOBN=0} {  if( (NR-1)%100==0 ){JOBN=JOBN+1; file=sprintf("ampin%d.fasta",JOBN);} printf(">%s\t%s\n%s\n",$name,$comment,$seq) >> file}'

cat cache/arab_proteome_standardaa.fasta | \
bioawk -c 'fastx' 'BEGIN{JOBN=000} {  if( (NR-1)%100==0 ){JOBN=JOBN+1; \
file=sprintf("ampin%d.fasta",JOBN);} \
printf(">%s\t%s\n%s\n",$name,$comment,$seq) \
>> file}'



