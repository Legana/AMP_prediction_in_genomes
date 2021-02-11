#!/bin/zsh

# split up Homo sapiens proteome into different fasta files containing 100 proteins each
#add 2 leading 0s with %03d to keep file numbering in order

# cat cache/homo_proteome_standardaa.fasta | \
# bioawk -c 'fastx' 'BEGIN{JOBN=0} {  if( (NR-1)%100==0 ){JOBN=JOBN+1; \
# file=sprintf("cache/homo_file_spam/""ampin%03d.fasta",JOBN);} \
# printf(">%s\t%s\n%s\n",$name,$comment,$seq) \
# >> file}'

# without leading zeros as HPC job array doesn't like it
cat cache/homo_proteome_standardaa.fasta | \
bioawk -c 'fastx' 'BEGIN{JOBN=0} {  if( (NR-1)%100==0 ){JOBN=JOBN+1; \
file=sprintf("cache/homo_file_spam/""ampin%d.fasta",JOBN);} \
printf(">%s\t%s\n%s\n",$name,$comment,$seq) \
>> file}'