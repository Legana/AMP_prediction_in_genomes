#!/bin/zsh

# get the full paths for ampgram files and modify print statement to add to build_data_package.sh

ls data/prediction_results/ampgram/homo/*.rds > homo1_names.txt 

filename='homo1_names.txt'
while read p; do 
    echo echo $p ">> data_list"
done < $filename > homo1_names_output.txt

rm homo1_names.txt

ls data/prediction_results/ampgram/homo/leftovers/*.rds > homo2_names.txt

filename='homo2_names.txt'
while read p; do 
    echo echo $p ">> data_list"
done < $filename > homo2_names_output.txt

rm homo2_names.txt

ls data/prediction_results/ampgram/arab/*.rds > arab_names.txt

filename='arab_names.txt'
while read p; do 
    echo echo $p ">> data_list"
done < $filename > arab_names_output.txt

rm arab_names.txt