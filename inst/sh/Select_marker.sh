#!/bin/bash

file1=./tmp.txt # cree par Extraction.sh
file2=$1 # argument donne en entree (path vers la liste de snp)
file_tmp_snp=./tmp_snp.txt # file tempo cree pour que la sortie soit tjr ./tmp.txt
file_res=./tmp.txt # file de sortie

snpList=$(cat "$file2")

awk 'NR==1 {print; next} {snp_name=substr($0, 1, index($0, "-A")-1); if (!snp_name) snp_name=substr($0, 1, index($0, "-B")-1); if (snp_list ~ "\\<" snp_name "\\>") print}' snp_list="$snpList" tmp.txt > tmp_snp.txt

cp $file_tmp_snp $file_res
rm $file_tmp_snp
