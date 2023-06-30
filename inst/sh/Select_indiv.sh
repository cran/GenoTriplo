#!/bin/bash

file1=./tmp.txt
file2=$1
file_tmp_indiv=./tmp_indiv.txt
file_tmp_list=./tmp_list.txt
file_res=./tmp.txt

# Prepare file2
pattern1="_[A-Z][0-9][0-9].CEL"
pattern2="_[A-Z][0-9].CEL"
pattern3="_[A-Z][0-9][0-9]"
pattern4="_[A-Z][0-9]"
  
sed "s/$pattern1//g; s/$pattern2//g; s/$pattern3//g; s/$pattern4//g" "$file2" > $file_tmp_list

# Read the column names from the first line of file1
read -r first_line < $file1

column_names=$first_line

indices=()
indices+=(1)
for word in $column_names;
do
  if [ "$(grep -c $word $file_tmp_list)" -gt "0" ];then
    index=$(echo $column_names | awk -v word="$word" '{for(i=1;i<=NF;i++) if ($i==word) print i}')
    indices+=($index)
  fi  
done

cut -f "${indices[*]}" $file1 > $file_tmp_indiv

cp $file_tmp_indiv $file_res
rm $file_tmp_indiv
rm $file_tmp_list

