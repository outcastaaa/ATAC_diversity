#!/bin/bash

echo "请输入文件名列表，以空格或换行符分隔："
read -r file_list

IFS=$' \n' read -ra file_array <<< "$file_list"

count=1

for ((i=0; i<${#file_array[@]}; i++)); do
    for ((j=i+1; j<${#file_array[@]}; j++)); do
        sample_a="${file_array[i]}"
        sample_b="${file_array[j]}"
        output_file="${count}$(basename "${sample_a%_pool_merge.bed}_${sample_b%_pool_merge.bed}_0.5.bed")"
        ((count++))
        bedtools intersect -a "$sample_a" -b "$sample_b" -sorted -f 0.5 -r > "$output_file"
    done
done




