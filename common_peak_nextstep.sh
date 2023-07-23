#!/bin/bash

echo "请输入文件名列表，以空格或换行符分隔："
read -r file_list

IFS=$' \n' read -ra file_array <<< "$file_list"

first_file="${file_array[0]}"

for ((i=1; i<${#file_array[@]}; i++)); do
    output_file="cross$((i)).bed"
    input_file="${file_array[i]}"
    bedtools intersect -a "$first_file" -b "$input_file" > "$output_file"
    first_file="$output_file"
done