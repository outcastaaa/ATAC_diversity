#!/bin/bash

if [ -z "$1" ]; then
  echo "请提供文件列表作为参数！"
  exit 1
fi

file_list_path="$1"

file_list=$(cat "$file_list_path")

IFS=$'\n' read -ra file_array <<< "$file_list"

count=1

for sample_a in "${file_array[@]}"; do
    for sample_b in "${file_array[@]}"; do
        if [ "$sample_a" != "$sample_b" ]; then
            output_file="${count}$(basename "$sample_a")_$(basename "$sample_b")_0.5.bed"
            ((count++))
            bedtools intersect -a "$sample_a" -b "$sample_b" -sorted -f 0.5 -r > "$output_file"
        fi
    done
done

