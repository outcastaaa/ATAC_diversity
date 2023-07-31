#!/bin/bash

input_bam="input.bam"
output_prefix="output"
num_parts=4

read_ids=$(samtools view -H "$input_bam" | grep -v '^@' | cut -f 1 | shuf)

num_reads_per_part=$(( ($(echo "$read_ids" | wc -l) + num_parts - 1) / num_parts ))

temp_dir=$(mktemp -d)

split -d -l $num_reads_per_part <(echo "$read_ids") "$temp_dir/part"

for i in $(seq -w 0 $((num_parts-1))); do
    output_bam="${output_prefix}_part${i}.bam"
    read_ids_file="$temp_dir/part$i"
    samtools view -b -L "$read_ids_file" "$input_bam" > "$output_bam"
done

rm -rf "$temp_dir"
