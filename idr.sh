#!/bin/bash

# 获取输入参数
input_files=("$@")  # 从命令行参数获取输入文件列表
output_dir="${input_files[-1]}"  # 最后一个参数作为输出目录
unset 'input_files[${#input_files[@]}-1]'  # 移除输出目录参数

# 创建输出目录（如果不存在）
mkdir -p "$output_dir"

# 函数：运行IDR分析
run_idr_analysis() {
    local input1="$1"
    local input2="$2"
    local output_name="$3"
    
    idr --samples "$input1" "$input2" \
        --input-file-type narrowPeak \
        --rank p.value \
        --soft-idr-threshold 0.05 \
        --use-best-multisummit-IDR \
        --output-file "${output_dir}/${output_name}_0.05.txt" \
        --log-output-file "${output_dir}/${output_name}_0.05.log" \
        --plot
    
    echo "IDR analysis completed for samples: $input1, $input2"
    echo "Output files: ${output_dir}/${output_name}_0.05.txt, ${output_dir}/${output_name}_0.05.log"
}

# 运行IDR分析
for ((i=0; i<${#input_files[@]}; i++)); do
    for ((j=i+1; j<${#input_files[@]}; j++)); do
        file1="${input_files[i]}"
        file2="${input_files[j]}"
        file1_name=$(basename "$file1" | sed 's/\..*//')  # 提取文件名，去除文件扩展名
        file2_name=$(basename "$file2" | sed 's/\..*//')  # 提取文件名，去除文件扩展名
        output_name="idr_${file1_name}_${file2_name}"
        
        run_idr_analysis "$file1" "$file2" "$output_name"
    done
done

