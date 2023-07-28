#!/bin/bash
# this script is to genrate total diff peaks.


# 从标准输入中读取文件列表
readarray -t files

# 获取脚本所在的文件夹路径
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 迭代处理每个文件
for ((i=0; i<${#files[@]}; i++))
do
    # 去除文件名中的换行符
    file="${files[$i]//$'\n'/}"

    # 构建输出文件名
    output_file="${file%%.*}_totaldiff.bed"
    output_path="$script_dir/$output_file"

    # 构建当前文件作为 -a 的输入
    a_file="$file"

    # 构建其余文件作为 -b 的输入
    b_files=("${files[@]:0:i}" "${files[@]:i+1}")

    # 构建命令
    command="bedtools intersect -a \"$a_file\" -b ${b_files[@]} -sorted -v > \"$output_path\""

    # 执行命令
    eval "$command"
done
