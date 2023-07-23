#!/bin/bash

# 检查是否提供了文件列表作为参数
if [ -z "$1" ]; then
  echo "请提供文件列表作为参数！"
  exit 1
fi

# 获取文件列表的路径
file_list_path="$1"

# 从文件列表中读取文件名，文件列表包含多个文件名，以空格或换行符分隔
file_list=$(cat "$file_list_path")

# 将文件列表转换为数组，以空格或换行符作为分隔符
IFS=$' \n' read -ra file_array <<< "$file_list"

# 定义计数器变量
count=1

# 外层循环遍历样本a
for sample_a in "${file_array[@]}"; do
    # 内层循环遍历样本b
    for sample_b in "${file_array[@]}"; do
        # 避免对同一样本进行两次比较
        if [ "$sample_a" != "$sample_b" ]; then
            # 生成输出文件名，计数器加1
            output_file="${count}$(basename "$sample_a")_$(basename "$sample_b")_0.5.bed"
            ((count++))
            # 执行bedtools intersect命令
            bedtools intersect -a "$sample_a" -b "$sample_b" -sorted -f 0.5 -r > "$output_file"
        fi
    done
done
