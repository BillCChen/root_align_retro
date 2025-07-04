#!/bin/bash

# 配置参数
model_dir="models/TMPtoR_fromUSPTO50K_P2R_high_dpt"
prefix="finetune_model.product-template_step_"
output_prefix="finetune_average_model"

# 用户输入参数
read -p "请输入起始检查点编号 (e.g. 110000): " start_step
read -p "请输入结束检查点编号 (e.g. 150000): " end_step
read -p "请输入间隔步数 (e.g. 10000): " interval
read -p "请输入每组平均的模型数量 (e.g. 5): " group_size

# 计算组数
total_steps=$((end_step - start_step))
groups=$(( (total_steps / interval) / group_size + 1 ))

echo "即将处理 $groups 组模型平均..."

# 主循环
for (( i=0; i<groups; i++ )); do
    # 计算当前组的起始和结束step
    group_start=$((start_step + i*group_size*interval))
    group_end=$((group_start + (group_size-1)*interval))
    
    # 生成模型列表
    models=()
    for (( step=group_start; step<=group_end; step+=interval )); do
        model_path="$model_dir/${prefix}${step}.pt"
        if [ -f "$model_path" ]; then
            models+=("$model_path")
        else
            echo "警告: 模型文件 $model_path 不存在"
        fi
    done
    
    # 执行平均操作
    if [ ${#models[@]} -ge 2 ]; then
        output_file="$model_dir/${output_prefix}_${group_start:0:2}-${group_end:0:2}.pt"
        echo "正在平均 ${#models[@]} 个模型: ${models[@]}"
        echo "输出到: $output_file"
        
        onmt_average_models -output "$output_file" -m "${models[@]}"
        
        if [ $? -eq 0 ]; then
            echo "成功完成第 $((i+1))/$groups 组平均"
        else
            echo "错误: 第 $((i+1)) 组平均失败"
        fi
    else
        echo "跳过第 $((i+1)) 组 - 可用模型少于2个"
    fi
done

echo "所有模型平均操作完成!"