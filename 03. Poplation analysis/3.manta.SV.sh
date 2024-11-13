#! /bin/bash
#
#
# manta检测SV
author: caominghao  date: 2024.7.5
****
## 1.软连接
```shell
ln -s /home/caominghao/07.HJB_SV/01.bam/*bam ./
ln -s /home/caominghao/07.HJB_SV/01.bam/*.bai ./
ln -s /home/caominghao/02.QCB/Reference_genome/GCF_000003025.6_Sscrofa11.1_genomic.fna ./
ln -s /home/caominghao/02.QCB/Reference_genome/GCF_000003025.6_Sscrofa11.1_genomic.fna.fai ./
```
## 2.激活manta环境
```shell
# 安装manta，源码
wget -c https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
tar -xjf manta-1.6.0.release_src.tar.bz2
# Ensure that CC and CXX are updated to target compiler if needed, e.g.:
#     export CC=/path/to/cc
#     export CXX=/path/to/c++
mkdir build && cd build
../manta-1.6.0.release_src/configure  --jobs=4 --prefix=../manta-1.6.0.release_src/
make -C /home/caominghao/software/manta/build
make -j256 install
conda activate manta


# 导出一个环境
conda env export --name manta --file manta.yaml
# 安装一个环境
conda env create --file manta.yaml

```
## 3.列出文件清单
```shell
for i in *.bam;do echo -e "`basename $i .dedup.bam`" >> list;done 
```
## 4.
```shell
for i in `cat list`;do 
time configManta.py \
--bam ${i}.dedup.bam \
--referenceFasta GCF_000003025.6_Sscrofa11.1_genomic.fna \
--runDir ./${i};
done

time /home/caominghao/software/manta/manta-1.6.0.release_src/bin/configManta.py \
--bam H54A.dedup.bam \
--referenceFasta GCF_000003025.6_Sscrofa11.1_genomic.fna \
--runDir ./H54A

```
## step two
```shell
for i in `cat list`;do /home/caominghao/02.HJB/y01.vcf_sv/01.merge.vcf/03.manta/${i}/runWorkflow.py;done
```


# 需要补充一个外群的manta检测
## 1. 下咋fq数据
```shell
# 安装aspera
conda install -c hcc aspera-cli
# 下载

ascp  -vQT -l 5000m -P33001 -k 1 -i /home/caominghao/anaconda3/envs/ascp/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR173/ERR173203/ERR173203_1.fastq.gz  ./
ascp  -vQT -l 5000m -P33001 -k 1 -i /home/caominghao/anaconda3/envs/ascp/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR173/ERR173203/ERR173203_2.fastq.gz  ./

```

## 合并，使用SURVIVOR
```shell
# （1）过滤
for i in `cat list`;do bcftools view -i 'FILTER=="PASS"' ${i}.manta.vcf > PASS_${i}.manta.vcf;done
# （2）删除不准确数据
for i in `cat list`;do cat PASS_${i}.manta.vcf | grep -v 'IMPRECISE' > PRECISE_PASS_${i}.manta.vcf;done
# (3).保留有SVLEN标签的SV
for i in `cat list`;do cat PRECISE_PASS_${i}.manta.vcf | grep -E '^#|SVLEN=' > SVLEN_PRECISE_PASS_${i}.manta.vcf;done
# (4).合并
# 耗时：6m23.205s 265907行
time SURVIVOR merge vcf_list 1000 2 1 0 1 50 492.vcf
# 耗时：7m15.420s 266255行
time SURVIVOR merge vcf_list 1000 2 1 1 1 50 492.vcf
# 耗时：6m8.754s 183233行
time SURVIVOR merge vcf_list 1000 5 1 1 1 50 492.vcf
# 耗时：6m4.899s 137245行
time SURVIVOR merge vcf_list 1000 10 1 1 1 50 492.vcf
# 耗时：5m53.321s 60431行
time SURVIVOR merge vcf_list 1000 50 1 1 1 50 492.vcf
# 耗时：8m49.984s 36762行
time SURVIVOR merge vcf_list 1000 100 1 1 1 50 492.vcf
# 耗时：9m6.755s 17391行
time SURVIVOR merge vcf_list 1000 200 1 1 1 50 492.vcf
# 耗时：8m52.944s 1415行
time SURVIVOR merge vcf_list 1000 400 1 1 1 50 492.vcf
# 耗时：8m56.156s 647行
time SURVIVOR merge vcf_list 1000 492 1 1 1 50 492.vcf

# 保留10%检测到的sv
# 耗时：1m50.079s 49636行
time SURVIVOR merge vcf_list 1000 30 1 0 1 50 292.vcf
```


## 统计
```shell
MantaBND
MantaDEL
MantaDUP
MantaINS
## 统计每个个体每一种变异的数量
for i in `cat 10-25fold.list`;do for j in MantaDEL MantaDUP MantaINS;do echo -en "${i}\t${j} number is : ";cat SVLEN_PRECISE_PASS_${i}.manta.vcf | grep ${j} | wc -l ;done; done > sv_number.10-25fold.txt
## 路径
E:\04.文件\02.西北农林科技大学\02.文件\6_2023-2024学年第二学期\01.HJB项目\2024.7.5gse\HJB\sv-version
## 准备数据
cat data.csv| awk -F "\t" 'BEGIN {OFS=","} {print $1,$2,$3,$4,$5,$6,$7,$8}' > data.csv2
## 画图的python代码
#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Author @caominghao
# Description @ plot_percentage_stacked_bar_chart.py
# CreateTime @ 2024-07-12 17:55:25


import pandas as pd
import matplotlib.pyplot as plt

def plot_percentage_stacked_bar_chart_with_gaps(file_path, output_path):
    # 读取数据
    df = pd.read_csv(file_path)
    
    # 将相关列转换为数值类型
    columns_to_convert = ["BND", "DEL", "DUP", "INS", "sum"]
    for column in columns_to_convert:
        df[column] = pd.to_numeric(df[column], errors='coerce')
    
    # 移除sum为0或NaN的行
    df = df[df["sum"] > 0].dropna(subset=["sum"])
    
    # 按region列和DEL列排序
    df.sort_values(by=["region", "DEL"], ascending=[True, False], inplace=True)
    
    # 计算百分比
    df["BND_pct"] = df["BND"] / df["sum"] * 100
    df["DEL_pct"] = df["DEL"] / df["sum"] * 100
    df["DUP_pct"] = df["DUP"] / df["sum"] * 100
    df["INS_pct"] = df["INS"] / df["sum"] * 100
    
    # 创建堆积柱状图
    fig, ax = plt.subplots(figsize=(14, 8))
    bar_width = 0.8
    
    # 获取唯一的 region 列表
    regions = df["region"].unique()
    
    # 初始化位置
    current_position = 0
    
    # 定义颜色
    colors = {
        "BND": "blue",
        "DEL": "orange",
        "DUP": "green",
        "INS": "red"
    }
    
    for region in regions:
        region_df = df[df["region"] == region]
        positions = range(current_position, current_position + len(region_df))
        
        # 绘制每个类别的柱状图
        ax.bar(positions, region_df["BND_pct"], bar_width, label="BND" if current_position == 0 else "", color=colors["BND"])
        ax.bar(positions, region_df["DEL_pct"], bar_width, bottom=region_df["BND_pct"], label="DEL" if current_position == 0 else "", color=colors["DEL"])
        ax.bar(positions, region_df["DUP_pct"], bar_width, bottom=region_df["BND_pct"] + region_df["DEL_pct"], label="DUP" if current_position == 0 else "", color=colors["DUP"])
        ax.bar(positions, region_df["INS_pct"], bar_width, bottom=region_df["BND_pct"] + region_df["DEL_pct"] + region_df["DUP_pct"], label="INS" if current_position == 0 else "", color=colors["INS"])
        
        current_position += len(region_df) + 1  # 添加空隙
    
    # 添加标签和标题
    ax.set_xlabel("Sample")
    ax.set_ylabel("Percentage")
    ax.set_title("Percentage Stacked Bar Chart with Gaps Between Regions")
    ax.legend()
    
    # 设置x轴刻度和标签
    ax.set_xticks(range(current_position))
    ax.set_xticklabels(df["FID"].tolist() + [''] * (current_position - len(df)), rotation=90)
    
    # 保存图表为PDF
    plt.tight_layout()
    plt.savefig(output_path, format='pdf')

# 示例调用
file_path = 'data.csv2'  # 请替换为你的数据文件路径
output_path = 'percentage_stacked_bar_chart_with_gaps.pdf'  # 图表输出路径
plot_percentage_stacked_bar_chart_with_gaps(file_path, output_path)


```
## 

