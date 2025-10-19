
`genelist` 文件是从 GTF（基因组注释文件）和 FASTA（参考序列文件）中提取出的**基因及转录本信息的表格化文件**，通常包括 `trans_name`、`gene_name`、`trans_id`、`gene_id`、`gene_biotype`、染色体、起止坐标、链方向和长度等字段。生成 `genelist` 的主要目的有三个：

1. **方便后续分析**：通过把 GTF 中复杂的文本格式转为结构化表格，可以快速筛选、统计或绘图，例如按染色体统计基因数、筛选特定类型的基因、分析转录本长度分布等。
2. **作为过滤和注释的基础**：例如在后续步骤中，使用 `awk` 或 Python 过滤掉伪染色体或无效染色体，或者用于 RNA-seq、转座子注释、差异表达分析等。
3. **提高可读性与可复现性**：GTF 文件原始格式冗长且属性列难以直接使用，生成 `genelist` 之后，每行就对应一个转录本或基因，信息清晰且便于存档。

生成 `genelist` 的方法是通过 Python 脚本读取 GTF 文件：

* 先解析 GTF 的前 8 列（染色体、起止、方向、类型等）和第9列的 attributes（如 `gene_id`、`gene_name`、`gene_biotype` 等）。
* 根据类型（`gene` 或 `transcript`）提取对应信息，并计算基因或转录本长度。
* 将每条记录保存到列表中，按 `trans_name` 和 `gene_name` 排序，最后写入一个制表符分隔的文本文件，这个文件就是 `genelist`。
* 如果需要，还可以结合 FASTA 文件，将每条序列的长度或 ID 信息整合到 `genelist` 中。


```

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Cong Liu (Yi_lab, Peking University)
Revised & Annotated: ChatGPT (GPT-5)
Date: Oct, 2025
Email: liucong-1112@pku.edu.cn
Usage: extract gene & transcript annotation info from GTF + FASTA
"""

### ----------------------- 导入模块 -----------------------
import argparse              ### argparse 代替 optparse，更现代、可靠
from Bio import SeqIO        ### Biopython 的序列解析器
import re                    ### 正则表达式模块
import sys                   ### 用于异常退出

### ----------------------- 辅助函数 -----------------------
def parse_attributes(attr_field):
    """
    将 GTF 第 9 列的 attributes 字段解析成字典。
    例如：gene_id "AT1G01010"; gene_name "NAC001"; gene_biotype "protein_coding";
    """
    attrs = {}
    for part in attr_field.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r'(\S+)\s+"([^"]+)"', part)
        if m:
            key, value = m.groups()
            attrs[key] = value
    return attrs

### ----------------------- 核心函数 -----------------------
def get_gene_and_transcript_info(gtf_path, fasta_path, output_path):
    ### 初始化
    gene_seen = set()        ### 用 set() 存已见基因，快速去重
    records = []             ### 存放最终结果
    seq_length = {}          ### 记录 FASTA 序列长度

    ### 读取 FASTA
    for seq in SeqIO.parse(fasta_path, "fasta"):
        seq_length[seq.id] = len(seq.seq)
    ### 目的：为每个序列保存长度（备用，或用于比对 transcript_id 与序列长度）

    ### 读取 GTF
    with open(gtf_path, "r") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue       ### 跳过注释行
            row = line.strip().split("\t")
            if len(row) < 9:
                continue       ### 非标准行跳过

            feature_type = row[2]           ### 第 3 列：gene / transcript / exon 等
            chrom = row[0]                  ### 染色体
            start, end = int(row[3]), int(row[4])  ### 起止坐标
            strand = row[6]                 ### 正负链
            attr_field = row[8]             ### 属性字段
            attrs = parse_attributes(attr_field)

            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue                    ### 没有 gene_id 的行直接跳过

            gene_name = attrs.get("gene_name", gene_id)
            gene_biotype = attrs.get("gene_biotype", "unknown")
            length = str(end - start + 1)   ### 按 GTF 闭区间计算长度

            ### ---- 如果是 transcript 级别 ----
            if feature_type == "transcript":
                trans_id = attrs.get("transcript_id", gene_id)
                trans_name = attrs.get("transcript_name", trans_id)
                records.append([trans_name, gene_name, trans_id, gene_id,
                                gene_biotype, chrom, start, end, strand, length])
                gene_seen.add(gene_id)

            ### ---- 如果是 gene 级别 ----
            elif feature_type == "gene" and gene_id not in gene_seen:
                trans_id = trans_name = gene_id   ### 基因级别没有 transcript 信息
                records.append([trans_name, gene_name, trans_id, gene_id,
                                gene_biotype, chrom, start, end, strand, length])
                gene_seen.add(gene_id)

    ### 排序输出
    records_sorted = sorted(records, key=lambda x: (x[0], x[1]))

    ### 写入文件
    with open(output_path, "w") as out:
        header = ["trans_name", "gene_name", "trans_id", "gene_id",
                  "gene_biotype", "chr", "start", "end", "strand", "length"]
        out.write("\t".join(header) + "\n")
        for rec in records_sorted:
            out.write("\t".join(map(str, rec)) + "\n")

    print(f"✅ Finished! Parsed {len(records_sorted)} records and wrote to {output_path}")

### ----------------------- 主程序入口 -----------------------
if __name__ == "__main__":
    ### 设置命令行参数
    parser = argparse.ArgumentParser(
        description="Extract gene/transcript annotation from GTF + FASTA."
    )
    parser.add_argument("-i", "--input", required=True, help="Input GTF file")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    args = parser.parse_args()

    ### 调用主函数
    try:
        get_gene_and_transcript_info(args.input, args.fasta, args.output)
    except Exception as e:
        sys.stderr.write(f"❌ Error: {e}\n")
        sys.exit(1)

```



```bash
awk '$6!~/_/&&$6!="na"' {GRCh38_prefix}_genomic.gtf_change2Ens.genelist > {GRCh38_prefix}_genomic.gtf_change2Ens.genelist2
```

的作用是对上一步 Python 脚本生成的 `genelist` 文件进行筛选，具体逻辑是：保留第6列（也就是染色体列）**不包含下划线 `_`** 且**不等于 `na`** 的行，这样可以去掉伪染色体、未定位序列或者无效标记，从而只保留标准染色体（如 1–22、X、Y、MT）对应的基因或转录本信息。这一步在基因组分析中很常用，因为伪染色体和未定位序列可能干扰后续分析，例如基因注释统计、表达量计算或转座子注释。等价的 Python 处理方式是在写入输出文件前对记录列表进行过滤，例如可以用列表推导式 `records_filtered = [r for r in records_sorted if "_" not in str(r[5]) and str(r[5]).lower() != "na"]`，然后再把 `records_filtered` 写入文件，这样就不需要单独执行 `awk`，过滤逻辑直接整合在 Python 脚本中，既简化了操作，也保证了处理的一致性。
