```
"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to build index for GLORI-tools"""
"""Input: [.fasta]"""
# 头部注释：作者、实验室、时间和用途。
# 该脚本用于为 GLORI-tools 构建索引文件。

import sys
import os
import argparse
from Bio import SeqIO                # Biopython 模块，用于读取/写入 fasta 序列
from Bio.Seq import Seq              # 序列对象
from Bio.Seq import reverse_complement
import time
import numpy as np                   # 用于数学计算（STAR 索引参数中用到）
import pysam                         # 用于读取 fasta 文件长度（pysam.FastaFile）
import re                            # 正则表达式替换碱基
import subprocess                    # 执行外部命令（如 bowtie-build, STAR）

# -------------------- 参数解析 --------------------
parser = argparse.ArgumentParser(description = "building index for reference with three bases")

parser.add_argument("-r", "--reads", nargs="?", type=str, default='read2', 
                    help="read1,read2,paired")
# 指定读段类型。GLORI 实验有不同的化学转化方式：
# read1 和 read2 代表不同链的测序方向，决定 T→C 或 A→G 转换规则。

parser.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, 
                    help="reference file")
# 输入参考序列 fasta 文件

parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', 
                    help = "number of alignment threads to launch")
# 并行线程数（传递给比对软件）

parser.add_argument("-t", "--tools", nargs="?", type=str, default='bowtie', 
                    help="bowtie,bowtie2,bwa,hisat2,STAR")
# 选择构建索引的比对工具类型（目前脚本支持 bowtie / bowtie2 / STAR）

parser.add_argument("-mate_length", "--mate_length", nargs="?", type=int, default=100, 
                    help="for STAR building index")
# STAR 特定参数，用于计算索引参数时（未在此版本使用）

parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',
                    help = "--outname_prefix")
# 输出文件名前缀

parser.add_argument("-o", "--outputdir", nargs="?", type=str, default='./', 
                    help="outputdir")
# 输出目录

args = parser.parse_args()

# 将命令行参数赋给变量，便于后续使用
reads = args.reads
reference = args.reference
Threads = args.Threads
tools = args.tools
mate_length = args.mate_length
outputdir = args.outputdir
outname_prx = args.outname_prefix


# -------------------- 函数 1: change_reference --------------------
def change_reference(reads, reference, outputdir, outname_prx):
    """
    根据 read 类型（read1 或 read2）对参考序列进行碱基替换转换。
    read1：T -> C  (即把所有 T 转换成 C)
    read2：A -> G  (即把所有 A 转换成 G)
    这是 GLORI 技术的核心步骤之一，用于模拟实验中实际发生的化学碱基转换。
    """

    # 确定输出文件的前缀名
    if outname_prx != 'default':
        refer_name = outname_prx
    else:
        # 若未指定前缀，则用输入 fasta 文件名（去掉扩展名）
        refer_name = "_".join(os.path.basename(reference).split(".")[:-1])

    # 创建输出目录（若不存在）
    if os.path.exists(outputdir):
        pass
    else:
        os.mkdir(outputdir)

    # -------------------- 如果是 read1 --------------------
    if reads == "read1":
        changed_refer = outputdir + "/" + refer_name + ".TC_conversion.fa"
        # 删除已有文件（防止旧文件残留）
        subprocess.call("rm -f " + changed_refer + " 2>/dev/null", shell=True)
        changed_refer2 = open(changed_refer, 'a')
        # 遍历 fasta 文件中的每条序列
        for record in SeqIO.parse(reference, "fasta"):
            # 将所有 T 替换为 C
            record.seq = Seq(re.sub('T', 'C', str(record.seq).upper()))
            # 在序列 ID 上加上标识后缀
            record.id = record.id + "_TC_converted"
            # 写入新 fasta 文件
            SeqIO.write(record, changed_refer2, "fasta")

    # -------------------- 如果是 read2 --------------------
    elif reads == "read2":
        changed_refer = outputdir + "/" + outname_prx + ".AG_conversion.fa"
        subprocess.call("rm -f " + changed_refer + " 2>/dev/null", shell=True)
        changed_refer2 = open(changed_refer, 'a')
        for record in SeqIO.parse(reference, "fasta"):
            # 将所有 A 替换为 G
            record.seq = Seq(re.sub('A', 'G', str(record.seq).upper()))
            record.id = record.id + "_AG_converted"
            SeqIO.write(record, changed_refer2, "fasta")

    return changed_refer
# ➤ 该函数最终输出一个“碱基转换后的参考 fasta 文件”。
#    这种文件对应 GLORI 实验设计中“实验链的化学转化模型”。
#    不同 reads 需要不同的转换方式，以保持与实际测序数据匹配。


# -------------------- 函数 2: build_index --------------------
def build_index(changed_refer, tool, Threads):
    """
    根据指定的比对软件（bowtie, bowtie2, STAR）对转换后的参考序列建索引。
    """

    # ① Bowtie1 索引
    if tool == "bowtie":
        command2 = 'bowtie-build -q ' + changed_refer + ' ' + changed_refer
        print(command2)
        subprocess.call(command2, shell=True)

    # ② Bowtie2 索引
    elif tool == "bowtie2":
        command2 = 'bowtie2-build -q ' + changed_refer + ' ' + changed_refer
        print(command2)
        subprocess.call(command2, shell=True)

    # ③ STAR 索引
    elif tool == "STAR":
        filedir_STAR = changed_refer[:-3]  # 建索引的输出目录（取去掉 .fa）
        fh = pysam.FastaFile(changed_refer)  # 打开 fasta
        transcriptomesize = sum(fh.lengths)  # 计算参考总长度
        # 计算 STAR 索引参数 Nbases（根据 STAR 官方推荐公式）
        Nbases = int(round(min(14, np.log2(transcriptomesize)/2 - 1)))

        command2 = (
            'STAR --runMode transcriptomeGenerate '
            '-runThreadN ' + str(Threads) +
            ' --transcriptomeDir ' + filedir_STAR +
            ' --transcriptomeFastaFiles ' + changed_refer +
            ' --transcriptomeSAindexNbases ' + str(Nbases) +
            ' --limittranscriptomeGenerateRAM 84807429045'
        )
        print(command2)
        subprocess.call(command2, shell=True)
# ➤ 该函数负责调用对应比对工具的索引构建命令，
#   使 GLORI 能快速在转化后参考序列上进行比对。


# -------------------- 主程序入口 --------------------
if __name__ == "__main__":
    refername2 = os.path.basename(reference)
    print("**********changing transcriptome ************")
    # Step 1: 对参考序列进行碱基转换
    changed_refer = change_reference(reads, reference, outputdir, outname_prx)

    print("**********Building transcriptome index for " + changed_refer + " with " + tools + "************")
    # Step 2: 为转换后的序列构建索引
    build_index(changed_refer, tools, Threads)

```