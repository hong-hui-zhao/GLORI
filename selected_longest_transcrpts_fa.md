```


import argparse            # 用来解析命令行参数
import sys                 # 提供对 sys.stdin 等系统级对象的访问
import time                # 计时用（计算脚本运行时间）
import pandas as pd        # 用于表格操作（读取/排序/去重等）
from Bio import SeqIO      # Biopython，用于读取/写入 fasta 序列
import subprocess          # 用于执行外部 shell 命令（这里用于删除文件）
import re                  # 正则表达式（脚本中未大量使用，但保留）
from collections import defaultdict, OrderedDict
# defaultdict：便于自动创建字典的默认值
# OrderedDict：保证插入顺序（用于保存 exons/introns 时保持顺序信息）

# 命令行参数解析器：定义输入输出参数
parser = argparse.ArgumentParser(description="parse the mpileup file")
parser.add_argument("-anno", "--annofile", nargs="?", type=str, default=sys.stdin, help="annofile")
parser.add_argument("-fafile", "--fafile", nargs="?", type=str, default=sys.stdin, help="fafile")
parser.add_argument("-outname_prx", "--outname_prx", nargs="?", type=str, default=sys.stdin, help="outname_prx")
args = parser.parse_args()

# 将解析后的参数绑定到变量
annofile = args.annofile
fafile = args.fafile
outname_prx = args.outname_prx

# -------- 函数：read_anno(fn) -----------
def read_anno(fn):
    output = defaultdict(dict)
    # 以默认字典初始化：output[trans_id] 会自动变成一个 dict
    with open(fn, 'r') as input:
        line = input.readline()
        while (line):
            line = line.strip().split("\t")
            # 假定传入的是脚本前面生成的注释表，每行用制表符分割
            trans_id = line[1]         # 第2列是转录本 ID（脚本前面产物格式）
            chr = line[2]              # 第3列：染色体
            dir = line[3]              # 第4列：链方向 ('+' 或 '-')
            # All to 0-based
            exonStarts = map(lambda x: int(x), line[9].split(",")[:-1])  # 第10列：exonStarts, 以逗号结尾，取除最后空项，转换为 int（0-based）
            exonEnds = map(lambda x: int(x) - 1, line[10].split(",")[:-1])  # 第11列：exonEnds, 原为 1-based，减1 转成 0-based closed end
            # 注：此处把 exonEnds 的坐标减 1，变成闭区间结束坐标
            gene_id = line[12]         # 第13列：gene id（ensg）
            bins = list(zip(exonStarts, exonEnds))  # 将每个 exon 的 start/end 配对成区间列表
            # 初始化输出结构（每个 transcript 保存 chr/strand/ensg/exons/introns）
            enst_start = 0
            output[trans_id]['dir'] = dir
            output[trans_id]['ensg'] = gene_id
            output[trans_id]['chr'] = chr
            output[trans_id]['introns'] = OrderedDict()
            output[trans_id]['exons'] = OrderedDict()
            last_end = None
            last_start = None
            enst_start = -1  # 以 -1 起始，后面每遇到一个 exon 先 +1 再计算
            if dir == "+":
                # 正链：按 bins 的原始顺序处理（5' 到 3'）
                for start, end in bins:
                    enst_start += 1
                    enst_end = enst_start + end - start
                    # 用 transcript 内坐标 (enst_start, enst_end) 作为 key，基因组坐标 (start,end) 作为 value
                    output[trans_id]['exons'][(enst_start, enst_end)] = (start, end)
                    enst_start = enst_end
                    last_end = end
            elif dir == "-":
                # 负链：反转 bins（因为转录本方向是反的），并在 key 顺序上反向记录
                bins = bins[::-1]
                for start, end in bins:
                    enst_start += 1
                    enst_end = enst_start + end - start
                    # 注意这里 key 使用 (enst_end, enst_start) 来表明负链的方向性（脚本作者的实现细节）
                    output[trans_id]['exons'][(enst_end, enst_start)] = (start, end)  # '-' strand is reverse
                    enst_start = enst_end
                    last_start = start
            line = input.readline()
    return output
# read_anno 的作用是把注释表 (.tbl2 格式) 解析成一个嵌套字典：
# output[transcript]['exons'] = OrderedDict(( (tstart,tend) -> (genomic_start, genomic_end) , ... ))

# -------- 函数：get_length(annofile) -----------
def get_length(annofile):
    annotation = read_anno(annofile)   # 先解析注释文件，得到结构化 annotation
    list_output = []
    for i in annotation:
        genome_info = annotation.get(i)
        dir = genome_info['dir']
        gene_name = genome_info['ensg']
        Chr = genome_info['chr']
        genome_info_iter = list(genome_info["exons"].items())  # 列表 [( (tstart,tend), (gstart,gend) ), ...]
        list_maxend = []
        list_sites = []
        for key, values in genome_info_iter:
            # key 是 transcript 内坐标 (tstart,tend)，values 是 genome 坐标 (gstart,gend)
            list_maxend += [key[0], key[1]]   # 收集 transcript 内的起点和终点
            list_sites += [values[0], values[1]]  # 收集 genome 坐标的最小/最大
        len_transcript = max(list_maxend)   # transcript 长度（使用 transcript 内坐标的最大值）
        min_sites = min(list_sites)         # 基因组上最小坐标（exon 左端）
        max_sites = max(list_sites)         # 基因组上最大坐标（exon 右端）
        # 拼一行输出：Chr, min_sites, max_sites, len_transcript, gene_name, dir, transcript_id
        list_output.append("\t".join(map(str, [Chr, min_sites, max_sites, len_transcript, gene_name, dir, i])))
    # 生成输出文件名：把 annofile 的后缀 _change2Ens.tbl2 替换成 .translength
    genefile = annofile.split("_change2Ens.tbl2")[0] + '.translength'
    OUTPUT = open(genefile, 'w+')
    OUTPUT.writelines("\n".join(list_output) + "\n")
    return genefile
# get_length 的作用是：为每个 transcript 计算转录本在基因组上的 min/max 坐标与 transcript 长度，
# 并写出一个 .translength 文件，该文件后续用于选出最长转录本等操作。

# -------- 函数：select_anno(fafile, genefile, outname_prx) -----------
def select_anno(fafile,genefile,outname_prx):
    # 读取 get_length 生成的 genefile，指定列名
    pd_gene = pd.read_csv(genefile, sep="\t", names=['Chr', 'min_sites', 'max_sites', 'Trans_length', 'Gene', 'Strand', 'Transcript'])
    total_gene_list = list(set(pd_gene['Gene'].values.tolist()))
    # 过滤掉以 XR/XM 开头的转录本（RefSeq 的预测或非优选模型，通常 XR/XM 表示预测转录本）
    pd_gene2 = pd_gene[(~pd_gene['Transcript'].str.contains('XR')) & (~pd_gene['Transcript'].str.contains('XM'))]
    # 按 Chr, Gene, Trans_length 降序排序（目的是把最长的 transcript 排到前面）
    pd_data_partial = pd_gene2.sort_values(by=['Chr', 'Gene', 'Trans_length'], ascending=False)
    # 对每个 (Chr,Gene) 保留第一条（即最长的）
    xx_p = pd_data_partial.drop_duplicates(['Chr', 'Gene'], keep='first')
    trans_list_p = xx_p['Transcript'].values.tolist()
    gene_list_p = xx_p['Gene'].values.tolist()
    # 对剩下没被选中的 gene（即没有在 pd_gene2 中被保留的），再按相同逻辑选取
    pd_gene3 = pd_gene[~pd_gene['Gene'].isin(gene_list_p)]
    pd_data_other = pd_gene3.sort_values(by=['Chr', 'Gene', 'Trans_length'], ascending=False)
    xx_o = pd_data_other.drop_duplicates(['Chr', 'Gene'], keep='first')
    xx_t = pd.concat([xx_p,xx_o])  # 合并两部分：优选集 + 备选集
    trans_list_o = xx_o['Transcript'].values.tolist()
    trans_list = trans_list_o + trans_list_p
    # 把最终选出的最长转录本信息写入文件（方便检查）
    xx_t.to_csv(genefile + ".longest.trans", sep="\t", index=False)
    index = 0
    subprocess.call('rm ' + outname_prx, shell=True)  # 删除已有的输出 fasta（如果存在），注意：直接用 rm 可能有风险
    changed_refer2 = open(outname_prx, 'a')  # 以追加模式打开新的输出 fasta 文件
    for record in SeqIO

```

### 输入是什么、要得到什么

脚本的输入是：一个注释表（annofile，通常是前面 gtf2anno.py/类似脚本生成的 .tbl2 文件）和一个包含所有转录本序列的 fasta（fafile）。

目标是：为每个基因选择一个“代表转录本”（通常选择最长的可靠转录本），并从 fasta 中提取这些代表转录本的序列写成新的 fasta（outname_prx）。此外，脚本会输出一个 .translength 文件，包含每个转录本的基因组范围与转录本长度，和一个 .longest.trans 表记录每个基因被选中的代表转录本信息，方便可追溯和检查。

### 为什么要计算 transcript 长度并选最长的？

在许多下游分析（例如基因表达量汇总、功能注释、基因集分析或对齐参考序列选择）中，需要为每个基因选择一个统一的代表转录本来避免重复统计或歧义。常见简单策略是：选择最长转录本（因为它通常包含最多的外显子/功能域，覆盖最完整的编码或转录空间）。

这里脚本还优先排除 XR / XM（RefSeq 的预测/模型号，通常代表低置信度或未验证的预测模型），以优先保留高质量的已注释转录本（例如以 NM/NR 或 Ensembl IDs 为主）。

### 为什么要先解析注释再从 fasta 提取？

注释文件告诉哪条转录本属于哪个基因、转录本在基因组上的范围和转录本长度。根据这些信息可以决定“代表转录本”的名单。接着从 fasta 中提取这些转录本序列，这样得到的输出 fasta 就是每个基因唯一、可代表的序列集合，便于下游分析（比如构建参考转录本集合、做比对索引、提取 CDS/UTR 等）。

实现细节上的考虑

read_anno 把 exon 的 genome 坐标与 transcript 内的坐标对应起来，便于计算 transcript 长度与坐标范围。

get_length 基于这些信息生成 .translength，既可供脚本内部使用，也方便人工审查（例如你可以检查某个基因被选到哪个转录本以及它的长度）。

select_anno 在过滤与排序时先优先非预测转录本（去掉 XR/XM），再在每个基因内按长度选最长；对于没有合规转录本的基因，脚本会在原始全集中再选最长的作为备选（保证每个基因至少有一个代表）。

为什么删除已有输出并重写？

用 subprocess.call('rm ' + outname_prx, shell=True) 删除旧文件确保输出是干净的，不会把新序列追加到旧文件中。但这种做法需要注意安全（若 outname_prx 为空或不受信任会有风险）；更稳妥方式是先判断文件是否存在再删除，或用 with open(outname_prx, 'w') 直接覆写。

总体用途（常见场景）

生成一个用于后续分析的“代表转录本 fasta”集，例如需要构建转录本参考库、做转录本水平的比对或注释、或对每个基因做单一序列代表的序列比对/注释工作。

该脚本是将注释表（GTF → .tbl2）和对应的 fasta（所有转录本序列）整合、筛选、精简的自动化工具。