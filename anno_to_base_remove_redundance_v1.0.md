
这段代码的核心目的是**从注释文件中提取每个基因位置的最佳转录本信息，并输出成 bed 文件格式**，在转录本存在多个重叠或不同选择的情况下，优先保留最长或者生物学意义更高的转录本。这在基因组注释分析中非常常见，因为一个基因往往有多个转录本，而不同转录本长度、类型、功能可能不同，如果分析没有规则地选择转录本，会导致统计或下游分析结果混乱。

整个程序可以分成几个模块来看：

1. **参数解析**
   通过 `argparse` 模块，程序接收四个主要参数：

   * `-i` 输入注释文件（anno file），通常是基因组注释表，包含染色体、位置、方向、基因ID、转录本ID等信息；
   * `-o` 输出 bed 文件路径；
   * `-g` 基因列表文件（genelist），用于记录每个转录本长度、基因类型（biotype）、名称等信息；
   * `--score` 可选的 biotype 打分表，用于自定义不同类型基因的优先级。

   解析后，程序会把这些参数存入 `options` 对象，方便后续调用。

2. **初始化数据结构**
   程序定义了几个字典：

   * `biotype` 存储每个基因的类型，比如 `protein_coding`、`lncRNA` 等；
   * `transLen` 存储转录本长度，用于后续选择最长或最短转录本；
   * `geneName`、`isoformName` 分别存储基因名和转录本名；
   * `GENES` 用于临时存储某个基因位置已经选择的最佳转录本信息。

   另外，`order` 字典是**生物类型优先级表**，数值越高表示优先级越高，比如蛋白编码 (`protein_coding`) 优先级最高，伪基因 (`pseudogene`) 优先级最低，这保证在同一位置有多个基因时程序能自动挑选最重要的那个。

3. **读取可选的自定义打分表**
   如果用户提供了 `--score` 文件，就会用文件里的优先级覆盖默认的 `order` 字典。这样可以灵活调整分析规则，比如把某些 lncRNA 优先级调高。

4. **读取基因列表文件**
   这里的 `genelist` 文件类似索引表，每一行包含：转录本名、基因名、转录本ID、基因ID、biotype、长度等信息。程序把这些信息读进对应字典里：

   * `transLen[trans]` 保存长度
   * `biotype[gene]` 保存类型
   * `geneName[gene]` 和 `isoformName[trans]` 保存名称

   这样就可以在处理注释文件时快速查询某个转录本或基因的属性，而不必每次都从文件中重新搜索。

5. **处理注释文件，挑选最佳转录本**
   核心逻辑在读取输入注释文件时完成：

   * 每一行代表一个基因或转录本的注释，包含染色体、起止位置、方向、基因ID、转录本ID等信息；
   * 程序用 `(chr, pos_0, dir)` 作为 key，记录该位置已经选择的最佳转录本；
   * 当同一位置出现新的转录本时，程序会根据以下规则决定是否替换原来的转录本：

     1. **同一基因不同转录本**：优先选择**长度更短或更长**的转录本，具体取决于 biotype 类型和代码逻辑。
     2. **不同基因**：

        * 如果新基因是蛋白编码且原基因不是，则替换；
        * 如果都是蛋白编码，选择**最长转录本**；
        * 如果都是非蛋白编码，按照 `order` 字典的优先级选择，如果优先级相同，再选择更长转录本。

   这个逻辑确保每个位置最终保留的转录本是**最重要、最有代表性的**，避免重复或无意义的转录本干扰分析。

6. **输出结果**
   程序把选出的转录本信息写入输出文件，每行包含：染色体、起始、结束、方向、基因ID、转录本ID、基因名、转录本名、biotype。
   生成的文件可以直接用于下游分析，例如生成 bed 文件用于绘制基因覆盖图、做转录本表达分析等。

---

**为什么要这样做？**

在基因组分析中，一个位置可能有多个重叠转录本（同一基因不同 isoform 或不同基因重叠），如果不选择规则的转录本，很多分析会出现问题：

* **重复计数**：同一位点算多次表达
* **注释混乱**：下游工具无法判断哪条转录本代表这个基因
* **生物学意义丢失**：比如伪基因、非编码 RNA 会覆盖真正的蛋白编码基因

这段代码的目标就是建立一套自动选择“最佳转录本”的规则，综合考虑：

* **biotype 优先级**（蛋白编码 > 非编码 > 伪基因）
* **转录本长度**（最长或最短，根据类型规则）
* **同一位置只保留一个转录本**


```Python
import os, sys  # 导入操作系统相关模块和系统相关模块
import argparse  # 导入参数解析模块，用于命令行参数处理
from collections import defaultdict  # 导入字典默认值功能
import time  # 导入时间模块

description = """
"""  # 脚本描述，可在命令行帮助中显示

parser = argparse.ArgumentParser(prog="", fromfile_prefix_chars='@', description=description, formatter_class=argparse.RawTextHelpFormatter)  
# 创建命令行解析器，支持从文件读取参数

# Required 参数组
group_require = parser.add_argument_group("Required")  # 创建必填参数组
group_require.add_argument("-i","--input",dest="input",required=True,help="anno file")  
# 必填参数：输入注释文件
group_require.add_argument("-o","--output",dest="output",required=True,help="output bed file")  
# 必填参数：输出 bed 文件
group_require.add_argument("-g","--genelist",dest="genelist",required=True,help="genelist, mark the longest transcript")  
# 必填参数：基因列表文件，用于标记最长转录本
group_require.add_argument("--score",dest="biotype_score",required=False,help="Score table for biotypes: biotype[tab]score. If not given, use preset")  
# 可选参数：biotype 打分表文件，如果未提供则使用预设顺序
group_other = parser.add_argument_group("Other")  # 其他可选参数组
options = parser.parse_args()  # 解析命令行参数

# 定义全局字典
biotype = {}       # 存储每个基因的 biotype 类型
transLen = {}      # 存储每个转录本的长度
geneName = {}      # 存储基因名
isoformName = {}   # 存储转录本名
GENES = {}         # 用于存储同一位置的基因信息

# 定义函数：更新 GENES 字典
def changeDict(key, gene_id, trans_id):
    GENES[key] = {  # 根据位置 key 保存基因信息
        'gene': gene_id,  # 基因ID
        'trans': trans_id,  # 转录本ID
        'biotype': biotype.get(gene_id),  # 基因biotype类型
        'length': transLen.get(trans_id),  # 转录本长度
        'geneName': geneName.get(gene_id),  # 基因名称
        'isoformName': isoformName.get(trans_id),  # 转录本名称
        'order': order[biotype.get(gene_id)]  # biotype 优先级顺序
    }

# biotype 优先级字典，数值越大优先级越高
order = {
    '3prime_overlapping_ncrna': 1,
    '3prime_overlapping_ncRNA': 1,
    'antisense': 0,
    'IG_C_gene': 9,
    'IG_C_pseudogene': -1,
    'IG_D_gene': 2,
    'IG_J_gene': 2,
    'IG_J_pseudogene': -1,
    'IG_V_gene': 2,
    'IG_V_pseudogene': -1,
    'lincRNA': 1,
    'lncRNA': 1,
    'miRNA': 1,
    'misc_RNA': -2,
    'Mt_rRNA': 9,
    'Mt_tRNA': 9,
    'polymorphic_pseudogene': -1,
    'processed_pseudogene': -1,
    'processed_transcript': 5,
    'protein_coding': 10,
    'pseudogene': -1,
    'rRNA': 9,
    'sense_intronic': 5,
    'sense_overlapping': -2,
    'snoRNA': 9,
    'snRNA': 9,
    'Spike_in': 10,
    'TR_C_gene': 2,
    'TR_D_gene': 2,
    'IG_D_pseudogene': -1,
    'TR_J_gene': 2,
    'TR_J_pseudogene': -1,
    'TR_V_gene': 2,
    'TR_V_pseudogene': -1,
    'IG_LV_gene': 2,
    'IG_pseudogene': -1,
    'TEC': -2,
    'unprocessed_pseudogene': -1,
    'transcribed_processed_pseudogene': -1,
    'transcribed_unprocessed_pseudogene': -1,
    'transcribed_pseudogene': -1,
    'transcribed_unitary_pseudogene': -1,
    'scaRNA': 9,
    'ribozyme': 0,
    'scRNA': 9,
    'bidirectional_promoter_lncRNA': 0,
    'unitary_pseudogene': -1,
    'macro_lncRNA': 0,
    'sRNA': 0,
    'guide_RNA': 0,
    'pre_miRNA': 9,
    'tRNA': 9,
    'SRP_RNA': 8,
    'ncRNA': 8,
    'nontranslating_CDS': 8,
    'RNase_MRP_RNA': 8,
    'antisense_RNA': 7,
    'C_region': 2,
    'C_region_pseudogene': -1,
    'D_segment': 2,
    'D_segment_pseudogene': -1,
    'J_segment': 2,
    'J_segment_pseudogene': -1,
    'ncRNA_pseudogene': -1,
    'other': 10,
    'RNase_P_RNA': 8,
    'telomerase_RNA': 8,
    'vault_RNA': 2,
    'V_segment': 2,
    'V_segment_pseudogene': -1,
    'Y_RNA': 2
}

# 如果提供了自定义 biotype 打分表，则读取覆盖默认优先级
if options.biotype_score:
    with open(options.biotype_score, 'r') as input:
        for line in input.readlines():
            line = line.strip().split("\t")  # 按制表符拆分
            order[line[0]] = float(line[1])  # 更新对应 biotype 的顺序值

# 读取基因列表文件，保存基因和转录本信息
with open(options.genelist, 'r') as input:
    line = input.readline()  # 读取第一行
    while line:
        line = line.strip().split("\t")  # 按制表符拆分
        if line[-1] != "None":
            length = int(line[-1])  # 转录本长度
        else:
            length = 0  # 如果没有长度，设为0
        trans = line[2]  # 转录本ID
        gene = line[3]   # 基因ID
        transLen[trans] = length  # 保存转录本长度
        biotype[gene] = line[4]  # 保存基因biotype
        geneName[gene] = line[1]  # 保存基因名
        isoformName[trans] = line[0]  # 保存转录本名
        line = input.readline()  # 读取下一行

# 处理注释文件，生成输出
with open(options.input, 'r') as input, open(options.output, 'w') as output:
    line = input.readline()  # 读取第一行
    line = line.strip()  # 去掉首尾空格
    row = line.split("\t")  # 按制表符拆分
    chr = row[0]       # 染色体
    pos_0 = row[1]     # 起始位置
    pos_1 = row[2]     # 结束位置
    gene_id = row[4]   # 基因ID
    type = biotype.get(gene_id)  # 获取biotype类型
    dir = row[3]       # 链方向
    trans_id = row[5]  # 转录本ID
    key = (chr, pos_0, dir)  # 以染色体-起始-方向为key
    print(line)  # 打印当前行
    changeDict(key, gene_id, trans_id)  # 更新GENES字典
    LINE = "\t".join([chr, pos_0, pos_1, dir, gene_id, trans_id, geneName.get(gene_id), isoformName.get(trans_id), type])
    tmp = LINE  # 保存当前行
    line = input.readline()  # 读取下一行
    while line:
        line = line.strip()
        row = line.split("\t")
        chr = row[0]
        pos_0 = row[1]
        pos_1 = row[2]
        gene_id = row[4]
        type = biotype.get(gene_id)
        dir = row[3]
        trans_id = row[5]
        key = (chr, pos_0, dir)
        LINE = "\t".join([chr, pos_0, pos_1, dir, gene_id, trans_id, geneName.get(gene_id), isoformName.get(trans_id), type])
        
        if key not in GENES:  # 如果该位置还没保存过
            output.write(tmp)  # 写入之前保存的行
            output.write("\n")
            tmp = LINE  # 更新tmp
            GENES = {}  # 重置GENES字典
            changeDict(key, gene_id, trans_id)  # 保存新行信息
        else:
            if GENES[key]['gene'] == gene_id:  # 相同基因不同转录本
                if GENES[key]['length'] > transLen.get(trans_id):  # 选择最短转录本
                    tmp = LINE
                    changeDict(key, gene_id, trans_id)
            else:  # 不同基因
                if type == "protein_coding":  # 当前基因是蛋白编码
                    if GENES[key]['biotype'] != "protein_coding":  # 替换非蛋白编码
                        tmp = LINE
                        changeDict(key, gene_id, trans_id)
                    else:  # 都是蛋白编码
                        print(line, GENES[key]['length'], transLen.get(trans_id))
                        if GENES[key]['length'] < transLen.get(trans_id):  # 选择最长转录本
                            tmp = LINE
                            changeDict(key, gene_id, trans_id)
                else:  # 非蛋白编码基因
                    if order.get(type) > GENES[key]['order']:  # 按biotype优先级选择
                        tmp = LINE
                        changeDict(key, gene_id, trans_id)
                    elif order.get(type) == GENES[key]['order']:  # 相同优先级
                        if GENES[key]['length'] < transLen.get(trans_id):  # 选择最长转录本
                            tmp = LINE
                            changeDict(key, gene_id, trans_id)
        line = input.readline()  # 读取下一行
    output.write(tmp)  # 写入最后一行
    output.write("\n")

```
