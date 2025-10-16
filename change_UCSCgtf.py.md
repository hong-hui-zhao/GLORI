GLORI 在对参考数据进行了转化，其首先使用NCBI 数据提供的 GTF 文件进行转化，转化为 UCSC 风格的 GTF 文件。为什么不直接从UCSC 网站直接下载的原因如下

虽然 UCSC 网站也提供了 GTF 注释文件，但在实际分析中，很多人仍然选择从 NCBI 下载官方基因组 GTF 文件并自行进行染色体名转换，主要有几个原因。

首先，NCBI 的 GTF 注释对应官方的参考基因组版本（例如 GRCh38.p13），其基因和转录本信息通常是最新、最完整的，而 UCSC 的 GTF 注释版本可能略有延迟或经过裁剪，某些基因或转录本可能缺失或被调整。

其次，NCBI 的 GTF 文件中染色体名称采用原始命名方式（如 `1, 2, MT`），而 UCSC 风格的染色体名称是加了 `chr` 前缀的（如 `chr1, chr2, chrM`），如果直接使用 UCSC 的 GTF，而下游分析使用的是 NCBI 的基因组序列或比对结果，染色体名不匹配就会导致工具无法正确识别或定位基因。通过自行转换，可以在保持 NCBI 最新基因组和注释完整性的同时，将染色体名统一为 UCSC 风格，总的来说，自行转换的目的就是**在保证基因组版本和注释一致性的前提下，实现工具兼容性和数据完整性**，而直接使用 UCSC 文件可能在版本、坐标或注释完整性上存在差异。

### 小鼠的下载途径：https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/GCF_000001635.27-RS_2024_02/
### 人类的下载途径：https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/


转化分为三个步骤
 
step 1： 建立NCBI 染色体名和UCSC 染色体名的映射关系
```
# -----------------------------
# step1_read_assembly_report.py
# -----------------------------

# 导入 defaultdict 数据结构，方便后续存储映射字典
from collections import defaultdict

# 定义函数：读取 assembly_report.txt，生成 NCBI -> UCSC 染色体名映射
def get_change(changefile):
    """
    读取 NCBI 提供的 assembly_report.txt 文件，生成染色体名转换字典
    参数：
        changefile: assembly_report.txt 文件路径
    返回：
        dict_change: 字典，key = NCBI 染色体名, value = UCSC 染色体名
    """
    dict_change = {}  # 初始化空字典，用于存储映射关系
    with open(changefile) as f:  # 打开文件
        for line in f:  # 逐行读取
            if line.startswith("#"):  # 忽略注释行
                continue
            cols = line.strip().split("\t")  # 按制表符拆分每一列
            ncbi = cols[6]  # 第7列是 NCBI 风格染色体名（如 "NC_000001.11"）
            ucsc = cols[9]  # 第10列是 UCSC 风格染色体名（如 "chr1", "chrM"）
            dict_change[ncbi] = ucsc  # 将 NCBI -> UCSC 对应关系加入字典
    return dict_change  # 返回字典

# -----------------------------
# 测试函数是否正常
# -----------------------------
assembly_file = "/path/to/GCF_000001405.39_GRCh38.p13_assembly_report.txt"  # 替换为你的 assembly_report.txt 文件路径
dict_change = get_change(assembly_file)  # 调用函数生成映射字典

# 打印前10条映射，方便检查是否正确
print("前10个映射：", list(dict_change.items())[:10])
```
![[Pasted image 20251013183434.png]]
建立如图的映射关系


step 2 :从原始 GTF 中提取 gene_id、transcript_id、gene_biotype，以及每个 transcript 的 exon 坐标，为后续生成 UCSC 风格注释文件做准备
```
# -----------------------------
# step2_parse_gtf.py
# -----------------------------

import re  # 导入正则表达式模块，用于匹配 GTF 属性列
from collections import defaultdict  # 导入 defaultdict，方便存储嵌套字典

# 初始化全局字典
dict_trans = defaultdict(dict)       # 存储 gene_id -> transcript_id -> exon 区间列表
dict_gene_biotype = defaultdict(str) # 存储 gene_id -> gene_biotype

# 输入 GTF 文件路径
gtf_file = "/path/to/GCF_000001405.39_GRCh38.p13_genomic.gtf"

# 打开 GTF 文件逐行解析
with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"):  # 忽略注释行（以 # 开头的行）
            continue
        
        cols = line.strip().split("\t")  # 按制表符分割每一列
        feature_type = cols[2]           # 第3列表示 feature 类型：gene, transcript, exon 等

        # 从属性列中提取 gene_id 等信息
        # GTF 属性列通常在第9列（cols[8]），格式类似：
        # gene_id "ENSG000001"; transcript_id "ENST000001"; gene_biotype "protein_coding";
        attr = [i.split(" ") for i in cols if "gene_id" in i][0]  # 找到包含 gene_id 的列并拆分成列表
        gene_id = attr[attr.index("gene_id")+1].replace('"','').strip(";")  # 提取 gene_id 并去掉引号和分号

        # -----------------------------
        # 处理 gene 行
        # -----------------------------
        if feature_type == "gene":
            # 提取 gene_biotype 信息
            gene_biotype = attr[attr.index("gene_biotype")+1].replace('"','').strip(";")
            dict_gene_biotype[gene_id] = gene_biotype  # 保存到字典中

        # -----------------------------
        # 处理 exon 行
        # -----------------------------
        elif feature_type == "exon":
            # 提取 transcript_id
            transcript_id = attr[attr.index("transcript_id")+1].replace('"','').strip(";")
            start, end = int(cols[3]), int(cols[4])  # 获取 exon 的起始和结束位置（GTF 是 1-based）
            
            # 将 exon 区间加入字典，如果 transcript 已存在则追加，否则新建列表
            if transcript_id in dict_trans[gene_id]:
                dict_trans[gene_id][transcript_id] += [start, end]
            else:
                dict_trans[gene_id][transcript_id] = [start, end]

# -----------------------------
# 测试输出
# -----------------------------
# 查看前10个 gene 的 gene_biotype
print("前10个 gene_biotype:", list(dict_gene_biotype.items())[:10])

# 查看前1个 gene 的 transcript 区间
print("前1个 gene 的 transcript 区间:", list(dict_trans.items())[:1])

```

![[Pasted image 20251013183511.png]]
提取如图的内容


```
# -----------------------------
# step3_replace_chr_and_output.py
# -----------------------------

# 输入文件路径和输出文件路径
gtf_file = "/path/to/GCF_000001405.39_GRCh38.p13_genomic.gtf"  # 原始 GTF 文件
output_file = "/path/to/output.anno"  # 最终输出文件路径
dict_change = ...  # 从 step1 得到的 NCBI -> UCSC 染色体名映射字典
dict_trans = ...  # 从 step2 得到的 exon 区间字典
dict_gene_biotype = ...  # 从 step2 得到的 gene_biotype 字典

trans_old = ''  # 用于记录上一行处理的 transcript_id，避免重复输出 transcript 行

# 打开 GTF 文件读取，同时打开输出文件写入
with open(gtf_file) as f, open(output_file, 'w') as out:
    for line in f:
        if line.startswith("#"):  # 忽略注释行
            continue

        cols = line.strip().split("\t")  # 按制表符分列
        chr_name = cols[0]  # 原始染色体名
        feature_type = cols[2]  # feature 类型（gene, exon, transcript）
        
        # 从属性列提取 gene_id
        attr = [i.split(" ") for i in cols if "gene_id" in i][0]
        gene_id = attr[attr.index("gene_id")+1].replace('"','').strip(";")

        # 如果染色体名在字典里，则替换为 UCSC 风格
        if chr_name in dict_change:
            cols[0] = dict_change[chr_name]

        # -----------------------------
        # 输出 gene 行
        # -----------------------------
        if feature_type == "gene":
            out.write("\t".join(cols)+"\n")  # gene 行直接写入输出文件

        # -----------------------------
        # 输出 transcript 和 exon 行
        # -----------------------------
        else:
            transcript_id = attr[attr.index("transcript_id")+1].replace('"','').strip(";")
            
            # -----------------------------
            # 生成 transcript 行（只生成一次，避免重复）
            # -----------------------------
            if transcript_id != trans_old and 'unknown_transcript' not in transcript_id:
                # 获取该 transcript 的最小起点和最大终点
                trans_line = cols[:2] + ['transcript', min(dict_trans[gene_id][transcript_id]),
                                          max(dict_trans[gene_id][transcript_id])] + \
                             cols[5:8] + [" ".join(attr[:-2]) + ' gene_biotype "' + dict_gene_biotype[gene_id] + '"']
                # 写入 transcript 行
                out.write("\t".join(map(str, trans_line)) + "\n")
                trans_old = transcript_id  # 更新已处理的 transcript_id

            # -----------------------------
            # 输出 exon 行，并附加 gene_biotype
            # -----------------------------
            cols[-1] += ' gene_biotype "' + dict_gene_biotype[gene_id] + '"'
            out.write("\t".join(map(str, cols)) + "\n")

```


转化后的 GTF 文件与原始 GTF 的主要差别并不在 alignment 坐标本身，而是在**注释信息的完整性和标准化上**。首先，染色体名称被统一为 UCSC 风格，保证下游工具和可视化软件能够正确识别。其次，原本可能缺失的 transcript 行被生成，每个 transcript 的起止坐标由其所有 exon 的最小起点和最大终点确定，使得 gene → transcript → exon 的层级关系完整清晰。此外，每个行的属性列都附加了 gene_biotype 信息，方便后续根据基因类型筛选或统计。总结来看，转化后的文件虽然不会改变比对结果，但在**兼容性、注释完整性和数据结构规范化**方面有明显提升，使得下游 RNA-seq 分析、功能注释和可视化工作更加可靠和方便。

