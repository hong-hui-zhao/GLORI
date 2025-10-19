

---



```python

import os, sys
import argparse
from collections import defaultdict, Counter
import copy
import sqlite3
import time
from time import strftime
import numpy as np


#=========================#
#      数据库相关函数     #
#=========================#

def create_table(cursor):
    """在内存数据库中创建表，用来存储注释信息"""
    sql = """ CREATE TABLE IF NOT EXISTS 
			locations (
				chr text NOT NULL,       -- 染色体名
				pos integer NOT NULL,    -- 位点坐标
				dir text NOT NULL,       -- 链方向（+ 或 -）
				gene text NOT NULL,      -- 基因名
				trans text NOT NULL,     -- 转录本ID
				name text NOT NULL,      -- 碱基名或特征名
				isoform text NOT NULL,   -- 异构体
				biotype text NOT NULL    -- 基因类型
			); 
		   """
    cursor.execute(sql)


def insert_var(cursor, chr, pos, dir, gene, trans, name, isoform, biotype):
    """插入单条注释记录"""
    cursor.execute(
        "insert into locations (chr, pos, dir, gene, trans, name, isoform, biotype) values ( ?, ?, ?, ?, ?, ?, ?, ?)",
        (chr, pos, dir, gene, trans, name, isoform, biotype))


def insert_many_var(cursor, rows):
    """批量插入注释记录，加快速度"""
    cursor.executemany("""
	insert into locations (chr, pos, dir, gene, trans, name, isoform, biotype) values ( ?, ?, ?, ?, ?, ?, ?, ?)
	""", rows)


def create_index(cursor):
    """为 chr, dir, pos 创建索引，加速查询"""
    cursor.execute("CREATE INDEX locations_index ON locations (chr,dir,pos)")


def build_database():
    """从注释文件构建内存数据库或字典"""
    global database, cursor, conn
    with open(options.database, 'r') as input:
        line = input.readline()
        if options.method == "sql":
            # 使用 SQLite 内存数据库（速度快）
            conn = sqlite3.connect(":memory:")
            cursor = conn.cursor()
            create_table(cursor)
            # 调整数据库参数，提高内存性能
            cursor.execute("PRAGMA cache_size=65536")
            cursor.execute("PRAGMA page_size=65536")
            cursor.execute('PRAGMA temp_store=MEMORY')
            cursor.execute("PRAGMA synchronous=OFF")
            cursor.execute('PRAGMA journal_mode=MEMORY')
        elif options.method == "dict":
            pass  # 使用字典存储模式（备选）

        n = 0
        rows = []
        # 按行读取注释文件
        while (line):
            line = line.strip().split("\t")
            chr = line[0]
            pos_1 = line[2]
            dir = line[3]
            gene = line[4]
            trans = line[5]
            name = line[6]
            isoform = line[7]
            biotype = line[8]

            # 选择不同的存储方式
            if options.method == "sql":
                rows.append([chr, int(pos_1), dir, gene, trans, name, isoform, biotype])
            elif options.method == "dict":
                database[(chr, int(pos_1), dir)] = [gene, name, trans, isoform, biotype]

            n += 1
            # 每读 100 万条就写入一次数据库（防止占太多内存）
            if n % 1000000 == 0:
                if options.method == "sql":
                    insert_many_var(cursor, rows)
                    rows = []
                sys.stderr.write("[%s] %d items processed...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()), n))
            line = input.readline()

        # 把剩余的行写入数据库
        if options.method == "sql" and rows:
            insert_many_var(cursor, rows)
            rows = []

    # 建立索引
    if options.method == "sql":
        sys.stderr.write("[%s] All loaded. Creating index...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        create_index(cursor)
        sys.stderr.write("[%s] All finished!\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


#=========================#
#     统计与输出函数      #
#=========================#

def write_CR():
    """计算并输出每个基因的转换率 (Conversion Rate, CR) 文件"""
    overall = []      # 存储所有位点的non-AG比例
    lines = []

    with open(options.conversion, 'w') as CR_output:
        if all_Aov != 0:   # 确保有统计数据
            # 全局Non-A-to-G比例 = A / (A+G)
            ratio_nonCR = str(all_A / (all_Aov + 0.0))
            CR_output.write("\t".join(["#ALL", str(all_Aov), str(all_A), ratio_nonCR, str(all_mappinglength)]))
            CR_output.write("\n")

            # 遍历每个基因，写入统计值
            for key in sorted(result_cov.keys()):
                Acounts = result_A.get(key)
                cov = result_cov.get(key)
                gene_length = result_genelength.get(key)
                if cov is not None and cov > 0:
                    if Acounts is None:
                        Acounts = 0
                    nonCR = Acounts / (cov + 0.0)  # Non-A-to-G 转换率
                    line = "\t".join([key, str(cov), str(Acounts), str(nonCR), str(gene_length)]) + "\n"
                    lines.append(line)
                    overall.append(nonCR)

            # 计算所有基因的统计量（中位数、均值、分位数）
            if len(overall) > 0:
                CR_output.write("#Median\t%.8f\n" % np.median(overall))
                CR_output.write("#Mean\t%.8f\n" % np.mean(overall))
                CR_output.write("#90%%\t%.8f\n" % np.percentile(overall, 90))
                CR_output.write("#75%%\t%.8f\n" % np.percentile(overall, 75))
                CR_output.write("#50%%\t%.8f\n" % np.percentile(overall, 50))
                CR_output.write("#25%%\t%.8f\n" % np.percentile(overall, 25))
                CR_output.write("#10%%\t%.8f\n" % np.percentile(overall, 10))

            # 写入每个基因详细数据
            for line in lines:
                CR_output.write("#")
                CR_output.write(line)


def search(chr, dir, pos, database=None, cursor=None):
    """查询某个位点属于哪个基因"""
    if options.method == "sql":
        cursor.execute("SELECT gene,name,trans,isoform,biotype FROM locations WHERE chr=? AND dir=? AND pos=?",
                       (chr, dir, pos))
        rows = cursor.fetchall()
        if rows:
            return [str(i) for i in rows[0]]
        else:
            return None
    else:
        return database.get((chr, dir, pos))


def Tideup(chr, pos_1, dir, PIPLEUPS, SURROUNDINGS, database=None, cursor=None):
    """核心函数：统计每个位点的碱基分布并写出结果"""
    global output, all_A, all_Aov, all_mappinglength

    #---------------- 注释关联 ----------------#
    if database != 'Notusingdatabase':
        if options.method == "sql":
            query = search(chr, dir, pos_1, cursor=cursor)
        elif options.method == "dict":
            query = search(chr, dir, pos_1, database=database)
        if query is not None:
            gene, name, trans, isoform, biotype = query
            GENE = gene
        else:
            gene, name, trans, isoform, biotype = ["NA"] * 5
            GENE = "ELSE"
    else:
        gene, name, trans, isoform, biotype = ["NA"] * 5
        GENE = "ELSE"

    #---------------- pileup统计 ----------------#
    bases = sorted(Counter(zip(PIPLEUPS, SURROUNDINGS)).items(), key=lambda x: x[0][1])
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'total': 0, 'AG': 0}  # 各碱基计数
    bins = {i: {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'total': 0, 'AG': 0} for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]}
    limits = [99999, 20, 15, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]  # 用于分区统计

    # 逐个碱基统计
    for (base, surr), count in bases:
        if surr > limits[-1]:
            while surr > limits[-1]:
                index = limits.pop()
                bins[index] = copy.copy(counts)
        if base != "N":              # 排除不确定碱基
            counts[base] += count
            counts["total"] += count
        if base == "A" or base == "G":
            counts["AG"] += count     # 统计与 A/G 相关的 reads

    for item in limits:
        bins[item] = copy.copy(counts)

    #---------------- 写出结果 ----------------#
    output.write("\t".join([chr, str(pos_1), dir, gene, name, trans, isoform, biotype]))  # 基本信息
    output.write("\t")
    output.write("\t".join([
        str(counts['total']), str(counts['AG']),
        str(counts['A']), str(counts['T']),
        str(counts['C']), str(counts['G'])
    ]))
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 99999][:-1]:
        output.write("\t")
        output.write(str(i))
        output.write(";")
        output.write(",".join([
            str(bins[i]['total']),
            str(bins[i]['AG']),
            str(bins[i]['A'])
        ]))
    output.write("\n")

    #---------------- 汇总基因层统计 ----------------#
    if counts['AG'] >= 15:  # 覆盖度足够高才统计
        result_A[GENE] += counts['A']
        result_cov[GENE] += counts['AG']
        result_genelength[GENE] += 1
        all_A += counts['A']
        all_Aov += counts['AG']
        all_mappinglength += 1


def check(chr, pos_1, dir, pileup, surrounding, database, cursor):
    """简单封装，调用 Tideup 处理每个位点"""
    PIPLEUPS = pileup
    SURROUNDINGS = surrounding
    Tideup(chr, pos_1, dir, PIPLEUPS, SURROUNDINGS, database=database, cursor=cursor)


#=========================#
#        主程序入口       #
#=========================#

if __name__ == "__main__":

    description = """输入应为无冗余的mpileup结果"""
    parser = argparse.ArgumentParser(prog="", fromfile_prefix_chars='@', description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # 必需参数
    group_require = parser.add_argument_group("Required")
    group_require.add_argument("-i", "--input", dest="input", required=True, help="输入pileup文件，按位置排序")
    group_require.add_argument("-o", "--output", dest="output", required=True, help="输出文件名前缀")
    group_require.add_argument("--CR", dest="conversion", required=True, help="基因转换率输出文件名")
    group_require.add_argument("--method", dest="method", default="sql", choices=["sql", "dict"],
                               help="数据库存储类型（sql或dict）")

    # 注释文件
    group_anno = parser.add_argument_group("m6A anno")
    group_anno.add_argument("--db", dest="database", default='Notusingdatabase', help="碱基注释文件")

    options = parser.parse_args()

    sys.stderr.write("[%s] Reading input...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    # 如果提供了注释文件，加载数据库
    if options.database != 'Notusingdatabase':
        if options.method == "sql":
            database = None
            cursor = None
            conn = None
        elif options.method == "dict":
            database = {}
            cursor = None
            conn = None
        build_database()
        sys.stderr.write("[%s] In-memory database connection setup\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    else:
        database = options.database
        cursor = None

    # 初始化统计变量
    result_cov = defaultdict(int)
    result_A = defaultdict(int)
    result_genelength = defaultdict(int)
    all_Aov = 0
    all_A = 0
    all_mappinglength = 0

    #=========================#
    #       逐位点处理       #
    #=========================#
    with open(options.input, 'r') as input, open(options.output + "_tmp", 'w') as output:
        line = input.readline()
        while (line):
            row = line.strip().split("\t")
            chr = row[0].split("_AG_converted")[0]  # 去掉后缀
            pos_1 = row[1]
            dir = row[2]
            base = row[3]  # 参考碱基
            pileup = row[4].split(",")              # pileup的碱基序列
            surrounding = [int(i) for i in row[5].split(",")]  # 在read中的相对位置

            # 仅分析正链A或负链T（即参考A）
            if dir == "+" and base == "A":
                check(chr, pos_1, dir, pileup, surrounding, database, cursor)
            elif dir == "-" and base == "T":
                check(chr, pos_1, dir, pileup, surrounding, database, cursor)

            line = input.readline()

    sys.stderr.write("[%s] Finished reading input.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    # 计算转换率并输出
    write_CR()
    sys.stderr.write("[%s] Calculating conversion rates...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    # 合并 CR 文件与格式化结果，生成最终输出
    os.system("cat %s %s > %s" % (options.conversion, options.output + "_tmp", options.output))
    os.remove(options.output + "_tmp")

    # 清理内存
    if options.database == "sql":
        database.close()
    else:
        del database
```

---


> 这个脚本把每个位点的 pileup 信息按染色体分开处理，
> 对参考碱基为 A/T 的位置统计 A→G 转换情况，
> 并结合注释信息计算每个基因及全局的 A→G 转换率。
> 最后输出两个结果文件：
>
> * **formatted.txt**：逐位点统计结果
> * **CR.txt**：每个基因的转换率（Conversion Rate）

---

这部分工作的核心目标是把前面产生的 mpileup/`referbase` 数据变成可以用于识别候选 m6A 位点的、结构化且可统计分析的输入：
**为什么要做**

因为 m6A 修饰会在实验或测序步骤中导致参考位点的碱基替换模式（特别是 A → G 的观测增加）与未修饰位点不同，直接在整段原始 pileup 文本上做统计既难以按基因或位点汇总，也难以控制覆盖度／读长偏倚、难以联合注释信息判断位点生物学背景，所以必须把每个位置的 pileup 信息格式化（把每个位点的碱基计数、读内相对位置、参考链信息、所属基因/转录本等都提取并统一列出），再按基因做汇总统计以得到可靠的“转换率”（Conversion Rate，通常用 A / (A+G) 或其补数来表示），从而为后续的显著性检验和候选 m6A 位点调用提供标准化、可复现的输入。

**做了什么（技术步骤）**

Shell 脚本按染色体逐条拆分全基因组的 `referbase.mpi`，把每条染色体单独写成文件并（若给出）同时按染色体切出对应的单碱基注释子集；然后对每个染色体调用 `m6A_pileup_formatter.py`：脚本首先（可选）把注释文件读入内存 SQLite（也可选字典模式）并建立针对 (chr, dir, pos) 的索引以实现快速位点到基因的映射；随后逐行读取该染色体的 pileup 条目，解析出参考染色体、1-based 位点、链方向、参考碱基、观测到的 pileup 碱基数组及这些碱基在读中的相对位置；脚本只对“参考为 A 的位点”进行分析（对于负链上的等价情况，参考为 T 时也要处理），因为 m6A 导致的 A→G 观测只在参考 A 上有生物学意义。对于每个位点，脚本用 `Counter` 统计各碱基（A、T、C、G、以及总数 total）和 A/G 相关计数（AG），同时按照 read 内相对位置把统计值分箱（代码中用了一系列阈值 1,2,3,...,10,15,20,99999 来记录不同位置累积的计数），这既能检测测序端点偏好也能评估是否存在读内位置相关的偏倚；把位点级别的信息写进格式化输出行（包括染色体、坐标、链、注释字段、total、AG、各碱基计数以及各分箱的累计统计），并且只有当位点上 A+G 覆盖度达到阈值（脚本里用 `counts['AG'] >= 15`）时才会把该位点的 A 与 A+G 计入基因层面的汇总统计，以避免低覆盖位点带来的噪声。脚本同时维护按基因的三个累积量：A 的总数（result_A）、A+G 的覆盖总数（result_cov）以及参与统计的位点数量（result_genelength）；在全部位点处理完后，`write_CR()` 会把全局统计（例如全局 A/(A+G) 比率）、以及每个基因的覆盖、A 数、非转换率（A/(A+G)）和基因长度（被统计的位点数）写入 CR 文件，并计算基因层面的中位数、均值和若干分位数（90/75/50/25/10%），这些汇总指标用于衡量样本或基因组范围内的背景转换分布，便于后续确定哪些基因/位点的转换率显著偏离全局或局部背景；最后脚本把 CR 文件与位点级的临时输出合并成最终输出文件并清理临时文件，整个过程通过日志记录（`5_call_sites.log`）把每一步的标准输出和错误输出保存下来以便重现与排查。

**为什么这么设计（参数与权衡）**

使用内存 SQLite 与索引是为了在处理包含百万级注释条目的文件时仍能快速查询；按染色体分割输入可以并行化/减小单次内存压力并更容易定位问题；只分析参考 A（和负链对应的 T）避免把无关碱基引入统计；设置覆盖阈值（如 15）是典型的折中：既要保证统计稳定性、避免低覆盖导致的假阳性，又不能把真实但低表达的信号完全丢弃；按读内位置分箱可以揭示测序或文库制备引入的系统偏倚（例如若大多数 G 来自读端，很可能是测序错误或比对问题而不是生物学修饰）；把位点与注释关联能让后继的显著性检验或富集分析在基因/转录本层面上进行，从而更容易解释生物学意义，也便于过滤掉位于低复杂区域或重复区域的位点。

**潜在注意事项与局限**

结果高度依赖于输入 mpileup 的质量与比对准确性（错配、多重比对、并入的 PCR duplicate 会扭曲计数），注释文件的坐标系统必须与参考一致（否则注释映射会错位），覆盖阈值和分箱策略可以根据样本的测序深度与实验设计调整；此外，A→G 增加并非 m6A 的唯一原因（例如 RNA 编辑、比对错误或序列上下文特异的技术偏差都可能造成类似模式），因此这个格式化与统计只是候选位点识别的必要但不充分步骤，后续还需要做对照比较、统计检验（例如 Fisher/GLM/混合模型）、复现性验证或实验验证来确认真正的甲基化位点。总之，这一步的价值在于把杂乱的 pileup 原始记录转成按位点与基因可比较的、带质量控制（覆盖阈、读内位置分箱）和可追溯的统计表格，从而为后续严格的候选位点调用与生物学解释奠定了坚实的数据层基础。
