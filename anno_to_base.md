```
import os  # 导入操作系统接口，用于执行系统命令
import sys # 导入系统模块，用于访问命令行参数（这里主要备用）
import argparse  # 导入 argparse，用于解析命令行参数

# ---------------------------------------------
# Step 1: Argument parsing
# ---------------------------------------------
description = """
Convert a gene annotation file to a BED-like file where each base of each exon
is represented as a separate line. Output is sorted by strand and genomic coordinates.
"""  # 对脚本功能的描述

# 创建 ArgumentParser 对象，用于从命令行获取输入输出文件
parser = argparse.ArgumentParser(
    prog="",  # 程序名称，空字符串表示使用默认
    fromfile_prefix_chars='@',  # 支持从文件读取参数，每行一个参数
    description=description,  # 使用上面定义的描述
    formatter_class=argparse.RawTextHelpFormatter  # 保留换行格式
)

# Required arguments（必需参数）
group_require = parser.add_argument_group("Required")  # 创建必需参数组
group_require.add_argument("-i", "--input", dest="input", required=True, help="Input annotation file")  # 输入文件
group_require.add_argument("-o", "--output", dest="output", required=True, help="Output BED file")  # 输出文件

# Optional arguments（可选参数，目前未使用）
group_other = parser.add_argument_group("Other")  # 创建可选参数组

# 解析命令行参数，获取输入输出文件路径
options = parser.parse_args()

# ---------------------------------------------
# Step 2: Process input annotation file
# ---------------------------------------------
with open(options.input, 'r') as input, open(options.output, 'w') as output:  # 打开输入文件和输出文件
    line = input.readline()  # 读取输入文件的第一行
    
    while line:  # 当行不为空时循环处理
        line = line.strip().split("\t")  # 去掉首尾空格并按制表符分割
        
        trans_id = line[1]  # 获取转录本ID（第2列）
        chr = line[2]       # 获取染色体（第3列）
        dir = line[3]       # 获取链方向（第4列，+ 或 -）
        exonCounts = int(line[8])       # 获取外显子数量（第9列）
        exonStarts = line[9].split(",") # 获取外显子起始位置列表（第10列，逗号分割）
        exonEnds = line[10].split(",")  # 获取外显子结束位置列表（第11列，逗号分割）
        gene_id = line[12]              # 获取基因ID（第13列）
        
        for i in range(exonCounts):  # 遍历每个外显子
            start = int(exonStarts[i])  # 外显子起始位置，0-based
            end = int(exonEnds[i])      # 外显子结束位置，1-based
            
            for pos in range(start, end):  # 遍历外显子范围内的每个碱基
                out = "\t".join([chr, str(pos), str(pos+1), dir, gene_id, trans_id])  # 拼接输出行
                output.write(out + "\n")  # 写入输出文件
        
        line = input.readline()  # 读取下一行

# ---------------------------------------------
# Step 3: Sort output file
# ---------------------------------------------
sorted_out = options.output + ".sorted"  # 排序后的输出文件名

# 使用 Linux 系统命令 sort 对输出文件进行排序
# 按列排序规则：
# -k 4,4 -> 按第4列（链方向）排序
# -k 1,2 -> 再按第1列（染色体）和第2列（起始位置）排序
# -S 10G -> 使用10GB内存进行排序
# -T ./   -> 临时文件存储在当前目录
# --parallel 10 -> 使用10线程加速排序
os.system(
    "sort -k 4,4 -k 1,2 -S 10G -T ./ --parallel 10 %s > %s" % (options.output, sorted_out)
)


```

---

### 1️⃣ **输入文件是什么**

脚本处理的输入文件类似 **UCSC genePred 格式**，每一行是一条转录本的注释信息，关键列包括（以0-based索引计）：

| 列号 | 内容                       |
| -- | ------------------------ |
| 1  | 转录本 ID (`transcript_id`) |
| 2  | 染色体 (`chr`)              |
| 3  | 链方向 (`+` 或 `-`)          |
| 8  | 外显子数量 (`exonCounts`)     |
| 9  | 外显子起始位置（用逗号分隔，0-based）   |
| 10 | 外显子结束位置（用逗号分隔，1-based）   |
| 12 | 基因 ID (`gene_id`)        |

**特点**：

* `exonStarts` 是 **0-based**，表示碱基计数从 0 开始。
* `exonEnds` 是 **1-based**，表示外显子结束位置不包含这个坐标本身（类似 BED 文件标准）。

---

### 2️⃣ **目标输出是什么**

输出文件是 **单碱基 BED-like 文件**，每行内容如下：

```
chr   start   end   strand   gene_id   transcript_id
```

* `start` 是碱基的起始位置（0-based）
* `end` 是 `start+1`，也就是每行只覆盖 **一个碱基**
* 每个外显子的每个碱基都会生成一行
* 最终文件可以直接用于覆盖度统计、绘制基因图谱或其他下游分析

---

### 3️⃣ **如何生成单碱基的**

核心逻辑在这几行：

```python
for i in range(exonCounts):  # 遍历每个外显子
    start = int(exonStarts[i])  # 外显子起始位置，0-based
    end = int(exonEnds[i])      # 外显子结束位置，1-based

    for pos in range(start, end):  # 遍历外显子范围内的每个碱基
        out = "\t".join([chr, str(pos), str(pos+1), dir, gene_id, trans_id])
        output.write(out + "\n")  # 写入输出文件
```

**详细解释**：

1. `for i in range(exonCounts)` → 遍历这个转录本的每一个外显子
2. `start = int(exonStarts[i])` → 取得外显子起始位置（0-based）
3. `end = int(exonEnds[i])` → 取得外显子结束位置（1-based）
4. `for pos in range(start, end)` → 生成从 `start` 到 `end-1` 的每个碱基位置

   * 例如一个外显子 `start=100, end=105`，循环会生成 100,101,102,103,104 五个碱基
5. `out = "\t".join([chr, str(pos), str(pos+1), dir, gene_id, trans_id])` → 拼接每行内容

   * `pos` → 起始碱基位置
   * `pos+1` → 结束位置，保证每行只覆盖 **1 个碱基**
6. `output.write(out + "\n")` → 写入输出文件，每个碱基单独占一行

---

### 4️⃣ **为什么叫单碱基级别**

* 因为 **每行只对应一个碱基**，而不是一个外显子区间
* 所有外显子碱基都被拆开成单独的行
* 优点：可以进行单碱基覆盖统计或可视化
* 类似做法在 RNA-seq 或 ChIP-seq 分析中，用于生成 **per-base coverage tracks**

---

### 5️⃣ **排序步骤**

最后，生成的文件需要排序：

```python
os.system(
    "sort -k 4,4 -k 1,2 -S 10G -T ./ --parallel 10 %s > %s" % (options.output, sorted_out)
)
```

* `-k 4,4` → 按链方向排序（+/-）
* `-k 1,2` → 再按染色体和起始位置排序
* 结果是文件按照基因组顺序和链方向排列，方便后续分析
* 排序后文件名为 `原输出名.sorted`

---

### 6️⃣ **总结单碱基处理流程**

1. 读取每条转录本信息
2. 遍历每个外显子
3. 对外显子内每个碱基生成一行 BED-like 数据
4. 写入输出文件
5. 文件生成完成后按链方向和基因组坐标排序

```
geneX   tx1   chr1   +   0   1000   1000   2000   1   1000,   2000,   0   geneX
```

### 脚本生成的输出行格式：

每一行有 **6 列**，对应：

| 列号 | 内容              | 解释                                                   |
| -- | --------------- | ---------------------------------------------------- |
| 1  | `chr`           | 染色体名，例如 `chr1`                                       |
| 2  | `start`         | 外显子碱基的起始位置（0-based）                                  |
| 3  | `end`           | 外显子碱基的结束位置（1-based），保证每行只覆盖一个碱基，所以 `end = start + 1` |
| 4  | `strand`        | 链方向 `+` 或 `-`                                        |
| 5  | `gene_id`       | 基因ID，例如 `geneX`                                      |
| 6  | `transcript_id` | 转录本ID，例如 `tx1`                                       |

---

### 举例说明

假设外显子从 1000 到 1002（0-based start，1-based end），链方向 `+`，基因ID `geneX`，转录本ID `tx1`：

* 脚本会生成如下三行（每行代表一个碱基）：

```
chr1    1000    1001    +    geneX    tx1
chr1    1001    1002    +    geneX    tx1
chr1    1002    1003    +    geneX    tx1
```

解释：

* 第1行 → 碱基在位置 1000
* 第2行 → 碱基在位置 1001
* 第3行 → 碱基在位置 1002
* 每行的 `end = start + 1` 保证单碱基覆盖

---

### 总结

**每一行对应一个外显子上的单碱基**，字段说明如下：

```
chr     start    end    strand    gene_id    transcript_id
```

* **chr** → 碱基所在染色体
* **start** → 碱基起始坐标（0-based）
* **end** → 碱基结束坐标（1-based）
* **strand** → 正链 `+` 或负链 `-`
* **gene_id** → 基因ID
* **transcript_id** → 转录本ID

这样做的好处是可以 **对每个碱基独立分析覆盖度、绘制基因图谱或做单碱基统计**。

