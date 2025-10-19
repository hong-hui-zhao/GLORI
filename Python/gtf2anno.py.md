
这个 Python 脚本的核心目的，是将标准的 GTF（Gene Transfer Format）基因注释文件转换成一种更易于下游分析使用的结构化注释文件，特别是针对每个转录本的外显子结构信息。
GTF 文件本身是一种广泛使用的基因注释格式，它记录了基因、转录本、外显子、CDS 等信息，但每条信息分散在文件的不同列和不同行中，并且直接使用 GTF 文件做下游分析（如 RNA-seq 定量分析、可视化或者外显子可变剪接分析）不够方便，因此需要对 GTF 文件进行解析和整理。

脚本首先通过导入 `re` 和 `argparse` 模块实现基本功能，其中 `re` 用于在字符串中匹配和提取 `gene_id`、`transcript_id` 等关键信息，而 `argparse` 则用于在命令行中传入输入文件路径和输出文件路径。

脚本定义了 `get_geneanno` 函数，该函数的作用是根据前面整理好的基因列表生成注释文件，首先打开输出文件，并初始化一个列表用于存储每条注释记录。

然后打开输入的 GTF 文件，逐行读取，按制表符分割每行数据，如果当前行是注释行（以 `#` 开头）或者不是 `gene` 类型，则跳过，因为这些行不包含直接需要的注释信息。

在解析每条 `gene` 行时，脚本提取染色体信息（`chr`）、链方向（`dir`），以及基因的起止位置（`start`, `end`），这里特别注意 GTF 文件的 `start` 是 1-based 坐标，因此需要减 1 转换为 0-based 坐标以兼容下游工具，如 BED 文件格式。

随后，脚本使用正则表达式从属性列中提取 `gene_id` 和 `gene_biotype`，并对字典 `dictInfo` 进行初始化，如果该基因或转录本尚未记录，则创建一个包含外显子起止位置列表、链方向、染色体、基因 ID 的字典条目，然后将当前外显子起止位置加入到对应列表中。整个过程使用了 `try-except` 语句块来跳过缺失字段或无法解析的行，从而提高脚本的鲁棒性，避免在面对不规范的 GTF 文件时报错中断。

在主程序中，脚本再次遍历 GTF 文件，但这次关注的是 `transcript_id` 相关的行，跳过注释和 `gene` 类型行，提取每条转录本的外显子位置信息，并将其存入同样的 `dictInfo` 字典中，每个转录本包含其链方向、染色体、基因 ID 以及所有外显子的起止位置，外显子类型包括 `exon` 和 `tRNAscan`。脚本还尝试提取 `gene_name`，如果没有则用 `gene_id` 代替。完成遍历后，脚本通过收集所有转录本对应的 `gene_id` 来生成一个去重的基因列表 `list_gene`，确保在生成最终注释文件时不会重复记录同一基因。随后，`get_geneanno` 函数被调用，对每个转录本进行处理：首先计算外显子数量，将外显子起止位置按升序排序并以逗号分隔，拼接成一行标准化的注释信息，包括转录本 ID、染色体、链方向、外显子数量、外显子起止坐标、基因 ID 等字段，并写入输出文件。

整个脚本的设计逻辑体现了几个关键目的：首先，它实现了 **GTF 文件到结构化注释文件的转换**，将分散在不同列和行的信息整合到以转录本为单位的字典中，便于后续分析。其次，它精确记录了每个转录本的外显子结构，这是 RNA-seq 定量分析、可变剪接分析和基因结构可视化的基础数据。第三，脚本处理了 **1-based 到 0-based 坐标的转换**，保证生成的注释文件与大多数下游分析工具（如 BED 文件处理工具、可视化工具）兼容。第四，它对异常情况进行了容错处理，例如缺少 `gene_id` 或 `transcript_id` 的行会被自动跳过，从而增强了脚本的鲁棒性。最后，通过去重基因 ID，保证输出文件不会重复记录相同基因，便于筛选特定基因集合进行分析。总而言之，作者通过这一系列步骤实现了从原始 GTF 文件提取、整理、去重并格式化输出转录本和外显子信息的功能，生成的注释文件不仅结构清晰，而且直接可用于下游分析流程，极大地提高了基因注释数据的可操作性和兼容性。


``` Python
import re      # 导入正则表达式模块，用于在字符串中搜索 gene_id 或 transcript_id
import argparse # 导入命令行参数解析模块，用于传入输入输出文件路径

# 定义函数 get_geneanno，用于根据 list_gene 生成注释文件
def get_geneanno(list_gene):
    lw = open(options.output, 'w+')  # 打开输出文件，以写入模式
    list_tt = []                     # 初始化列表，用于存储每条基因注释信息

    # 打开输入的 GTF 文件
    with open(options.input, 'r') as gtf:
        line = gtf.readline()        # 读取文件的第一行
        while line:                  # 循环读取文件的每一行
            row = line.strip().split("\t")  # 按制表符分割每行，得到各列信息
            if line.startswith("#") or row[2] != 'gene':
                # 如果是注释行（#开头）或者不是 gene 类型，跳过
                line = gtf.readline()
                continue

            # 提取染色体信息
            chr = row[0]            
            # 提取链方向（+ 或 -）
            dir = row[6]            

            # 提取基因起止位置
            # GTF 文件 start 是 1-based，需要减 1 转为 0-based
            start, end = int(row[3]) - 1, int(row[4])

            try:
                # 用正则在每列里找包含 gene_id 的字符串，并分割成列表
                row2 = [i.split(" ") for i in row if re.search('gene_id', i)][0]
                # 获取 gene_id 并去掉引号和分号
                gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")

                if gene_id not in list_gene:  # 如果基因不在目标列表里
                    trans_id = gene_id        # 使用 gene_id 作为转录本 ID
                    # 获取基因类型（gene_biotype）
                    gene_biotype = row2[row2.index("gene_biotype") + 1].replace('\"', '').strip(";")

                    # 如果该转录本 ID 尚未记录，则初始化字典
                    if trans_id not in dictInfo:
                        dictInfo[trans_id] = {
                            'start': [],   # 外显子起始位置列表
                            'end': [],     # 外显子结束位置列表
                            'dir': "",     # 链方向
                            'chr': "",     # 染色体
                            'gene_id': ""  # 基因ID
                        }
                        dictInfo[trans_id]['dir'] = dir
                        dictInfo[trans_id]['chr'] = chr
                        dictInfo[trans_id]['gene_id'] = gene_id
                        # 保存外显子起止位置
                        dictInfo[trans_id]['start'].append(start)
                        dictInfo[trans_id]['end'].append(end)
                # 读取下一行
                line = gtf.readline()
            except ValueError:
                # 如果行里没有 gene_id，则跳过
                line = gtf.readline()

        # 遍历 dictInfo 中所有转录本，生成注释行
        for trans_id in dictInfo.keys():
            exonCount = str(len(dictInfo[trans_id]['start']))  # 外显子数量
            # 外显子起始位置，按升序排序并用逗号分隔，末尾加逗号
            exonStarts = ','.join([str(i) for i in sorted(dictInfo[trans_id]['start'])]) + ","
            # 外显子结束位置，按升序排序并用逗号分隔，末尾加逗号
            exonEnds = ','.join([str(i) for i in sorted(dictInfo[trans_id]['end'])]) + ","
            # 拼接注释行
            list_tt.append("\t".join([
                '.', trans_id, dictInfo[trans_id]['chr'], dictInfo[trans_id]['dir'],
                'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
                exonCount, exonStarts, exonEnds,
                ".", dictInfo[trans_id]['gene_id'], '.', '.', '.'
            ]))
    # 将所有注释行写入输出文件
    lw.writelines("\n".join(list_tt) + "\n")
    lw.close()  # 关闭输出文件


if __name__ == "__main__":
    # 命令行参数解析
    usage = "Usage: python %prog -i <gtf> > output.anno"
    parser = argparse.ArgumentParser(usage=usage)
    # 输入 GTF 文件路径
    parser.add_argument("-i", dest="input", help="Input GTF file")
    # 输出注释文件路径
    parser.add_argument("-o", dest="output", help="Output annotation file")
    options = parser.parse_args()

    # 初始化字典和基因列表
    dictInfo = {}  # 用于存储转录本信息
    list_gene = [] # 用于存储基因 ID

    # 第一遍读取 GTF 文件，获取转录本和外显子信息
    with open(options.input, 'r') as gtf:
        line = gtf.readline()  # 读取第一行
        while line:            # 循环每行
            row = line.strip().split("\t")  # 分列
            if line.startswith("#") or row[2] == 'gene':
                # 如果是注释行或者 gene 类型行，跳过
                line = gtf.readline()
                continue

            # 提取染色体和链方向
            chr = row[0]
            dir = row[6]

            # 提取外显子起止位置（0-based）
            start, end = int(row[3]) - 1, int(row[4])

            try:
                # 找到 transcript_id 所在列
                row2 = [i.split(" ") for i in row if re.search('transcript_id', i)][0]
                trans_id = row2[row2.index("transcript_id") + 1].replace('\"', '').strip(";")
                gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
                try:
                    # 尝试提取 gene 名称，如果没有则用 gene_id
                    gene_name = row2[row2.index("gene") + 1].replace('\"', '').strip(";")
                except ValueError:
                    gene_name = gene_id

                type = row[2]  # 获取类型：exon, CDS 等

                # 如果该转录本 ID 尚未记录，则初始化
                if trans_id not in dictInfo:
                    dictInfo[trans_id] = {
                        'start': [], 'end': [], 'dir': "", 'chr': "", 'gene_id': ""
                    }
                    dictInfo[trans_id]['dir'] = dir
                    dictInfo[trans_id]['chr'] = chr
                    dictInfo[trans_id]['gene_id'] = gene_id

                # 如果是 exon 或 tRNAscan，则保存起止位置
                if type == "exon" or type == "tRNAscan":
                    dictInfo[trans_id]['start'].append(start)
                    dictInfo[trans_id]['end'].append(end)
                # 读取下一行
                line = gtf.readline()
            except ValueError:
                # 如果没有 transcript_id 列，跳过
                line = gtf.readline()

    # 收集所有基因 ID 到 list_gene
    for trans_id in dictInfo.keys():
        gene_name2 = dictInfo[trans_id]['gene_id']
        if gene_name2 not in list_gene:
            list_gene.append(gene_name2)

    # 调用函数生成最终注释文件
    get_geneanno(list_gene)


```


## 🔹 一、先看 GTF 文件结构

GTF 文件的典型一行大致如下（以 Ensembl 为例）：

```
NC_000001.11    BestRefSeq      gene    11874   14409   .       +       .       gene_id "DDX11L1"; db_xref "GeneID:100287102"; db_xref "HGNC:HGNC:37102"; description "DEAD/H-box helicase 11 like 1"; gbkey "Gene"; gene "DDX11L1"; gene_biotype "transcribed_pseudogene"; pseudo "true"; 
```

GTF 文件每行分为 **9 列**（用 `\t` 分隔）：

| 列号  | 名称         | 示例内容                                                                                                                                                                                                       | 说明                                  |
| --- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------- |
| 1   | seqname    | NC_000001.11                                                                                                                                                                                               | 染色体(需要跟NCBI对应)                      |
| 2   | source     | BestRefSeq                                                                                                                                                                                                 | 数据来源                                |
| 3   | feature    | gene                                                                                                                                                                                                       | 特征类型（gene, transcript, exon, CDS 等） |
| 4   | start      | 11874                                                                                                                                                                                                      | 起始位置                                |
| 5   | end        | 14409                                                                                                                                                                                                      | 终止位置                                |
| 6   | score      | .                                                                                                                                                                                                          | 分值或空                                |
| 7   | strand     | +                                                                                                                                                                                                          | 链方向                                 |
| 8   | frame      | .                                                                                                                                                                                                          | 读码框                                 |
| 9   | attributes | gene_id "DDX11L1"; db_xref "GeneID:100287102"; db_xref "HGNC:HGNC:37102"; description "DEAD/H-box helicase 11 like 1"; gbkey "Gene"; gene "DDX11L1"; gene_biotype "transcribed_pseudogene"; pseudo "true"; | 额外属性字段                              |

前 8 列结构固定，第 9 列 `"attributes"` 存放了多个键值对，用分号 `;` 分隔，用引号包裹值。

---

## 🔹 二、代码逐行解释


```python
row2 = [i.split(" ") for i in row if re.search('transcript_id', i)][0]
trans_id = row2[row2.index("transcript_id") + 1].replace('\"', '').strip(";")
gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
try:
    gene_name = row2[row2.index("gene") + 1].replace('\"', '').strip(";")
except ValueError:
    gene_name = gene_id
```

假设这时：

```python
row = [
 'NC_000001.11'  ,  'BestRefSeq'   ,   'gene' ,   '11874' ,  '14409'  , '.'   ,    '+'   ,    '.'   ,    'gene_id' ,'DDX11L1;' ,'db_xref' ,'GeneID:100287102;' ,'db_xref' ,'HGNC:HGNC:37102;' 'description "DEAD/H-box helicase 11 like 1";' 'gbkey' ,'"Gene";' 'gene', '"DDX11L1";' ,'gene_biotype' ,'"transcribed_pseudogene";' 'pseudo', "true"
]
```

---

### 第一步：

```python
row2 = [i.split(" ") for i in row if re.search('transcript_id', i)][0]
```

解释：

* `for i in row`：遍历当前行的每个字段。
  因为 `row` 是通过 `row = line.strip().split("\t")` 得到的，所以第 9 列是整个属性字段。
* `if re.search('transcript_id', i)`：用正则表达式查找包含 “transcript_id” 的那一列。
  一般情况下，只有第 9 列（attributes）会匹配,取exon。
* `i.split(" ")`：把这一整列按照空格分割成单词列表。

举例：

```python
i = 'gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_biotype "pseudogene";'
i.split(" ") = [
  'gene_id', '"ENSG00000223972";', 'transcript_id', '"ENST00000456328";',
  'gene_name', '"DDX11L1";', 'gene_biotype', '"pseudogene";'
]
```

* `[ ... ][0]` 表示取列表的第一个元素（因为只会找到一个匹配）。

此时：

```python
row2 = [
  'gene_id', '"ENSG00000223972";', 'transcript_id', '"ENST00000456328";',
  'gene_name', '"DDX11L1";', 'gene_biotype', '"pseudogene";'
]
```

---

### 第二步：提取 transcript_id

```python
trans_id = row2[row2.index("transcript_id") + 1].replace('\"', '').strip(";")
```

解释：

* `row2.index("transcript_id")` 找到 `"transcript_id"` 所在的位置索引。
  在上例中，它是索引 2。
* `+ 1` 取其后一个元素：`'"ENST00000456328";'`。
* `.replace('\"', '')` 删除引号 `" "`。
* `.strip(";")` 去掉末尾的分号 `;`。

结果：

```python
trans_id = 'ENST00000456328'
```

---

### 第三步：提取 gene_id

```python
gene_id = row2[row2.index("gene_id") + 1].replace('\"', '').strip(";")
```

逻辑完全相同：

* 找到 `"gene_id"` 的位置（索引 0）。
* 取下一个元素 `"ENSG00000223972";`。
* 去掉引号和分号。

结果：

```python
gene_id = 'ENSG00000223972'
```

---

### 第四步：尝试提取 gene_name

```python
try:
    gene_name = row2[row2.index("gene") + 1].replace('\"', '').strip(";")
except ValueError:
    gene_name = gene_id
```

解释：

* `row2.index("gene")` 尝试在列表中查找 `"gene"`。
  注意，这里有点模糊，因为正确的键应该是 `"gene_name"`，但代码写的是 `"gene"`。
  因此，如果列表中确实有 `"gene_name"`，`index("gene")` 仍然能匹配成功（因为 "gene" 是 "gene_name" 的前缀）。
* 取它后一个元素，并去掉引号和分号，得到基因名称。
* 如果找不到（即出现 `ValueError` 异常），说明这一行没有基因名称字段，则默认用 `gene_id` 代替。

结果：

```python
gene_name = 'DDX11L1'
```

---

## 🔹 三、总结逻辑流程图

| 步骤 | 操作                                        | 结果示例                                                                                |
| -- | ----------------------------------------- | ----------------------------------------------------------------------------------- |
| 1  | 找到包含 `"transcript_id"` 的那一列（即 attributes） | `'gene_id "ENSG..."; transcript_id "ENST..."; gene_name "DDX11L1";'`                |
| 2  | 用空格分割成词列表                                 | `['gene_id', '"ENSG..."', 'transcript_id', '"ENST..."', 'gene_name', '"DDX11L1";']` |
| 3  | 找到 `"transcript_id"` 索引并取下一个值             | `'ENST00000456328'`                                                                 |
| 4  | 找到 `"gene_id"` 索引并取下一个值                   | `'ENSG00000223972'`                                                                 |
| 5  | 尝试提取 `"gene"` 或 `"gene_name"`             | `'DDX11L1'`                                                                         |
| 6  | 如果没有 gene_name，就用 gene_id 替代              | `'ENSG00000223972'`                                                                 |

---

## 🔹 四、为什么作者要这样做

1. **GTF 第 9 列不是固定结构**：
   不同版本的 GTF 文件中属性顺序可能不同，有的可能写在前，有的在后，所以不能用固定列索引。
   → 因此作者使用正则搜索 `re.search('transcript_id', i)` 来灵活匹配包含目标关键词的那一列。

2. **用空格分割提取键值**：
   GTF 的属性列是由空格分隔的键值对形式，直接用 `split(" ")` 可快速获得键和值。

3. **用字符串替换清理格式**：
   每个值都带有引号 `"` 和分号 `;`，必须删除才能得到干净的 ID。

4. **容错设计**：
   并非所有 GTF 文件都包含 `gene_name` 字段，例如一些自动注释的文件只包含 `gene_id`。
   → 所以 `try-except` 确保即使缺少 gene_name，也不会报错，而是自动使用 gene_id 代替。

