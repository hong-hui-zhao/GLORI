

---

```python
"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to merge multiple BAM files to one, then sort and index it"""
"""Input: [.bam]"""

import sys, os                    # 系统与文件操作库，用于文件删除、错误输出等
import argparse                   # 命令行参数解析模块，让脚本支持 -i -o 等参数
import pysam                      # pysam 库用于读写、排序、索引 BAM 文件
import time                       # 时间模块，用于打印日志时间戳
from time import strftime          # 格式化时间输出，如 "2025-10-19 17:00:00"

# -------------------------------------------------------------------
# 函数：合并单个 BAM 到输出文件
# -------------------------------------------------------------------
def merge_bam(fin, fout):                                             # fin 是输入 bam 路径，fout 是输出句柄
    with pysam.AlignmentFile(fin, 'rb') as INPUT:                     # 以只读二进制方式打开输入 BAM
        for read in INPUT:                                            # 遍历输入 BAM 中的每条 read
            read.reference_id = lift_over[fin][read.reference_id]     # 将原 reference_id 映射到新 header 的 id
            read.next_reference_id = -1                               # 去掉配对 read 的 reference 信息（防止错配）
            read.next_reference_start = 0                             # 同上，将 mate 的起始位置清零
            fout.write(read)                                          # 写入输出 BAM 文件

# -------------------------------------------------------------------
# 函数：读取 header 并建立映射（lift_over）
# -------------------------------------------------------------------
def read_headers(fn, hid, new_header, hid_dict, lift_over):           
    lift_over[fn] = {}                                                # 初始化当前输入文件的映射表
    with pysam.AlignmentFile(fn, 'rb') as INPUT:                      # 打开输入 BAM 文件
        n = 0                                                         # n 是输入 BAM 内部的 contig 索引计数
        for header in INPUT.header["SQ"]:                             # 遍历 header 中的每个序列 (SQ) 项
            if header['SN'] not in hid_dict:                          # 如果这个染色体名还没在全局字典中
                hid_dict[header['SN']] = hid                          # 分配一个新的全局序列 ID
                new_header['SQ'].append(header)                       # 把该 header 信息加入新的 header 中
                hid += 1                                              # 全局序列计数加一
            lift_over[fn][n] = hid_dict[header['SN']]                 # 建立当前输入的索引 n -> 新索引的映射
            n += 1                                                    # 输入文件中的索引增加
    return hid, new_header, hid_dict, lift_over                       # 返回更新后的结构

# -------------------------------------------------------------------
# 主程序部分
# -------------------------------------------------------------------
if __name__ == "__main__":                                            # 主入口（仅当直接运行此脚本时执行）
    description = ""                                                  # 描述文本（可填帮助说明）
    parser = argparse.ArgumentParser(                                 # 创建参数解析器
        prog="concat_bam", 
        fromfile_prefix_chars='@', 
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # ------------------ 必需参数组 ------------------
    group_required = parser.add_argument_group("Required")            
    group_required.add_argument("-i", "--input", dest="input", nargs='*', required=True,
                                help="Input BAM files")               # 输入 BAM 文件，可多个
    group_required.add_argument("-o", "--output", dest="output", required=True,
                                help="Output BAM file")               # 输出 BAM 文件名

    # ------------------ 可选参数组 ------------------
    group_optional = parser.add_argument_group("Optional")            
    group_optional.add_argument("--sort", dest="sort", default=False, action="store_true",
                                help="Sort BAM (and delete unsorted file)")  # 是否排序
    group_optional.add_argument("--no-del-bam", dest="no_del_bam", default=False, action="store_true",
                                help="Do not delete unsorted BAM file")      # 排序后是否保留原 BAM
    group_optional.add_argument("--index", dest="index", default=False, action="store_true",
                                help="Index the sorted BAM file")            # 是否建立索引
    group_optional.add_argument("-t", "--threads", dest="threads", default=1, type=int,
                                help="Threads for samtools sort (default=1)")# 排序线程数
    group_optional.add_argument("-m", "--memory", dest="memory", default="4G",
                                help="Memory per thread for sorting (default=4G)") # 每线程内存限制
    options = parser.parse_args()                                    # 解析命令行参数

    # ------------------ 初始化结构 ------------------
    hid = 0                                                          # 全局 contig 索引计数器
    hid_dict = {}                                                    # 记录 contig 名称 -> 全局索引 的字典
    lift_over = {}                                                   # 记录每个输入文件的索引映射
    new_header = {}                                                  # 合并后的 BAM 头部
    new_header['HD'] = {'SO': 'unsorted', 'VN': '1.0'}               # 设置 header 版本与排序状态
    new_header['SQ'] = []                                            # 初始化序列信息列表

    # ------------------ 构建新 header ------------------
    for fn in options.input:                                         # 遍历所有输入 BAM
        hid, new_header, hid_dict, lift_over = read_headers(         # 读取 header 并建立索引映射
            fn, hid, new_header, hid_dict, lift_over
        )

    # ------------------ 合并 BAM 内容 ------------------
    with pysam.AlignmentFile(options.output, 'wb', header=new_header) as OUTPUT:  # 打开输出 BAM
        for fn in options.input:                                     # 对每个输入文件执行合并
            merge_bam(fn, OUTPUT)                                    # 写入转换后的 reads

    # ------------------ 可选：排序与索引 ------------------
    if options.sort:                                                 # 如果指定了 --sort 参数
        sys.stderr.write("[%s] Sorting BAM...\n" %                   # 打印时间日志
                         strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        base, ext = os.path.splitext(options.output)                 # 拆分输出文件名 (去掉 .bam)
        sorted_out = f"{base}.sorted{ext}"                           # 新的排序后输出名 (xxx.sorted.bam)

        try:
            if options.threads > 1:                                  # 多线程排序
                pysam.sort("-@", str(options.threads), 
                           "-m", options.memory, 
                           "-o", sorted_out, options.output)
            else:                                                    # 单线程排序
                pysam.sort("-m", options.memory, 
                           "-o", sorted_out, options.output)
        except Exception as e:                                       # 捕获错误
            sys.stderr.write("Error during sorting: %s\n" % str(e))
            sys.exit(1)

        # ------------------ 可选：建立索引 ------------------
        if options.index:                                            # 如果指定了 --index
            sys.stderr.write("[%s] Indexing BAM...\n" %
                             strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            try:
                pysam.index(sorted_out)                              # 调用 pysam.index 生成 .bai 文件
            except Exception as e:
                sys.stderr.write("Error during indexing: %s\n" % str(e))
                sys.exit(1)

        # ------------------ 可选：删除未排序文件 ------------------
        if not options.no_del_bam:                                   # 如果用户没有指定保留原 BAM
            os.remove(options.output)                                # 删除未排序的中间文件节省空间
```

---

### ✅ 功能总结（含解释）

| 功能步骤               | 作用                  | 为什么这样做                                |
| ------------------ | ------------------- | ------------------------------------- |
| **读取 header**      | 获取每个 BAM 的染色体信息（SQ） | 不同文件可能有不同的 contig，需要统一                |
| **建立 `lift_over`** | 记录原始索引 → 新索引 的对应关系  | BAM 内部用数字存 contig，需要重映射               |
| **合并 reads**       | 写入统一 header 的新 BAM  | 确保所有 reads 的 reference_id 与 header 对应 |
| **清除 mate 信息**     | 避免跨文件配对错乱           | 各输入文件来源不同，mate 信息无效                   |
| **排序**             | 按坐标排序 BAM           | 为后续分析（如 IGV、bcftools）做准备              |
| **索引**             | 生成 `.bai` 文件        | 支持快速随机访问                              |
| **删除未排序文件**        | 节省空间                | 合并后排序才是最终产物                           |

这个脚本把 **多个 BAM 文件合并为一个 BAM**，在合并时统一并重写参考序列（header 中的 SQ），保证所有输入 BAM 的 reference ID 在新输出文件中是一致的。合并完成后可选地对输出做 **排序（samtools/pysam sort）**、**建立索引（.bai）**，并可选择是否删除未排序的中间文件。主要用到 `pysam` 来读写 BAM。



# merge_bam 函数

```python
def merge_bam(fin,fout):
	with pysam.AlignmentFile(fn, 'rb') as INPUT: 
		for read in INPUT:
			# name = read.reference_name
			# new_id = hid_dict.get(name)
			# read.reference_id = new_id
			# print(fin,read.reference_id,lift_over[fin][read.reference_id],read)
			# time.sleep(1)
			read.reference_id = lift_over[fin][read.reference_id]
			read.next_reference_id = -1
			read.next_reference_start = 0
			# read.set_tag("TS", "-")
			# if read.next_reference_id is not None:
				# read.next_reference_id = lift_over[fin][read.next_reference_id]
			fout.write(read)
```

1. `def merge_bam(fin,fout):`
   定义一个函数，用来把单个输入 BAM (`fin`) 的所有记录写入到已打开的输出 BAM (`fout`)。`fin` 预期是输入文件路径（字符串），`fout` 是一个已经打开的 `pysam.AlignmentFile` 写句柄。

2. `with pysam.AlignmentFile(fn, 'rb') as INPUT:`
   **意图**：打开输入 BAM 以只读二进制方式遍历 reads。
   **为什么用 with**：`with` 确保文件正确关闭（上下文管理器）。

3. `for read in INPUT:`
   遍历输入 BAM 文件中的每条对齐记录（`pysam.AlignedSegment` 对象）。

4. `read.reference_id = lift_over[fin][read.reference_id]`
   **做什么**：把当前 `read` 的 `reference_id`（一个数字索引，指向输入 BAM header 的 SQ 列表）转换为新 BAM 中对应的 reference id。`lift_over` 是在脚本别处构建的字典（见 `read_headers`），用于映射每个输入 BAM 中的 reference 索引到合并后新 header 中的索引。
   **为什么要这么做**：不同输入 BAM 的 header 中 reference 的顺序或存在的 contig 可能不同。合并时要创建一个包含所有 contig 的统一 header（`new_header['SQ']`），因此原来的 `reference_id`（数字）不能直接沿用，否则读者/工具会指向错误的 contig。必须把数字重映射到新 header 的索引（这就是 `lift_over` 的用途）。

5. `read.next_reference_id = -1`
   **做什么**：把 fragment/paired-end read 的 mate（next）参考 ID 清除/设为 -1（表示没有 mate 的引用）。
   **为什么**：如果不调整 `next_reference_id`，在合并不同文件时 mate 信息可能指向旧的 reference id 或不一致的 contig 索引，导致错误或不一致。脚本作者选择简化处理：移除 mate 的引用信息（把它设置为未指向任何 reference）。这会丢失 mate 位置信息，但通常合并来自独立文件的 read 时，保持 mate 信息的一致性很难；有时把 mate 信息移除比插入错误映射更安全。注：如果输入 BAM 中包含正确的 mate 信息且你想保留 mate 信息，应当相应地也做 `next_reference_id = lift_over[...]` 并设置 `next_reference_start`，但前提是 mate 仍在同一合并后文件中并且可映射。

6. `read.next_reference_start = 0`
   同上，把 mate 的起始位置设置为 0（与设为 -1 一起代表没有 mate 信息）。这避免了不一致的 mate start 导致 downstream 工具出错。

7. 注释：`# read.set_tag("TS", "-")`
   可能作者曾考虑给每个 read 设置一个 tag（如 `TS`）来标注来源或时间，但被注释掉了。

8. 注释掉的另外一段：

   ```
   # if read.next_reference_id is not None:
   #     read.next_reference_id = lift_over[fin][read.next_reference_id]
   ```

   表示如果想保留 mate 的 reference id，可以按同样方式映射 `next_reference_id`。但是必须小心：如果 mate 在别的文件里或者映射后仍然存在一致性问题，这样做可能会产生错误。脚本选择更保守的做法（删除 mate 指向）。

9. `fout.write(read)`
    将修改后的 `read` 写入输出 BAM（`fout` 是写句柄）。`pysam.AlignmentFile` 负责正确写入二进制 BAM 格式。

**总结：** 这个函数的核心功能是遍历输入 BAM，把每条 read 的 reference id 映射到合并后新的 header（通过 `lift_over`），并移除 mate 指向，最后写入输出。主要目的是保证合并后的 BAM header 与每条记录的 reference id 一致。

---

# read_headers 函数

```python
def read_headers(fn,hid,new_header,hid_dict,lift_over):
	lift_over[fn] = {}
	with pysam.AlignmentFile(fn, 'rb') as INPUT:
		n = 0
		for header in INPUT.header["SQ"]:
			if header['SN'] not in hid_dict:
				# print(header['SN'])
				hid_dict[header['SN']] = hid
				new_header['SQ'].append(header)
				hid += 1
			lift_over[fn][n] = hid_dict[header['SN']]
			n += 1
	return hid,new_header,hid_dict,lift_over
```

逐行解释：

1. `def read_headers(fn,hid,new_header,hid_dict,lift_over):`
   读取指定输入 BAM (`fn`) 的 header 中的 SQ（sequence/contig）列表，把新的 contig 添加到 `new_header` 并构建 `lift_over` 映射。`hid` 是下一个可用的新 header 索引（整数），`hid_dict` 是从 contig 名称（`SN`）到全局 ID 的字典。

2. `lift_over[fn] = {}`
   为当前输入文件初始化一个映射表；之后会把输入 BAM 中每个 contig 的原始 index（按遍历顺序）映射到新的全局 `hid`。

3. `with pysam.AlignmentFile(fn, 'rb') as INPUT:`
   打开输入 BAM。`INPUT.header["SQ"]` 是 header 中 contig 列表（每个是 dict，至少含 `SN` 字段）。

4. `n = 0`
   `n` 用来表示输入 BAM header 中 contig 的原始索引（0、1、2...），这是 read 中 `read.reference_id` 使用的索引（数值）。

5. `for header in INPUT.header["SQ"]:`
   遍历所有 contig header（字典，例如 `{'SN': 'chr1', 'LN': 248956422}`）。

6. `if header['SN'] not in hid_dict:`
   如果这个 contig 名称还没有被加入全局 `hid_dict`（即还没见过），就把它加入全局 header 列表中。

7. `hid_dict[header['SN']] = hid`
   把 contig 名称映射到当前全局 id（`hid`）。

8. `new_header['SQ'].append(header)`
   把这个 contig 的 header 条目追加到 `new_header` 的 `SQ` 列表中。这样 `new_header` 最终包含所有输入文件的所有 contig（去重）。

9. `hid += 1`
   增加下一个可用 ID。

10. `lift_over[fn][n] = hid_dict[header['SN']]`
    把输入文件中第 `n` 个 contig（原始索引 `n`）映射到全局 id（`hid_dict[SN]`）。这就是后面在 `merge_bam()` 中使用的映射：`lift_over[fn][read.reference_id]` 会得到新 header 的 id。

11. `n += 1`
    增加输入文件的 contig 原始索引计数。

12. `return hid,new_header,hid_dict,lift_over`
    返回更新后的计数和结构，便于在主流程中累积多个输入文件的 header 信息。

**为什么这样做**：合并 BAM 必须确保输出 BAM 的 header 包含所有可能的 contig，并且每个 read 中的 `reference_id`（数值）要对应到输出 header 的正确位置。因为 `read.reference_id` 在文件内部是整数索引而不是名字，必须通过这种“lift over”把不同文件的整数索引统一到一个全局索引体系上。

---

# 主流程：参数解析与准备

```python
if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="concat_bam",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	# Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="input", nargs='*',required=True,help="Input bam files")
	group_required.add_argument("-o","--output",dest="output",required=True,help="Output bam")
	group_optional = parser.add_argument_group("Optional")
	group_optional.add_argument("--sort",dest="sort",default=False,action="store_true",help="Sort bam (and delete unsort)")
	group_optional.add_argument("--no-del-bam",dest="no_del_bam",default=False,action="store_true",help="Do not del bam file after sorting")
	group_optional.add_argument("--index",dest="index",default=False,action="store_true",help="Index sorted bam")
	group_optional.add_argument("-t","--threads",dest="threads",default=1,type=int,help="Threads for samtools sort, default=1")
	group_optional.add_argument("-m","--memory",dest="memory",default="1G",help="Memory for samtools sort, default=4G")
	options = parser.parse_args()
```

逐行解释：

* 这是标准 Python 脚本入口。`argparse` 配置了必要参数和可选参数：

  * `-i/--input`：一个或多个输入 BAM（`nargs='*'` 允许多个）。`required=True` 强制必须提供。
  * `-o/--output`：输出 BAM 路径（字符串）。
  * `--sort`：布尔开关。如果给了 `--sort`，脚本会对输出 BAM 做排序并产生 `.sorted.bam`。
  * `--no-del-bam`：若给出该开关，则在排序后不要删除原始（未排序的）输出 BAM。
  * `--index`：是否对排序后的 BAM 建立索引（`.bai`）。
  * `-t/--threads`：排序时使用的线程数（传给 `pysam.sort`）。
  * `-m/--memory`：`pysam.sort` 的内存参数（如 `"4G"`）。**注意**：help 字符串写了 `default=4G`，但实际 `default="1G"`，二者不一致（可能是笔误）。应该统一。

**为什么要这样**：用参数让脚本更灵活、可在管线中被调用。排序和索引通常是下游分析（例如使用 `samtools index`、`bcftools`、IGV 查看）常见需求，提供开关而不是强制可以适应不同场景。

---

# 构建全局 header 与合并写入

```python
hid = 0
hid_dict = {}
lift_over = {}
new_header = {}
new_header['HD'] = {'SO': 'unsorted', 'VN': '1.0'}
new_header['SQ'] = []
	
for fn in options.input:
	hid,new_header,hid_dict,lift_over = read_headers(fn,hid,new_header,hid_dict,lift_over)
	
with pysam.AlignmentFile(options.output, 'wb', header = new_header) as OUTPUT: 
	for fn in options.input:
		merge_bam(fn,OUTPUT)
```

逐行解释 & 目的：

* `hid = 0`：全局 contig 索引从 0 开始（与 pysam 的内部索引一致）。

* `hid_dict = {}`：字典把 contig 名称 (`SN`) 映射到全局 id。

* `lift_over = {}`：字典用来为每个输入文件存放索引映射（`lift_over[fn][old_index] = new_index`）。

* `new_header = {}` 初始化合并后 header。`new_header['HD']` 设置 header 数据，`SO='unsorted'` 表明输出未排序（如果后面做了排序，应该改为 `SO: 'coordinate'`）。

* `new_header['SQ'] = []`：保存合并后所有 contig 的列表（去重）。

* `for fn in options.input:` 调用 `read_headers` 逐个输入文件读取 header 并更新 `new_header`、`hid_dict`、`lift_over`。这样 `new_header` 在写入合并输出前就包含所有 contig 信息。

* `with pysam.AlignmentFile(options.output, 'wb', header = new_header) as OUTPUT:`
  以写模式（`'wb'`）打开输出 BAM，并用 `new_header` 作为 header。之后把每个输入文件的 reads 写入该文件。

* `for fn in options.input: merge_bam(fn,OUTPUT)`
  对每个输入 BAM 调用 `merge_bam`，把处理后的 reads 写入 `OUTPUT`。

**为什么要这样**：先预建好包含所有 contig 的 `new_header`，这样在写输出 BAM 时，读写器知道输出 BAM 的 header 并能正确把 `reference_id`（数值）映射到 header 的 SQ 中。必须要先把 header 构造好再写入 reads。

---

# 可选排序、索引、删除原 BAM

```python
if options.sort == True:
	sys.stderr.write("[%s]Sorting bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	if options.threads > 1:
		pysam.sort("-@",str(options.threads),"-m",options.memory,"-o", options.output.replace(".bam",".sorted.bam"),options.output)
	else:
		pysam.sort("-m",options.memory,"-o", options.output.replace(".bam",".sorted.bam"), options.output)
	if options.index == True:
		sys.stderr.write("[%s]Indexing bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		pysam.index(options.output.replace(".bam",".sorted.bam"))
	if options.no_del_bam == False:
		os.remove(options.output)
```

逐行解释与理由：

1. `if options.sort == True:`
   如果用户要求排序，则进行以下操作。

2. `sys.stderr.write("[%s]Sorting bam...\n" % strftime(...))`
   把带时间戳的日志写到标准错误，便于脚本在管线中运行时区分日志输出和正常输出。使用 `stderr` 是良好习惯，尤其当脚本被其他程序捕获 stdout。

3. `if options.threads > 1: pysam.sort("-@",str(options.threads),"-m",options.memory,"-o", output.sorted.bam,output)` else ...
   使用 `pysam.sort`（实际上是调用 `samtools sort` 的包装）来对未排序的 `options.output` 进行排序，输出 `options.output.replace(".bam",".sorted.bam")`。- `-@` 指定线程数，`-m` 指定内存每线程或每临时文件的内存（取决于 samtools 版本）。
   **注意**：

   * `pysam.sort` 接受命令行风格的参数。这里脚本区分是否指定 `-@`（线程）来避免在单线程时传 `-@ 1`（两种写法都可，但作者选择仅在多线程时加 `-@`）。
   * `options.memory` 默认是 `"1G"`（但 help 误写为 4G）。实际可根据数据量调整（大型 BAM 需更大内存或更改临时目录）。

4. `if options.index == True:`
   如果用户要索引，调用 `pysam.index()` 对排序后的 BAM 生成 `.bai` 文件。索引对于随机访问、可视化（IGV）或某些分析工具是必须的。

5. `if options.no_del_bam == False: os.remove(options.output)`
   如果没有指定 `--no-del-bam`（默认会删除），则删除原始未排序的输出 BAM，从而只保留 `.sorted.bam`（节省磁盘空间）。
   **为何这样做**：排序会生成一个新的文件（`.sorted.bam`），原始的合并文件通常只是中间产物。删除它能节省磁盘空间，但也失去未排序的备份数据。脚本提供 `--no-del-bam` 让用户保留原始文件。

**重要注意点**：

* 在删除之前最好确认排序成功并索引完成，否则可能误删有效数据。脚本没有检查排序命令的返回值或异常捕获，这是一个改进点（建议在生产环境加上错误处理）。
* 排序后 `new_header['HD']['SO']` 没有更新为 `coordinate`，理论上可以在写入排序后修改 header 或在排序时 samtools 会在输出 header 中写入正确的 `SO`。

---

这个脚本的设计目的是解决“**不同 BAM 文件 header/contig 不一致** 时如何安全合并”的问题：通过合并所有输入的 SQ 到一个 `new_header`，并用 `lift_over` 把每条 read 的 `reference_id` 重写成新 header 的索引，从而保证输出 BAM 的一致性。排序/索引/删除中间文件是常见的后处理步骤，脚本把它们做成可选项以增加灵活性。
