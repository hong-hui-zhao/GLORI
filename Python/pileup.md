
```python

# 导入系统库
import sys  # 系统相关操作
import os   # 文件路径操作
import argparse  # 命令行参数解析
import pysam  # 处理 BAM/SAM 文件
from Bio.Seq import reverse_complement  # 获取互补序列
from Bio import SeqIO  # 读取 fasta 文件
import multiprocessing  # 多进程处理
from multiprocessing import Process, Pool  # 创建进程池
import signal  # 捕获 Ctrl+C 信号
import time  # 时间相关操作
from time import strftime  # 格式化时间

# 定义信号处理函数，按 Ctrl+C 时结束进程
def signal_handler(sig, frame):
    pool.terminate()  # 终止所有进程
    sys.exit()  # 退出程序

# worker 函数，每个进程执行
def worker(contig, start, stop):
    global parent_pid, dictRefSeq, options  # 全局变量
    pid = os.getpid()  # 获取当前进程 PID
    tmp_file_name = 'tmp_' + str(parent_pid) + '_' + str(pid)  # 临时文件名，防止冲突

    # 根据 contig 判断最大深度
    max_depth = options.max_depth  # 默认最大深度
    if "GL" in contig:  # 如果是 rRNA 或特殊 contig
        max_depth = options.rRNA_max_depth  # 使用 rRNA 最大深度
    else:
        max_depth = options.max_depth

    # 打开 BAM 文件和临时文件
    with pysam.AlignmentFile(options.input, 'rb') as samfile, open(tmp_file_name, 'a') as tmp_file:
        try:
            # 对指定区域做 pileup
            for pileupcolumn in samfile.pileup(
                contig,
                start=start,
                stop=stop,
                max_depth=max_depth,
                ignore_orphans=False,
                ignore_overlaps=False,
                min_base_quality=0,
                truncate=True
            ):
                chr = pileupcolumn.reference_name  # 当前染色体/contig
                ref_seq = dictRefSeq.get(chr)  # 获取参考序列
                ref_base = ref_seq[pileupcolumn.pos].upper()  # 当前位点参考碱基
                pos = pileupcolumn.pos + 1  # 1-based 位置

                # 初始化正链和负链的字典
                PF_positive_base = {}  # 正链碱基
                PF_negative_base = {}  # 负链碱基
                PF_positive_qual = {}  # 正链质量
                PF_negative_qual = {}  # 负链质量
                PF_positive_A = {}  # 正链 A 数量
                PF_negative_A = {}  # 负链 A 数量

                # 遍历覆盖该位点的所有 reads
                for pileupread in pileupcolumn.pileups:
                    if pileupread.query_position is not None and pileupread.alignment.query_qualities[pileupread.query_position] >= options.qual:
                        not_at_end = True  # 默认不在 read 头尾
                        # 判断是否在 read 的头尾需要修剪
                        if pileupread.alignment.is_reverse:  # 负链
                            if options.trim_tail and pileupread.query_position < options.trim_tail:
                                not_at_end = False
                            if not_at_end and options.trim_head and pileupread.alignment.query_length - options.trim_head <= pileupread.query_position:
                                not_at_end = False
                        else:  # 正链
                            if options.trim_tail and pileupread.alignment.query_length - options.trim_tail < pileupread.query_position:
                                not_at_end = False
                            if not_at_end and pileupread.query_position < options.trim_head:
                                not_at_end = False

                        # 如果通过头尾过滤
                        if not_at_end:
                            # 获取 read 信息
                            query_name = pileupread.alignment.query_name  # read 名称
                            query_base = pileupread.alignment.query_sequence[pileupread.query_position]  # 当前碱基
                            query_qual = pileupread.alignment.query_qualities[pileupread.query_position]  # 当前碱基质量

                            # 处理负链
                            if pileupread.alignment.is_reverse:
                                A_counts = pileupread.alignment.query_sequence.count('T')  # 负链用 T 统计 A 数量
                                query_base = reverse_complement(query_base)  # 取互补
                                if PF_negative_base.get(query_name) is None:  # 第一次出现
                                    PF_negative_base[query_name] = query_base
                                    PF_negative_qual[query_name] = query_qual
                                    PF_negative_A[query_name] = A_counts
                                else:  # 已存在 read
                                    lastRead = PF_negative_base.get(query_name)
                                    lastAcount = PF_negative_A.get(query_name)
                                    if lastRead != query_base:  # 碱基冲突
                                        if options.omit:  # 冲突直接丢弃
                                            PF_negative_base.pop(query_name)
                                            PF_negative_qual.pop(query_name)
                                            PF_negative_A.pop(query_name)
                                        else:  # 取质量更高的
                                            lastQual = PF_negative_qual.get(query_name)
                                            if lastQual < query_qual:
                                                PF_negative_base[query_name] = query_base
                                                PF_negative_qual[query_name] = query_qual
                                                PF_negative_A[query_name] = A_counts
                                    elif lastRead == query_base:  # 碱基相同，比较 A 数量
                                        if lastAcount < PF_negative_A.get(query_name):
                                            PF_negative_base[query_name] = query_base
                                            PF_negative_qual[query_name] = query_qual
                                            PF_negative_A[query_name] = A_counts
                            else:  # 处理正链
                                A_counts = pileupread.alignment.query_sequence.count('A')  # 正链统计 A
                                if PF_positive_base.get(query_name) is None:  # 第一次出现
                                    PF_positive_base[query_name] = query_base
                                    PF_positive_qual[query_name] = query_qual
                                    PF_positive_A[query_name] = A_counts
                                else:  # 已存在 read
                                    lastRead = PF_positive_base.get(query_name)
                                    lastAcount = PF_positive_A.get(query_name)
                                    if lastRead != query_base:  # 碱基冲突
                                        if options.omit:
                                            PF_positive_base.pop(query_name)
                                            PF_positive_qual.pop(query_name)
                                            PF_positive_A.pop(query_name)
                                        else:
                                            lastQual = PF_positive_qual.get(query_name)
                                            if lastQual < query_qual:
                                                PF_positive_base[query_name] = query_base
                                                PF_positive_qual[query_name] = query_qual
                                                PF_positive_A[query_name] = A_counts
                                    elif lastRead == query_base:  # 相同碱基比较 A 数量
                                        if lastAcount <= PF_positive_A.get(query_name):
                                            PF_positive_base[query_name] = query_base
                                            PF_positive_qual[query_name] = query_qual
                                            PF_positive_A[query_name] = A_counts

                # 输出正链 pileup
                if len(PF_positive_base) > 0:
                    positive_seq = PF_positive_base.values()
                    positive_A_seg = [str(i) for i in PF_positive_A.values()]
                    tmp_file.write("\t".join([chr, str(pos), '+', ref_base, ','.join(positive_seq), ','.join(positive_A_seg)]))
                    tmp_file.write("\n")

                # 输出负链 pileup
                if len(PF_negative_base) > 0:
                    negative_seq = PF_negative_base.values()
                    negative_A_seg = [str(i) for i in PF_negative_A.values()]
                    ref_base_rev = reverse_complement(ref_base)  # 负链参考碱基取互补
                    tmp_file.write("\t".join([chr, str(pos), '-', ref_base_rev, ','.join(negative_seq), ','.join(negative_A_seg)]))
                    tmp_file.write("\n")

        except ValueError:
            print("Contig [%s] does not exist in @SQ header, pass" % contig)


# 主程序入口
if __name__ == "__main__":
    description = "Pileup genome using multiprocessing"
    parser = argparse.ArgumentParser(
        prog="pileup_genome_multiprocessing",
        fromfile_prefix_chars='@',
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
    )

    # 必填参数
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i", dest="input", required=True, help="Sorted bam and indexed bam input")
    group_required.add_argument("-o", dest="output", default="genome_merge.pileup.txt", help="Output text")
    group_required.add_argument("-f", dest="fasta", nargs="*", required=True, help="Fasta files, can be multiple (-f a.fa b.fa ...)")

    # 可选参数
    group_optional = parser.add_argument_group("Optional")
    group_optional.add_argument("-P", dest="process", type=int, default=1, help="Process number, default=1")
    group_optional.add_argument("-q", dest="qual", type=int, default=10, help="Quality filter, default=10")
    group_optional.add_argument("-s", dest="step", type=int, default=100000, help="Steps for pileup, default=100000")
    group_optional.add_argument("-m", "--max-depth", dest="max_depth", type=int, default=10000000, help="Max depth for pileup, default=10,000,000")
    group_optional.add_argument("-M", "--max-depth-rRNA", dest="rRNA_max_depth", type=int, default=10000, help="Max depth for rRNA pileup, default=10,000")
    group_optional.add_argument("--omit-confilct", dest="omit", action="store_true", default=False, help="Omit conflict base, default is using higher sequencing quality")
    group_optional.add_argument("--trim-head", dest="trim_head", type=int, default=0, help="Trim the head of reads, default=0")
    group_optional.add_argument("--trim-tail", dest="trim_tail", type=int, default=0, help="Trim the tail of reads, default=0")

    # 其他
    group_other = parser.add_argument_group("Other")
    group_other.add_argument("--version", action="version", version="%(prog)s 1.3")

    # 解析参数
    options = parser.parse_args()

    parent_pid = os.getpid()  # 主进程 PID
    sys.stderr.write("[%s] Pileup genome, processes: %d, pid: %d\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()), options.process, parent_pid))
    t1 = time.time()  # 记录开始时间

    dictRefSeq = {}  # 保存参考序列
    RefBins = []  # 保存分块信息

    # 读取 fasta 序列并拆分 bin
    for genome in options.fasta:
        for seq in SeqIO.parse(genome, "fasta"):
            if seq.id not in dictRefSeq:
                dictRefSeq[str(seq.id)] = seq.seq
                if "GL" in seq.id:  # rRNA 或特殊 contig
                    for start in range(0, len(seq.seq), 100):
                        RefBins.append((seq.id, start, start+100))
                else:  # 普通染色体
                    for start in range(0, len(seq.seq), options.step):
                        RefBins.append((seq.id, start, start+options.step))

    # 捕获 Ctrl+C 信号
    signal.signal(signal.SIGINT, signal_handler)
    multiprocessing.freeze_support()  # Windows 支持多进程

    # 创建进程池
    pool = multiprocessing.Pool(options.process)
    try:
        # 为每个 bin 分配 worker
        for item in RefBins:
            chr, start, stop = item
            pool.apply_async(func=worker, args=(chr, start, stop,))

        pool.close()  # 不再添加新任务
        print("Merging tmp files222...")
        pool.join()  # 等待所有进程完成
        print("Merging tmp files333...")

        # 合并临时文件
        tmp_names = "tmp_" + str(parent_pid) + "*"
        tmp_merge_name = options.output
        print("Merging tmp files...")
        sys.stderr.write("[%s] Merging TEMPs\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        os.system("cat %s > %s" % (tmp_names, tmp_merge_name))
        print("Removing pileup tmp files...")
        os.system('rm %s' % tmp_names)
    finally:
        pool.terminate()  # 确保终止进程

    t2 = time.time()  # 记录结束时间
    print("time spent:", t2-t1)
    sys.stderr.write("[%s] Genome pileup finished.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
```


> 给定一个或多个已比对好的 BAM 文件（即 reads 已经映射到参考基因组上），它会对参考基因组每个位置统计：
>
> * 该位点覆盖的所有 reads；
> * 每条 reads 在该位点的碱基（A/T/C/G）；
> * 该碱基的测序质量；
> * 区分正链（+）和反链（-）；
> * 统计每条 reads 内部 A 的数量（或者 T，如果是反向互补）；
> * 把结果输出成一个文本 pileup 文件。

它相当于是 `samtools mpileup` 的一个 Python 并行实现，但输出格式更为定制化，支持 trim、质量过滤、冲突解决、以及多进程。

---
通过 argparse 设置运行参数：

| 参数                           | 含义                               | 默认值                       |
| ---------------------------- | -------------------------------- | ------------------------- |
| `-i`                         | 输入 BAM 文件（必须是 sorted 并且带有 index） | 无                         |
| `-o`                         | 输出文件名                            | `genome_merge.pileup.txt` |
| `-f`                         | 输入 fasta 序列，可多文件                 | 无                         |
| `-P`                         | 进程数                              | 1                         |
| `-q`                         | 最低碱基质量                           | 10                        |
| `-s`                         | 每个进程处理的分块大小（bin 步长）              | 100,000                   |
| `-m`                         | 最大 pileup 深度                     | 10,000,000                |
| `-M`                         | rRNA 区域 pileup 深度上限              | 10,000                    |
| `--omit-confilct`            | 当同一 read 不同位置出现冲突时是否丢弃           | 否                         |
| `--trim-head`, `--trim-tail` | 修剪 read 前后端                      | 0                         |

---

### 2️⃣ 读取参考序列

```python
for genome in options.fasta:
	for seq in SeqIO.parse(genome, "fasta"):
		dictRefSeq[str(seq.id)] = seq.seq
		if "GL" in seq.id:
			for start in range(0, len(seq.seq), 100):
				RefBins.append((seq.id, start, start + 100))
		else:
			for start in range(0, len(seq.seq), options.step):
				RefBins.append((seq.id, start, start + options.step))
```

* `dictRefSeq`：保存每条染色体或 contig 的完整序列。
* `RefBins`：把每条染色体分块（bin），每块由 `(contig, start, stop)` 表示。
  比如人类染色体 1 的 2.5 亿个碱基，就会被拆成大约 2500 个块（如果 step=100,000）。
* `"GL"` 的 contig（rRNA 或短 scaffold）会用更小的 bin（100 bp）。

这些 bin 将被分配给多个进程并行 pileup。

---

### 3️⃣ 多进程分发

主进程创建一个 `multiprocessing.Pool`：

```python
for item in RefBins:
	chr, start, stop = item
	pool.apply_async(func=worker, args=(chr, start, stop,))
```

每个 `worker` 处理一个分块，对该区域做 pileup。

---

### 4️⃣ worker函数核心逻辑

每个进程独立执行 `worker(contig, start, stop)`：

```python
with pysam.AlignmentFile(options.input, 'rb') as samfile:
	for pileupcolumn in samfile.pileup(contig, start, stop, ...):
```

这里使用 `pysam.pileup` 方法从 BAM 读取 pileup 数据。
它会返回一个 `pileupcolumn` 对象，表示参考基因组上一个特定位置（`pileupcolumn.pos`）。

对于每个位置，程序：

1. 取出该位点参考碱基；
2. 对所有覆盖该位点的 reads 进行遍历（`for pileupread in pileupcolumn.pileups:`）；
3. 判断：

   * 该碱基质量是否 ≥ `--qual`
   * 是否位于 read 的头尾（根据 `--trim-head/tail` 过滤）
4. 如果通过过滤，则记录：

   * 碱基（正向用原序列，反向互补）
   * 质量
   * 该 read 中 A（或 T）的数量（作为额外特征）

---

### 5️⃣ 冲突解决策略

脚本中一段很长的判断：

```python
if lastRead != query_base:
	if options.omit:
		# 冲突就删掉
	else:
		# 取质量更高的那个
elif lastRead == query_base:
	# 如果相同碱基但 A_count 不同，则取 A_count 更多的
```

含义是：

* 同一条 read 对同一位置出现多个 pileup 记录（比如 indel 附近），需要决定保留哪一个；
* 默认策略是“取质量更高的”；如果 `--omit-conflict` 启动，则直接丢弃该 read。

---

### 6️⃣ 输出 pileup 格式

每个位点写两行（若有正负链覆盖）：

```
chr   pos   strand   ref_base   read_bases   A_counts
```

例如：

```
chr1	12345	+	A	A,A,A,G	10,12,9,11
chr1	12345	-	T	T,C	8,6
```

其中：

* `+/-` 表示正链或反链；
* `ref_base` 是参考基因组该位点碱基；
* `read_bases` 是所有覆盖 reads 的碱基；
* 最后一列表示每条 reads 内部 A 的数量（一个额外特征，用于 GLORI 的甲基化分析）。

---

### 7️⃣ 临时文件合并与清理

每个进程输出一个临时文件（如 `tmp_1234_5678`），最后主进程执行：

```bash
cat tmp_1234_* > genome_merge.pileup.txt
rm tmp_1234_*
```

合并成最终结果文件，再清理中间文件。

---

## 🧾 输出结果总结

输出文件（`genome_merge.pileup.txt`）的每一行结构如下：

| 列名       | 内容                          | 示例        |
| -------- | --------------------------- | --------- |
| 染色体      | pileupcolumn.reference_name | chr1      |
| 位置       | pileupcolumn.pos + 1        | 12345     |
| 链方向      | + / -                       | +         |
| 参考碱基     | ref_base                    | A         |
| 碱基列表     | 所有 reads 在该位点的碱基            | A,A,G,C   |
| A_counts | 每条 read 内 A 的数量             | 10,12,9,8 |

---

> 这个脚本是 GLORI 流程里生成碱基统计矩阵的关键环节，它用 `pysam.pileup` 并行扫描 BAM，输出每个位点正负链上的 reads 覆盖碱基和 A 含量，为后续 RNA 修饰检测（比如 m6A/m1A）提供基础数据。


1. **读取参考基因组序列**

   * 脚本先读取一个或多个 fasta 文件，把每条染色体或 scaffold 的序列存储在字典中。
   * 然后把每条染色体拆分成多个小段（bin），rRNA 或特殊 scaffold 用小 bin（100 bp），普通染色体用默认 bin（100,000 bp）。
   * **为什么要拆分？**

     * 因为整个染色体可能非常大，如果一次性处理会占用大量内存并且处理缓慢。
     * 拆成 bin 可以让每个 bin 分配给不同的进程，**并行处理，提高速度**。

2. **设置多进程池处理每个 bin**

   * 使用 Python 的 `multiprocessing.Pool`，每个 bin 分配一个 worker 进程。
   * 每个 worker 都独立处理自己负责的区域，生成临时文件（文件名带进程 PID，防止冲突）。
   * **为什么要多进程？**

     * BAM 文件可能包含亿级 reads，顺序处理非常慢。
     * 多进程利用 CPU 多核能力，同时处理不同区域，大幅缩短运行时间。

3. **对每个位点做 pileup**

   * 每个 worker 调用 `pysam.AlignmentFile.pileup` 遍历指定区间的每个位点。
   * 对每个位点，它会获取：

     * 当前参考碱基
     * 覆盖该位点的所有 reads
     * 每条 read 在该位置的碱基和质量
   * **为什么要在每个位点做 pileup？**

     * pileup 可以统计每个位点被多少 reads 覆盖、各碱基的分布，这是 RNA 修饰分析和 SNP 检测的基础。

4. **过滤低质量和端点碱基**

   * 每条 read 的碱基只有在质量大于等于设定阈值、并且不在 read 头尾被 trim 的情况下才被统计。
   * **为什么要做这一步？**

     * 测序末端往往错误率高，直接统计可能产生假阳性。
     * trim 也可以去掉 adapter 或化学处理引入的偏差。

5. **区分正链和负链**

   * 对正链统计 A 的数量，对负链统计 T 的数量，并且对负链取互补碱基。
   * **为什么区分链？**

     * RNA 修饰在不同链上可能影响逆转录效率不同。
     * 正负链分开统计，有助于判断修饰信号来源。

6. **处理冲突碱基**

   * 如果一条 read 对同一位置出现多个不同碱基，脚本可以选择：

     * `omit=True`：丢弃该 read
     * `omit=False`：取质量最高的碱基
   * **为什么要这样做？**

     * 高覆盖数据中，有些 read 可能在不同位置出现碱基冲突或测序错误，直接统计会产生噪声。
     * 通过这个策略可以提高 pileup 的准确性。

7. **记录 read 内碱基分布（A 或 T 数量）**

   * 统计每条 read 内的 A/T 数量，并与该位置的碱基一起输出。
   * **为什么不仅看 pileup 中位点分布，还要看 read 内分布？**

     * RNA 修饰影响的是逆转录产物在 read 内的碱基分布模式，不仅仅是单个位点。
     * 统计 read 内的 A/T 数量可以过滤噪声，增强修饰信号可靠性。

8. **输出临时文件**

   * 每个 worker 输出自己负责的 bin 的 pileup 数据到临时文件，格式如下：

     ```
     chr  pos  strand  ref_base  read_bases  A_counts/T_counts
     ```
   * 正链和负链分开输出。

9. **合并临时文件生成最终 pileup 文件**

   * 所有 worker 完成后，主进程用 `cat` 命令合并临时文件，并删除临时文件。
   * **为什么要分两步输出再合并？**

     * 避免多进程同时写入同一个文件造成冲突。
     * 保证输出的完整性和顺序。

10. **记录耗时和完成提示**

    * 计算整个 pileup 过程耗时，输出日志信息。

---

###  总结

这段 Python 脚本**不仅仅是统计每个位点碱基分布**，而是做了**深度优化和定制化处理**：

1. **多进程加速** → 处理大 BAM 文件。
2. **分段处理（bin）** → 控制内存，方便并行。
3. **质量和端点过滤** → 避免低质量测序错误。
4. **冲突处理策略** → 保证统计碱基可靠。
5. **统计 read 内碱基分布** → 为 RNA 修饰分析提供上下文信息。
6. **正负链分开** → 增强对链特异性修饰信号的检测。

> 简单一句话理解：
>
> **它是在从 BAM 文件中提取最可靠、最详细的碱基覆盖信息，并保留每条 read 内部的碱基分布，为 RNA 修饰分析提供高质量输入。**

---


好的，我来给你做一个**事无巨细、带完整注释和逻辑解释**的讲解，帮你理解这段针对 mpileup 文件的后续处理脚本。

---

## 🧬 代码作用概览

这个 Python 脚本的作用是：

> 对之前生成的 mpileup 文件进行处理，**将每个位点的参考碱基替换成原始参考基因组的碱基**，只处理那些覆盖碱基中存在多个碱基的位点（即真正有信息的位点）。

这个步骤通常是 GLORI 或其他 RNA 修饰检测流程中的 **后处理步骤**，目的是生成一个“带有参考碱基的标准化 pileup 文件”，方便下游分析。

---
## 核心函数 FilterAll

```python
def FilterAll(referfa, seqment, OUTPUT):
```

* 功能：处理 mpileup 文件中的单行数据。
* 参数：

  * `referfa`：pysam 读取的参考 fasta 文件对象。
  * `seqment`：mpileup 文件的一行字符串。
  * `OUTPUT`：输出文件对象，用于写入处理后的行。

---
### 3.1 分割行内容

```python
line = seqment.strip().split("\t")
```

* 去掉换行符后，用 `\t` 分割 mpileup 的字段。
* mpileup 格式通常是：

  ```
  chr  pos  ref_base  read_bases  quality  ...
  ```
* **目的**：方便后续修改参考碱基字段（line[3]）。

---

### 3.2 判断是否需要处理

```python
if len(line[5].split(",")) > 1:
```

* `line[5]` 是 pileup 中存储的覆盖碱基或其他信息（可能是 read 内碱基列表）。
* `line[5].split(",")` 将覆盖碱基按逗号分开。
* 条件 `>1`：

  * 只处理 **多条覆盖 read 的位点**。
  * **为什么要这样做？**

    * 单 read 位点信息不可靠，噪声大。
    * 多 read 覆盖的位点才有统计意义，适合后续修饰或变异分析。

---

### 3.3 获取真实参考碱基

```python
site_chr = line[0].split('_AG_converted')[0]
site_loci = int(line[1])
sites_base = referfa.fetch(reference=site_chr, start=site_loci - 1, end=site_loci).upper()
```

* `line[0]`：染色体名，可能带 `_AG_converted` 后缀（因为实验中有化学转换或基因组修改）。
* `.split('_AG_converted')[0]`：恢复原始染色体名，确保可以从原始 fasta 获取碱基。
* `line[1]`：位点位置（1-based）。
* `referfa.fetch(start=end)`：

  * pysam 读取 fasta 时是 **0-based**，所以 `start=site_loci-1`。
  * `end=site_loci`，取单个碱基。
* `.upper()`：确保碱基为大写。

**目的**：

* mpileup 文件可能是经过 AG/T/C 转换后的“实验参考”，原始参考碱基可能丢失。
* 下游分析需要知道该位置在原始基因组中的真实碱基。

---

### 3.4 替换行中的参考碱基

```python
line[3] = sites_base
```

* mpileup 第 4 列（line[3]）是参考碱基字段。
* 将其替换成从原始参考基因组中获取的碱基。

**目的**：

* 统一参考碱基信息，避免因 AG 转换或实验处理导致的错配。
* 下游分析（如 RNA 修饰检测）需要对照原始参考碱基。

---

### 3.5 写入输出文件

```python
OUTPUT.write("\t".join(map(str,line)) + "\n")
```

* 将修改后的行写入输出文件。
* 使用 `map(str, line)` 确保所有字段都是字符串。
* 加上换行符。

---

## 4️⃣ 主程序逻辑

```python
if __name__ == "__main__":
```

* 标准 Python 主入口，确保脚本可作为模块导入时不会执行。

---

### 4.1 参数解析

```python
parser = argparse.ArgumentParser(description="filered mpilup")
group_required = parser.add_argument_group("Required")
group_required.add_argument("-input", "--input", nargs="?", type=str, default='default', help="mpileup file")
group_required.add_argument("-referFa", "--referFa", nargs="?", type=str, default='default', help="reference fa")
group_required.add_argument("-outname_prx", "--outname_prx", nargs="?", type=str, default=sys.stdin, help="outname_prx")
options = parser.parse_args()
```

* `-input`：输入 mpileup 文件。
* `-referFa`：原始参考 fasta 文件。
* `-outname_prx`：输出文件前缀，最终文件名是 `<outname_prx>.referbase.mpi`。

**设计目的**：

* 灵活指定输入输出，便于批量处理或管道连接。

---

### 4.2 打开参考 fasta

```python
referfa = pysam.FastaFile(options.referFa)
```

* pysam 以 **高效随机访问**方式打开 fasta 文件。
* 下游可以直接通过 `fetch` 获取任意染色体、任意位置碱基。
* **为什么不用普通文件读取？**

  * fasta 文件可能很大（几百 MB 到几十 GB）。
  * pysam 支持随机访问，比读取整个序列到内存更高效。

---

### 4.3 打开输入输出文件并处理

```python
with open(options.input,'r') as INPUT, open(options.outname_prx + ".referbase.mpi",'w+') as OUTPUT:
    for segment in INPUT:
        FilterAll(referfa, segment, OUTPUT)
```

* 遍历 mpileup 文件每行，调用 `FilterAll` 处理：

  1. 判断是否覆盖多个碱基
  2. 获取原始参考碱基
  3. 替换行中的参考碱基
  4. 写入输出文件

**设计原因**：

* **按行处理**：减少内存占用，适合大文件。
* **只处理多 read 位点**：过滤低信息量或单 read 噪声位点。
* **生成标准化参考碱基文件**：方便 RNA 修饰检测或下游统计。

---

## 5️⃣ 输出文件结构

* 输出文件 `<outname_prx>.referbase.mpi` 与原 mpileup 文件结构类似，但第 4 列参考碱基是 **原始参考基因组碱基**。

示例：

| chr          | pos  | ref (旧) | ref (替换) | read_bases | ... |
| ------------ | ---- | ------- | -------- | ---------- | --- |
| chr1_AG_conv | 1234 | G       | A        | A,G,G,A    | ... |
| chr2_AG_conv | 5678 | T       | T        | T,T,C,T    | ... |

* `_AG_conv` 后缀被去掉，用原始染色体名访问参考 fasta。
* ref 列被替换为真实参考碱基。

---

## 6️⃣ 为什么要这样做（逻辑总结）

1. **实验序列可能经过 AG/T/C 转换** → 参考碱基在 mpileup 中可能不准确。
2. **下游 RNA 修饰分析需要真实参考** → m6A/m5C 判断需要知道原始 A/C/T/G。
3. **只处理多 read 位点** → 避免单 read 噪声导致假阳性。
4. **按行处理** → 大文件也能高效处理，节省内存。
5. **输出标准化 pileup 文件** → 方便与其他分析流程对接。

> 总结一句话：
>
> 这段脚本的核心目的就是“还原 mpileup 中的参考碱基”，把经过实验转换或化学处理的碱基替换成原始参考基因组碱基，同时只保留有统计意义的位点，为下游 RNA 修饰分析提供可靠输入。

---

这段 Python 脚本是整个 GLORI 或类似 RNA 修饰检测流程中，在 `mpileup` 文件生成之后的一个“修正参考碱基”的后处理步骤。因为在前面的比对或 pileup 阶段，研究者使用的是经过 **AG 转换（A-to-G conversion）** 或其他化学改造的“人工参考基因组”（所谓 ternary genome），这种参考序列的染色体名带有 `_AG_converted` 等后缀，并且参考碱基字段（pileup 文件中第 4 列）也可能已经不是原始基因组的真实碱基了。例如，在检测 RNA 上的 m⁶A 修饰时，实验通过化学反应把 A 转成 G，再比对回修改过的参考基因组，结果导致 mpileup 文件里显示的参考碱基可能是 G，而实际上在原始基因组中该位点的真实参考碱基是 A。为了保证后续修饰识别或突变分析能在“真实的碱基背景”下进行，脚本必须重新读取原始参考基因组（用 `pysam.FastaFile` 打开，可随机访问 fasta 文件的任意位置），并把每个位点的参考碱基重新取出来，替换掉 mpileup 文件里错误或转化过的那一个。

脚本的主程序部分首先解析命令行参数，要求用户输入三个关键文件：① 一个 mpileup 文件（由上一步 pileup 脚本生成，记录每个位点的 reads 碱基分布），② 原始参考基因组的 fasta 文件，③ 输出文件的前缀。然后程序逐行读取 mpileup 文件，对于每一行（代表一个基因组位点），调用 `FilterAll()` 函数进行处理。函数首先按制表符分割这一行的各个字段，并检查第 6 列（`line[5]`）里记录的覆盖碱基信息。如果该列中用逗号分隔后有多个元素，说明这个位点有多个 reads 覆盖，信息可靠；反之如果只有一个 reads 覆盖，则被认为不具统计意义或可能是噪声，不处理也不输出。对于需要处理的位点，脚本会先从染色体名中去掉 `_AG_converted` 后缀（因为原始 fasta 里没有这个后缀），然后取出该行记录的位点位置（mpileup 是 1-based，而 pysam 是 0-based，所以要减 1），用 `referfa.fetch(reference=site_chr, start=site_loci-1, end=site_loci)` 精确读取该位置的单个碱基。取出的碱基转成大写后，被替换进 mpileup 行的第 4 列（`line[3]`），这样这一行的“参考碱基”字段就恢复成原始基因组的真实碱基。最后，程序把修正后的整行重新写入输出文件（`outname_prx.referbase.mpi`）。

整个过程的意义在于把一个受实验处理影响的 mpileup 文件“还原”回一个能反映原始基因组真实碱基信息的标准 pileup 文件。只有经过这样的还原，下游的 RNA 修饰检测算法（例如统计 A/G 转换比例判断 m⁶A 修饰）才能基于真实的参考 A 位点进行分析，而不会被 AG 转换参考序列误导。此外，通过只保留那些有多个 reads 覆盖的位点，这个脚本还在同时进行了噪声过滤，避免单 read 位点带来假阳性。最终输出的 `.referbase.mpi` 文件既保留了每个位点的所有 reads 覆盖信息，又补充了正确的参考碱基，是后续修饰识别、变异统计、或可视化分析的标准输入文件。换句话说，这个脚本的任务是“为 pileup 文件重新找回真实的参考碱基”，确保后续计算建立在真实、可靠的序列基础上。
