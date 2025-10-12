**GLORI（Glyoxal and Nitrite-mediated deamination for quantitative m⁶A mapping）
##  GLORI 全流程概览

GLORI 的分析流程可以分为四个主要阶段：

| 阶段  | 名称                                          | 主要内容                      |
| --- | ------------------------------------------- | ------------------------- |
| ①   | **原始测序数据预处理（Preprocessing）**                | 去接头、质量控制、去重、去掉UMI等        |
| ②   | **基因组与注释文件准备（Annotation & Index Building）** | 下载、统一染色体名、构建转换后的基因组/转录组索引 |
| ③   | **m⁶A位点定量分析（Read Alignment & m⁶A Calling）** | 比对序列并计算每个位点的A-to-G转换率     |
| ④   | **对照样本分析（Untreated Sample Alignment）**      | 未处理样本的比对，用于基因表达分析或背景比较    |

---

##  ① 原始测序数据预处理

###  选择不同类型的文库：

GLORI 支持两种文库类型的预处理方式：

| 文库类型                   | 选项           | 特点                             | 主要差别            |
| ---------------------- | ------------ | ------------------------------ | --------------- |
| eCLIP library          | **Option A** | 带有 **5′端UMI标签**，需要去除UMI并去重     | 先去接头、再去重、再去掉UMI |
| NEB #E7300 kit library | **Option B** | 标准文库，无UMI，只需去接头、质量控制、剪除3′端1个碱基 | 简单清洗流程          |
https://github.com/FelixKrueger/TrimGalore
https://journal.embnet.org/index.php/embnetjournal/article/view/200/458
---

### **(A) eCLIP 文库预处理步骤**

| 步骤             | 命令                                                                                                                                                | 工具                                      | 功能说明               |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------- | ------------------ |
| 1️⃣ 去接头 & 质量控制 | `$trim_galore-q20--stringency1-e0.3--length{length(5'UMI)+25nt}\<br>--path_to_cutadapt {cutadapter} --dont_gzip -o {output_dir1} <br>{inputfile}` | **Trim Galore (v0.6.6)** + **Cutadapt** | 去除低质量碱基（Q<20）与接头序列 |
| 2️⃣ 去重复        | `$seqkit rmdup -j {Thread} -s -D {output_dupname} {filter_file1} \> {filter_file2}`                                                               | **SeqKit**                              | 基于序列和UMI去除PCR重复    |
| 3️⃣ 去掉5′端UMI   | `$fastx_trimmer -Q 33 -f {length(5'UMI) +1nt} -i {filter_file2} \-o {filter_file3}`                                                               | **FASTX toolkit**                       | 删除UMI部分，保留真实序列     |
![[Pasted image 20251012143554.png]]

----
### Option A  解释
| 参数                                | 含义                                                        | 在 GLORI 中的作用                                                |
| --------------------------------- | --------------------------------------------------------- | ----------------------------------------------------------- |
| `trim_galore`                     | 主程序，基于 **Cutadapt** 和 **FastQC** 的自动化封装，用于清洗测序reads。      | 自动完成去接头 + 质量剪切。                                             |
| `-q 20`                           | **质量阈值**。修剪掉 Phred 质量值 < 20 的碱基（错误率 ≥1%）。                 | 去除低质量碱基，减少假阳性 A/G 转化。                                       |
| `--stringency 1`                  | **接头匹配严格度**。要求至少1个碱基重叠即可认为是adapter。                       | 放宽检测，避免GLORI短reads中接头无法识别。                                  |
| `-e 0.3`                          | **最大错配率** (error rate ≤0.3)。                              | 要求至少1个碱基重叠即可认为是adapter。                                     |
| `--length {length(5'UMI)+25nt}`   | **最短保留长度**。小于此长度的reads会被丢弃。                               | 确保去除UMI后仍保留 ≥25nt 真正序列用于比对。<br>例如：UMI=10bp → `--length 35`。 |
| `--path_to_cutadapt {cutadapter}` | 指定 **Cutadapt 程序路径**。                                     | 告诉 Trim Galore 使用哪个 cutadapt 程序执行核心剪切。                      |
| `--dont_gzip`                     | 不压缩输出文件（默认是 `.fq.gz`）。                                    | 保持输出为 `.fastq` 格式，方便后续 UMI 处理脚本直接读取。                        |
| `-o {output_dir1}`                | 输出目录。                                                     | 指定结果输出文件夹。                                                  |
| `{inputfile}`                     | 输入文件。                                                     | 原始测序reads（如 `.fastq.gz`）。                                   |
| 参数                                | 作用                                                        |                                                             |
| `-j {Thread}`                     | 指定使用的线程数，加快处理速度。`{Thread}` 是你希望使用的 CPU 核心数                |                                                             |
| `-s`                              | 按序列本身去重（sequence-based deduplication），不考虑序列名称             | 如果不加 -s，默认可能会 按照序列内容 + 序列名称一起去重也就是说，即使序列相同，但名字不同，也会被保留      |
| `-D {output_dupname}`             | 将重复序列写入 `{output_dupname}` 文件，方便检查或统计                     |                                                             |
| `{filter_file1}`                  | 输入文件，通常是去重前的 FASTA/FASTQ 文件                               |                                                             |
| `{filter_file2}`                  | 输出文件，去重后的序列写入这个文件                                         |                                                             |
| 参数                                | 作用                                                        |                                                             |
| `-Q 33`                           | 指定 FASTQ 文件的质量编码为 **Phred+33**（ASCII 33 对应质量值 0）          | 在官方文档中，该命令并没有这个参数                                           |
| `-f {length(5'UMI)+1nt}`          | 从 reads **第 {length(5'UMI)+1} 个碱基开始保留**（即剪掉前端 UMI + 1 nt） |                                                             |
| `-i {filter_file2}`               | 输入文件，通常是去重后的 FASTQ 文件                                     |                                                             |
| `-o {filter_file3}`               | 输出文件，保存剪切后的序列                                             |                                                             |
![[Pasted image 20251012145356.png]]
### **(B) NEB #E7300 文库预处理步骤**

| 步骤             | 命令                                                                                                                                          | 工具                         | 功能说明            |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------- | --------------- |
| 1️⃣ 去接头 & 质量控制 | `$trim_galore -q 20 -a AGATCGGAAGAGCA --stringency 1 --length 20 \--path_to_cutadapt {cutadapter} --dont_gzip -o {output_dir1} {inputfile}` | **Trim Galore + Cutadapt** | 去除adapter和低质量片段 |
| 2️⃣ 去除3′端1个碱基  | `$fastx_trimmer -Q 33 -t 1 -i {filter_file1} -o {filter_file2}`                                                                             | **FASTX toolkit**          | 去掉3′端可能的测序偏差    |
-t 指 去掉3‘端的一个碱基

NEB 文库构建中有 末端修复和加 A，容易在测序读段末端出现人工碱基。
eCLIP 文库构建则采用 直接连接接头，RNA 片段末端保持原始状态，不会多出碱基。



---

##  ② 准备注释与基因组索引

这一步是 GLORI 的基础。你需要准备 **基因组序列（.fa）**、**注释文件（.gtf）** 并构建 **A→G转换版本** 的索引文件。

---

### **Step 124. 下载并统一注释文件**

| 步骤          | 工具                               | 功能说明                              |
| ----------- | -------------------------------- | --------------------------------- |
| 下载 GTF 注释文件 | `wget`                           | 从 NCBI RefSeq 下载人类基因组注释（GRCh38）   |
| 统一染色体名      | `change_UCSCgtf.py`（GLORI-tools） | 统一基因组和注释文件的染色体命名（例如将“chr1”↔“1”统一） |

---

### **Step 125. 构建基因组索引**

| 步骤                    | 工具                                   | 功能说明                                   |
| --------------------- | ------------------------------------ | -------------------------------------- |
| 下载基因组                 | `wget`                               | 从 UCSC 下载 `hg38.fa`                    |
| 构建A→G转换版基因组 + 反向互补基因组 | `build_genome_index.py`（GLORI-tools） | 创建两个版本：forward & reverse strand A→G转换版 |
| 构建 STAR 索引            | **STAR (v2.7.5c)**                   | 为两个基因组版本构建索引用于后续比对                     |

---

### **Step 126. 构建转录组索引**

| 步骤                 | 工具                                                  | 功能说明                                    |
| ------------------ | --------------------------------------------------- | --------------------------------------- |
| 下载转录组文件 (.fna)     | `wget`                                              | 从 NCBI RefSeq 下载                        |
| 选取最长转录本            | `gtf2anno.py` + `selected_longest_transcrpts_fa.py` | 保留每个基因的最长isoform，构建最终转录组                |
| 构建 A→G 转换版转录组 & 索引 | `build_transcriptome_index.py`                      | 创建 A→G 转换的转录组并构建 **Bowtie (v1.3.0)** 索引 |

---

### **Step 127 (可选). 构建单碱基分辨率注释**

| 步骤       | 工具                                       | 功能说明                              |
| -------- | ---------------------------------------- | --------------------------------- |
| 基因组碱基级注释 | `anno_to_base.py`                        | 把 GTF 注释转化为单碱基分辨率注释               |
| 基因列表文件   | `gtf2genelist.py`                        | 生成包含基因注释的列表文件                     |
| 去除冗余注释   | `anno_to_base_remove_redundance_v1.0.py` | 去掉重复注释，得到非冗余版本（用于gene-specific分析） |

---

## 🧠 ③ 读取比对与 m⁶A 位点定量

主程序：**`run_GLORI.py`**
GLORI 将比对 + m⁶A位点定量整合为一步。

有两种模式：

| 模式                    | 背景校正方法              | 使用的文件                  | 推荐用途          |
| --------------------- | ------------------- | ---------------------- | ------------- |
| **(A)** 基因特异性 A→G 转换率 | 基于单基因转换率文件（step127） | `{.noredundance.base}` | 精准分析单基因位点     |
| **(B)** 转录组整体转换率      | 基于全转录组平均转换率         | 无需单基因文件                | 常规m⁶A calling |

---

### **(A) 基因特异性分析命令**

```bash
python run_GLORI.py -i {GLORI-tools} -q {clean_reads} -T {threads} \
-f {hg38}.AG_conversion.fa -f2 {hg38.fa} \
-rvs {hg38}.rvsCom.fa -Tf {Transcriptome_Index} \
-a {GRCh38_prefix}_genomic.gtf_change2Ens.tbl2 \
-b {GRCh38_prefix}_genomic.gtf_change2Ens.tbl2.noredundance.base \
-pre {prefix} -o {output_dir} --combine --rvs_fac
```

### **(B) 转录组整体分析命令**

```bash
python run_GLORI.py -i {GLORI-tools} -q {clean_reads} -T {threads} \
-f {hg38}.AG_conversion.fa -f2 {hg38.fa} \
-rvs {hg38}.rvsCom.fa -Tf {Transcriptome_Index} \
-a {GRCh38_prefix}_genomic.gtf_change2Ens.tbl2 \
-pre {prefix} -o {output_dir} --combine --rvs_fac
```

---

## 🧩 ④ 对照样本（未GLORI处理）比对分析

这一步用于**背景或表达量分析**，不会调用m⁶A位点：

```bash
python run_GLORI.py -i {GLORI-tools} -q {file} -T {threads} \
-f {genome.fa} -Tf {Transcriptome_Index} \
-a {GRCh38_prefix}_genomic.gtf_change2Ens.tbl2 \
-pre {prefix} -o {outputdir} --combine --untreated
```

输出 `.bam` 文件可用于表达分析或与GLORI样本比较。

---

## 🧩 总结：GLORI分析完整流程表

| 阶段                | 主要命令 / 脚本                                                                | 工具                            | 主要功能                     |
| ----------------- | ------------------------------------------------------------------------ | ----------------------------- | ------------------------ |
| **① 数据预处理**       | trim_galore / seqkit / fastx_trimmer                                     | Trim Galore, SeqKit, FASTX    | 去adapter、质量控制、去重复、去UMI   |
| **② 下载注释与构建索引**   | change_UCSCgtf.py / build_genome_index.py / build_transcriptome_index.py | GLORI-tools, STAR, Bowtie     | 统一注释、构建A→G转化索引           |
| **③ 构建单碱基注释（可选）** | anno_to_base.py / anno_to_base_remove_redundance_v1.0.py                 | GLORI-tools                   | 生成单碱基注释用于gene-specific分析 |
| **④ m⁶A位点定量分析**   | run_GLORI.py                                                             | GLORI-tools (整合STAR + Bowtie) | 读取比对 + m⁶A位点识别与定量        |
| **⑤ 对照样本分析**      | run_GLORI.py --untreated                                                 | GLORI-tools                   | 比对未处理样本以供表达分析            |

