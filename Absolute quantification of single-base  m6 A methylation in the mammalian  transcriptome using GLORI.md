### 背景
m^6A是最常见的RNA内源性修饰之一，存在于0.2-0.6%的腺苷残基中。它由METTL3、METTL14等甲基转移酶复合体添加，并可通过FTO和ALKBH5去除。m^6A影响RNA的稳定性、翻译效率、可变剪接和转运等多种生物学功能，进而影响生理和病理过程。

### 现有的m^6A检测方法：

1. 抗体依赖方法（如m^6A-seq、MeRIP-seq等）：能广泛检测m^6A，但分辨率差，且无法提供特定位点的定量信息，还可能因抗体的非特异性结合而降低准确性。

2. 抗体无关方法（如MAZTER-seq、m^6A-REF-seq等）：通过酶切法（如MazF核糖核酸酶）识别m^6A，但这些方法的灵敏度不足，且只能检测ACA序列中的m^6A位点，无法做到全转录组的高分辨率和高定量。

现有方法在分辨率、定量能力和特异性方面仍存在局限性。

GLORI方法：

为了克服现有方法的缺点，研究者开发了GLORI（Glyoxal and nitrite-mediated deamination of unmethylated adenosines）方法，这是一种基于测序的高通量、单碱基分辨率、绝对定量的RNA m^6A检测方法，类似于DNA 5mC的重亚硫酸氢盐测序。

GLORI的工作原理：通过乙二醛（glyoxal）和亚硝酸盐（nitrite）诱导未甲基化的腺苷（A）转化为肌苷（I），而m^6A保持不变，从而区分已甲基化与未甲基化的A。通过这种方法，能够定量每个m^6A位点的甲基化水平。

### 结果
亚硝酸脱氨反应：早在半个世纪前，已有研究报道了亚硝酸对腺苷（A）的脱氨作用，即将A转化为肌苷（I），这一反应被认为可以用于RNA中的m^6A检测 。但亚硝酸不仅对A有效，鸟嘌呤（G）和胞嘧啶（C）也会受到脱氨作用，导致分别转化为黄嘌呤（X）和尿嘧啶（U） 。此外，G的脱氨速率比A更快，这使得该方法在m^6A检测中的应用受到限制。

GLORI方法的创新：

为了解决上述问题，研究者们提出通过 **乙二醛（glyoxal** 来保护G，进而实现特异性地加速A的脱氨反应。具体过程如下：

G保护：

乙二醛能够与G的氨基反应，形成可逆的N1,N2-二羟基鸟嘌呤加合物，这种加合物在后续的脱氨反应中能够稳定存在，避免G过早脱氨 。这一保护作用也被广泛应用于其他测序方法，如iSeq和EndoVIPER-seq 。

A脱氨：

在乙二醛保护G的同时，研究者发现，使用乙二醛保护的RNA进行脱氨反应时，A的脱氨反应效率显著提高，达到约70%【Supplementary Fig. 1b】。而如果不使用乙二醛保护，A的脱氨效率几乎为零（<1%）【Supplementary Fig. 1b】。

进一步优化反应条件（包括亚硝酸钠浓度、反应温度和时间）后，成功实现了几乎100%的A→I转化，G→X转化约为3%，C→U转化约为4%【Fig. 1e】【Supplementary Fig. 2b】。重要的是，m^6A在此反应条件下并不会发生脱氨【Fig. 1e】。

反应机制与催化作用：

乙二醛通过反应生成的半缩醛结构可能作为催化剂，显著提高A→I的脱氨效率【Supplementary Fig. 1d】。这一反应的特异性和效率使得m^6A能够被保留，而未修饰的A则被去除，从而为m^6A的定量提供了一个高效的化学手段。

结果：

在优化条件下，使用**液相色谱-串联质谱（LC–MS/MS）进行定量分析，研究者成功实现了99%**的A→I转化，且G→X和C→U的转化率低于5%【Fig. 1e】【Supplementary Fig. 2b】。这一结果表明，GLORI方法能够高效且特异性地识别并定量RNA中的m^6A修饰。
![[Pasted image 20251011110416.png]]

实验原理：

将mRNA打断为小片段；用乙二醛（glyoxal）+亚硝酸盐（nitrite）脱氨：在文库构建中加入 UMI（Unique Molecular Identifier）；测序时，I被读作G；通过统计每个位点中A的比例（A reads / total reads）计算m⁶A水平。

好的 — 我会**逐句逐句**把你贴的那段 GLORI 生信方法文本翻译成中文，并对每一句**逐句解读（为什么要这样做 / 设计目的与生物信息学原因）**。我**不会**做总结，只按原文句子逐条处理。开始：

---

1. 原句：**Detection and quantification of m6A Filtering of original data of the present invention.**
   翻译：**m6A 的检测与定量 — 原始数据的过滤。**
   解读 / 为什么：这是小节标题，说明接下来要讲两件事：如何从测序数据中检测并量化 m6A，以及在此之前对“原始（raw）测序数据”进行哪些预处理/过滤以保证后续结果可靠。数据预处理是所有高通量测序分析的第一步，目的是去除接头、低质量序列与 PCR 伪重复，从而减少下游假阳性/噪声。

---

2. 原句：**Strand orientation of the original RNA was preserved on the process of library construction and reads R2 yields sequences sense to the original RNA.**
   翻译：**在文库构建过程中保留了原始 RNA 的链定向（strand orientation），并且测序的 R2 读段对应于原始 RNA 的正义（sense）序列。**
   解读 / 为什么：说明该文库是有链特性的（strand-specific），并且 R2 是对应转录本本义链的那端。保留链信息可以区分正义/反义转录本，避免将来自基因反义链的信号误认为来自目标转录本。这一步决定了后续只需使用哪一端（R2）来做可靠检测，从而简化分析并避免双向映射带来的混淆。

---

3. 原句：**Thus, only reads R2 was used in our study.**
   翻译：**因此，在本研究中仅使用了 R2 读段。**
   解读 / 为什么：既然 R2 已经给出与原始 RNA 同向的序列，使用单端（R2）就足够且更安全：避免使用 R1（可能是反向或包含引物/UMI信息）造成比对冲突或重复计数问题；同时降低计算量和复杂性。

---

4. 原句：**Illumina sequencing reads were firstly treated with trim_galore (v.0.6.6) for adapter removal with Cutadapt (v.2.10) and quality trimming.**
   翻译：**首先用 trim_galore（v0.6.6，内部调用 Cutadapt v2.10）对 Illumina 测序读段进行接头去除和质量剪切。**
   解读 / 为什么：接头和低质量尾端会导致错误比对与假阳性碱基计数（例如把接头的碱基算进 A/G 计数），因此必须先剪切掉。trim_galore 是常用的封装工具，自动调用 Cutadapt 去除接头并做质量过滤，保证后续比对和碱基计数的准确性。

---

5. 原句：**Working command are as follows: ‘trim_galore -q 20 --stringency 1 -e 0.3 --length 35’.**
   翻译：**实际运行命令为：`trim_galore -q 20 --stringency 1 -e 0.3 --length 35`。**
   解读 / 为什么（逐参数）：

   * `-q 20`：保留碱基的最小质量阈值为 Phred 20（对应错误率 1%），以去除低质量碱基，减少假阳性碱基替换。
   * `--stringency 1`：接头识别时最小重叠长度设为 1（较宽松），确保少量接头残留也能被识别并去掉。
   * `-e 0.3`：允许的最大接头比对错误率为 0.3（30%），平衡去接头的敏感性与特异性。
   * `--length 35`：剪切后保留的最短序列长度为 35 nt，低于此长度的 read 被丢弃，以保证比对质量和碱基统计可靠性。
     这些设置是为了在尽量保留有用数据的前提下去掉可能引入误差的片段。

---

6. 原句：**We then used Seqkit (v.0.13.2) to deduplicate PCR on the basis of 10-base-pair UMI at the 5′ end of reads R2.**
   翻译：**随后使用 Seqkit（v0.13.2）基于位于 R2 读段 5′ 端的 10bp UMI（唯一分子标识符）来去除 PCR 复制（去重）。**
   解读 / 为什么：实验中加入 UMI 是为了解决 PCR 扩增偏倚：同一分子经 PCR 多次扩增会产生大量相同序列，若不去重会导致覆盖与 A/G 计数虚高。用 UMI 去重可以把 PCR 复制的影响去掉，使每个唯一分子只计一次，从而更准确估算 A/G 比例（A rate）。

---

7. 原句：**Key process parameters are as follows: ‘seqkit rmdup -s’.**
   翻译：**关键运行参数如下：`seqkit rmdup -s`。**
   解读 / 为什么：`seqkit rmdup` 用于基于序列内容或 UID 去重，`-s` 表示按序列（sequence）或指定策略进行去重。这里强调用该命令把带 UMI 的 reads 基于 UMI（以及序列）进行去重，以保证统计使用的是去重后的唯一分子，而非 PCR 产物。*（注：不同工具参数含义可略有差异，但目的均为去除 PCR 假重复。）*

---

8. 原句：**Finally, FASTX-Toolkit (v.0.0.13) was used to remove the 10-base-pair UMI in the deduplication reads whose key parameters are defined as follows: ‘fastx_trimmer -f 11’.**
   翻译：**最后使用 FASTX-Toolkit（v0.0.13）从去重后的读段中移除 10bp 的 UMI，运行参数为：`fastx_trimmer -f 11`。**
   解读 / 为什么：去重操作需要读取 UMI 信息，但进行比对时 UMI 并非真实生物序列的一部分，会导致比对错误或多余的 mismatches。因此去重后需把 UMI 从 read 序列中剔除。`-f 11` 表示从第 11 个碱基开始保留（即去掉前 10bp），保证比对只用真实的转录本序列部分，减少伪匹配。

---

9. 原句：**GLORI-tools and code availability.**
   翻译：**GLORI-tools 与代码可用性。**
   解读 / 为什么：这是小节标题，说明后续会介绍用于自动化分析的工具、功能与获取方式（代码托管位置）；开放代码便于他人复现和修改。

---

10. 原句：**After the read-cleaning steps, analysis of GLORI data mainly includes two major parts: alignment of reads and identification of m6A sites with quantified methylation level.**
    翻译：**在完成读段清洗后，GLORI 数据分析主要包括两大部分：读段比对（alignment）以及基于定量甲基化水平的 m6A 位点识别。**
    解读 / 为什么：明确分析流程结构。先把 reads 比对到参考（获得每个位点的覆盖与碱基组成），再基于 A/G 比例识别并量化 m6A。这是从原始测序数据到生物结论的必要两步。

---

11. 原句：**To make the analysis pipeline easy to implement, we developed GLORI-tools, a computational pipeline that quantifies the A rate of each A site and follows multiple filters for false positives.**
    翻译：**为便于实现分析流程，我们开发了 GLORI-tools —— 一个计算管线，用于量化每个 A 位点的 A rate 并执行多重过滤以去除假阳性。**
    解读 / 为什么：自动化管线能标准化处理流程，保证不同样本、不同研究者能复现结果。多重过滤是必要的，因为化学脱氨并非 100% 完美，存在未完全转化、测序错误、结构阻挡等，会产生假阳性；故需要一系列经验或统计过滤以提高置信度。

---

12. 原句：**GLORI-tools takes cleaned reads as input and finally reports the conversion rate (A-to-G) of GLORI and single-base loci of m6A sites with corresponding A rate as modification level.**
    翻译：**GLORI-tools 以清洗后的 reads 为输入，最终输出 GLORI 的总体转化率（A→G）以及逐个位点的 m6A 单碱基位点列表，并给出每个位点对应的 A rate（作为其修饰水平）。**
    解读 / 为什么：两类输出都重要：总体 A→G 转化率用于 QC（判断化学反应是否成功），而单碱基位点与其 A rate 是最终的生物学产物（哪些位点有 m6A 以及它们的甲基化比例）。这使得研究既有实验层面的 QC，又有精确的位点信息。

---

13. 原句：**The code of the GLORI-tools is available in the following GitHub repositories: [https://github.com/liucongcas/GLORI-tools](https://github.com/liucongcas/GLORI-tools).**
    翻译：**GLORI-tools 的代码托管在以下 GitHub 仓库：`https://github.com/liucongcas/GLORI-tools`。**
    解读 / 为什么：公开源码便于社区审查、复现、优化与二次开发，增加方法透明性和可重复性。

---

14. 原句：**Alignment of reads.**
    翻译：**读段比对。**
    解读 / 为什么：小节标题，接下来详细说明比对策略。GLORI 的比对需要特殊处理（因为大量 A 被转换为 G），所以要单独强调。

---

15. 原句：**As the reads R2 of GLORI were A-to-G-converted reads, GLORI-tools firstly transformed the GLORI reads to the A-to-G version (transformed unconverted A to G) and then aligned them to equivalently preconverted forms of the reference genome (GRCh38) using STAR (v.2.7.5.c) with key parameters as follows: ‘--outFilterMismatchNmax 2 --outSAMprimaryFlag AllBestScore --outFilterMultimapNmax 1’.**
    翻译：**由于 GLORI 的 R2 读段中大量 A 被转化为 G，GLORI-tools 首先把读段处理成“全部 A→G” 的版本（把未转化的 A 也临时替换为 G），然后将其比对到同样事先“预先转换”（preconverted）的参考基因组（GRCh38）的对应版本上，使用 STAR（v2.7.5.c）进行比对，关键参数为：`--outFilterMismatchNmax 2 --outSAMprimaryFlag AllBestScore --outFilterMultimapNmax 1`。**
    解读 / 为什么（分点详解）：

    * **为什么把 reads 转成 A→G 版本并对“预转换”参考比对？** 因为化学处理会使绝大多数原始 A 变为 I（测序上读作 G），这会导致常规比对出现大量 A↔G 不匹配，从而降低比对率或引入错误比对。把 reads 和参考都做相同的“预转换”可在比对阶段把 A→G 替换当作“匹配”，极大提高比对效率与准确度。
    * **为什么用 STAR？** STAR 是对 RNA-seq 处理常用的快速比对工具，支持高效的大规模比对与 splice-aware（跨剪接）映射。
    * **参数含义与目的：**

      * `--outFilterMismatchNmax 2`：允许的最大 mismatch 数为 2，限制过多不匹配以保证比对质量；
      * `--outSAMprimaryFlag AllBestScore`：在多重最佳比对时指定主比对标记策略，保证后续只使用评分最高的比对或以一致方式报告；
      * `--outFilterMultimapNmax 1`：只保留唯一比对（multimap 最大为 1），减少重复或多位点来源造成歧义。
        这些设置旨在在“预转换”策略与严格比对约束下得到高质量、唯一的比对结果，便于后续逐位点碱基统计。

---

16. 原句：**Unmapped and multiple mapped reads were then realigned to preconverted forms of the reverse complementary sequences of the reference genome (GRCh38).**
    翻译：**对于未比对上的以及发生多重比对的读段，接着将其重新比对到参考基因组（GRCh38）反向互补序列的“预转换”版本上。**
    解读 / 为什么：有些 reads 可能来自反义链或由于复杂基因结构在初步比对时未被匹配。把反向互补序列也做“预转换”并重比对，是为了补救那些在正向参考上漏掉的 reads，增加比对回收率并确保链方向问题不会丢失真实信号。

---

17. 原句：**Unmapped and multiple mapped reads were realigned to preconverted forms of the transcriptome (GCF_000001405.39) by bowtie (v.1.3.0) with default parameters.**
    翻译：**未比对或多重比对的读段还会被重新比对到“预转换”的转录组参考（GCF_000001405.39），使用 bowtie（v1.3.0，默认参数）。**
    解读 / 为什么：将 reads 映射到转录组（转录本序列）可以救回跨剪接或 splice junction 导致在基因组层面未比对成功的片段。bowtie 对短序列或转录本比对速度快且稳定，用作补救映射以提高覆盖度和位点检测灵敏度。

---

18. 原句：**Considering a certain proportion of C-to-U in GLORI data, the process of mapping was mismatch tolerance (alignment with ≤2 mismatch).**
    翻译：**鉴于 GLORI 数据中存在一定比例的 C→U（C 被转化为 U）的情况，比对过程允许一定的错配（即允许 ≤2 个 mismatches）。**
    解读 / 为什么：除了 A→I（读作 G）之外，化学处理也会引入少量 C→U、G→X 等非目标转换（副反应）。因此比对时放宽错配限制（但仍很严格，≤2）可以兼容这些真实的化学副产物，避免把真实读段误判为不匹配从而丢失数据；但同时仍限制太多错配以防止误比对。

---

19. 原句：**The unconverted As were then transformed back.**
    翻译：**随后将那些在预转换过程中被临时替换的未转化 A 还原回去（即恢复为原始 A）。**
    解读 / 为什么：预转换只是为了比对方便；在完成比对后必须把序列还原回原始碱基，以便准确统计哪些读段在该位置呈现 A（即可能为 m6A）或 G（即原本未甲基化的 A 经转化）。还原是为了得到精确的碱基计数与 A rate。

---

20. 原句：**Uniquely mapped reads with reconverted As in bam format was then merged together with genome coordinates and sorted with samtools (v.1.10.2) for the detection of m6A.**
    翻译：**具有唯一比对且已恢复原始 A 的 reads（以 BAM 格式）被合并、按基因组坐标排序（使用 samtools v1.10.2），以便后续进行 m6A 检测。**
    解读 / 为什么：把所有高质量、唯一比对的 BAM 文件合并并按坐标排序是做 pileup / 位点统计的前提。排序后的 BAM 能有效地按位置检索读段并统计在每个位点上的碱基组成（A/G 覆盖等），这是计算 A rate 的基础。

---

21. 原句：**Call m6A sites and filter false positives.**
    翻译：**调用 m6A 位点并过滤假阳性。**
    解读 / 为什么：这是本小节标题，表明接下来要介绍如何从比对结果中统计位点 A/G，并通过一系列过滤/统计方法剔除由不完全转换、测序错误或结构阻挡产生的假阳性位点。

---

22. 原句：**The mapping output were used to detect putative m6A sites based on A rate (percentage of A in total coverage).**
    翻译：**利用比对结果，通过计算 A rate（即在总覆盖中 A 的百分比）来检测候选 m6A 位点。**
    解读 / 为什么：这正是 GLORI 的核心量化策略：在化学处理后，未甲基化的 A 多被转为 G，甲基化的 A 仍然读作 A；因此某个位点中读到 A 的比例直接反映了该位点上 m6A 的甲基化比例，称为 A rate。

---

23. 原句：**To check all A positions throughout the genome, pileup_genome_multiprocessing command of GLORI-tools was used to convert the mapping output from bam format to mpileup format, and then count the number of unconverted As in each mapped read.**
    翻译：**为了遍历基因组上所有的 A 位点，使用 GLORI-tools 的 `pileup_genome_multiprocessing` 命令将比对输出（BAM）转换为 mpileup 格式，然后统计每个比对到该位点的 read 中未被转化的 A 数目。**
    解读 / 为什么：mpileup（pileup）是一种按基因组位置汇总所有对该位置覆盖的读段及其碱基的常用格式。用并行化命令（multiprocessing）可以高效地在全基因组范围内统计每个 A 位点上读到 A 的次数与 G 的次数，为后续 A rate 计算提供基础数据。

---

24. 原句：**Next, the m6A_pileup_formatter command could format the mpileup file into a specific FMAT format file.**
    翻译：**接着，使用 `m6A_pileup_formatter` 命令将 mpileup 文件格式化为特定的 FMAT 格式文件。**
    解读 / 为什么：将原始 mpileup 数据转换为内置的 FMAT（格式化）文件方便后续统一处理、索引、与注释信息结合，也便于导出表格结果或用于后续的过滤与统计步骤。

---

25. 原句：**Combined with gene annotation information, this step could mark the genomic coordinates of each A site and the corresponding gene, coverage of total reads, and A rate as methylation level, simultaneously.**
    翻译：**结合基因注释信息，该步骤能够同时标注每个 A 位点的基因组坐标与所属基因、总覆盖度，以及将 A rate 记录为其甲基化水平。**
    解读 / 为什么：把位点与基因注释连接可以实现按基因或转录本下游分析（比如基因层面的甲基化汇总、位点在 3′UTR/CDS/5′UTR 的分布等）。同时输出 coverage 与 A rate 便于后续过滤（低覆盖剔除）与生物学解释。

---

26. 原句：**FMAT format file was used for subsequent analysis.**
    翻译：**FMAT 格式文件被用于后续分析。**
    解读 / 为什么：这是流程说明：FMAT 是后续过滤、统计与显著性检测的标准中间文件格式，有利于管线模块化。

---

27. 原句：**To further eliminate the influence of the incomplete conversions and false positives in nitrite-treatment-resistant regions, a series of filters was introduced in m6A_caller: (1) filtered mapped reads containing >3 unconverted As; (2) filtered A sites with low coverage and A rate (coverage ＜15 and A rate ≤ 0.1); and (3) filtered A sites with high ‘signal-to-noise ratio’ (≥20%).**
    翻译：**为了进一步消除不完全转化以及对亚硝酸处理耐受的区域带来的假阳性影响，m6A_caller 中引入了一系列过滤规则：
    (1) 过滤包含 >3 个未转化 A 的映射 read；
    (2) 过滤覆盖度和 A rate 较低的 A 位点（coverage < 15 且 A rate ≤ 0.1）；
    (3) 过滤具有高“信号比（signal-to-noise ratio）”的 A 位点（≥ 20%）。**
    解读 / 为什么（逐项）：

    * (1) **过滤包含 >3 个未转化 A 的 reads**：一条 read 上如果有多个（>3）未转化的 A，说明该条 read 可能来自化学处理失败的分子或来自强二级结构/修饰导致的脱氨阻断，这会带来系统性假阳性，将这些 read 去除可以减少“多数未转化 A”造成的位置偏倚。
    * (2) **过滤低覆盖且低 A rate 的位点（coverage<15 且 A rate ≤0.1）**：低覆盖位点本身统计不稳定，且 A rate≤0.1 接近背景噪声门限（研究中用 0.1 作为初步阈值），这样的位点证据不足，容易是测序或随机误差导致的伪信号；剔除它们可以提高整体结果置信度。
    * (3) **过滤具有高“信号比”的位点（≥20%）**：文中将“signal ratio”定义为“reads 中有 ≥2 个未转化 A 的数量 / 总覆盖”，如果该比值很高（≥20%），说明该位点所在区域存在大量携带多个未转化 A 的 reads，可能是“亚硝酸处理耐受区域”或结构复杂区，容易产生系统性偏差。因此把这类位点过滤掉可以减少由区域特性导致的伪阳性。

---

28. 原句：**The signal ratio is defined as the signal (the number of reads with ≥2 unconverted As) divided by the total covered reads48.**
    翻译：**所谓的“信号比”定义为：携带 ≥2 个未转化 A 的 reads 数（即信号）除以该位点的总覆盖 reads 数。**
    解读 / 为什么：这是对第 (3) 条过滤的量化定义。衡量一个位点是否位于“耐处理区域”或受局部结构/修饰影响的关键指标：若很多 reads 在该位点附近都显示多个未转化 A（≥2），这表明化学反应在该区域普遍失败，产生的 A rate 不能被信赖，故需过滤。

---

29. 原句：**The candidate sites were firstly selected by m6A_caller and should be supported by ≥5 variant nucleotides, ≥15 coverage of A + G and ≥0.1 A rate.**
    翻译：**m6A_caller 首先选出候选位点，候选位点须满足：由 ≥5 个变异核苷（variant nucleotides）支持、A+G 覆盖 ≥15 次，以及 A rate ≥0.1。**
    解读 / 为什么：这是候选位点的硬阈值：

    * **≥5 个“变异核苷”支持**（原文措辞略含糊）：通常应理解为至少有 ≥5 个读段显示与参考不同的碱基（即有变体证据），或至少 5 次支持信号（例如 5 个读段显示 A 而不是 G），以避免单个读段的偶发错误造成假阳性。
    * **A+G 覆盖 ≥15**：要求足够的覆盖度，以保证 A rate 的估计统计稳定。
    * **A rate ≥0.1**：初步筛除低于 10% 的位点（接近背景噪声），把灵敏性与特异性平衡。
      这些阈值组合保证候选位点在测序证据上有基本的可靠性。

---

30. 原句：**To call m6A sites from these candidate sites, we next applied statistical tests, which were based on binomial distribution or Poisson distribution, to each of them on the basis of the gene-specific A-to-G conversion rate.**
    翻译：**为了从这些候选位点中最终判定 m6A 位点，接着对每个候选位点进行统计检验（基于二项分布或泊松分布），检验参照量是对应基因的 A→G 转化率（gene-specific conversion rate）。**
    解读 / 为什么：候选位点的 A rate 需要与该基因背景（全基因或染色体层面的 A→G 转化基线）比较，以判断该位点的 A 富集是否显著高于背景（即不是整体转换效率低或高导致的假象）。使用二项或泊松模型可在统计上评估读到当前 A 个数的罕见性（P 值），从而为是否是真正 m6A 提供显著性证据。基于基因特异的转换率能控制基因内部化学效率差异。

---

31. 原句：**The calculated significance P value are required to be <0.005.**
    翻译：**计算得到的显著性 P 值要求小于 0.005。**
    解读 / 为什么：设定较严格的单点显著性阈值（0.005 而非常见的 0.05）以减少假阳性。在全基因组多重检验的背景下，单点阈值更严格可以保留高置信度位点。

---

32. 原句：**Finally, sites with false-discovery rate adjusted P < 0.005 were determined as the finally detected high-confidence m6A sites and their corresponding A rates are determined as m6A modification level.**
    翻译：**最后，对 P 值进行假发现率（FDR）校正后，只有 FDR 调整后的 P < 0.005 的位点被确定为最终的高置信度 m6A 位点；这些位点各自的 A rate 即被用作其 m6A 的修饰水平。**
    解读 / 为什么：通过 FDR 校正控制多重比较带来的假阳性率（即在大量位点检验中仍能保证预计的假发现比例）。把最终阈值设为 FDR < 0.005 非常严格，目的是只保留高度可信的 m6A 位点，并用 A rate 直接报告其甲基化比例，形成最终定量结果。

---

33. 原句：**All these analyses could be carried out by the m6A_caller command and the default parameter ‘-c 15 -C 5 -r 0.1 -s 0.8 ’ is set to obtain the candidate m6A sites.**
    翻译：**以上所有分析可以通过 `m6A_caller` 命令执行，默认参数为 `-c 15 -C 5 -r 0.1 -s 0.8`，用于获取候选的 m6A 位点。**
    解读 / 为什么：把复杂的过滤与统计逻辑封装为单个命令便于重复执行与参数调整。默认参数对应的含义（结合前文）大致为：`-c 15` → 最小 coverage 15；`-C 5` → 最小支持 reads（或“变异碱基数”）5；`-r 0.1` → 最小 A rate 0.1；`-s 0.8` → 这里文中未直接说明 `-s` 对应哪个阈值（前文 signal ratio 的阈值为 20%），但可理解为一个可调的信号相关阈值参数。工具暴露这些参数供用户按需放宽或收紧过滤标准，以适配不同的实验条件或研究目的。

---

34. 原句：**Removal of m1A from the detected sites.**
    翻译：**去除 m1A 对检测位点的干扰。**
    解读 / 为什么：m1A（N1-甲基腺苷）也是一种 RNA 修饰，会干扰逆转录并产生特有的信号（比如截断）。需要额外步骤识别并剔除由于 m1A 引起的假阳性，以免误把 m1A 位置当作 m6A。

---

35. 原句：**We calculated the sequence-truncation signals (calculated as stop ratio) at the detected sites; sites with a >5% stop ratio was removed from the final list.**
    翻译：**我们计算了检测到的位点处的序列截断信号（以 stop ratio 表示）；stop ratio 大于 5% 的位点会从最终列表中移除。**
    解读 / 为什么：m1A 常导致逆转录在修饰位点发生停转（reverse transcription stops），从而在该位置观测到较高比例的 read 中断（截断）。如果某个位点的 stop ratio（被截断的 reads / 总覆盖）超过 5%，很可能该位点受 m1A 或其它导致终止的修饰/结构影响，从而不可靠或不是 m6A，故排除以避免混淆。

---

36. 原句：**Calculated A-to-G conversion rate.**
    翻译：**A→G 转化率的计算。**
    解读 / 为什么：小节标题，接下来说明如何在染色体/基因层面估算转换效率（QC 指标），用于评估化学处理在不同区域的总体效果。

---

37. 原句：**The m6A_pileup_formatter command could report chromosome-specific and gene-specific A-to-G conversion rates, that is, the median A-to-G conversion rate of all A loci with covered reads on each chromosome or gene.**
    翻译：**`m6A_pileup_formatter` 命令可以报告按染色体和按基因的 A→G 转化率（即在每条染色体或每个基因上所有有覆盖的 A 位点的 A→G 转化率的中位数）。**
    解读 / 为什么：按染色体/基因汇总的转化率是重要 QC 指标，可以用来检查化学反应是否在不同染色体或不同基因上均一（若某些染色体/基因转化率异常，可能是文库制备、富集或结构导致的偏差）。采用中位数可以减少极端位点对整体估计的影响。

---

38. 原句：**All of these operations are wrapped in GLORI-tools. The cutoff of all filters can be customised by the user.**
    翻译：**所有这些操作都封装在 GLORI-tools 中，并且所有过滤阈值均可由用户自定义。**
    解读 / 为什么：提供可定制参数使方法具有灵活性：不同实验/细胞类型/测序深度下可能需要不同的阈值来平衡灵敏度与特异性，开放参数便于下游用户根据数据特性调整策略。

---

39. 原句：**Calculation of methylation level for gene**
    翻译：**基因层面的甲基化水平计算。**
    解读 / 为什么：除了逐位点报告，研究者常希望获得基因级别（或转录本级别）的甲基化概览，因此需要把位点-level 的信息汇总成基因-level 的度量。

---

40. 原句：**The number of m6A reads (reads containing unconverted A that can support the locus) for each detected m6A site was calculated as A rate × coverage.**
    翻译：**对于每个检测到的 m6A 位点，其“m6A reads”（即含有未转化 A、能够支持该位点的 reads 数）被计算为：**`A rate × coverage`**。**
    解读 / 为什么：A rate × coverage 给出该位点中预估的绝对支持读数（例如 A rate 0.3 × coverage 100 → 30 个支持的 m6A reads），这是从比例回到绝对支持数的简单方法，有利于后续汇总与统计计算（比如基因加权）。

---

41. 原句：**The methylation level for one gene was then calculated as: the sum of m6A reads of all m6A sites divided by the sum of coverage of all sites in the longest isoform.**
    翻译：**某一基因的甲基化水平计算为：该基因（取最长转录本）的所有 m6A 位点的 m6A reads 之和除以该转录本所有位点的总覆盖（coverage）之和。**
    解读 / 为什么：这是按基因（或按最长 isoform）计算总体甲基化负载的一种归一化方法：分子为所有位点上的“支持 m6A 的绝对读数之和”，分母为相应位点的总测序覆盖，得到一个介于 0–1 的比例，代表该基因在测序覆盖范围内整体被 m6A 修饰的“平均比例”。选择最长 isoform 作为参考可以统一位点集合（避免不同转录本长度差异带来的比较偏差），这种加权方法兼顾了位点的读数支持度与测序覆盖，从而比简单的位点数统计更能反映真实的转录本修饰水平。

---

如果你希望，我可以按你给出的这些步骤把**GLORI-tools 的关键命令行调用顺序**（trim → dedup → trim UMI → convert/read-preconvert → STAR 比对 → bowtie 补救 → pileup → m6A_caller）写成一个**可直接执行的 shell 脚本示例**，并注释每一步的参数含义与可修改项。要我帮你生成吗？
