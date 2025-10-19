为 GLORI 生成专门的“转换参考序列”并建立索引，以便后续把测序得到的、经过 glyoxal+nitrite 脱氨处理的 RNA 读段准确比对回参考、从而识别 m6A 位点.

具体地说，GLORI 会将未甲基化的 A 脱氨为 I（在反转录/测序时读作 G），而 m6A 位点抵抗脱氨仍读作 A，所以生信上要用把参考中 A 统一“转换”为 G 的参考（即 A→G 的 ternary reference）来实现对比：

未甲基化位点的读段（测为 G）能匹配上转换参考的 G，而真正的 m6A 位点在比对到转换参考时会表现为读到 A（即参考是 G、读端为 A 的 G→A 信号），从而被识别。

代码实现上首先通过 argparse 解析命令行（包括 read1/read2、参考 fasta、线程数、工具类型、输出前缀和输出目录等），然后 `change_reference` 函数根据 `reads` 参数对参考做碱基替换并写出新的 fasta：当 `reads=='read2'` 时（本 pipeline 的关键路径），它将参考序列中的所有 A 替换为 G（`A->G`），并在序列 ID 上加后缀 `_AG_converted`，生成供比对未甲基化 A（被读作 G）的读段使用的参考；当 `reads=='read1'` 时代码会执行 `T->C` 的替换并写出 `.TC_conversion.fa`（这一分支在 GLORI 的 A→I 原理上并非核心，通常是为了兼容其他三碱基转换的流程或特定文库/比对策略而保留的兼容选项）。

实现细节上每次写入前会尝试删除旧文件以避免残留，若输出目录不存在则创建；

`get_reversecom` 函数用于生成反向互补的转换参考（先对原序列取反向互补、在该方向上做 A→G 转换，然后再反向互补回去并写成 `.rvsCom.fa`），这样做的目的是保证来自不同链向或文库构建方向的读段都有对应的转换参考可比对，避免链偏倚导致的漏检。

为什么进行两次反向互补的原因：
首先hg38.fa文件只有正义链的信息，因此反义链的信息需要将其进行互补后得到。但是在习惯上，我们喜欢从左至右从5‘端到3’端。因此进行反向。此时得到的就是在习惯上的反义链基因组，从5‘端至3’端。但是我们下一步需要进行比对。但是star比对的时候只接受正义链的输入。此时我们就需要对反义链基因组转化为我们想要的正义链基因组。这个时候比对就可以得到了。
比如一个正义链是
```ATGC```
其反义链基因组就是 
```CGAT（5'-3'）```
在GLORI中，需要将A转化为G，此时得到的是
```CGGT（5'-3')```
进行反向互补后，得到的序列是
```ACGC（5'-3'```
但antisense 比对的时候，就会从3‘-5’进行互补配对。因此能够匹配到正义链的的需要就匹配不上反义链的序列了。此时如果有一个M6A 就没有T进行匹配了（因此T已经转化为C了）而如果是没有修饰的A 此时就会因为被乙二醇以及NaNO2转化为G 而匹配上。

`build_index` 函数根据用户选择调用不同索引构建程序：对 bowtie/bowtie2 简单调用它们的 build 命令；对 STAR 则先以去掉 `.fa` 后缀的路径作为 genomeDir，新建或重建该目录，使用 pysam 计算总基因组长度并据此用 `Nbases = int(round(min(14, np.log2(genomesize)/2 - 1)))` 自动设置 `--genomeSAindexNbases`（这是为在不同基因组规模下优化 STAR SA 索引性能与内存占用而计算的经验公式），然后以指定线程数和内存限制调用 `STAR --runMode genomeGenerate` 来生成索引。

脚本主流程依次运行：先生成（或再次生成）转换后的参考，接着生成其反向互补转换参考，然后对 AG 转换参考与反向互补参考分别构建索引并把结果放到指定输出目录；同时通过一系列打印语句记录进度与生成的路径。

总体上，这段代码的设计目标是把化学脱氨造成的 A→I（测序中读作 G）效应“预先”映射到参考序列上，通过产生 A→G 的转换参考和对应的反向互补版本来保证不同读向和链的读段都能被正确比对，进而在比对结果中以“参考为 G、读端为 A（G→A）”的突变信号来定位耐受脱氨的 m6A 位点；同时通过索引构建、目录/文件重建与序列 ID 标记等措施保证输出有序且可追溯。

FASTA 中序列默认按正链给出，基因的链向须参考 GTF/GFF 注释中的“+”/“–”字段。
```
import sys  ## 导入 sys 模块，用于处理命令行参数和标准输入
import os  ## 导入 os 模块，用于文件和目录操作
import argparse  ## 导入 argparse 模块，用于解析命令行参数
from Bio import SeqIO  ## 导入 BioPython 的 SeqIO，用于读取和写入 FASTA 文件
from Bio.Seq import Seq  ## 导入 Seq 类，用于表示 DNA 序列
from Bio.Seq import reverse_complement  ## 导入生成反向互补序列的函数
import time  ## 导入 time 模块，用于计时或日志
import numpy as np  ## 导入 numpy 模块，用于数值计算
import pysam  ## 导入 pysam 模块，用于读取基因组 FASTA 文件信息
import re  ## 导入正则模块，用于碱基替换
import subprocess  ## 导入 subprocess 模块，用于调用外部程序（如 STAR, bowtie）

parser = argparse.ArgumentParser(description = "building index for reference with three bases")  ## 创建命令行参数解析器，描述程序用途
parser.add_argument("-r", "--reads", nargs="?", type=str, default='read2', help="read1,read2,paired")  ## 参数：指定测序类型（read1/read2/paired），决定碱基转换方式
parser.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, help="reference file")  ## 参数：参考基因组 FASTA 文件路径
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "number of alignment threads to launch")  ## 参数：对齐线程数
parser.add_argument("-t", "--tools", nargs="?", type=str, default='STAR', help="bowtie,bowtie2,bwa,hisat2,STAR")  ## 参数：选择索引构建工具
parser.add_argument("-mate_length", "--mate_length", nargs="?", type=int, default=100, help="for STAR building index")  ## 参数：STAR 建索引时的 mate length
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',help = "--outname_prefix")  ## 参数：输出文件前缀
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default='./', help="outputdir")  ## 参数：输出目录

args = parser.parse_args()  ## 解析命令行参数，得到 args 对象

reads = args.reads  ## 获取测序类型
reference = args.reference  ## 获取参考基因组路径
Threads = args.Threads  ## 获取线程数
tools = args.tools  ## 获取工具类型
mate_length = args.mate_length  ## 获取 STAR mate length 参数
outputdir = args.outputdir  ## 获取输出目录
outname_prx = args.outname_prefix  ## 获取输出文件前缀

def change_reference(reads,reference,outputdir,outname_prx):  ## 定义函数：对参考基因组进行碱基转换
    if outname_prx !='default':  ## 如果用户提供了输出前缀
        refer_name = outname_prx  ## 使用用户指定前缀
    else:
        refer_name = "_".join(os.path.basename(reference).split(".")[:-1])  ## 否则使用参考文件名去掉扩展名
    if os.path.exists(outputdir):  ## 如果输出目录存在
        pass  ## 什么都不做
    else:
        os.mkdir(outputdir)  ## 否则创建目录

    if reads == "read1":  ## 如果测序类型是 read1（正向链）
        changed_refer = outputdir + "/" + refer_name + ".TC_conversion.fa"  ## 构建输出文件路径，标记为 TC 转换
        subprocess.call("rm -f " + changed_refer + " 2>/dev/null", shell=True)  ## 删除旧文件，避免写入冲突
        changed_refer2 = open(changed_refer, 'a')  ## 打开文件以追加模式写入
        for record in SeqIO.parse(reference, "fasta"):  ## 遍历参考基因组的每一条序列
            record.seq = Seq(re.sub('T', 'C', str(record.seq).upper()))  ## 将 T 替换为 C，这是 bisulfite 测序对应 read1 的转换
            record.id = record.id + "_TC_converted"  ## 修改序列 ID，标记已转换
            SeqIO.write(record, changed_refer2, "fasta")  ## 写入转换后的序列到文件

    elif reads == "read2":  ## 如果测序类型是 read2（反向链）
        changed_refer = outputdir+"/" + outname_prx + ".AG_conversion.fa"  ## 输出文件路径，标记为 AG 转换
        subprocess.call("rm -f " + changed_refer + " 2>/dev/null", shell=True)  ## 删除旧文件
        changed_refer2 = open(changed_refer,'a')  ## 打开文件写入
        for record in SeqIO.parse(reference, "fasta"):  ## 遍历每条序列
            record.seq = Seq(re.sub('A','G',str(record.seq).upper()))  ## 将 A 替换为 G，对应 bisulfite read2
            print('Genome',record.id)  ## 打印当前序列 ID，方便调试
            record.id = record.id + "_AG_converted"  ## 修改序列 ID，标记已转换
            SeqIO.write(record, changed_refer2, "fasta")  ## 写入文件
    return changed_refer  ## 返回转换后的文件路径

def get_reversecom(filename,output):  ## 定义函数：生成反向互补序列
    subprocess.call("rm -f " + output + ".rvsCom.fa",shell=True)  ## 删除旧文件
    changed_refer = open(output+".rvsCom.fa", 'a')  ## 打开输出文件
    for record in SeqIO.parse(filename, "fasta"):  ## 遍历每条参考序列
        reversed_seq = reverse_complement(record.seq)  ## 得到序列的反向互补
        converted_seq = re.sub('A', 'G', str(reversed_seq).upper())  ## 对反向互补序列做 A->G 转换
        reversed_seq2 = reverse_complement(converted_seq)  ## 再次反向互补，恢复方向，便于索引
        record.seq = Seq(reversed_seq2)  ## 更新序列
        record.id = record.id  ## 保持原序列 ID
        print('reversecomplement Genome',record.id)  ## 打印序列信息
        record.id = record.id + "_AG_converted"  ## 修改序列 ID
        SeqIO.write(record, changed_refer, "fasta")  ## 写入文件
    return output+".rvsCom.fa"  ## 返回反向互补文件路径

def build_index(changed_refer,tool,Threads):  ## 定义函数：构建索引
    if tool == "bowtie":  ## 如果工具是 bowtie
        subprocess.call('bowtie-build -q ' + changed_refer + ' ' + changed_refer,shell=True)  ## 构建 bowtie 索引
    elif tool == "bowtie2":  ## 如果工具是 bowtie2
        subprocess.call('bowtie2-build -q ' + changed_refer + ' ' + changed_refer,shell=True)  ## 构建 bowtie2 索引
    elif tool == "STAR":  ## 如果工具是 STAR
        filedir_STAR = changed_refer[:-3]  ## 去掉 .fa 作为 STAR genomeDir
        if not os.path.exists(filedir_STAR):  ## 如果目录不存在
            os.makedirs(filedir_STAR)  ## 创建目录
        else:  ## 如果目录存在
            print("Path is still exist,deleting the original filedir")  ## 打印信息
            os.rmdir(filedir_STAR)  ## 删除原目录
            os.makedirs(filedir_STAR)  ## 重建目录
        fh = pysam.FastaFile(changed_refer)  ## 打开 FASTA 文件
        genomesize = sum(fh.lengths)  ## 计算基因组总长度
        Nbases = int(round(min(14, np.log2(genomesize)/2 - 1)))  ## STAR 参数公式，优化 SA index
        command2 = 'STAR --runMode genomeGenerate -runThreadN ' + str(Threads) + ' --genomeDir ' + filedir_STAR + \
                 ' --genomeFastaFiles ' + changed_refer + ' --genomeSAindexNbases '+ str(Nbases) + ' --limitGenomeGenerateRAM 84807429045'  ## 构建 STAR 命令
        print(command2)  ## 打印命令
        subprocess.call(command2,shell=True)  ## 执行命令

if __name__ == "__main__":  ## 主程序入口
    refername2 = os.path.basename(reference)  ## 获取参考基因组文件名
    print("**********changing genome ************")  ## 打印日志
    changed_refer = change_reference(reads,reference,outputdir,outname_prx)  ## 调用碱基转换函数
    print("**********reverse complementary genome ************")  ## 打印日志
    rvs_refer = get_reversecom(reference,outputdir+outname_prx)  ## 调用反向互补函数
    changed_refer = outputdir+"/" + outname_prx + ".AG_conversion.fa"  ## 设置 AG 转换文件路径
    print("**********Building genome index for " + changed_refer + " with " + tools + "************")  ## 打印日志
    build_index(changed_refer,tools,Threads)  ## 构建索引
    print("**********Building reverse complementary genome index for " + rvs_refer + " with " + tools + "************")  ## 打印日志
    build_index(rvs_refer,tools,Threads)  ## 构建反向互补索引
    print("**********Results will be found in "+outputdir + "************")  ## 完成提示


```