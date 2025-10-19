
```python
# import mxnet
import sys  ## 导入 sys 模块，用于访问与 Python 解释器有关的变量和函数（例如 sys.stdin、sys.argv）
import os  ## 导入 os 模块，用于文件路径处理、文件/目录操作等
import re  ## 导入正则表达式模块，用于字符串模式匹配与处理
import argparse  ## 导入 argparse 模块，用于解析命令行参数
import subprocess  ## 导入 subprocess 模块，用于在 Python 中调用外部 shell 命令
import itertools  ## 导入 itertools 模块，提供高效的迭代器工具（如 islice）
import time  ## 导入 time 模块，用于获取时间戳与时间格式化
from heapq import merge  ## 从 heapq 导入 merge，用于合并多个已排序序列并保持排序
import glob  ## 导入 glob，用于文件路径模式匹配（例如查找分片文件）
from time import strftime  ## 从 time 模块导入 strftime，用于时间格式化输出
from Bio.Seq import reverse_complement  ## 从 Biopython 导入 reverse_complement，用于获取 DNA/RNA 序列的反向互补

## -------------------------
## 命令行参数解析
## -------------------------
parser = argparse.ArgumentParser(description = "reads alignment")  ## 创建 ArgumentParser 对象，description 为帮助信息
parser.add_argument("-q", "--fastq", nargs="?", type=str, default=sys.stdin, help = "fastqfiles with surfix as _1.fq;_1.fastq;_2.fq;_2.fastq")  ## 添加参数 -q/--fastq：输入 fastq 文件路径，默认 stdin（可改为文件）
parser.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, help = "Index file for the genome")  ## 添加参数 -f/--reference：基因组索引（比对索引）目录或文件
parser.add_argument("-rvs", "--rvsref", nargs="?", type=str, default=sys.stdin, help = "Index file for the minus strand of the genome")  ## 添加参数 -rvs/--rvsref：反向基因组索引（可选）
parser.add_argument("-Tf", "--transref", nargs="?", type=str, default=sys.stdin, help = "Index file for the minus strand of the transcriptome")  ## 添加参数 -Tf/--transref：转录组索引（用于转录组比对）
parser.add_argument("-t", "--tools", nargs="?", type=str, default=sys.stdin,
                    help="We recommend using STAR for genome alignment and Bowtie for transcriptome alignment")  ## 添加参数 -t/--tools：比对工具名（例如 STAR 或 bowtie）
parser.add_argument("-m", "--mismatch", nargs="?", type=int, default=2, help="Permitted mapping mismatches")  ## 添加参数 -m/--mismatch：允许的错配数，默认 2
parser.add_argument("-F", "--FilterN", nargs="?", type=str, default=0.5, help="The setting for the STAR parameter --outFilterScoreMinOverLread")  ## 添加参数 -F/--FilterN：STAR 特定过滤参数（字符串形式）
parser.add_argument("-mulMax", "--mulMax", nargs="?", type=int, default=1, help="Suppress all alignments if > <int> exist")  ## 添加参数 -mulMax/--mulMax：多重比对阈值，默认 1
parser.add_argument("--combine", "--combine", help="Whether mapping to transcriptome",action="store_true")  ## 添加 flag --combine：若存在则同时比对到转录组
parser.add_argument("--untreated", "--untreated", help="If the input is untreated",action="store_true")  ## 添加 flag --untreated：表示输入未做 A->G 替换等处理
parser.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default',help = "--outname_prefix")  ## 添加参数 -pre/--outname_prefix：输出文件前缀，默认 'default'
parser.add_argument("-o", "--outputdir", nargs="?", type=str, default=sys.stdin, help="outputdir")  ## 添加参数 -o/--outputdir：输出目录
parser.add_argument("--rvs_fac", "--rvs_fac",help="Whether to map to the reverse strand of the transcriptome", action="store_true")  ## 添加 flag --rvs_fac：是否比对到反向转录组/基因组
parser.add_argument("-p", "--Threads", nargs="?", type=str, default='1', help = "Used threads")  ## 添加参数 -p/--Threads：线程数，作为字符串传递给外部命令
args = parser.parse_args()  ## 解析命令行参数，结果存放在 args 对象中

## 将解析好的参数赋值到局部变量，方便后续使用
fastq = args.fastq  ## 输入 fastq 文件路径
Threads = args.Threads  ## 使用的线程数（字符串）
reference = args.reference  ## 基因组索引路径
rvsref = args.rvsref  ## 反向基因组索引路径
transref = args.transref  ## 转录组索引路径
tools = args.tools  ## 比对工具名称
global FilterN  ## 声明全局变量 FilterN（在某些函数中会直接使用）
FilterN = args.FilterN  ## 将命令行提供的 FilterN 赋值给全局变量

mismatch = args.mismatch  ## 允许的错配数
mulMax = args.mulMax  ## 多重比对上限
outputdir = args.outputdir  ## 输出目录
outname_prx = args.outname_prefix  ## 输出文件名前缀

re_digits = re.compile(r'(\d+)')  ## 编译正则表达式，用于在字符串中匹配数字序列（用于自然数排序）

## 定义一个用于“嵌入数字的自然排序”键函数：把字符串切分为文本和整数的序列，便于按数字顺序排序
def embedded_numbers(s):  ## 定义函数 embedded_numbers(s) 作为排序键
    s2=s.strip().split("\t")  ## 先按制表符分割行，取第一列（通常为染色体或名称）
    pieces = re_digits.split(s2[0])  ## 使用正则将字符串按数字分割为片段列表，如 ['chr', 10, 'part']
    pieces[1::2] = map(int, pieces[1::2])  ## 将分割出的数字片段（在奇数位）转换为 int，便于数值比较
    return pieces  ## 返回混合了字符串与整数的片段列表，作为排序键

## 将大的 bed 文件分块排序并合并，解决一次性内存排序问题
def sort_bedfiles(bedfiles,outputfiles):  ## 定义 sort_bedfiles：输入一个 bed 文件路径，输出排序后的文件路径
    prx2 = bedfiles[:-4]  ## 去掉文件名末尾的后缀（假设为 .bed），得到前缀
    path = prx2+"_chunk_*.bed"  ## 分块文件的通配路径，用来收集分块
    chunksize = 5000000  ## 每个分块的行数阈值（这里设置为 5,000,000 行），可调整以适配内存
    fid = 1  ## 分块文件编号（从 1 开始）
    lines = []  ## 临时行缓冲区
    with open(bedfiles, 'r') as f_in:  ## 打开输入 bed 文件进行逐行读取
        f_out = open(prx2+'_chunk_{}.bed'.format(fid), 'w')  ## 打开第一个分块文件用于写入
        for line_num, line in enumerate(f_in, 1):  ## 遍历输入文件，line_num 从 1 开始
            lines.append(line)  ## 将当前行加入缓冲区
            if not line_num % chunksize:  ## 如果行号达到 chunksize 的倍数（即缓冲区满）
                lines = sorted(lines, key=embedded_numbers)  ## 对缓冲区按嵌入数字排序
                f_out.writelines(lines)  ## 将排序后的块写入分块文件
                f_out.close()  ## 关闭当前分块文件
                lines = []  ## 清空缓冲区
                fid += 1  ## 分块编号加 1
                f_out = open(prx2+'_chunk_{}.bed'.format(fid), 'w')  ## 打开下一个分块文件
        # last chunk
        if lines:  ## 处理最后一个未满的缓冲区
            lines = sorted(lines, key=embedded_numbers)  ## 排序
            f_out.writelines(lines)  ## 写入
            f_out.close()  ## 关闭
            lines = []  ## 清空

    chunks = []  ## 用于打开并收集所有分块的文件句柄列表
    for filename in glob.glob(path):  ## 使用 glob 查找所有分块文件
        chunks += [open(filename, 'r')]  ## 将每个分块文件以读取模式打开并加入列表

    with open(outputfiles, 'w') as f_out:  ## 打开最终输出文件用于写入合并结果
        print('merging bedfiles')  ## 打印提示，表示正在合并分块
        f_out.writelines(merge(*chunks, key=embedded_numbers))  ## 使用 heapq.merge 合并多个有序文件流（保持排序）
    subprocess.call("rm -f " + prx2 + "_chunk_*.bed", shell=True)  ## 删除临时分块文件，清理工作目录

## 对原始 fastq 做 A->G（或其他）替换并生成“changed” fastq 与一个记录 A 位点的 bed 文件
def change_reads(fastq,changename,output_bed,outputdir,change_fac):  ## 定义 change_reads：输入原 fastq，输出替换后的 fastq（changename）与 A 位点 bed（output_bed）
    file = open(fastq,'r')  ## 打开输入 fastq 文件用于读取
    fac_t=change_fac[1]  ## change_fac（例如 'AG'）的第二个字符（目标替换字符），本脚本默认 'G'
    fac_q=change_fac[0]  ## change_fac 的第一个字符（被替换的字符），本脚本默认 'A'
    if os.path.exists(outputdir):  ## 如果输出目录已存在则跳过创建
        pass
    else:
        os.makedirs(outputdir)  ## 否则创建输出目录
    subprocess.call("rm -f " + changename + " 2>/dev/null",shell=True)  ## 删除已有的 changename 文件（容错），避免追加冲突
    subprocess.call("rm -f " + output_bed + " 2>/dev/null",shell=True)  ## 删除已有的 bed 文件（容错）
    file_change = open(changename,'a+')  ## 以追加模式打开 changename 文件，用于写入替换后的 fastq
    file_bed = open(output_bed,'a+')  ## 以追加模式打开 output_bed，用于写入每条 read 的 A 位点信息
    fac = True  ## 控制循环的标志
    while fac:
        fr = list(itertools.islice(file, step))  ## 一次读取 step（全局定义，后面会设置 step=10000）行，作为一个批次
        list_change = []  ## 保存替换后要写入 fastq 的行（临时）
        list_bed = []  ## 保存要写入 bed 的行（临时）
        if len(fr) != 0:  ## 如果此次批次读取到内容
            list_fr = [lines.strip().split("\t") for lines in fr]  ## 去掉行尾换行并按制表符分割（通常 fastq 每行没有 tab，但这里沿用原作者处理）
            for x in range(0,len(list_fr),4):  ## fastq 文件每 4 行为一条 read（header, seq, +, qual）
                yr = list_fr[x+1][0].upper()  ## 取得序列行并转换为大写
                reads_name=list_fr[x][0].split(" ")[0]  ## 取得 read 名称（去除 header 中可能的额外注释）
                A_sites = [m.start() for m in re.finditer('A', yr)]  ## 找出序列中所有 'A' 的位置（索引基于 0）
                if len(A_sites)>=1:  ## 如果找到至少一个 A
                    A_sites2 = "_".join(map(str,A_sites))  ## 将所有位置以下划线连接成字符串（如 "0_5_10"）
                    list_bed.append([reads_name[1:],A_sites2])  ## 将 read 名称（去掉首字符 '>' 或 '@'）和 A 位点信息加入 bed 列表
                    list_change += [reads_name,yr.replace(fac_q,fac_t), \
                            list_fr[x+2][0],list_fr[x+3][0]]  ## 构建替换后的 fastq 四行（注意 header 保留原样，seq 做替换）
                else:
                    list_bed.append([reads_name[1:],'NA'])  ## 若没有 A，则在 bed 中记录 'NA'
                    list_change += [reads_name, yr, list_fr[x + 2][0], list_fr[x + 3][0]]  ## 序列不变，直接加入 list_change
            file_change.writelines("\n".join(list_change) + "\n")  ## 将本批次所有替换后行写入 changefastq（以换行连接）
            list_bed1 = ['\t'.join(map(str, it)) for it in list_bed]  ## 将 bed 列表按制表符拼接成字符串行
            file_bed.writelines("\n".join(list_bed1) + "\n")  ## 将 bed 信息写入 output_bed 文件
        else:
            file_change.close()  ## 如果没有读取到更多行，关闭文件句柄
            file_bed.close()  ## 关闭 bed 文件句柄
            fac = False  ## 结束循环
    sort_bedfiles(output_bed, output_bed + "_sorted")  ## 对生成的 bed 文件进行排序并输出为 _sorted
    subprocess.call("rm -f " + output_bed, shell=True)  ## 删除未排序的临时 bed 文件，保留排序后的 _sorted 文件

## 调用外部比对软件（bowtie 或 STAR），并处理各自需要的参数与输出
def mapping_files(tool,fastq,reference,Threads,muta_N,fqname,outputdir,mulMax,flag):  ## 定义 mapping_files：调用比对命令并返回输出 sam 路径与未比对 fastq（若有）
    outputfile = outputdir +fqname+".sam"  ## 输出 SAM 文件路径（在 STAR 情况下为前缀，会在内部产生 Aligned.out.bam）
    unmapfastq = outputdir +fqname+"_un_2.fq"  ## 未比对（unmapped）reads 输出 fastq 文件路径
    if tool == "bowtie":  ## 如果使用 bowtie（通常用于转录组比对）
        para_0 = 'bowtie -k 1 -m '+ str(mulMax)  ## -k 1: 报告 1 个最佳比对；-m <int>: 如果比对位置 > int 则抑制输出
        para_A = ' -v '+ str(muta_N)  ## -v N: 允许 N 个错配
        para_B = ' --best --strata -p ' + Threads  ## 使用 --best --strata 策略并设置线程数
        para_C = ' -x '+ reference +" "+ fastq +' -S ' + outputfile  ## 指定索引、输入 fastq 与输出 SAM
        para_unmap = ' --un ' + unmapfastq  ## 指定未比对的 reads 输出文件
        para_end = ' 2>' + outputfile +'.output'  ## 将标准错误重定向到 .output 文件（以捕获 bowtie 日志）
        command = para_0+para_A+para_B+para_C+para_unmap+para_end  ## 拼接完整命令
        print(command)  ## 打印命令（调试信息）
        subprocess.call(command,shell=True)  ## 以 shell 模式执行命令（同步阻塞）
    elif tool == "STAR":  ## 如果使用 STAR（通常用于基因组比对）
        para_0 = "STAR --runThreadN "+ Threads  ## 指定 STAR 使用的线程数
        para_g = " --genomeDir "+ reference[:-3]  ## 注意：作者通过 reference[:-3] 推断 genomeDir（可能因为 reference 参数传入带后缀），这里依赖输入格式
        para_A = " --limitOutSJcollapsed 5000000 "  ## STAR 参数，允许更大量的剪接位点合并上限
        para_B = " --outFilterMismatchNmax " + str(muta_N)  ## 指定允许的最大错配数
        # para_B_2 = " --outFilterMismatchNoverLmax 0.3"
        # para_B_3 = " --outFilterMismatchNoverReadLmax 1"
        para_B_2=''  ## 预留参数（当前为空）
        para_B_3=' --outFilterScoreMinOverLread '+FilterN+' --outFilterMatchNminOverLread '+FilterN+' --seedSearchStartLmax 30 '# increase overall mapping sensitivity  ## 一系列过滤与灵敏度参数（使用全局 FilterN）
        para_C = " --outSAMattributes All --outSAMprimaryFlag AllBestScore --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMtype BAM Unsorted"  ## STAR 输出格式与属性
        para_D = " --outFilterMultimapNmax " + str(mulMax)  ## 多重比对上限
        para_E = " --outFileNamePrefix " + outputfile[:-3] + " --readFilesIn " + fastq  ## STAR 输出文件名前缀（取 outputfile[:-3]）与输入 fastq
        para_unmap = " --outSAMunmapped Within --outReadsUnmapped Fastx"  ## 指示 STAR 将 unmapped reads 以 fastx 格式输出
        line_command = para_0+para_g+para_A+para_B+para_B_2+para_B_3+para_C+para_D+para_E + para_unmap  ## 拼接完整 STAR 命令行
        print(line_command)  ## 打印命令（调试信息）
        subprocess.call(line_command, shell=True)  ## 执行 STAR 命令（同步）
        print("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile[:-3] + 'Aligned.out.bam | samtools sort -n -O SAM > ' + outputfile)  ## 打印后续 samtools 命令（调试）
        subprocess.call("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile[:-3] + 'Aligned.out.bam | samtools sort -n -O SAM > ' + outputfile, shell=True)  ## 将 STAR 输出的 Aligned.out.bam 按 read name 排序并输出为 SAM（-F flag 过滤某些 reads）
        subprocess.call("mv " + outputfile[:-3] + 'Unmapped.out.mate1 ' + unmapfastq, shell=True)  ## 将 STAR 产生的未比对文件重命名为 unmapfastq
        subprocess.call("rm -f " + outputfile[:-3] + 'Aligned.out.sam', shell=True)  ## 清理中间文件
        subprocess.call("rm -f " + outputfile[:-3] + 'Aligned.out.bam', shell=True)  ## 清理中间文件
    return outputfile,unmapfastq  ## 返回输出 SAM 文件路径与未比对 fastq 路径

## 将 SAM 转成按指定 flag 过滤并排序的 BAM，并建立索引
def getbamfiles(outputfile,fac,Threads,flag):  ## 定义 getbamfiles：将 sam 转成 bam 并索引，fac 为输出后缀（如 '_s.bam'）
    output_bam = outputfile[:-4] + fac  ## 通过替换后缀来构建输出 bam 路径
    print("samtools view -F " + flag + " -bS -@ " + Threads+" -h "+outputfile + \
                    " | samtools sort > " + output_bam)  ## 打印将要执行的 samtools 命令（调试）
    subprocess.call("samtools view -F " + flag + " -bS -@ " + Threads+" -h "+outputfile + \
                    " | samtools sort > " + output_bam,shell=True)  ## 使用 samtools 将 SAM 转为 BAM 并排序
    subprocess.call("samtools index " + output_bam, shell=True)  ## 为产生的 BAM 创建索引（.bai）
    subprocess.call("rm -f " + outputfile, shell=True)  ## 删除原始 SAM 文件，节省空间
    return output_bam  ## 返回输出 bam 路径

## 将指定位置替换为给定字符（用于在反向或线性替换时把 A 位点改回或修改）
def multi_sub(string,sitesA,repl):  ## 定义 multi_sub：在字符串 string 的指定位置 sitesA（如 "1_3_5"）替换为字符 repl
    string_change = ''  ## 初始化结果字符串
    if sitesA!="NA":  ## 如果 sitesA 有效（不是 'NA'）
        A_list = map(int, sitesA.split("_"))  ## 将位点字符串拆分并转换为整数列表（map 对象）
        new = []  ## 新的字符列表
        for s in string:  ## 将原始字符串每个字符拆成列表
            new.append(s)
        for index in A_list:  ## 遍历所有要替换的位置
            new[index] = repl  ## 在指定位置替换为 repl
        string_change=''.join(new)  ## 将字符列表重新拼接为字符串
    else:
        string_change=string  ## 若 sitesA 为 'NA'，直接返回原始字符串
    return string_change  ## 返回修改后的字符串

## 将 mapping 后的 SAM 中的 reads 根据 bed 文件所记录的位点做“反向替换/修正”，用于恢复 A 或生成反向互补修复
def reverseReads2(outputfile_change,output_bed,reverse_fac,Threads,flag):  ## 定义 reverseReads2：为某些特殊 mapping 情况（反向互补）使用
    sorted_sam=outputfile_change[:-4] +"_sorted.sam"  ## 构建临时排序后的 sam 名（用于按 name 排序）
    print("samtools view -F " + flag + " -q 255 -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam)  ## 打印 samtools 命令（调试）
    subprocess.call("samtools view -F " + flag + " -q 255 -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam,shell=True)  ## 用 samtools 按 read name 排序并输出为临时 SAM（过滤与质量阈值）
    # print("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > "+sorted_sam )
    # subprocess.call("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam,shell=True)
    subprocess.call("rm -f "+outputfile_change,shell=True)  ## 删除原始未排序的 SAM，保留排序后的临时 SAM
    reverse_sam = outputfile_change[:-4] + "_r.sam"  ## 最终的“反向修正” SAM 文件路径
    f1 = open(sorted_sam,'r')  ## 打开排序后的 SAM 以读取
    f2 = open(output_bed+"_sorted",'r')  ## 打开之前生成并排序的 bed 文件（记录 A 位点）
    print(sorted_sam,output_bed+"_sorted")  ## 打印正在处理的文件名（调试信息）
    subprocess.call("rm -f " + reverse_sam + " 2>/dev/null",shell=True)  ## 删除已有的 reverse_sam（容错）
    file_reverse = open(reverse_sam,'a+')  ## 以追加模式打开 reverse_sam，用于写入修正后的 SAM
    fac=True  ## 循环标志
    index=0  ## 计数器：匹配到的 reads 数目（修正计数）
    index2=0  ## 计数器：遍历的 reads 总数
    old_items='la'  ## 用于缓存上一条处理过的 read（避免重复处理）
    while fac:
        fr1 = list(itertools.islice(f1, step))  ## 批量读取排序后 SAM 的若干行（batch）
        list_reverse = []  ## 存储本批次要写入的行
        if len(fr1) != 0:  ## 如果读到内容
            list_fr = [lines.strip().split("\t") for lines in fr1]  ## 将每行按制表符分割
            for items in list_fr:  ## 遍历本批次每行（每个 items 是一行拆分后的字段列表）
                index2 += 1  ## 总遍历计数器加 1
                reads_A = items[0]  ## 取得 SAM 行的第一个字段（read 名称）
                if reads_A[0] != '@' and reads_A != old_items[0]:  ## 若不是 header 且不是与上一个相同的 read
                    for row in f2:  ## 在 bed 文件中逐行查找与当前 read 匹配的条目（注意：此处使用文件游标继续读，会前进）
                        its = row.strip().split("\t")  ## 分割 bed 行
                        reads_S = its[0]  ## bed 中记录的 read 名
                        if reads_A == reads_S:  ## 找到匹配的 read
                            index += 1  ## 匹配计数加一
                            reverse_1 = reverse_complement(items[9])  ## 取得原始 mapping read 的序列（假设 items[9] 为 SEQ），并取反向互补
                            reverse_reads = multi_sub(reverse_1, its[1], reverse_fac)  ## 在反向互补序列上将 bed 中记录的位置替换为 reverse_fac（如 'A'）
                            items[9] = reverse_complement(reverse_reads)  ## 将替换后的序列再取反向互补以得到正确方向的序列并赋回 items[9]
                            list_reverse.append(items)  ## 将修改后的 SAM 行加入 list_reverse
                            old_items = items  ## 缓存当前 items
                            break  ## 退出对 bed 文件的遍历（找到后不必继续）
                elif reads_A[0] == '@':  ## 如果是 header 行（@SQ 等）
                    list_reverse.append(items)  ## 直接加入输出（保持 header）
                    index += 1  ## 计数器加一
                elif reads_A == old_items[0]:  ## 如果当前 read 与上一个相同（重复），直接复用 old_items
                    index += 1
                    list_reverse.append(old_items)
            list_reverse1 = ['\t'.join(map(str, it)) for it in list_reverse]  ## 将列表重新拼接为文本行
            file_reverse.writelines("\n".join(list_reverse1) + "\n")  ## 写入文件
        else:
            f1.close()  ## 没有读取到内容则关闭文件句柄
            f2.close()  ## 关闭 bed 文件
            file_reverse.close()  ## 关闭输出文件
            fac = False  ## 结束循环
    print("************reversed_reads == mapped reads***************",index,index2)  ## 打印修正后匹配的统计信息
    subprocess.call("rm -f " + sorted_sam, shell=True)  ## 删除临时排序 SAM
    return reverse_sam  ## 返回修正后的 SAM 路径

## 类似 reverseReads2，但更简单的修正流程（不做反向互补）
def reverseReads(outputfile_change,output_bed,reverse_fac,Threads,flag):  ## 定义 reverseReads：用于常规的按位点替换修正（非反向互补）
    sorted_sam=outputfile_change[:-4] +"_sorted.sam"  ## 构建临时排序 SAM 路径
    # print("samtools view -F " + flag + " -q 255 -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam)
    # subprocess.call("samtools view -F " + flag + " -q 255 -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam,shell=True)
    print("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam)  ## 打印命令（调试）
    subprocess.call("samtools view -F " + flag + " -@ " + Threads+" -h " + outputfile_change + " | samtools sort -n -O SAM > " +sorted_sam,shell=True)  ## 用 samtools 对 SAM 按 read name 排序

    subprocess.call("rm -f " + outputfile_change, shell=True)  ## 删除原始未排序的 SAM（清理）
    reverse_sam = outputfile_change[:-4] + "_r.sam"  ## 最终修正后的 SAM 路径
    f1 = open(sorted_sam,'r')  ## 打开排序后的 SAM
    f2 = open(output_bed+"_sorted",'r')  ## 打开 bed_sorted（A 位点记录）
    print(sorted_sam,output_bed+"_sorted")  ## 打印文件名（调试）
    subprocess.call("rm -f " + reverse_sam + " 2>/dev/null",shell=True)  ## 删除已有的 reverse_sam（容错）
    file_reverse = open(reverse_sam,'a+')  ## 以追加模式打开输出文件
    fac=True  ## 循环标志
    index=0  ## 已处理计数
    index2=0  ## 总计数
    old_items='la'  ## 缓存上一个处理项
    while fac:
        fr1 = list(itertools.islice(f1, step))  ## 批读若干行
        list_reverse = []  ## 存储本批次写入项
        if len(fr1) != 0:
            list_fr = [lines.strip().split("\t") for lines in fr1]  ## 分割每行
            for items in list_fr:
                index2 += 1  ## 总计数加一
                reads_A = items[0]  ## 取得 read 名称
                if reads_A[0] != '@' and reads_A != old_items[0]:  ## 非 header 且不是重复行
                    for row in f2:  ## 遍历 bed_sorted 找匹配
                        its = row.strip().split("\t")  ## 分割 bed 行
                        reads_S = its[0]  ## bed 中记载的 read 名称
                        if reads_A == reads_S:  ## 找到匹配
                            index += 1  ## 匹配计数 +1
                            reverse_reads = multi_sub(items[9], its[1], reverse_fac)  ## 在 items[9]（序列）指定位置替换为 reverse_fac（通常 'A'）
                            items[9] = reverse_reads  ## 将修正后的序列放回
                            list_reverse.append(items)  ## 加入输出队列
                            old_items = items  ## 缓存当前条目
                            break  ## 跳出对 bed 文件的遍历
                elif reads_A[0] == '@':  ## 如果是 header
                    list_reverse.append(items)  ## 保留 header
                    index += 1  ## 计数
                elif reads_A == old_items[0]:  ## 如果与上一个相同（重复）
                    index += 1
                    list_reverse.append(old_items)  ## 复用缓存的上一个条目
            list_reverse1 = ['\t'.join(map(str, it)) for it in list_reverse]  ## 拼接为文本行
            file_reverse.writelines("\n".join(list_reverse1) + "\n")  ## 写入文件
        else:
            f1.close()  ## 关闭排序 SAM
            f2.close()  ## 关闭 bed_sorted
            file_reverse.close()  ## 关闭输出文件
            fac = False  ## 退出循环
    print("************reversed_reads == mapped reads***************",index,index2)  ## 打印统计（匹配数 vs 遍历数）
    subprocess.call("rm -f " + sorted_sam, shell=True)  ## 删除临时排序 SAM
    return reverse_sam  ## 返回修正后的 SAM 文件路径

## -------------------------
## 主程序入口
## -------------------------
if __name__ == "__main__":  ## 当以脚本方式运行时，执行以下主程序逻辑
    global step  ## 声明 step 为全局变量（前面函数中也使用）
    step = 10000  ## 批量处理时每次读取的行数（用于 itertools.islice）
    global change_fac,fqname2  ## 声明全局变量 change_fac 与 fqname2
    change_fac = 'AG'  ## 默认替换因子，表示将 A->G（change_fac[0] 为原字符 A，change_fac[1] 为目标字符 G）
    if outname_prx != 'default':  ## 如果用户通过 -pre 指定了输出前缀
        fqname = outname_prx  ## 使用用户提供的前缀作为 fqname
    else:
        fqname = "_".join(os.path.basename(fastq).split(".")[:-1])  ## 否则从 fastq 文件名生成前缀（去掉扩展名）
    outputdir2 = outputdir+"/"  ## 确保输出目录末尾有斜杠，方便后续拼接
    if os.path.exists(outputdir2):  ## 如果输出目录存在则跳过创建
        pass
    else:
        os.makedirs(outputdir2)  ## 否则创建输出目录

    fqname2= outname_prx  ## 将 outname_prx 赋给 fqname2（用于后续命名）

    if args.untreated:  ## 若指定 --untreated，表示输入文件未经 A->G 处理，直接比对
        sys.stderr.write("[%s]untreated...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))  ## 向标准错误输出时间戳和提示
        outputfile_untreated, unmapfastq = mapping_files(tools, fastq, reference, Threads, mismatch,
                                                        fqname2, outputdir2, mulMax,'4')  ## 使用 mapping_files 对原始 fastq 进行比对（flag '4' 用于后续过滤）
        untreated_bam = getbamfiles(outputfile_untreated,"_s.bam",Threads,'4')  ## 将产生的 SAM 转为 BAM 并索引
        if args.combine:  ## 如果同时指定了 --combine，还会把未比对的 reads 比对到转录组
            sys.stderr.write("[%s]untreated map to transcriptome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))  ## 打印提示
            outputfile_untreated_unmap, _, = mapping_files('bowtie', unmapfastq,  transref, Threads,
                                                mismatch,fqname2 + "_un", outputdir2, mulMax,'4')  ## 使用 bowtie 将 unmapfastq 比对到转录组
            bamAG_unmap = getbamfiles(outputfile_untreated_unmap,"_s.bam", Threads,'4')  ## 将转录组比对结果转为 bam 并索引
    else:  ## 如果不是 untreated（即脚本默认流程：先做 A->G 替换再比对）
        changefastq = outputdir2 + "/" + fqname2 + "_" + change_fac + "changed_2.fq"  ## 生成替换后 fastq 的路径
        output_bed = outputdir2 + "/" + fqname2 + "_A.bed"  ## 记录 A 位点的 bed 文件路径
        if os.path.exists(changefastq):  ## 如果替换后的 fastq 已存在（避免重复处理）
            sys.stderr.write("[%s] Warning：the changed files already exists, please make sure the input file is correct\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))  ## 输出警告，提示检查输入文件
            pass
        else:
            sys.stderr.write("[%s] change to A>G...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))  ## 打印正在进行 A->G 替换的提示（含时间戳）
            change_reads(fastq, changefastq,output_bed,outputdir2, change_fac)  ## 调用 change_reads 执行替换并生成 bed 文件
        sys.stderr.write("[%s] map to genome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))  ## 打印正在比对到基因组的提示
        outputfile_changeAG, unmapfastq = mapping_files(tools, changefastq, reference, Threads, mismatch,
                                                      fqname2, outputdir2,mulMax,'20')  ## 使用 mapping_files 将替换后的 fastq 比对到基因组（flag '20'）
        outputfile_changeAG = outputdir2 + fqname2 + ".sam"  ## 明确 outputfile_changeAG 的路径（覆盖或规范化）
        reverse_samAG = reverseReads(outputfile_changeAG,output_bed, 'A',Threads,'20')  ## 对比对产生的 SAM 做按位点的修正（把替换位置恢复或按 A 标记）
        reversed_bamAG = getbamfiles(reverse_samAG,'s.bam', Threads,'20')  ## 将修正后的 SAM 转为 BAM 并索引
        if args.rvs_fac:  ## 如果指定了 --rvs_fac，脚本还会对未比对的 reads 进行反向基因组比对并做反向互补修正
            sys.stderr.write("[%s]map to reversegenome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))  ## 打印提示
            rvsmapfastq = outputdir2 + fqname2 + "_un_2.fq"  ## 取得未比对 reads 的 fastq 路径（来自 STAR 的输出或 mapping_files）
            outputfile_changeAG_rvsmap, _, = mapping_files(tools, rvsmapfastq, rvsref, Threads, mismatch,
                                                          fqname2 + "_rvs", outputdir2, mulMax, '4')  ## 将未比对 reads 比对到反向基因组
            reverse_samAG_rvsmap = reverseReads2(outputfile_changeAG_rvsmap, output_bed, 'A', Threads, '4')  ## 对反向比对结果做反向互补修正
            reversed_bamAG_rvsmap = getbamfiles(reverse_samAG_rvsmap, 's.bam', Threads, '4')  ## 将反向修正后的 SAM 转为 BAM 并索引
        else:
            pass  ## 若未请求反向基因组比对，则跳过

        if args.combine:  ## 若启用了 --combine，还会把 rvs（或特定未比对）reads 比对到转录组并做修正
            unmapfastq = outputdir2 + fqname2 + "_rvs_un_2.fq"  ## 指定用于转录组比对的未比对 fastq（来自上一阶段）
            sys.stderr.write("[%s]map to transcriptome...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))  ## 打印提示
            outputfile_changeAG_unmap2, _, = mapping_files('bowtie', unmapfastq, transref, Threads, mismatch,
                                                          fqname2 + "_tf", outputdir2, mulMax,'20')  ## 使用 bowtie 将这些 reads 比对到转录组（flag '20'）
            reverse_samAG_unmap2 = reverseReads(outputfile_changeAG_unmap2, output_bed, 'A', Threads,'20')  ## 对转录组比对结果做修正
            reversed_bamAG_unmap2 = getbamfiles(reverse_samAG_unmap2, 's.bam', Threads, '20')  ## 将修正后的 SAM 转为 BAM 并索引
```

---


   * 这是一个用于将测序 reads（fastq）与基因组/转录组比对的脚本，设计目的是支持一种 *A → G 模拟替换并据此判断 RNA 修饰位点（例如 m⁶A 等）* 的流程。
   * 脚本支持两种主要工作模式：`--untreated`（输入未做 A→G 替换，直接比对）和默认模式（先把序列中的 A 替换成 G，记录 A 位点，再比对，随后根据比对信息恢复/修正序列以确认 mapping 情况）。另外还有 `--combine`（同时比对到转录组） 和 `--rvs_fac`（进行反向基因组/互补处理）等可选行为。

2. **命令行参数解析**

   * 使用 `argparse` 定义了一组参数，包括输入 fastq、参考基因组/转录组索引、比对工具名称、比对参数（错配数、多重比对阈值）、输出目录、线程数等。
   * `FilterN` 用来传递给 STAR 的 `--outFilterScoreMinOverLread` 等参数以调整灵敏度。

3. **预处理：查找 A 位点并做替换**（函数 `change_reads`）

   * 将 fastq 按批次读取（`step` 行一组，脚本中默认 `step=10000`），每条 read 的第二行（即序列）中查找所有 'A' 位点（不区分大小写）。
   * 如果 read 中存在 A，则将这些位置记录到一个 bed 样式的文件（每行：read_name \t positions，以下划线分隔多个位置；没有 A 的记录为 'NA'）。
   * 同时生成一个 “changed” fastq：将序列中的所有 A 替换为 G（默认 `change_fac='AG'`），写入新的 fastq 用于后续 mapping。
   * 这个处理的目的是让原本含 A 的碱基被替换为 G，从而测试在比对后那些位置是否仍然能被 map 上（或用于判断 editing/修饰）。
   * 处理完成后会对 bed 文件做外部分块排序（`sort_bedfiles`），该函数避免一次性将超大 bed 文件全部加载到内存，先分块排序再合并（归并排序策略）。

4. **比对（mapping）步骤**（函数 `mapping_files`）

   * 支持两种工具：`bowtie`（通常用于短 read/转录组）和 `STAR`（推荐用于基因组比对）。
   * `bowtie`：构建相应参数如 `-v`（错配数）、`-m`（多重映射过滤）、`--un`（未比对输出）等，然后直接调用。
   * `STAR`：构建更复杂的参数（线程数、genomeDir、各种过滤阈值和输出选项），运行 STAR 之后再调用 `samtools` 将 STAR 的 `Aligned.out.bam` 转换并按 read name 排序输出为 SAM 文件，然后把 STAR 产生的 unmapped 输出改名为统一的 unmapfastq。
   * 函数返回的是 `outputfile`（SAM 文件路径）与 `unmapfastq`（未比对 reads 的 fastq），供后续处理使用。

5. **SAM -> BAM 并索引**（函数 `getbamfiles`）

   * 使用 `samtools view`、`samtools sort` 将 SAM 转为 BAM 并排序（这里的排序依据和参数取决于传入的 flag 与命令），然后 `samtools index` 为结果 BAM 建立索引以便后续查看与分析。
   * 完成后会删除中间 SAM，节省磁盘空间。

6. **修正 mapping 的 reads（reverseReads / reverseReads2）**

   * 在做 A->G 替换后，比对得到的 reads 里序列字段（SAM 中的第 10 列，脚本用 items[9]）可能是替换后的序列。为了恢复原始信息或根据 A 位点做特殊修正，需要把替换的位置恢复或按位替换成期望的碱基。
   * `reverseReads`：适用于不需要反向互补计算的情形（直接在 items[9] 上指定位置替换成 'A' 等）。
   * `reverseReads2`：适用于需要先取反向互补再替换的位置（如某些比对到反向基因组的情形），逻辑为：取 items[9] 的反向互补、在反向互补上做替换（multi_sub），然后再取反向互补回去赋值 items[9]，确保替换在正确的方向上进行。
   * 两个函数的实现方式都是：将 SAM 按 read name 排序（`samtools sort -n`），批量读入排序后的 SAM，按 read name 在预先生成并排序好的 bed 文件中查找匹配条目，找到后对 items[9] 做位置替换，然后写入新的修正 SAM。
   * 注意：在这两个函数中，bed 文件的读取是顺序进行的（使用文件游标逐行读取）；假如 SAM 与 bed_sorted 都按相同的 read 名称顺序，这样的线性遍历是高效的，但如果顺序不一致可能会导致查找效率问题或错误。总体设计依赖于 SAM 与 bed_sorted 都按 read 名排序或至少以可串行查找的方式匹配。

7. **主流程分支（untreated / 默认）**

   * `--untreated`：直接把输入 fastq 传给 mapping_files 做比对，得到 BAM 并索引；若 `--combine`，还会把 unmapped reads 比对到转录组。此分支跳过 A->G 替换步骤。
   * 默认流程（非 untreated）包含：

     1. 先运行 `change_reads`：从原始 fastq 生成 `changefastq`（A->G）与 `*_A.bed`（记录 A 位点并排序）。
     2. 将 `changefastq` 比对到基因组（mapping_files）；得到 SAM 与未比对 fastq。
     3. 用 `reverseReads` 把比对结果中受替换影响的序列按 bed 中记录的位置恢复或替换为目标碱基（例如恢复 A），输出修正后的 SAM（`_r.sam`）。
     4. 将修正后的 SAM 转为 BAM 并建立索引（getbamfiles）。
     5. 可选：如果开启 `--rvs_fac`，对未比对的 reads 还要比对到反向基因组并用 `reverseReads2` 做更复杂的反向互补修正，然后转为 BAM。
     6. 可选：如果开启 `--combine`，对 rvs 未比对 reads 比对到转录组并修正、转 BAM。

8. **关键参数与注意事项（实现细节）**

   * `change_fac='AG'`：默认替换 A -> G。若想改为其他替换（例如 T->C），可修改该变量或添加参数支持。
   * `step=10000`：批量读取行数，避免一次性内存占用。如果样本非常大，可能需增大或减小此值以平衡 I/O 与内存。
   * `reference[:-3]`：脚本中对 STAR 的 `--genomeDir` 使用了 `reference[:-3]` 的切片，这是基于作者对传入 `reference` 字符串格式（可能包含后缀）的假设。务必确保传入的 `reference` 与此约定一致，否则 STAR 的 genomeDir 可能被设置错误。建议在调用时传入正确的 genomeDir 或调整此行。
   * `flag` 在不同调用中传入 '4' 或 '20'，这些是 samtools 的 FLAG 过滤参数（`-F <int>`），代表过滤特定 bit 的 reads（例如 unmapped 等）。确保对这些 flag 的含义有清楚理解再进行调整。
   * `reverseReads` 与 `reverseReads2` 使用了文件游标在 bed 文件上的线性寻址策略（`for row in f2`）。这要求 SAM 排序和 bed_sorted 排序必须兼容（通常按 read 名排序）。若两者排序不一致，匹配可能失败或需要改为随机访问（例如把 bed_sorted 加载到 dict 中以支持随机查找，但会增加内存开销）。
   * 脚本大量使用外部工具（STAR、bowtie、samtools），需要在运行环境中安装这些工具并确保它们在 PATH 中可用。Biopython 也用于 `reverse_complement`，需安装 `biopython`。
   * 脚本对文件读写使用了简单的追加模式（`a+`），运行前最好清理旧输出或使用唯一的 outputdir，以免历史文件干扰本次结果。

9. **输出结果**

   * 主要输出为若干 BAM 文件（带索引 .bai），标识例如 `*_s.bam`、`*_r.s` 等，视不同分支而定（genome 比对、rvs 比对、转录组比对等）。
   * 还会在输出目录生成 `*_A.bed_sorted`（记录 A 位点并排序）以及 `*_AGchanged_2.fq`（替换后的 fastq）以及若干中间文件（部分会被脚本删除）。

#### `para_0 = "STAR --runThreadN "+ Threads`

* 启动 STAR 并设置使用的 CPU 线程数。
* 例：`--runThreadN 8` 表示使用 8 核并行处理。

---

#### `para_g = " --genomeDir "+ reference[:-3]`

* 指定 STAR 的参考基因组索引目录。
* `reference[:-3]` 是把输入的 `reference` 文件路径去掉 `.fa` 或 `.fna` 等后缀，用于匹配 STAR 索引路径。
* STAR 索引是通过 `STAR --runMode genomeGenerate` 预先构建好的目录。

---

#### `para_A = " --limitOutSJcollapsed 5000000 "`

* 控制输出的 **剪接位点（splice junctions）** 的上限。
* 默认值一般是 1,000,000；这里调大到 5,000,000，以防复杂转录本数据被截断。

---

#### `para_B = " --outFilterMismatchNmax " + str(muta_N)`

* 设置 **允许的最大错配碱基数（mismatch）**。
* 比如 `--outFilterMismatchNmax 4` 表示每条 read 最多允许 4 个碱基不匹配。
* 对于 m⁶A 或 A→G 转化分析，这个参数要略宽松（允许一定突变）。

---

#### `para_B_2` 与 `para_B_3`

```python
para_B_2=''
para_B_3=' --outFilterScoreMinOverLread '+FilterN+' --outFilterMatchNminOverLread '+FilterN+' --seedSearchStartLmax 30 '
```

* `para_B_2` 被禁用（空字符串），因为它原本控制相对错配比例（已不用）。
* `para_B_3` 用来增强 STAR 的比对灵敏度，参数含义如下：

  * `--outFilterScoreMinOverLread <float>`：最小得分与 read 长度比值阈值。

    * 比如 `0.5` 表示得分低于 0.5×read_length 的比对会被过滤。设成较低值（如 0.5）表示允许低得分的比对保留。因为在 A→G 转化的 reads 中，存在碱基错配，如果得分阈值太高，STAR 会把这些转化 reads 当作“低质量比对”丢弃

  * `--outFilterMatchNminOverLread <float>`：最小匹配碱基数与 read 长度比。默认值：0.66,这里同样设置成 FilterN（如 0.5），即要求至少 50% 的碱基必须匹配上参考序列。举例：100 bp 的 read → 至少有 50 bp 匹配上参考基因组才会保留。


    * 用同一个 `FilterN`（比如 0.5）来控制灵敏度。
  * `--seedSearchStartLmax 30`：限制 STAR 在比对时的最大种子搜索长度，值越小越容易匹配偏移 reads（提高敏感性）。

---

#### `para_C = " --outSAMattributes All --outSAMprimaryFlag AllBestScore --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMtype BAM Unsorted"`

这行控制 **输出 BAM 文件的内容和多重比对策略**：

| 参数                                 | 含义                             |
| ---------------------------------- | ------------------------------ |
| `--outSAMattributes All`           | 输出所有可用的 SAM 字段（更完整的比对信息）。      |
| `--outSAMprimaryFlag AllBestScore` | 对多重比对保留所有“最佳得分”的比对结果。          |
| `--outMultimapperOrder Random`     | 当 read 有多个比对位置时，随机选取一个位置输出。    |
| `--outSAMmultNmax 1`               | 每条 read 最多输出 1 个多重比对结果（多余的丢弃）。 |
| `--outSAMtype BAM Unsorted`        | 输出未排序的 BAM 文件。                 |

---

#### `para_D = " --outFilterMultimapNmax " + str(mulMax)`

* 控制每条 read 允许的最大多重比对次数。
* 例：`--outFilterMultimapNmax 1` 表示只保留唯一比对（严格），
  若设置为 10 则允许最多 10 个位置的比对。

---

#### `para_E = " --outFileNamePrefix " + outputfile[:-3] + " --readFilesIn " + fastq`

* 设置输出文件前缀（`--outFileNamePrefix`）和输入 reads 文件路径（`--readFilesIn`）。
* `outputfile[:-3]` 是去掉结尾 `.bam`，避免 STAR 自动重复加后缀。

---

#### `para_unmap = " --outSAMunmapped Within --outReadsUnmapped Fastx"`

* 控制 **未比对的 reads** 如何保存：

  * `--outSAMunmapped Within` → 未比对 reads 仍写入 SAM 文件中。
  * `--outReadsUnmapped Fastx` → 另存为 `.fastq` 文件，方便后续二次比对（例如用 Bowtie）。

---

### 🧭 总结：


```bash
STAR --runThreadN 8 \
--genomeDir /path/to/genome_index \
--limitOutSJcollapsed 5000000 \
--outFilterMismatchNmax 4 \
--outFilterScoreMinOverLread 0.5 \
--outFilterMatchNminOverLread 0.5 \
--seedSearchStartLmax 30 \
--outSAMattributes All \
--outSAMprimaryFlag AllBestScore \
--outMultimapperOrder Random \
--outSAMmultNmax 1 \
--outFilterMultimapNmax 1 \
--outSAMtype BAM Unsorted \
--outFileNamePrefix sample_output \
--readFilesIn sample.fastq \
--outSAMunmapped Within \
--outReadsUnmapped Fastx
```

这确保 STAR 输出既高灵敏度（适合 A→G 检测），又能保留未比对 reads 用于下游分析。
