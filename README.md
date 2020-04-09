### mtDNA_annotation_of_structure




  #### 1.) Functional gene annotation
###### 针对植物线粒体基因组注释: Mitofy；可注释几乎真核生物所有的线粒体基因组：AGORA；注释植物叶绿体和线粒体：Geneious。
###### 1.1) 植物线粒体基因组注释: Mitofy
###### 1.1.1） 使用在线网址：http://dogma.ccbb.utexas.edu/mitofy/
###### 线粒体基因组上一般分布着所有编码基因，概括来说包含编码呼吸链复合体 I、II、III、IV和 V，核糖体亚基，核糖体 RNA、转运 RNA 和细胞色素 c 等的基因。 这些基因具体为复合体 I 基因（nad1、 nad2、 nad3、 nad4、 nad4L、 nad5、 nad6、 nad7 和 nad9） 、 复合体II 基因（sdh3 和 sdh4） 、 复合体 III 基因（cob） 、 复合体 IV 基因（cox1、 cox2 和 cox3） 、复合体 V 基因（atp1、 atp4、 atp6、 atp8 和 atp9）、 细胞色素 c 生物合成基因（ccmB、 ccmC、ccmFc 和 ccmFn） 、 核糖体蛋白基因（rps1、 rps2、 rps3、 rps4、 rps7、 rps10、 rps11、 rps12、rps13、 rps14、 rps19、 rpl2、 rpl5、 rpl10 和 rpl16）、 核糖体 RNA 基因（rrn5、 rrnL 和 rrnS） 、tRNA 基因（trnN、 trnD、 trnC、 trnE、 trnQ、 trnH、 trnI、 trnK、 trnM、 trnfM、 trnF、 trnP、trnS、 trnW和 trnY） 以及编码类成熟酶的 matR 基因和编码转运子的 mttB 基因等。
###### 1.1.2） 使用mitofyX脚本
###### 下载脚本
###### 或 $ wget https://github.com/jianzuoyi/mitofyX/archive/master.zip
$ git clone https://github.com/jianzuoyi/mitofyX.git
$ cd mitofyX
###### 将此文件夹内的所有文件复制到文件夹/etc/perl中，否则会报错
$ sudo cp * /etc/perl
###### 运行脚本。/etc/perl/blast_dbs/mt_genes是自带的参考文件的文件夹，tig_4.fasta是要分析的序列，tig_4是前缀。
$ /etc/perl/mitofyX.pl --gene_db /etc/perl/blast_dbs/mt_genes tig_4.fasta tig_4
###### 会报错，没有tig_4_tRNAscan.out文件，可以提前分析下一个步骤，将文件重命名为tig_4_tRNAscan.out，放在/blast_output/tig_4文件夹内，然后重新运行上面命令
###### 结果文件夹：/blast_output/tig_4
$ ls /blast_output/tig_4
###### 1.2) tRNAscan-SE鉴定tRNA基因
###### -O 适合于线粒体和叶绿体,选择该参数，则仅使用 Cove 进行分析，搜索速度会很慢，同时也不能给出 pseudogenes 检测。-C 仅使用 Cove 进行 tRNA 分析,虽然从一定程度上提高了准确性，但是会极慢，当然不建议了。-o <file> 将结果保存到文件。-f <file> 将 tRNA 的二级结构结果保存到文件。-m <file> 将统计结果保存到文件。
$ conda create -n trnascan-se trnascan-se -y
$ conda activate trnascan-se
$ tRNAscan-SE -o tRNA_2.out -f rRNA_2.ss -m tRNA_2.stats tig_2.fasta
$ tRNAscan-SE -O -o tRNA.out -f rRNA.ss -m tRNA.stats tig_2.fasta
###### 主要结果在*.out中:
$ less .out
$ less .ss
$ less .stats



  #### 2.) Identification of open reading frames
###### 1.）使用在线ORF：https://indra.mullins.microbiol.washington.edu/sms2/orf_find.html
###### 2.）使用ORFfinder脚本
###### 下载脚本
$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
###### 解压
$ gunzip ORFfinder.gz
###### 获取执行权限
$ sudo chmod 777 ORFfinder  
$ ./ORFfinder -h
###### 先测试一下，什么都没有显示就表示一切配置OK
$ ./ORFfinder -dryrun
###### 运行
###### 参数：-b 接一个整数，指定序列上从哪个位点开始翻译（begin）。-e 接一个整数，指定翻译到哪一个位点结束（end）。-c 接t/f（或true/false），指明序列是否是环形，默认是false。-g 接一个整数，指定是哪一套密码子（可选0-25，0表示起始位点为ATG的ORF，默认是0），大部分植物线粒体均使用通用密码子。-ml 接一个整数指定ORF的最小长度是多少个nt，默认值75。-n 接t/f（或true/false）指明是否忽略ORF里边的ORF(完全被包含的ORF)，默认为false。-out 接输出的文件名。-outfmt 接一个整数，如下（默认是0）：0 = list of ORFs in FASTA format, 1 = CDS in FASTA format, 2 = Text ASN.1, 3 = Feature table
$ ./ORFfinder -in tig_4.fasta -out orf.out
$ . /ORFfinder -in tig_4.fasta -out orf_0.out -outfmt 0
$ ./ORFfinder -in tig_4.fasta -out orf_1.out -outfmt 1
$ ./ORFfinder -in tig_4.fasta -out orf_2.out -outfmt 2
$ ./ORFfinder -in tig_4.fasta -out orf_3.out -outfmt 3
###### 查看结果：
$ less orf_0.out
$ less orf_1.out
$ less orf_2.out
$ less orf_3.out
###### 拿到氨基酸序列后，直接做blastp，如果有匹配到，就是正确的ORF区了。另外也可以同时用几个工具来验证，如blastp,Pfam,CDD。



  #### 3.) Repeat sequence analysis
###### 3.1) 随机重复序列
###### 3.1.1） 使用在线软件 bibiserv2（https://bibiserv.cebitec.uni-bielefeld.de/reputer?id=reputer_view_submission） 进行识别和定位， 主要识别对象包括线粒体基因组的正向重复， 反向重复和回文重复， Minimal Repeat Size 设置为500， 即取500bp 为最小重复长度， 两条重复序列间的相似度要求达 99%以上。 可以视图。
###### 3.1.2）使用vmatch软件：分析结果同在线软件 bibiserv2
###### $ wget http://www.vmatch.de/distributions/vmatch-2.3.1-Linux_x86_64-64bit.tar.gz
$ conda create -n vmatch vmatch -y
$ conda activate vmatch
$ vmatch -help
###### 构建索引
$ mkvtree -db tig_4.fasta -v -pl -sti1 -bwt -dna -bck -suf -lcp -tis -ois -skp
###### 运行软件
###### -d 正向直接匹配。-p 反向互补/回文序列。 -l 500 最小长度500。-best  99 两条重复序列间的相似度要求达 99%以上。
$ vmatch -d -p -h 3 -l 500 -best 99 -noscore -noidentity -absolute tig_4.fasta > tig_4_repeat.txt
$ less tig_4_repeat.txt
###### 3.1.3）#通过BLASTN搜索以99％的可信度鉴定出大的重复序列（> 500-bp）
$ conda create -n blast blast -y
$ conda activate blast
$ makeblastdb -in tig_4.fasta -dbtype nucl -parse_seqids -hash_index
$ blastn -db tig_4.fasta -query tig_4.fasta -out tig_4_blast.out -outfmt 6 -perc_identity 99 -num_threads 8
###### 结果文件：tig_4_blast.out
###### 第一列为Query(递交序列)，第二列为数据库序列(目标序列subejct)，第三列为: identity，第四列为：比对长度，第五列为：错配数，第六列为：gap数，第七列和第八列为：Query开始碱基位置和结束碱基位置，第九列和第十列为：Subject开始碱基位置和结束碱基位置，第十一列为：期望值，第十二列为：比对得分
$ less tig_4_blast.out
###### 3.2) 串联重复序列
###### 3.2.1）使用在线软件 Tandem Repeats Finder（http://tandem.bu.edu/trf/trf.basic.submit.html） 对串联重复序列进行分析和定位， 参数设置为默认。
###### 3.2.2）使用trf软件
###### 安装trf软件
$ conda create -n trf trf -y
$ conda activate trf
$ trf -h
###### 用官网默认参数运行
$ trf tig_4.fasta 2 7 7 80 10 50 500 -f -d -m
###### 结果文件：tig_4.fasta.2.7.7.80.10.50.500.1.html       tig_4.fasta.2.7.7.80.10.50.500.1.txt.html       tig_4.fasta.2.7.7.80.10.50.500.dat      tig_4.fasta.2.7.7.80.10.50.500.mask
###### 3.2.3) 利用MISA软件对串联重复中的简单序列重复SSR（长度不小于 8bp 的 1-6 个碱基的重复） 做进一步分析，筛选单核苷酸，二核苷酸，三核苷酸，四核苷酸，五核苷酸和六核苷酸重复序列
$ wget https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_22092015.zip
$ unzip misa_sourcecode_22092015.zip
###### 解压后有两个文件：misa.ini   misa.pl，运行前放在一起
###### misa.ini最小重复次数阈值分别设置为10、5、4、3、3和3 ， 并通过命令行： perl misa.pl <FASTAfile>进行分析。
$ perl misa.pl tig_4.fasta
###### 结果文件：tig_4.fasta.statistics     mtDNA.gff      tig_4.fasta.misa



  #### 4.) Search for MtDNA sequences homologous to ctDNA
$ conda create -n blast blast -y
$ conda activate blast
###### 通过blastn搜索相应植物叶绿体基因组的叶绿体同源序列，具有70％的可信度，e值为1e-5。
###### 构建核酸BLAST数据库
###### -in参数后面接将要格式化的数据库，-parse_seqids, -hash_index两个参数一般都带上，主要是为blastdbcmd取子序列时使用，-dbtype nucl告诉程序这是核酸数据库
$ makeblastdb -in chloroplast.fasta -dbtype nucl -parse_seqids -hash_index
###### 运行blast
###### 参数：-db：为前面用于建库的基因组文件。-out：为输出文件的文件名。-evalue：为筛选标准(evalue越低，相似性越高)。-perc_identity 70 : 相似度大于70。-query： 输入文件路径及文件名；-out：输出文件路径及文件名。-outfmt：输出文件格式，总共有12种格式，6是tabular格式对应BLAST的m8格式。-num_threads：线程数。
###### 结果文件说明见：http://www.omicsclass.com/article/505
$ blastn -db chloroplast.fasta -query tig_4.fasta -out chloroplast_blast.out -outfmt 6 -evalue 1e-5 -perc_identity 70 -num_threads 8
###### 结果文件：chloroplast_blast.out
$ less chloroplast_blast.out



  #### 5.) Mitochondrial genome circle diagram
###### 在线网站Genome Vx：http://wolfe.ucd.ie/GenomeVx/
###### 根据基因/重复序列/编码区域/手动输入需要格式的数据，然后绘图。或者输入Genbank格式的文件，生成基因组图。
