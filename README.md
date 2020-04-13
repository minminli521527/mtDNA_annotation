# mtDNA_annotation_of_structure



* ## 1) Functional gene annotation
###### Annotation for plant mitochondrial genome: Mitofy; annotation of almost all mitochondrial genomes of eukaryotes: AGORA; annotation of plant chloroplasts and mitochondria: Geneious.
- * ### 1.1) Plant mitochondrial genome annotation: Mitofy
- * #### 1.1.1) Mitofy
###### Use online URL: http://dogma.ccbb.utexas.edu/mitofy/
###### The mitochondrial genome is generally distributed with all coding genes, and in general contains genes encoding respiratory chain complexes I, II, III, IV, and V, ribosomal subunits, ribosomal RNA, transport RNA, and cytochrome c. These genes are specifically complex I genes (nad1, nad2, nad3, nad4, nad4L, nad5, nad6, nad7 and nad9), complex II genes (sdh3 and sdh4), complex III genes (cob), complex IV genes (Cox1, cox2 and cox3), complex V genes (atp1, atp4, atp6, atp8 and atp9), cytochrome c biosynthesis genes (ccmB, ccmC, ccmFc and ccmFn), ribosomal protein genes (rps1, rps2, rps3 , Rps4, rps7, rps10, rps11, rps12, rps13, rps14, rps19, rpl2, rpl5, rpl10, and rpl16), ribosomal RNA genes (rrn5, rrnL and rrnS), tRNA genes (trnN, trnD, trnC, trnD, trnC , TrnH, trnI, trnK, trnM, trnfM, trnF, trnP, trnS, trnW, and trnY), matR genes encoding mature enzymes and mttB genes encoding transporters, etc.

- * #### 1.1.2) Using mitofyX script
###### Download script
	$ wget https://github.com/jianzuoyi/mitofyX/archive/master.zip
###### or
	$ git clone https://github.com/jianzuoyi/mitofyX.git
	$ cd mitofyX
###### Copy all the files in this folder to the folder /etc/perl, otherwise an error will be reported.
	$ sudo cp * /etc/perl
###### Run the script. 
###### /etc/perl/blast_dbs/mt_genes is the folder of reference files that comes with it, tig_4.fasta is the sequence to be analyzed, and tig_4 is the prefix.
	$ /etc/perl/mitofyX.pl --gene_db /etc/perl/blast_dbs/mt_genes tig_4.fasta tig_4
###### An error will be reported, there is no file tig_4_tRNAscan.out, could analyze the next step in advance, rename the file to tig_4_tRNAscan.out, put it in the /blast_output/tig_4 folder, and then rerun the above command.
###### Results folder：/blast_output/tig_4
	$ ls /blast_output/tig_4
- * ### 1.2) tRNAscan-SE identifies tRNA genes
###### Parameter analysis: -O is suitable for mitochondria and chloroplasts, when this parameter is selected, only Cove is used for analysis, the search speed will be very slow, and pseudogenes detection cannot be given. -C only uses Cove for tRNA analysis, although the accuracy is improved to a certain extent, it will be extremely slow, of course, it is not recommended. -o, <file> save the result to a file. -f <file>, save tRNA secondary structure results to a file. -m <file>, save statistical results to a file.
	$ conda create -n trnascan-se trnascan-se -y
	$ conda activate trnascan-se
	$ tRNAscan-SE -o tRNA_2.out -f rRNA_2.ss -m tRNA_2.stats tig_2.fasta
	$ tRNAscan-SE -O -o tRNA.out -f rRNA.ss -m tRNA.stats tig_2.fasta
###### The main result is in * .out:
	$ less .out
	$ less .ss
	$ less .stats



* ## 2) Identification of open reading frames
- * ### 2.1) ORF
###### Use online URL：https://indra.mullins.microbiol.washington.edu/sms2/orf_find.html
- * ### 2.2) Using ORFfinder script
- * #### 2.2.1) Download script
###### Download script
	$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
###### 解压
	$ gunzip ORFfinder.gz
###### 获取执行权限
	$ sudo chmod 777 ORFfinder  
	$ ./ORFfinder -h
###### 先测试一下，什么都没有显示就表示一切配置OK
	$ ./ORFfinder -dryrun
- * #### 2.2.2) Run
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
###### After getting the amino acid sequence, do blastp analysis directly, if there is a match, it is the correct ORF region. In addition, could also use several tools to verify, such as blastp, Pfam, CDD.



* ## 3) Repeat sequence analysis
- * ### 3.1) 随机重复序列
- * #### 3.1.1) 使用在线软件 bibiserv2（https://bibiserv.cebitec.uni-bielefeld.de/reputer?id=reputer_view_submission） 进行识别和定位， 主要识别对象包括线粒体基因组的正向重复， 反向重复和回文重复， Minimal Repeat Size 设置为500， 即取500bp 为最小重复长度， 两条重复序列间的相似度要求达 99%以上。 可以视图。
- * #### 3.1.2) 使用vmatch软件：分析结果同在线软件 bibiserv2
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
- * #### 3.1.3) 通过BLASTN搜索以99％的可信度鉴定出大的重复序列（> 500-bp）
	$ conda create -n blast blast -y
	$ conda activate blast
	$ makeblastdb -in tig_4.fasta -dbtype nucl -parse_seqids -hash_index
	$ blastn -db tig_4.fasta -query tig_4.fasta -out tig_4_blast.out -outfmt 6 -perc_identity 99 -num_threads 8
###### 结果文件：tig_4_blast.out
###### 第一列为Query(递交序列)，第二列为数据库序列(目标序列subejct)，第三列为: identity，第四列为：比对长度，第五列为：错配数，第六列为：gap数，第七列和第八列为：Query开始碱基位置和结束碱基位置，第九列和第十列为：Subject开始碱基位置和结束碱基位置，第十一列为：期望值，第十二列为：比对得分
	$ less tig_4_blast.out
- * ### 3.2) 串联重复序列
- * #### 3.2.1) 使用在线软件 Tandem Repeats Finder（http://tandem.bu.edu/trf/trf.basic.submit.html） 对串联重复序列进行分析和定位， 参数设置为默认。
- * #### 3.2.2) 使用trf软件
###### 安装trf软件
	$ conda create -n trf trf -y
	$ conda activate trf
	$ trf -h
###### 用官网默认参数运行
	$ trf tig_4.fasta 2 7 7 80 10 50 500 -f -d -m
###### 结果文件：tig_4.fasta.2.7.7.80.10.50.500.1.html       tig_4.fasta.2.7.7.80.10.50.500.1.txt.html       tig_4.fasta.2.7.7.80.10.50.500.dat      tig_4.fasta.2.7.7.80.10.50.500.mask
- * #### 3.2.3) 利用MISA软件对串联重复中的简单序列重复SSR（长度不小于 8bp 的 1-6 个碱基的重复） 做进一步分析，筛选单核苷酸，二核苷酸，三核苷酸，四核苷酸，五核苷酸和六核苷酸重复序列
	$ wget https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_22092015.zip
	$ unzip misa_sourcecode_22092015.zip
###### 解压后有两个文件：misa.ini   misa.pl，运行前放在一起
###### misa.ini最小重复次数阈值分别设置为10、5、4、3、3和3 ， 并通过命令行： perl misa.pl <FASTAfile>进行分析。
	$ perl misa.pl tig_4.fasta
###### 结果文件：tig_4.fasta.statistics     mtDNA.gff      tig_4.fasta.misa



* ## 4) Search for MtDNA sequences homologous to ctDNA
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



* ## 5) Mitochondrial genome circle diagram
###### 在线网站Genome Vx：http://wolfe.ucd.ie/GenomeVx/
###### 根据基因/重复序列/编码区域/手动输入需要格式的数据，然后绘图。或者输入Genbank格式的文件，生成基因组图。


