### Cut adaptor and remove rRNA reads
- 建库信息：  
  - mRNA: 
    非链特异性建库，公司未告知接头序列
  - Ribo-Seq:
    3'接头序列： AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
    Ribo-Seq建库测序策略为PE150，插入片段25-38bp，所以一端reads就已满足分析要求，只用read1 分析即可
   
- Cut adptor
  - mRNA:
    Trimmomatic
    ```
    java -jar /BioII/lulab_b/chenyinghui/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -summary $script/${i}_trim.summary.txt $dataPath/${i}_1.fq.gz $dataPath/${i}_2.fq.gz $outDir/$i/${i}_1.fq.gz $outDir/$i/${i}_1.unpaired.fq.gz $outDir/$i/${i}_2.fq.gz $outDir/$i/${i}_2.unpaired.fq.gz ILLUMINACLIP:/BioII/lulab_b/chenyinghui/software/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:4:true LEADING:10 TRAILING:10 SLIDINGWINDOW:4:10 MINLEN:15
    ```
  - Ribo-seq:
    使用公司提供的clean data数据

- Remove rRNA reads
  - mRNA:
    Bowtie1
    ```
    /BioII/lulab_b/chenyinghui/software/bowtie/bowtie-1.2.3/bowtie --norc -v 1 -M 1 -m 10000 --best -t -p 4 --un $output/"$i".rm_rRNA.fq /BioII/lulab_b/chenyinghui/database/Homo_sapiens/RefSeq/bowtie1_rRNA_index/human_rRNA_Refseq -1 $input1/$i/${i}_1.fq.gz -2 $input1/$i/${i}_2.fq.gz $output/"$i".alngned_rRNA.txt
    ```
  - Ribo-seq:
    使用公司提供的clean data数据
### Alignment

    
