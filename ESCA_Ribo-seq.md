### Cut adaptor and remove rRNA reads
- 建库信息：  
  - mRNA: 
    非链特异性建库，公司未告知接头序列，未提供clean data
  - Ribo-Seq:
    3'接头序列： AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC   
    Ribo-Seq建库测序策略为PE150，插入片段25-38bp，所以一端reads就已满足分析要求，只用read1 分析即可  
    公司提供了去除adaptor以及rRNA的clean data  
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
使用STAR软件比对
- mRNA
```
/BioII/lulab_b/chenyinghui/software/STAR/STAR-2.7.3a/bin/Linux_x86_64_static/STAR \\
--runThreadN 8 \\
--outFilterType BySJout \\
--outFilterMismatchNmax 2 \\
--outFilterMultimapNmax 1 \\
--outFilterMatchNmin 16 \\
--genomeDir $STAR_genome_index \\
--readFilesIn $dataPath/$i.rm_rRNA_1.fq.gz $dataPath/$i.rm_rRNA_2.fq.gz \\
--readFilesCommand 'zcat' \\
--outFileNamePrefix  $outDir/${i}/$i. \\
--outSAMtype BAM SortedByCoordinate \\
--quantMode TranscriptomeSAM GeneCounts \\
--outSAMattributes All \\
--outSAMattrRGline ID:1 LB:mRNA_seq PL:ILLUMINA SM:${i} \\
--outBAMcompression 6 \\
--outReadsUnmapped Fastx

/BioII/lulab_b/chenyinghui/software/conda3/bin/samtools sort -T $outDir/${i}/$i.Aligned.toTranscriptome.out.sorted -o $outDir/${i}/$i.Aligned.toTranscriptome.out.sorted.bam $outDir/${i}/$i.Aligned.toTranscriptome.out.bam

/BioII/lulab_b/chenyinghui/software/conda3/bin/samtools index  $outDir/${i}/$i.Aligned.toTranscriptome.out.sorted.bam
/BioII/lulab_b/chenyinghui/software/conda3/bin/samtools index $outDir/${i}/$i.Aligned.sortedByCoord.out.bam
```
- Ribo-Seq
```
/BioII/lulab_b/chenyinghui/software/STAR/STAR-2.7.3a/bin/Linux_x86_64_static/STAR \\
--runThreadN 8 \\
--outFilterType BySJout \\
--outFilterMismatchNmax 2 \\
--outFilterMultimapNmax 1 \\
--outFilterMatchNmin 16 \\
--genomeDir $STAR_genome_index \\
--readFilesIn $dataPath/$i.clean.fa.gz \\
--readFilesCommand 'zcat' \\
--outFileNamePrefix  $outDir/${i}/$i. \\
--outSAMtype BAM SortedByCoordinate \\
--quantMode TranscriptomeSAM GeneCounts \\
--outSAMattributes All \\
--outSAMattrRGline ID:1 LB:ribo_seq PL:ILLUMINA SM:${i} \\
--outBAMcompression 6 \\
--outReadsUnmapped Fastx

/BioII/lulab_b/chenyinghui/software/conda3/bin/samtools sort -T $outDir/${i}/$i.Aligned.toTranscriptome.out.sorted -o $outDir/${i}/$i.Aligned.toTranscriptome.out.sorted.bam $outDir/${i}/$i.Aligned.toTranscriptome.out.bam

/BioII/lulab_b/chenyinghui/software/conda3/bin/samtools index  $outDir/${i}/$i.Aligned.toTranscriptome.out.sorted.bam
/BioII/lulab_b/chenyinghui/software/conda3/bin/samtools index $outDir/${i}/$i.Aligned.sortedByCoord.out.bam
```

### Read length selection of Ribo-Seq
选择比对到基因组上后匹配碱基数量在26=< X =< 32区间的Ribo-seq reads做后续定量分析。
方法参考该文献  
> Dmitry E Andreev, et al. Translation of 5′ leaders is pervasive in genes resistant to eIF2 repression. eLife. 2015.  
> In order to maximize the genuine ribosome footprints aligning to the transcriptome, ribo-seq reads with a length typical for monosomes (29–35 inclusive) were used for further analysis.   
> 由于使用[Ribocode](https://github.com/xryanglab/RiboCode) 分析得到本次Ribo-seq测序的3nt周期性出现在28nt/29nt reads，所以调整为26~32 nt.

```
perl /BioII/lulab_b/chenyinghui/project/ESCA_riboSeq/Ribo_Seq/bin/get_fit_reads.pl -inBam $outDir/${i}/$i.Aligned.sortedByCoord.out.bam -outBam $outDir/${i}/$i.sorted.fit_length.bam
```
```perl
#get_fit_reads.pl
use strict;
use FindBin '$Bin';
use Getopt::Long;
use File::Basename;
my ($inBam,$outBam,$samtools,$Rscript,$help);

GetOptions(
	"inBam|i:s" => \$inBam,
	"outBam|o:s" => \$outBam,
	"samtools|s:s" =>\$samtools,
	"Rscript|r:s" =>\$Rscript,
	"help|?" => \$help
);

$samtools ||= "/BioII/lulab_b/chenyinghui/software/conda3/bin/samtools";
$Rscript  ||= "/BioII/lulab_b/chenyinghui/software/conda3/bin/Rscript";

my $outSam = "$outBam".".sam";
my $outBai = "$outBam".".bai";

open IN, "$samtools view -h $inBam |" or die "can't open $inBam: $!\n";
open(TEMP,">$outSam");
while(<IN>){
	if(/^@/){
		print TEMP $_;
	}else{
		my @bam_info = split /\s+/;
		my $CIGAR = $bam_info[5];
		my @CIGAR_M = ($CIGAR =~ /(\d+)M/g);
		my $match_length = 0;
		foreach my $num (@CIGAR_M){
			$match_length += $num ;
		}
		if($match_length >= 26 and $match_length <= 32){
			print TEMP $_;
		}
	}
}
close IN;
close TEMP;

system("$samtools view -hb -o $outBam $outSam");
system("rm $outSam");
system("$samtools index -b $outBam > $outBai");
```
### Reads Counting  
使用featureCount统计基因上的read count数量。
对于Ribo-Seq数据，我们只统计比对到CDS区域的reads数量，对于mRNA-seq数据，计数落在整个转录本的reads。参考该文献做法：
> Dmitry E Andreev, et al. Translation of 5′ leaders is pervasive in genes resistant to eIF2 repression. eLife. 2015.  
> The normalized read counts of ribo-seq reads aligning to the coding regions (as determined by inferred locations of the A-site codons) and of mRNA-seq reads aligning to the entire transcript were used for the differential expression analysis.

- mRNA
```
$featureCount -T 2 -s 0 -p -t exon -g gene_id -a $GTF -o $outDir/${i}.featurecounts.txt $dataPath/$i/${i}.Aligned.sortedByCoord.out.bam
```
- Ribo-seq
```
$featureCount -T 2 -s 0 -p -t CDS -g gene_id -a $GTF -o $outDir/${i}.featurecounts.txt $dataPath/$i/${i}.sorted.fit_length.bam
```

### Diff TE gene detection
使用[Xtail](https://www.nature.com/articles/ncomms11194) 检测TE发生差异性变化的基因。
> [Xtail Github](https://github.com/xryanglab/xtail)

```
library(xtail)
ribo <- read.table('/BioII/lulab_b/chenyinghui/project/ESCA_riboSeq/Ribo_Seq/02.read_count_featurecount_merge/ctrl_mut.featurecounts.txt',header=T, quote='',check.names=F, sep='\t',row.names=1)
mrna <- read.table('/BioII/lulab_b/chenyinghui/project/ESCA_riboSeq/mRNA_Seq/04.read_count_featurecount_merge/ctrl_mut.featurecounts.txt',header=T, quote='',check.names=F, sep='\t',row.names=1)

condition <- c("control","control","treat","treat")
results <- xtail(mrna,ribo,condition,minMeanCount=1,bins=10000)
results_tab <- resultsTable(results,sort.by="pvalue.adjust",log2FCs=TRUE, log2Rs=TRUE)
write.table(results_tab,"05.TE_Xtail/ctrl_mut.TE.xls",quote=F,sep="\t")

pdf("05.TE_Xtail/ctrl_mut.plotFCs.pdf")
plotFCs(results)
dev.off()

pdf("05.TE_Xtail/ctrl_mut.plotRs.pdf")
plotRs(results)
dev.off()

pdf("05.TE_Xtail/ctrl_mut.volcanoPlot.pdf")
volcanoPlot(results)
dev.off()
```

### GSEA pre-ranked analysis

使用标有log2(TE_Foldchage)的 TE差异基因列表进行[GSEA pre-ranked analysis](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_GSEAPreranked_Page)，可以知道哪些基因集富集了TE上调基因或富集了TE下调基因。

```
#filter by pvalue
python ./bin/filter.TE.Xtail.py pvalue $outdir/ctrl_mut.TE.xls $mRNA_readcount_Dir/ctrl_mut.featurecounts.txt $ribo_readcount_Dir/ctrl_mut.featurecounts.txt  /BioII/lulab_b/chenyinghui/database/Homo_sapiens/GRCh38/gencode.v32.annotation.gene_info.bed $outdir/ctrl_mut.TE.annot.xls $outdir/ctrl_mut.TE.up.pvalue_sig.xls $outdir/ctrl_mut.TE.down.pvalue_sig.xls $outdir/ctrl_mut.TE.pvalue.mapStat.xls

#filter by  FDR
python ./bin/filter.TE.Xtail.py FDR $outdir/ctrl_mut.TE.xls $mRNA_readcount_Dir/ctrl_mut.featurecounts.txt $ribo_readcount_Dir/ctrl_mut.featurecounts.txt /BioII/lulab_b/chenyinghui/database/Homo_sapiens/GRCh38/gencode.v32.annotation.gene_info.bed $outdir/ctrl_mut.TE.annot.xls $outdir/ctrl_mut.TE.up.FDR_sig.xls $outdir/ctrl_mut.TE.down.FDR_sig.xls $outdir/ctrl_mut.TE.FDR.mapStat.xls

cat $outdir/ctrl_mut.TE.up.pvalue_sig.xls | grep -v "gene_ID" |awk -F "\t" '{print $1"\t"$13}' > $outdir/ctrl_mut.TE.pvalue_sig.genelist.rnk
cat $outdir/ctrl_mut.TE.down.pvalue_sig.xls| grep -v "gene_ID" |awk -F "\t" '{print $1"\t"$13}' >> $outdir/ctrl_mut.TE.pvalue_sig.genelist.rnk

```
### RPF差异基因分析，以及GO富集
使用EdgeR分析Ribo-seq中RPF差异的基因(wt-vs-mut)，RPF显著上调的基因说明Ribosome显著富集在该基因。
取EdgeR得到的RPF显著上调的基因（p-value < 0.05 & log2FC > 0.45）进行GO富集分析（软件：[clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/)）

### 散点图 (mRNA_log2FC -- RPF_log2FC)
mRNA_log2FC与RPF_log2FC都是由Xtail得到的结果表中的值。
Ribosome散点图中的基因列表来自KEGG
MAPK 通路基因来自李扬从Review上收集得到

### 热图 (mRNA_log2FC -- RPF_log2FC -- TE_log2FC)
mRNA_log2FC、RPF_log2FC与TE_log2FC都是由Xtail得到的结果表中的值(log2FC_TE_final)。
作图软件: [pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap)
