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
- Dmitry E Andreev, et al. Translation of 5′ leaders is pervasive in genes resistant to eIF2 repression. eLife. 2015.
- In order to maximize the genuine ribosome footprints aligning to the transcriptome, ribo-seq reads with a length
typical for monosomes (29–35 inclusive) were used for further analysis. 

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
