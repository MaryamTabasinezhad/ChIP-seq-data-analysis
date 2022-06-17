# ChIP-seq-data-analysis

# Introduction

# Datasets

# Basic workflow

[1) QC fastq files](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#1-qc-fastq-files)

[2) Trim fastq files](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#2-trim-fastq-files)

[3) QC trimmed fastq files](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#3-qc-trimmed-files)

[4) Align trimmed fastq files](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#4-align-trimmed-fastq-files)

[5) QC BAM files](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#5-qc-bam-files)

[6) Filter BAM files](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#6-filter-bam-files)

[7) Run ChIP-seq QC](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#7-run-chip-seq-qc)

[8) Generate BigWig files for visualization](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#8-generate-bigwig-files-for-visualization)

[9) Peak calling](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#9-peak-calling)

[10) Peak QC](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#10-peak-qc)

[11) Peak annotation](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#11-peak-annotation)

[12) Differential peak binding analysis](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/README.md#12-differential-peak-binding-analysis)

Before running any of the scripts, make sure the programs are properly installed and paths are set in your environment
move to appropriate folder before running scripts

## 1) QC fastq files

`FASTQC` used for quality metrics, Run in folder with all `fastq.gz` files. 
```ruby
for file in "*fastq.gz"
do
 echo $file
 fastqc $file
done
```
##  2) Trim fastq files

  Use `Skewer` or `trimmomatic` to trim the adapter of these reads
  
##  A. SKEWER used for trimming
run in folder with fastq files, double check adapter sequences, single vs paired end sequencing, length of reads, and cores available:
```ruby
skewer-0.2.2-linux-x86_64 -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m pe -q 3 -l 100 -o output_folder -t R1.fastq.gz R2.fastq.gz
```

B. trimmomatic
```ruby
module load StdEnv/2020 trimmomatic/0.39 java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE Patient43_week_5_S6_R1_001.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
```

> _Note: for timming more than one file, you can use the [autotrim.sh file](https://github.com/MaryamTabasinezhad/ChIP-seq-data-analysis/blob/main/auto_trim.sh)_

## 3) QC trimmed files
`FASTQC` used for quality metrics

run in folder with  all trimmed fastq files in fastq.gz format 
```ruby
for file in "*fastq.gz"
do
 echo $file
 fastqc $file
done
```
> _Note:  For QC  more than one trimmed file, you can use:
```ruby
module load fastqc
mkdir fastqc_after
fastqc -q -t 20 -o fastqc_after *.fastqgz
```
## 4) Align trimmed fastq files
`STAR` used for alignment, Gencode annotations used for genome index creation, but you can also use Refseq, use the trimmed files from step (2)
first need to create genome indices

```ruby
STAR --runThreadN 10
--runMode genomeGenerate
--genomeDir /path/to/genome/indices/output
--genomeFastaFiles genome.fa
--sjdbGTFfile path/to/annotation.gtf
--sjdbOverhang 99
```

#make new directory for STAR output
#create one folder for each pair of paired end reads
```ruby
for file in *pair1.fastq
do
	mkdir /path/to/output/folders/"${file%pair1.fastq}"_Gencode
done
```
#run STAR
#10 threads, 70GB of AM limit, output file format BAM
```ruby
for file in *pair1.fastq
do
	#print names of input fastq read files
	echo $file
	echo "${file%pair1.fastq}"pair2.fastq
	STAR --runThreadN 10
		--genomeDir ~/path/to/genome/indices
		--readFilesIn $file "${file%pair1.fastq}"pair2.fastq
		--outFileNamePrefix /output/folder/path/"${file%pair1.fastq}"_Gencode/"${file%pair1.fastq}"_Gencode_
		--outSAMtype BAM SortedByCoordinate
		--limitBAMsortRAM 70000000000
done

#gzip all unzipped fastq files to save space on disk
#find all files that matches .fastq --> compress these files to maximum compression
#move into folder with all fastq files
#xargs -n 1 means one gzip process per file
#gzip -9 means maximum compression

find . -name '*.fastq' -print0 | xargs -0 -n 1 gzip -9
#move all fastq.gz files into new folder
find . -name '*fastq.gz' -exec mv -t '/new/destination/folder' {} +
```

## 5) QC BAM files
use `SAMSTAT` for `BAM QC`

produce samstat statistics for all bam files
```ruby
for file in *.bam
do
	#run samstat on all bam files
	echo $file
	samstat $file
done
```
## 6) Filter BAM files
`SAMTOOLS` 6nano for filtering

this is likely one of the most arbitrary step, prone to introducing genome wide bias

you can try multiple different types of filtering: filtering just based on `MAPQ cutoff`, removing duplicates, removing read mapping to mitochondrial genome and contigs, or any combinations

I've tried various combinations but filtering based on MAPQ values, removing reads mapping to mitochondria/contigs, and keeping duplicates gave the best RSC and NSC values (metrics for measuring signal to noise)

MAPQ cutoff of 3 --> "uniquely mapped reads cutoff, 50% error probability for reads with MAPQ=3"

filter bam files based on MAPQ values of 3

rename .bam files to filtered.bam

```ruby
for file in *.bam
do
	echo $file
	samtools view -b -q 3 $file > "${file%.bam}.filtered.bam"
done
```

track the number of reads you're losing at each filtering step

to get `total reads` in each BAM file

```ruby
for file in *.bam
do
	echo $file
	samtools view -c $file
done
```
remove reads aligning to mit and contigs

make sure the bam files in the folder are the filtered bam files from the earlier step

```ruby
for file in *filtered.bam
do
	echo $file
	#list total number of reads in mapq.bam file
	samtools view -c $file
	#convert bam file to sam file
	samtools view -h "${file%filtered.bam}.bam" > "${file%filtered.bam}.sam"

	#remove chromosomes outside of chr1-22,X,Y
	sed '/chrM/d;/random/d;/chrUn/d' < "${file%filtered.bam}.sam" > "${file%filtered.bam}simple.sam"

	#remove intermediate files
	rm "${file%filtered.bam}.sam"

	#convert sam back to bam
	samtools view -bS "${file%filtered.bam}simple.sam" > "${file%filtered.bam}simple.bam"

	#remove intermediate files
	rm "${file%filtered.bam}simple.sam"

	#print new file name
	echo "${file%filtered.bam}simple.bam"

	#list total number of reads in simple.bam file
	samtools view -c "${file%filtered.bam}simple.bam"

	#generate bam index file
	picard-tools BuildBamIndex I="${file%filtered.bam}simple.bam"

	#generate bed files for bam files
	bedtools bamtobed -i "${file%filtered.bam}simple.bam" > "${file%filtered.bam}simple.bed"

done
```
use `PICARDTOOLS` to generat BAM index

run picard tools

```ruby
for file in *simple.bam
do
	picard-tools BuildBamIndex I=$file
done
```

## 7) Run ChIP-seq QC
several metrics checked, fragment size, strand cross-correlation, BAM file correlation, BAM file PCA

use `PHANTOMPEAKQUALTOOLS` for strand cross-correlation analysis (not reliable for broad signals like H3K27me3)

run phantompeakqualtools

gives `fragment size`, `NSC`, `RSC`

```ruby
for file in *simple.bam
do
	echo $file
	Rscript /path/to/phantompeakqualtools/run_spp.R -c=$file -savp -out="${file%simple.bam}"
done
```
use `DEEPTOOLS` for average binding profiles, plot heatmaps, multibamsummary, correlation matrix, coverage, pca, etc.

multibamsummary: create comrpessed summary of bam files for downstream correlation and `PCA analysis`
```ruby
multiBamSummary bins \
--bamfiles /list/all/bam/files \
-out output.npz \
--outRawCounts output_readcounts.tab \
-bl ~/path/to/wgEncodeDacMapabilityConsensusExcludable.bed \
-p 8 -v --extendReads 147

#plot correlation

plotCorrelation --corData output.npz --plotFile output.png --corMethod spearman --whatToPlot heatmap --outFileCorMatrix output_cormatrix

plotCorrelation --corData output.npz --plotFile output.png --corMethod pearson --whatToPlot heatmap --outFileCorMatrix output_cormatrix

#plot PCA
plotPCA --corData output.npz --plotFile output.png --outFileNameData output_PCA_Data

#use NGSPLOT to create average profile plots
#can create profiles and heatmaps for tss, enhancers, or over a specified list of genes
ngs.plot.r -G hg19 -R tss -P 8 -C /path/to/config.txt -O output_file

#also see step (10) for additional QC steps using CHIPQC R package (both for BAM files pre-peak calling and BED files post peak calling)
```
## 8) Generate BigWig files for visualization

`DEEPTOOLS` used for bigwig file generation

bigwig files can be opened with `IGV`
```ruby
for file in *.bam
do
	echo $file
	bamCoverage --bam $file -o "${file%.bam}.bw" -of bigwig --binSize 10 --normalizeTo1x 2451960000 \
	-bl ~/path/to/wgEncodeDacMapabilityConsensusExcludable.bed -p 8 --extendReads 147

done
```
## 9) Peak calling

`SICER` used for peak calling, better for broad peaks

However, `MACS2` is the most popular peak caller. MACS2 also tested but called peaks separate better when called using SICER (according to PCA plots)

parameters used for `BROAD` (H3K27me3): SICER.sh . "ip file" "input file" "output directory" hg19 1 200 147 0.87 600 .01

parameters used for `NARROW` (H3K4me3): SICER.sh . "ip file" "input file" "output directory" hg19 1 200 147 0.87 200 .01

run SICER on BED files -H3K4me3
```ruby
for file in *.bed
do
	echo $file
	#create directory for results
	mkdir /path/to/output/"${file%.bed}_sicer"
	#run sicer
	SICER.sh . $file "${file%.bed}INPUT.bed" /path/to/output/"${file%.bed}_sicer" hg19 2147483647 200 147 0.87 200 .01
done
```
run SICER on BED files -H3K27me3
```ruby
for file in *.bed
do
	echo $file
	#create directory for results
	mkdir /path/to/output/"${file%.bed}_sicer"
	#run sicer
	SICER.sh . $file "${file%.bed}INPUT.bed" /path/to/output/"${file%.bed}_sicer" hg19 2147483647 200 147 0.87 600 .01
done
```
## 10) Peak QC
`ChIPQC` and `ChIPSeeker` R packages used for peak metrics

## 11) Peak annotation
`HOMER` used for gene annotation

run HOMER

for SICER PEAKS
```ruby
for file in *.bed
do
	echo $file
	annotatePeaks.pl $file hg19 > "${file%.bed}annotated.txt" -annStats output_annstats.txt -go path/to/output/folder
done
```
## 12) Differential peak binding analysis
`SICER-df` and `Diffbind` used for differential binding
```ruby
#SICER-df for H3K4me3
SICER-df.sh condition1.bed condition1_input.bed condition2.bed condition2_input.bed 200 200 0.01 0.01

#SICER-df for H3K27me3
SICER-df.sh condition1.bed condition1_input.bed condition2.bed condition2_input.bed 200 600 0.01 0.01
```
