# ChIP-seq-data-analysis

# Intrduction

# Datasets

# Basic workflow

1) QC fastq files

2) Trim fastq files

3) QC trimmed fastq files

4) Align trimmed fastq files

5) QC BAM files

6) Filter BAM files

7) Run ChIP-seq QC

8) Generate BigWig files for visualization

9) Peak calling

10) Peak QC

11) Peak annotation

12) Differential peak binding analysis

#Before running any of the scripts, make sure the programs are properly installed and paths are set in your environment
#move to appropriate folder before running scripts

#--------------------------------------------------------------------------------

1) QC fastq files

FASTQC used for quality metrics

Run in folder with all fastq.gz files

for file in "*fastq.gz"

do

	echo $file
  
	fastqc $file
  
done

#--------------------------------------------------------------------------------
