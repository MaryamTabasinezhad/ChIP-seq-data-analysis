#!/bin/bash
#SBATCH --account=def-hpcg1641
#SBATCH --qos=privileged
#SBATCH -c 32
#SBATCH --mem 80G
#SBATCH --time 6:00:00
#SBATCH --j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maryam.tabasinezhad@queensu.ac

#--------------------------------------------------------------------------------

cd ~/hpc4798/chip_seq/run1

module load StdEnv/2020 trimmomatic/0.39


# this is the list of file names for the R1 samples. The R2 filenames are generated below
sample="INPUT1_S10_R1_001.fastq.gz INPUT2_S11_R1_001.fastq.gz Patient38_R1_001.fastq.gz \
Patient43_R1_001.fastq.gz Patient46_R1_001.fastq.gz"
# the trimmomatic trimming command. it's just here for convenience
trimmer="ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36"
# my datafiles are in a different directory so I define a path and store it in the variable data
data='/global/home/hpc4798/chip_seq/run1'
# loop over all the file names in $sample
for r1 in $sample; do
     # copy the filename, r1, to a new file name, r2 and substitute R1 in the name with R2
     # this generates the name of the R2 file
     r2=$r1
     r2="${r1/R1/R2}"
     #echo file: $dir$r1 $dir$r2
     # generate the names for the four output files, R1.unpaired, R1.paired, R2.unpaired, and R2.paired
     # from the names of the R1 and R2 input files
     # notice I skipped copying the name of r1 into r1p, and just substituted .fastq to .fastq.paired and put the result in a new variable
     r1p="${r1/.fastq.gz/.paired.fastq.gz}"
     r1u="${r1/.fastq.gz/.unpaired.fastq.gz}"
     r2p="${r2/.fastq.gz/.paired.fastq.gz}"
     r2u="${r2/.fastq.gz/.unpaired.fastq.gz}"
    # run the trimmomatic command, note the path to the datafiles, $data
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE $data/$r1 $data/$r2 $r1p $r1u $r2p $r2u $trimmer
done
# this is the end of the loop
############################################################
# post-cleaning - fastqc
############################################################
module load fastqc
mkdir fastqc_after
fastqc -q -t 20 -o fastqc_after *.fastqgz
