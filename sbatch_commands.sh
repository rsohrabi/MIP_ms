#Prepare to separate script files: File1= sbatchPipeline.sh

!/bin/bash --login
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --job-name MulSub
#SBATCH --mail-type=ALL
#SBATCH --mail-user= "write email here"
#SBATCH --output=%x-%j.SLURMout

#point to file path for Trim_Star.sh
SB_script="../Trim_Star.sh"

ls -1 *R1*.fastq.gz | awk -F '_' '{print $1 "_" $2 "_" $3 "_" $4}' | sort | uniq > ID

for i in `cat ./ID`;
do echo $i
   sbatch $SB_script ${i};
done


#File2: Trim_star.sh

#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=10  
#SBATCH --mem=50G            
#SBATCH --job-name PipeTAC
#SBATCH --mail-type=ALL
#SBATCH --mail-user="write email here"
#SBATCH --output=%x-%j.SLURMout

# The path for directory and several files needs to be defined before runing the sbatch. 
i=$1   
module load Trimmomatic/
#for i in `cat ./ID`;
#do echo 
#merging Lane1/2 files 
cat $i\_L001_R1_001.fastq.gz $i\_L002_R1_001.fastq.gz > $i\_R1.fastq.gz;
cat $i\_L001_R2_001.fastq.gz $i\_L002_R2_001.fastq.gz > $i\_R2.fastq.gz;

#trimming step, point to adaptor directory 
adapterDir= ../ProjectData/Ref/adapters
java -jar $TRIMMOMATIC PE $i\_R1.fastq.gz $i\_R2.fastq.gz $i\_R1_trimmomatic_paired.fastq.gz $i\_R1_trimmomatic_UNpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
gunzip $i\_R1_trimmomatic_paired.fastq.gz


#STAR alignment
module load GCC #/9.3.0 
module load STAR #/2.7.10a

GENOMEDIR=../ProjectData/Ref
GTF_file=../ProjectData/Ref/TAIR10_GFF3_genes.gtf 
GFF_file=../ProjectData/Ref/TAIR10_GFF3_genes.gff 

nb_cpu=10
genomeDir=../ProjectData/Ref/STAR		 
genomeNCBI=$GENOMEDIR/TAIR10_chr_all.fas			 
#point to fastq Directory
fastqDir= #directory path 
GTF_file=$GTF_file
fastq_file_pe1=$i\_R1_trimmomatic_paired.fastq                                                                                  
inputReadfile1=$fastqDir$fastq_file_pe1

limit_nb_multireads=4
overhang_unannotated_junctions=8
overhang_annotated_junctions=1
min_intron_length=20
max_intron_length=1000000
quality_threshold=255
norm_wiggle=None # RPM or None
read_length=50 
strand_wig="Unstranded"
outFilter_mismatchNMax=(0.06 * $read_length)
BAM_dir=$i\_R1_trimmomatic_BAM
mkdir $fastqDir$BAM_dir

STAR --runThreadN $nb_cpu --genomeDir $genomeDir --readFilesIn $inputReadfile1 --outFileNamePrefix $fastqDir${BAM_dir}/${BAM_dir} --alignSJoverhangMin $overhang_unannotated_junctions --alignSJDBoverhangMin $overhang_annotated_junctions --outFilterMismatchNmax $outFilter_mismatchNMax --alignIntronMin $min_intron_length --alignIntronMax $max_intron_length --outSAMmapqUnique $quality_threshold --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outWigType wiggle --outWigNorm $norm_wiggle --outWigStrand $strand_wig --sjdbGTFfile $GTF_file output_BAM=Aligned.sortedByCoord.out.bam

