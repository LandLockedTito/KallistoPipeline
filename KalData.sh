# Pseudoaligning and generating counts. 
# Now that pre processing Quality Control is done, its time to generate abundance files.



# In order to  align or pseudoalign data, there must be a reference transcriptome. 
# References data is usually given as a Fasta file (.fa). For it to be useful we will need to have an index made from that file. An index is something you only----
# need to do once. 
# There should be a current kallisto human cDNA index linked in the kallisto document. However, should a new version of the reference transcriptome be made available----
# or you just want to make your own kallisto index you can do that by first downloading Homo sapien cDNA Fasta from http://useast.ensembl.org/info/data/ftp/index.html ----
# unzip the data so that it is a .fa file and then run.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
 kallisto index -i Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa 

# The next steps will depend on if youre working with Single or Paired end data, and Bulk vs Single Cell RNA-seq data.

# Optional bootstrapping. Bootstrapping can be used to increase technical replicates by rerunning alignments. 
# Any of the kallisto quant functions can have bootstraps added by adding the argument -b Number. 

# For Bulk RNA-seq make sure youre in the rnaseq environment that has kallisto, samtools and multiqc installed.
conda activate rnaseq

# Single Bulk RNA-seq. This loop will pseudoalign samples to the refenence index. It will create output folders for each sample along with a log folder of the run.
# Read length (-l) value and standard deviation (-s) should be changed for your read specifications. That data should be provided by the sequencing core.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
for FASTQ in *fastq*
do
	OUT=${FASTQ//.fastq.gz/_mapped}
	LOG=${FASTQ//.fastq.gz/_mapped.log}
	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o $OUT --single -l 250 -s 30 $FASTQ -t 8 &> $LOG
done 

# Paired Bulk RNA-seq. This loop will pseudoalign samples to the refenence index. It will create output folders for each sample along with a log folder of the run.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
for R1_SRC in *R1.fastq.gz*
do
	R2_SRC=${R1_SRC//R1.fastq.gz/R2.fastq.gz}
	OUT=${R1_SRC//R1.fastq.gz/_mapped}
	LOG=${R1_SRC//R1.fastq.gz/_mapped.log}
	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o $OUT $R1_SRC $R2_SRC -t 8 &> $LOG
done

# For Single Cell RNA-seq make sure youre in the kb environment that has kb, kallisto, samtools and multiqc installed.
conda activate kb 

# Paired Single Cell RNA-seq. scRNA-seq data will require a few more additional steps. 

# Making a t2g.txt file. kb will access a prebuilt index for human genes using -d human. The prebuilt DNA index will be mapped to our kallisto cDNA index---
# to produce a transcript to gene or t2g file. 
kb ref -d human -i Homo_sapiens.GRCh38.cdna.all.index -g t2g.txt 

# Creating Cellranger compatible files. Using kb count we can process our forward and reverse reads to create cellranger compatible data that will allow for singl cell data analysis.
# The technology being here is 10X. To know what other technologies are compatible and that could be used in this command you can use----
# kb --list
kb count file1.fastq.gz file2.fastq.gz -i Homo_sapiens.GRCh38.cdna.all.index -x 10XV3 -g t2g.txt -t 8 --cellranger &> log

# Loop kb count. 
for R1_SRC in *R1.fastq.gz*
do
	R2_SRC=${R1_SRC//R1.fastq.gz/R2.fastq.gz}
	LOG=${R1_SRC//R1.fastq.gz/_mapped.log}
	eval kb count $R1_SRC $R2_SRC -i Homo_sapiens.GRCh38.cdna.all.index -x 10XV3 -g t2g.txt -t 8 --cellranger &> $LOG
done

# Optional-For paired end reads (Single Cell and Bulk) you may also want to create a BAM file. 
# kallisto command --pseudobam allows kallisto data to be turned into a pseudobam file. Using "|" the pseudobam file can be piped into samtools to be turned into a BAM file.
for R1_SRC in *R1.fastq.gz*
do
	R2_SRC=${R1_SRC//R1.fastq.gz/R2.fastq.gz}
	OUTR1=${R1_SRC//R1.fastq.gz/_mapped}
	LOGR1=${R1_SRC//R1.fastq.gz/_mapped1.log}
	OUTR2=${R2_SRC//R2.fastq.gz/_mapped}
	LOGR2=${R2_SRC//R2.fastq.gz/_mapped2.log}
	BAMR1=${R1_SRC//R1.fastq.gz/_kallisto.fr.bam}
	BAMR2=${R2_SRC//R2.fastq.gz/_kallisto.rf.bam}
	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o $OUTR1 --fr-stranded --pseudobam $R1_SRC $R2_SRC | samtools view -Sb - > $BAMR1 &> $LOGR1
	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o $OUTR2 --rf-stranded --pseudobam $R1_SRC $R2_SRC | samtools view -Sb - > $BAMR2 &> $LOGR2
done

# Compiled Quality control. Using MultiQC you can compile all the fastqc reports and log folders from quantifications together. 
# If data was trimmed, move old fastqc reports out of directory before running.
# Running multiqc will generate an html file that can be opened in your browser to view your data.
multiqc -d .

echo "Pseudoalignment was been completed!"
