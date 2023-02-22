# Quality control reports. Before alignment happens we need to check out our data to make sure it looks okay.

# Drag and drop your sample file data folder into the terminal or copy the file path and use.
#cd /your/file/path/here/

# Make sure scripts are in the same directory and the correct permissions are given. 

for R1_SRC in *R1.fastq.gz*
do
    DIRECTORY=${R1_SRC//R1.fastq.gz/output}
    eval mkdir -v -m 777 $DIRECTORY
done 

# Make sure youre in the trim environment that has fastqc, cutadapt, trim-galore and fq installed.
source activate trim

# Validation of paired end data using fq. Toggle off with # if using single data.
# Exits cleanly if no error. R1 and R2 being first and second or left and right read files of paired data. Must be .gz file.
 #fq lint  R1_SRC R2_SRC 

 # Looped Pair validation. If pairs are named the same except for the ending of either R1.fastq.gz or R2.fastq.gz then it might be faster to use a loop. 
 # For this we will create a parameter R1_SRC that can be tied to a list, *R1.fastq.gz* which is comprised of all files with that ending in the directory.
 # From R1_SRC a second parameter R2_SRC will be created by searching for R1 in the file names and replacing them with R2.
 for R1_SRC in *R1.fastq.gz*
 do
     R2_SRC=${R1//R1.fastq.gz/R2.fastq.gz}
    eval fq lint $R1_SRC $R2_SRC
done

echo "Paired end files validated"

# For a one sample fastq file run (toggle command on an off by removing/adding before running.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
#fastq yoursample.fastq.gz -t 8

# Creating a loop for fastqc. Realistically, youre probably going to have multiple samples to run. Creating a loop will automate the process. 
# For this we will create a variable FASTQ that can be tied to a list, *fastq.gz* which is comprised of all files with that ending in the directory.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
for R1_SRC in *R1.fastq.gz*
do
    DIRECTORY=${R1_SRC//R1.fastq.gz/output}
    R2_SRC=${R1_SRC//R1.fastq.gz/R2.fastq.gz}
	eval fastqc $R1_SRC -t 8 --outdir $DIRECTORY
    eval fastqc $R2_SRC -t 8 --outdir $DIRECTORY
done 

# Optional-Trimming your reads. Sometimes, near the end of a long fragment you may get poor quality scores for bases. Most aligners do soft trimmings of these sequences.
# Kallisto is a pseudoaligner and doesnt soft trim. It does not need to trim adaptor sequences. 
#However, if you notice your base quality decreases near the end then trim. No more than 80%.

# trim-galore works with cutadapt and fastqc to sense adaptor sequences and poor quality bases to do trimmings.
# This loop will perform trimmings on all files ending with fastq.gz and then generate: a new trimmed fastq.gz file, a trimming report and will do another fastqc report for the trimmed data.
for R1_SRC in *R1.fastq.gz*
do
    DIRECTORY=${R1_SRC//R1.fastq.gz/output}
    R2_SRC=${R1_SRC//R1.fastq.gz/R2.fastq.gz}
	eval trim_galore --fastqc $R1_SRC -o $DIRECTORY
    eval trim_galore --fastqc $R2_SRC -o $DIRECTORY
done 

echo "QCReports have been generated"

conda deactivate
# For Bulk RNA-seq make sure youre in the rnaseq environment that has kallisto, samtools and multiqc installed.
source activate rnaseq

# Paired Bulk RNA-seq. This loop will pseudoalign samples to the refenence index. It will create output folders for each sample along with a log folder of the run.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
for R1_SRC in $(find . -type f -name '*R1_trimmed.fq.gz'); do
	R2_SRC=${R1_SRC//R1_trimmed.fq.gz/R2_trimmed.fq.gz}
	OUT="${R1_SRC%R1_trimmed.fq.gz}_mapped"
	LOG="${R1_SRC%R1_trimmed.fq.gz}_mapped.log"
	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o "$OUT" "$R1_SRC" "$R2_SRC" -t 8 &> "$LOG"
done

echo "Pseudoalignment has been completed"

# Optional-For paired end reads (Single Cell and Bulk) you may also want to create a BAM file. 
# kallisto command --pseudobam allows kallisto data to be turned into a pseudobam file. Using "|" the pseudobam file can be piped into samtools to be turned into a BAM file.
for R1_SRC in $(find . -type f -name '*R1_trimmed.fq.gz'); do
	R2_SRC=${R1_SRC//R1_trimmed.fq.gz/R2_trimmed.fq.gz}
	OUTR1="${R1_SRC%R1_trimmed.fq.gz}_mapped"
	LOGR1="${R1_SRC%R1_trimmed.fq.gz}_mapped1.log"
	OUTR2="${R2_SRC%R2_trimmed.fq.gz}_mapped}"
	LOGR2="${R2_SRC%R2_trimmed.fq.gz}_mapped2.log"
	BAMR1="${R1_SRC%R1_trimmed.fq.gz}_kallisto.fr.bam"
	BAMR2="${R2_SRC%R2_trimmed.fq.gz}_kallisto.rf.bam"
	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o "$OUTR1" --fr-stranded --pseudobam "$R1_SRC" "$R2_SRC" &> "$LOGR1" | samtools view -Sb - > "$BAMR1" &> "$LOGR1"
	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o "$OUTR2" --rf-stranded --pseudobam "$R1_SRC" "$R2_SRC" &> "$LOGR2" | samtools view -Sb - > "$BAMR2" &> "$LOGR2"
done

echo "BAM file has been created"

# Compiled Quality control. Using MultiQC you can compile all the fastqc reports and log folders from quantifications together. 
# If data was trimmed, move old fastqc reports out of directory before running.
# Running multiqc will generate an html file that can be opened in your browser to view your data.
multiqc --dirs .

echo "multiqc report has been generated"


