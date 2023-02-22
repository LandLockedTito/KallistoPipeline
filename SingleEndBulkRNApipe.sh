# Quality control reports. Before alignment happens we need to check out our data to make sure it looks okay.

# Drag and drop your sample file data folder into the terminal or copy the file path and use.
#cd /your/file/path/here/

# Make sure scripts are in the same directory and the correct permissions are given. 

for FASTQ in *fastq.gz*
do
    DIRECTORY=${FASTQ//fastq.gz/output}
    eval mkdir -v -m 777 $DIRECTORY
done  

# Make sure youre in the trim environment that has fastqc, cutadapt, trim-galore and fq installed.
source activate trim

# For a single fastq file (toggle command on an off by removing/adding before running.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
#fastq yoursample.fastq.gz -t 8

# Creating a loop for fastqc. Realistically, youre probably going to have multiple samples to run. Creating a loop will automate the process. 
# For this we will create a variable FASTQ that can be tied to a list, *fastq.gz* which is comprised of all files with that ending in the directory.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
for FASTQ in *fastq.gz*
do
	DIRECTORY=${FASTQ//fastq.gz/output}
	eval fastqc $FASTQ -t 8 --outdir $DIRECTORY
done 

# Optional-Trimming your reads. Sometimes, near the end of a long fragment you may get poor quality scores for bases. Most aligners do soft trimmings of these sequences.
# Kallisto is a pseudoaligner and doesnt soft trim. It does not need to trim adaptor sequences. 
#However, if you notice your base quality decreases near the end then trim. No more than 80%.

# trim-galore works with cutadapt and fastqc to sense adaptor sequences and poor quality bases to do trimmings.
# This loop will perform trimmings on all files ending with fastq.gz and then generate: a new trimmed fastq.gz file, a trimming report and will do another fastqc report for the trimmed data.
for FASTQ in *fastq.gz*
do
	DIRECTORY=${FASTQ//fastq.gz/output}
	eval trim_galore --fastqc $FASTQ  --output_dir $DIRECTORY
done 

echo "QCReport has finished running"

conda deactivate
# For Bulk RNA-seq make sure youre in the rnaseq environment that has kallisto, samtools and multiqc installed.
source activate rnaseq

# Single Bulk RNA-seq. This loop will pseudoalign samples to the refenence index. It will create output folders for each sample along with a log folder of the run.
# Read length (-l) value and standard deviation (-s) should be changed for your read specifications. That data should be provided by the sequencing core.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
 for FASTQ in $(find . -type f -name '*.fq.gz'); do
    OUT="${FASTQ%.fq.gz}_mapped"
    LOG="${FASTQ%.fq.gz}_log"
    eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o "$OUT" --single -l 250 -s 30 "$FASTQ" -t 8 &> "$LOG"
done

echo "Pseudoalignment has been completed"

# Compiled Quality control. Using MultiQC you can compile all the fastqc reports and log folders from quantifications together. 
# If data was trimmed, move old fastqc reports out of directory before running.
# Running multiqc will generate an html file that can be opened in your browser to view your data.
multiqc --dirs . 

echo "multiqc report has been generated"