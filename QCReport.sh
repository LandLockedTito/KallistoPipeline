# Quality control reports. Before alignment happens we need to check out our data to make sure it looks okay.

# Drag and drop your sample file data folder into the terminal or copy the file path and use.
cd /your/file/path/here/

# Make sure scripts are in the same directory and the correct permissions are given. 

# Make sure youre in the trim environment that has fastqc, cutadapt, trim-galore and fq installed.
conda activate trim

# Depending on your application you will need to toggle codes on and off before running script.
# Place # infront of code to turn off and removed to turn on.

# Validation of paired end data using fq. Toggle off with # if using single data.
# Exits cleanly if no error. R1 and R2 being first and second or left and right read files of paired data. Must be .gz file.
 fq lint  R1_SRC R2_SRC 

 # Looped Pair validation. If pairs are named the same except for the ending of either R1.fastq.gz or R2.fastq.gz then it might be faster to use a loop. 
 # For this we will create a parameter R1_SRC that can be tied to a list, *R1.fastq.gz* which is comprised of all files with that ending in the directory.
 # From R1_SRC a second parameter R2_SRC will be created by searching for R1 in the file names and replacing them with R2.
 for R1_SRC in *R1.fastq.gz*
 do
     R2_SRC=${R1//R1.fastq.gz/R2.fastq.gz}
    eval fq lint $R1_SRC $R2_SRC
done

# Fastqc for quality control reports. Running fastqc will generate an html file for each sample that can be opened in your browser to view your data.

# For a single fastq file (toggle command on an off by removing/adding before running.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
fastq yoursample.fastq.gz -t 8

# Creating a loop for fastqc. Realistically, youre probably going to have multiple samples to run. Creating a loop will automate the process. 
# For this we will create a variable FASTQ that can be tied to a list, *fastq.gz* which is comprised of all files with that ending in the directory.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
for FASTQ in *fastq.gz*
do
	eval fastqc $FASTQ -t 8
done 

# Optional-Trimming your reads. Sometimes, near the end of a long fragment you may get poor quality scores for bases. Most aligners do soft trimmings of these sequences.
# Kallisto is a pseudoaligner and doesnt soft trim. It does not need to trim adaptor sequences. 
#However, if you notice your base quality decreases near the end then trim. No more than 80%.

# trim-galore works with cutadapt and fastqc to sense adaptor sequences and poor quality bases to do trimmings.
# This loop will perform trimmings on all files ending with fastq.gz and then generate: a new trimmed fastq.gz file, a trimming report and will do another fastqc report for the trimmed data.
for FASTQ in *fastq.gz*
do
	eval trim_galore --fastqc $FASTQ
done 

echo "QCReport has finished running!"