# Pseudoaligning and generating counts. 
# Now that pre processing Quality Control is done, its time to generate abundance files.

# In order to  align or pseudoalign data, there must be a reference transcriptome. 
# References data is usually given as a Fasta file (.fa). For it to be useful we will need to have an index made from that file. An index is something you only----
# need to do once. 
# There should be a current kallisto human cDNA index linked in the kallisto document. However, should a new version of the reference transcriptome be made available----
# or you just want to make your own kallisto index you can do that by first downloading Homo sapien cDNA Fasta from http://useast.ensembl.org/info/data/ftp/index.html ----
# unzip the data so that it is a .fa file and then run.

source activate rnaseq
kallisto index -i Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa 

# Paired Single Cell RNA-seq. scRNA-seq data will require a few more additional steps. 
# For Single Cell RNA-seq make sure youre in the kb environment that has kb, kallisto, samtools and multiqc installed.
source activate singleseq 

# Making a t2g.txt file. kb will access a prebuilt index for human genes using -d human. The prebuilt DNA index will be mapped to our kallisto cDNA index---
# to produce a transcript to gene or t2g file. 
kb ref -d human -i Homo_sapiens.GRCh38.cdna.all.index -g t2g.txt 

echo "Kallisto index and t2g.txt have been created"