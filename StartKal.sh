# Getting started with kallisto!

# Getting conda packagage manager installed and configured.
# Make sure miniconda script is downloaded and put into home folder before running.
bash ~/Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda 
source $HOME/miniconda/bin/activate

# Latest update for mac uses init zsh in terminal, if you prefer to use bash dont init for zsh 
conda init zsh

# Configure conda to allow channels for bioinformatics programmes.
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set offline false

# Creating a conda environment to house your tools for bulkRNA-seq. You will be given a base conda environment as soon as conda is installed.
# It is best to create separate environments for each type of analysis. Different programmes may have dependency conflicts, need different versions of python etc.
conda create --name rnaseq 

# Activate the new environment and install your programmes and dependencies. 
source activate rnaseq

# Kallisto is a pseudoaligner that will pseudoalign out data to a reference. 
conda install -c bioconda kallisto

# Fastqc can create quality control reports about our samples.
conda install -c bioconda fastqc

# Multiqc can compile all the data that we generate into one report.
conda install -c bioconda multiqc 

# Samtools can take kallisto data to generate BAM files.
conda install -c "bioconda/label/cf201901" samtools

# Final check to make sure everything is installed. This will display whats installed in your current environment.
conda list 

# Close that conda environment so we can make another one.
conda deactivate 

# Creating a new environemnt to house quality control and trimming programmes.
conda create --name trim

# Activate the new environment and install your programmes and dependencies.
source activate trim

# Fastqc can create quality control reports about our samples.
conda install -c bioconda fastqc

# Cutadapt can trim data files.
conda install cutadapt 

# Trim-galore wraps cutadapt and fastqc together to autodetect adaptor sequences and low quality bases and trim them.
conda install trim-galore 

# Multiqc can compile all the data that we generate into one report.
conda install -c bioconda multiqc

# fq can do many things for validating strand pairs to subsampling
conda install fq

# Final check to make sure everything is installed. This will display whats installed in your current environment.
conda list 

# Close that conda environment so we can make another one.
conda deactivate 

# Creating a conda environment to house tools for Single Cell RNA-seq. 
# Must specify this python version for this environmet for tools to work. 
conda create -y --name singleseq python=3.8 

# Activate the new environment and install your programmes and dependencies.
source activate kb

# KB, Kallisto Bustools allows for scRNA-seq analysis.
pip install kb-python

# Final check to make sure everything is installed. This will display whats installed in your current environment.
conda list 

# View some of kb capabilities for scRNA-seq data.
kb --list

# Should you ever need to delete a programme package from an environment, due to conflicts or whatever reason, 
# conda remove -n NameOfEnvironment NameOfProgramme

echo "Set up was successful successful!"

