# Getting Kallisto data into Rstudio 

#Start by creating a new R project.

# Accessing repositories to download programmes.
setRepositories()
# it will prompt you to choose which repositories to activate. select them via listing out the numbers. 1 2 3 4 5 6 7 
# Install or activate the packages listed below. 

# Packages to activate or install----
library(rhdf5) # For hdf5 files from bootstraps.
library(tidyverse) # Provides helpful R packages for data science.
library(tximport) # How Kallisto results are brought into R.
library(ensembldb) # Helps interface with ensembl.
library(EnsDb.Hsapiens.v86) # Human-specific database package.
library(beepr) # Surprise function, very necessary.
library(datapasta) # Adds copy and past function for data.

# A study design txt file shouldve been created containing atleast three columns: a column called Sample that has the same sample name used in kallisto,
# and then columns information relevant to the study. It will be used to provide paths to our data. ----

targets <- read_tsv("studydesign.txt") # Using the read_tsv function we can read data into R and turn them into an object.

path <- file.path(targets$sample, "abundance.tsv") # Each sample now has a filepath to its abundance file from the kallisto analysis. These filepaths are in the object path.

all(file.exists(path)) # This will output a TRUE response if all files have been correctly mapped and exist in path. 

# An annotated reference genome object needs to be created in the R environment. ----
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) # Data can be taken directly from Ensembl, specifying just the transcript IDs and gene names.

Tx <- as_tibble(Tx) # This data can then be transformed into tidyverse table or tibble. 

Tx <- dplyr::rename(Tx, target_id = tx_id) # Using dplyr, a programme in the tidyverse, we can rename the tibble column name from tx_id to target_id. 
# The column title is import for it to be recognised by the programmes later. Also note that for the commands rename and select, dplyr is specified since
# those functions are present in other packages in the project.

Tx <- dplyr::select(Tx, "target_id", "gene_name") # This step just insures that the column target_id is first in the table, order matters. 

# Using tximport to bring the data into Rstudio. ----
Txi_gene <- tximport(path, 
                          type = "kallisto", 
                          tx2gene = Tx, 
                          txOut = FALSE, # FALSE gives gene data, TRUE will give transcript data. 
                          countsFromAbundance = "lengthScaledTPM",
                          ignoreTxVersion = TRUE)

beep(sound = 3)
print("Transfer Complete, Goodbye kallisto and Hello Rstudio")