##1. Set up your Project Directory AND DOWNLOAD DATA
# Navigate to your project directory

cd /media/ekramah/DISK2/RNASEQ2

# Create subdirectories for raw data, reference files, results, etc.
mkdir 01_raw_data
mkdir 02_trimmed_reads
mkdir 03_alignment
mkdir 04_quantification
mkdir 05_DE_analysis
mkdir 06_variant_calling
mkdir 07_scripts
mkdir 08_logs

#. Install Miniconda (if you haven't already)

# Check if conda is installed
conda --version

# If not installed, download and install Miniconda (choose the appropriate version for your system)

# Go to https://docs.conda.io/en/latest/miniconda.html to get the latest installer link

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
rm Miniconda3-latest-Linux-x86_64.sh

# Initialize conda
source $HOME/miniconda3/bin/activate
conda init

# Close and reopen your terminal or run:
source ~/.bashrc # or ~/.zshrc if you use zsh
Part 1: Initial Setup, Data & Reference Download

    Navigate to your project base directory:

        This is the main folder where all your analysis subdirectories will reside.
    Bash

cd /media/ekramah/DISK2/RNASEQ2

Create Project Subdirectories:

    #This organizes your data, results, and logs neatly.

Bash

mkdir -p 00_reference 01_raw_data 02_qc_reports 02_trimmed_reads 03_alignment 04_quantification 05_DE_analysis 06_variant_calling 07_scripts 08_logs
echo "Project directories created."

#Create and Activate Conda Environment and Install Tools:

    This step sets up a dedicated environment named bioenv with all the necessary bioinformatics tools and R packages.

Bash

echo "Creating Conda environment 'bioenv' and installing tools..."
conda create -n bioenv multiqc fastqc trimmomatic hisat2 samtools subread sra-tools r-essentials r-biocmanager gatk4 -y
echo "Conda environment 'bioenv' created."

echo "Activating 'bioenv' environment and installing R Bioconductor packages..."
conda activate bioenv

R --vanilla <<EOF
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
BiocManager::install("ggplot2", update = FALSE, ask = FALSE)
BiocManager::install("dplyr", update = FALSE, ask = FALSE)
BiocManager::install("EnhancedVolcano", update = FALSE, ask = FALSE)
BiocManager::install("pheatmap", update = FALSE, ask = FALSE)
BiocManager::install("RColorBrewer", update = FALSE, ask = FALSE)
EOF
echo "R Bioconductor packages installed in 'bioenv'."
echo "Conda environment setup complete. 'bioenv' is now active."

Create sra_accessions.txt for Selected Samples:

    This file lists the SRA accession numbers for the 10 samples you will download.

Bash

cd /media/ekramah/DISK2/RNASEQ2/01_raw_data/
cat << EOF > sra_accessions.txt
SRR7498002
SRR7498004
SRR7498006
SRR7498008
SRR7498010
SRR7498016
SRR7498018
SRR7498020
SRR7498022
SRR7498024
EOF
echo "sra_accessions.txt created."

Download Raw RNA-Seq Data (FASTQ.gz files):

    Uses fasterq-dump from sra-tools to download and gzip the single-end reads.

Bash

# Ensure bioenv is active
conda activate bioenv
cd /media/ekramah/DISK2/RNASEQ2/01_raw_data/

echo "Downloading raw RNA-Seq data from SRA..."
while IFS= read -r accession; do
    echo "Downloading and gzipping $accession..."
    # fasterq-dump outputs uncompressed .fastq. Pipe it to gzip.
    fasterq-dump "$accession" -e 8 -O . --temp-dir /tmp/fasterqdump_temp_"$accession" 2>> ../08_logs/fasterq_dump.log
    if [ -f "$accession".fastq ]; then # Check if fastq was created before zipping
        gzip "$accession".fastq 2>> ../08_logs/gzip.log
        echo "$accession.fastq.gz created."
    else
        echo "Warning: $accession.fastq was not created by fasterq-dump. Check log for errors." | tee -a ../08_logs/fasterq_dump.log
    fi
done < sra_accessions.txt
echo "Raw data download and compression complete."
ls -lh *.fastq.gz # Verify downloaded files