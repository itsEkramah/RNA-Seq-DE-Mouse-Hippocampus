#Adapter and Quality Trimming (Trimmomatic):

    #Removes low-quality bases and adapter sequences from the compressed raw reads.

Bash

cd /media/ekramah/DISK2/RNASEQ2/02_trimmed_reads/

# Automatically find the adapter file path within your conda environment
ADAPTERS=$(find $CONDA_PREFIX/share/trimmomatic/adapters/ -name "TruSeq3-SE.fa" | head -n 1)
if [ -z "$ADAPTERS" ]; then
    echo "Error: TruSeq3-SE.fa not found. Ensure Trimmomatic is installed correctly." | tee -a ../08_logs/trimmomatic.log
    exit 1
fi
echo "Using Trimmomatic adapters from: $ADAPTERS"

echo "Starting Adapter and Quality Trimming with Trimmomatic..."
for f1_gz in ../01_raw_data/*.fastq.gz; do
    base=$(basename "$f1_gz" .fastq.gz)
    output_f1_trimmed="${base}.trimmed.fastq" # Output will be unzipped .fastq
    echo "Trimming sample: $base"
    trimmomatic SE -phred33 \
                   "$f1_gz" \
                   "$output_f1_trimmed" \
                   ILLUMINACLIP:"$ADAPTERS":2:30:10 \
                   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> ../08_logs/trimmomatic.log 2>&1
done
echo "Trimming complete. Trimmed reads are in 02_trimmed_reads."
ls -lh *.trimmed.fastq # Verify trimmed files
