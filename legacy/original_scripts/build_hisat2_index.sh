#Build HISAT2 Index for Mouse Genome:

 #   Creates the necessary index files from your reference FASTA for efficient alignment.

 cd /media/ekramah/DISK2/RNASEQ2/00_reference/

REF_FASTA="Mus_musculus.GRCm39.dna.primary_assembly.fa"
HISAT2_INDEX_NAME="mouse_grcm39_index"

echo "Building HISAT2 index for mouse genome..."
if [ ! -f "${HISAT2_INDEX_NAME}.1.ht2" ]; then
    hisat2-build -p 8 "$REF_FASTA" "$HISAT2_INDEX_NAME" 2>> ../08_logs/hisat2_build.log
    echo "HISAT2 index building complete."
else
    echo "HISAT2 index already exists. Skipping build."
fi