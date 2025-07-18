#Align Trimmed Reads to Mouse Reference Genome (HISAT2):

 #   Maps your cleaned reads to the indexed genome, generating BAM files.
 cd /media/ekramah/DISK2/RNASEQ2/03_alignment/

REF_FASTA="Mus_musculus.GRCm39.dna.primary_assembly.fa" # Redefine if not in scope
HISAT2_INDEX_NAME="mouse_grcm39_index" # Redefine if not in scope

echo "Starting HISAT2 alignment..."
for f_trimmed in ../02_trimmed_reads/*.trimmed.fastq; do
    base=$(basename "$f_trimmed" .trimmed.fastq)
    output_sam="${base}.sam"
    output_bam="${base}.bam"
    echo "Aligning sample: $base"
    hisat2 -p 8 \
           --dta \
           -x ../00_reference/"$HISAT2_INDEX_NAME" \
           -U "$f_trimmed" \
           -S "$output_sam" 2>> ../08_logs/"${base}_hisat2.log"
    echo "Converting SAM to BAM and sorting for sample: $base"
    samtools view -bS "$output_sam" | samtools sort -o "$output_bam" 2>> ../08_logs/"${base}_samtools_sort.log"
    echo "Indexing BAM file for sample: $base"
    samtools index "$output_bam" 2>> ../08_logs/"${base}_samtools_index.log"
    rm "$output_sam" # Remove large intermediate SAM file
done
echo "Alignment complete. Sorted and indexed BAM files are in 03_alignment."
ls -lh *.bam *.bai # Verify BAM and BAI files