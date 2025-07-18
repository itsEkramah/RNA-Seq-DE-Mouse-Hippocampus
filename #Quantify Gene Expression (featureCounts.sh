#Quantify Gene Expression (featureCounts):

 #   Counts reads that map to each gene based on your GTF annotation.

 cd /media/ekramah/DISK2/RNASEQ2/04_quantification/

REF_GTF="Mus_musculus.GRCm39.114.gtf" # Redefine if not in scope

echo "Starting gene expression quantification with featureCounts..."
BAM_FILES=$(find ../03_alignment/ -maxdepth 1 -name "*.bam" | tr '\n' ' ')
OUTPUT_COUNTS="gene_counts.txt"
featureCounts -T 8 \
              -a ../00_reference/"$REF_GTF" \
              -o "$OUTPUT_COUNTS" \
              $BAM_FILES 2>> ../08_logs/featureCounts.log
echo "Gene quantification complete. Counts are in 04_quantification/$OUTPUT_COUNTS."
ls -lh "$OUTPUT_COUNTS" # Verify gene counts file
head "$OUTPUT_COUNTS" # Peek at the file