#Quality Control (FastQC) & MultiQC Report:

    #Assess raw read quality and generate an aggregated report.
cd /media/ekramah/DISK2/RNASEQ2/01_raw_data/

echo "Running FastQC on raw data..."
for fastq_gz_file in *.fastq.gz; do
    echo "Processing $fastq_gz_file with FastQC..."
    fastqc "$fastq_gz_file" -o ../02_qc_reports/ 2>> ../08_logs/fastqc_errors.log
done

cd /media/ekramah/DISK2/RNASEQ2/02_qc_reports/
echo "Running MultiQC to aggregate FastQC reports..."
multiqc . 2>> ../08_logs/multiqc_errors.log
echo "FastQC and MultiQC complete. Check 02_qc_reports/multiqc_report.html."
