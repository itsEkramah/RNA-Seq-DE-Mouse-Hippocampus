#Downloads the Mus musculus GRCm39 primary assembly and GTF from Ensembl.

cd /media/ekramah/DISK2/RNASEQ2/00_reference/

# Define reference file names
REF_FASTA="Mus_musculus.GRCm39.dna.primary_assembly.fa"
REF_GTF="Mus_musculus.GRCm39.114.gtf"

echo "Downloading mouse reference genome and GTF (Ensembl GRCm39, Release 114)..."
if [ ! -f "$REF_FASTA" ]; then
    wget http://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/dna/"$REF_FASTA".gz 2>> ../08_logs/wget_ref.log
    gunzip "$REF_FASTA".gz 2>> ../08_logs/gunzip_ref.log
else
    echo "$REF_FASTA already exists. Skipping download."
fi
if [ ! -f "$REF_GTF" ]; then
    wget http://ftp.ensembl.org/pub/release-114/gtf/mus_musculus/"$REF_GTF".gz 2>> ../08_logs/wget_gtf.log
    gunzip "$REF_GTF".gz 2>> ../08_logs/gunzip_gtf.log
else
    echo "$REF_GTF already exists. Skipping download."
fi
echo "Reference files setup complete."