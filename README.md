

# 🧬 RNA-Seq and Differential Expression Analysis of Mouse Hippocampal Stem Cells

This repository contains a complete bioinformatics pipeline for RNA sequencing (RNA-Seq) data, including quality control, read trimming, alignment, gene quantification, differential expression (DE) analysis, and variant calling. The analysis focuses on understanding gene expression changes and genetic variations in mouse hippocampal stem cells under different conditions.

---

## 📚 Project Overview

This project analyzes publicly available RNA-Seq data from *Mus musculus* (mouse) hippocampal stem cells. The primary goal is to identify differentially expressed genes between Wild Type and Notch2 Conditional Knockout (CKO) samples and to perform variant calling to identify genetic differences.

---

## 📊 Data Source

The data used in this analysis is publicly available from the NCBI Gene Expression Omnibus (GEO) and Sequence Read Archive (SRA):

* **GEO Accession**: [GSE116773](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116773)
* **SRA Study Accession**: [PRJNA480214](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA480214)
* **ENA Browser Link**: [PRJNA480214](https://www.ebi.ac.uk/ena/browser/view/PRJNA480214)

**Specific Samples Used (10 total):**

* **Wild Type (5 samples)**: SRR7498002, SRR7498004, SRR7498006, SRR7498008, SRR7498010
* **Notch2 CKO (5 samples)**: SRR7498016, SRR7498018, SRR7498020, SRR7498022, SRR7498024

---

## 🧠 Biological Context

The study investigates the role of Notch2 in hippocampal stem cells. By comparing Wild Type (control) and Notch2 CKO (experimental) samples, we aim to uncover genes and pathways affected by Notch2 deletion, providing insights into its function in neurogenesis or related processes.

---

## 🔁 Workflow Overview

1. **Project Setup**
2. **Software Installation**
3. **Data Download**
4. **Reference Data Preparation**
5. **Quality Control (FastQC & MultiQC)**
6. **Read Trimming (Trimmomatic)**
7. **Read Alignment (HISAT2)**
8. **Gene Quantification (featureCounts)**
9. **Differential Expression Analysis (DESeq2)**

---

## 🧰 Software Requirements

The pipeline is designed for **Linux-based systems** using **conda** for environment management.

```bash
conda create -n bioenv multiqc fastqc trimmomatic hisat2 samtools subread sra-tools r-essentials r-biocmanager gatk4 -y
conda activate bioenv
```

Install R/Bioconductor packages:

```r
R --vanilla <<EOF
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ggplot2", "dplyr", "EnhancedVolcano", "pheatmap", "RColorBrewer"), update = FALSE, ask = FALSE)
EOF
```

---

## 📁 Project Structure

```
.
├── README.md
├── run_analysis.sh
├── 00_reference/
├── 01_raw_data/
├── 02_qc_reports/
├── 02_trimmed_reads/
├── 03_alignment/
├── 04_quantification/
├── 05_DE_analysis/
├── 06_variant_calling/
├── 07_scripts/
│   ├── DE_analysis.R
│   └── DE_analysis_Rstudio.R
└── 08_logs/
```

---

## 🚀 How to Run the Analysis

```bash
git clone https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus.git
cd RNA-Seq-DE-Mouse-Hippocampus

# Set up environment
conda activate bioenv

# Run the pipeline
bash run_analysis.sh
```

---

## 🧾 Key Outputs & Visualizations

---

### 🧪 MultiQC Report

Summarizes quality metrics across all samples:
`02_qc_reports/multiqc_report.html`

---

### 🧮 Gene Counts

Raw gene count matrix for DESeq2:
`04_quantification/gene_counts.txt`

---

### 📉 DESeq2 Results Table

Log2 fold change, p-values, padj:
`05_DE_analysis/deseq2_results.csv`

---

### 📈 Dispersion Plot

Visualizes gene-wise dispersion estimates.
![Dispersion Plot](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/dispersion_plot/dispersion_plot-1.png)

---

### 🧬 MA Plot

Log2 fold change vs. mean expression.
![MA Plot](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/MA_plot/MA_plot-1.png)

---

### 🧮 PCA Plot

Principal component analysis shows sample clustering.
![PCA Plot](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/PCA_plot/PCA_plot-1.png)

---

### 🔥 Heatmap of Top DE Genes

Heatmap of the top 50 differentially expressed genes.
![Heatmap](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/heatmap_top_DE_genes/heatmap_top_DE_genes-1.png)

---

### 🌋 Volcano Plot

Highlights statistically significant and high fold-change genes.
![Volcano Plot](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/volcano_plot/volcano_plot-1.png)

---

### 📊 Normalized Counts for Top DE Genes

![Top Gene Counts 1](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/top_gene_counts_plots/top_gene_counts_plots-1.png)
![Top Gene Counts 2](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/top_gene_counts_plots/top_gene_counts_plots-2.png)
![Top Gene Counts 3](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/top_gene_counts_plots/top_gene_counts_plots-3.png)
![Top Gene Counts 4](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/top_gene_counts_plots/top_gene_counts_plots-4.png)
![Top Gene Counts 5](https://github.com/itsEkramah/RNA-Seq-DE-Mouse-Hippocampus/blob/main/Converted_PDF_Images/top_gene_counts_plots/top_gene_counts_plots-5.png)

---

### 🧱 Sample Distance Heatmap

Euclidean distance between samples (VST transformed).
`05_DE_analysis/sample_distance_heatmap.pdf`

---

### 🧪 P-Value Histogram

Distribution of raw p-values.
`05_DE_analysis/pvalue_histogram.pdf`

---

### 🧬 Filtered Variants

Filtered VCF file with SNPs and Indels.
`06_variant_calling/filtered_variants.vcf.gz`

---

## 🔍 Interpretation Notes

* **multiqc\_report.html**: Quality Control (multiqc_report.html): Always check this first. Look for adapter contamination, low quality scores, and GC content biases.

* **PCA\_plot.pdf**: Clear clustering = Samples from the same condition (Wild Type vs. Notch2 CKO) should cluster together. If they don't, it might indicate batch effects or unexpected biological variation.
  
* **volcano\_plot.pdf**: Upper corners = high-confidence DE genes Genes that are both highly differentially expressed (large absolute log2FC) and statistically significant (low p-value) will appear in the top-left and top-right corners.
  
* **deseq2\_results.csv**: Sort by `padj`; focus on large log2FC.  Sort by padj (adjusted p-value) to find the most significant genes. log2FoldChange indicates the direction and magnitude of change.

* **filtered\_variants.vcf.gz**: Use tools like SnpEff/VEP for annotation It is recommended to understand the functional impact of these variants.

---

## ⚖️ License

This project is open-source and available under the [MIT License](LICENSE).

---
