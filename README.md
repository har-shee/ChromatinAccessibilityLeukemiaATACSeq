INSTALLATION REQUIREMENTS


Software dependencies (install via conda or modules):
- Nextflow >= 22.10.0
- fastqc
- trim-galore
- bowtie2
- samtools
- bedtools
- macs2
- subread (featureCounts)
- R with packages:
    * DESeq2
    * dplyr
    * ggplot2
    * pheatmap
    * ggrepel
    * gprofiler2
    * tidyr
Conda environment setup:
conda create -n atac_env -c bioconda nextflow fastqc trim-galore bowtie2 samtools bedtools macs2 subread
conda activate atac_env
conda install -c conda-forge r-base r-deseq2 r-ggplot2 r-pheatmap r-ggrepel r-dplyr r-tidyr
R -e 'install.packages("gprofiler2", repos="https://cloud.r-project.org")'






DATA DOWNLOAD INSTRUCTIONS


All the input fastq files, GrCh38, GTF files are uploaded in zip format in the Input Files folder.
Steps:
Download ENCODE ATAC-seq FASTQ files (example for K562 and GM12878):
Use ENCODE portal: https://www.encodeproject.org
Example IDs: ENCFF391BFJ, ENCFF186CQZ, ENCFF736ZWY, ENCFF411MVZ
Download via wget:
wget https://www.encodeproject.org/files/ENCFF391BFJ/@@download/ENCFF391BFJ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF186CQZ/@@download/ENCFF186CQZ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF736ZWY/@@download/ENCFF736ZWY.fastq.gz
wget https://www.encodeproject.org/files/ENCFF411MVZ/@@download/ENCFF411MVZ.fastq.gz


Download human genome (GRCh38) FASTA:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
Build Bowtie2 index:
bowtie2-build GRCh38.primary_assembly.genome.fa GRCh38_noalt_as
Download GTF annotation (Gencode v43):
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
gunzip gencode.v43.annotation.gtf.gz




USAGE INSTRUCTIONS


To run the full pipeline:
nextflow run atac_seq.nf \
  --input "path/to/*.fastq" \
  --outdir "path/to/output" \
  --genome_index "path/to/bowtie2/index/prefix" \
  --sample_info "path/to/sample_info.csv" \
  --gtf_file "path/to/annotation.gtf"






INPUT FILE DETAILS AND HOW THEY ARE USED




1. FASTQ files (*.fastq or *.fq.gz):
   - Used in: FASTQC_BEFORE → TRIM_GALORE → FASTQC_AFTER → ALIGN_BOWTIE2
   - Description: Raw sequencing reads to be processed and aligned.
2. Bowtie2 Genome Index:
   - Used in: ALIGN_BOWTIE2
   - Description: Prebuilt index files (.bt2) used by Bowtie2 to align reads.
3. Sample Information CSV (samples.txt):
   - Used in: READORDER_COUNTS → DESEQ2_ANALYSIS
   - Description: Maps each BAM file to a sample group (e.g. Leukemia or Healthy).
4. GTF Annotation File (Gencode):
   - Used in: SORT_GTF → EXTRACT_GENES_FROM_GTF → BEDTOOLS_CLOSEST_ANNOTATION
   - Description: Contains gene annotations for identifying closest genes to peaks.
5. MACS2 narrowPeak files (*.narrowPeak):
   - Used in: MERGE_NARROWPEAK_FILES
   - Description: Represent peaks of accessible chromatin from each BAM.
6. Merged Peaks BED (merged_peaks.bed):
   - Used in: CONVERT_BED_TO_SAF
   - Description: Union of all narrowPeak files, sorted and merged using bedtools.
7. SAF file (generated from merged BED):
   - Used in: FEATURECOUNTS
   - Description: Tabular format describing genomic intervals for read counting.
8. BAM files (from alignment):
   - Used in: FEATURECOUNTS
   - Description: Aligned read files used to quantify accessibility over merged peaks.
9. read_counts.txt:
   - Used in: REORDER_COUNTS
   - Description: Raw count matrix produced by featureCounts; columns = samples, rows = peaks.
10. processed_counts.csv:
    - Used in: DESEQ2_ANALYSIS → HEATMAP_VST_PLOT → PCA_PLOT
    - Description: Reordered count matrix used for downstream statistical analysis and plots.
11. deseq2_results.csv:
    - Used in: FILTER_SIGNIFICANT_PEAKS → FINAL_MERGE_SIGPEAKS_GENES → VOLCANO_PLOT
    - Description: DESeq2 output of differential peak accessibility (log2FC, padj, etc).
12. significant_peaks.csv:
    - Used in: CREATE_SIGNIFICANT_BED
    - Description: Filtered peaks with padj < 0.001 and |log2FC| > 2.
13. significant_peaks.bed:
    - Used in: SORT_PEAKS_BED → ADD_CHR_TO_PEAKS
    - Description: BED-formatted version of significant peaks.
14. significant_peaks.chr.sorted.bed:
    - Used in: BEDTOOLS_CLOSEST_ANNOTATION
    - Description: Final BED-formatted peaks with chromosome prefix for annotation.
15. genes.sorted.gtf:
    - Used in: BEDTOOLS_CLOSEST_ANNOTATION
    - Description: GTF containing only gene entries, needed for peak-to-gene mapping.
16. peaks_with_nearest_genes.csv:
    - Used in: FINAL_MERGE_SIGPEAKS_GENES
    - Description: Output of bedtools closest; links each peak to its nearest gene.
17. final_cleaned_peaks.csv:
    - Used in: VOLCANO_PLOT, ENRICHMENT_ANALYSIS
    - Description: Merged table with peak statistics and gene annotations.








PROCESS DESCRIPTIONS


FASTQC_BEFORE: Quality control check before trimming
TRIM_GALORE: Adapter trimming and quality filtering
FASTQC_AFTER: Quality control check after trimming
ALIGN_BOWTIE2: Align reads to genome using Bowtie2
MACS2_PEAKCALLING: Call peaks from aligned BAM files using MACS2
MERGE_NARROWPEAK_FILES: Combine all narrowPeak files into one unified peak set
CONVERT_BED_TO_SAF: Convert merged peaks to SAF format for featureCounts
FEATURECOUNTS: Count reads overlapping peaks
REORDER_COUNTS: Format and reorder the count matrix to match sample info
DESEQ2_ANALYSIS: Perform differential accessibility analysis using DESeq2
FILTER_SIGNIFICANT_PEAKS: Select significantly differentially accessible peaks
CREATE_SIGNIFICANT_BED: Convert significant peaks to BED format for annotation
SORT_PEAKS_BED: Sort significant peaks
ADD_CHR_TO_PEAKS: Add 'chr' prefix to peaks if not present
SORT_GTF: Sort GTF file for annotation
EXTRACT_GENES_FROM_GTF: Filter GTF to retain only gene entries
BEDTOOLS_CLOSEST_ANNOTATION: Annotate peaks with nearest genes using bedtools
FINAL_MERGE_SIGPEAKS_GENES: Merge DESeq2 and gene annotation results
VOLCANO_PLOT: Plot volcano plot of DESeq2 results with gene labels
HEATMAP_VST_PLOT: Create heatmap of top variable peaks using VST transformation
PCA_PLOT: Create PCA plot based on chromatin accessibility
ENRICHMENT_ANALYSIS: Perform GO/KEGG/Reactome enrichment analysis on significant genes
PLOT_ENRICHMENT_BAR_PLOTS: Generate bar plots and dotplots of top enriched pathways
IGV_SNAPSHOT: Generate IGV snapshot of peak regions using a batch script
