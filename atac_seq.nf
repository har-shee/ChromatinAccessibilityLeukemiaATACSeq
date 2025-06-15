// ATAC-seq Analysis Pipeline
// Authors: Ramya Harshitha Bolla, Hemika Amilineni, George Steven Muvva

// ------------------------------------------------------------------------------------------------------------------------------------ 
// INSTALLATION REQUIREMENTS
// Software dependencies (install via conda or modules):
// - Nextflow >= 22.10.0
// - fastqc
// - trim-galore
// - bowtie2
// - samtools
// - bedtools
// - macs2
// - subread (featureCounts)
// - R with packages:
//     * DESeq2
//     * dplyr
//     * ggplot2
//     * pheatmap
//     * ggrepel
//     * gprofiler2
//     * tidyr

// Conda environment setup:
// conda create -n atac_env -c bioconda nextflow fastqc trim-galore bowtie2 samtools bedtools macs2 subread
// conda activate atac_env
// conda install -c conda-forge r-base r-deseq2 r-ggplot2 r-pheatmap r-ggrepel r-dplyr r-tidyr
// R -e 'install.packages("gprofiler2", repos="https://cloud.r-project.org")'

// -----------------------------------------------------------------------------------------------------------------------------------------------
// DATA DOWNLOAD INSTRUCTIONS

// Download human genome (GRCh38) FASTA:
// wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
// gunzip GRCh38.primary_assembly.genome.fa.gz

// Build Bowtie2 index:
// bowtie2-build GRCh38.primary_assembly.genome.fa GRCh38_noalt_as

// Download GTF annotation (Gencode v43):
// wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
// gunzip gencode.v43.annotation.gtf.gz

// Download ENCODE ATAC-seq FASTQ files (example for K562 and GM12878):
// Use ENCODE portal: https://www.encodeproject.org
// Example IDs: ENCFF391BFJ, ENCFF186CQZ, ENCFF736ZWY, ENCFF411MVZ
// Download via wget:
// wget https://www.encodeproject.org/files/ENCFF391BFJ/@@download/ENCFF391BFJ.fastq.gz
// wget https://www.encodeproject.org/files/ENCFF186CQZ/@@download/ENCFF186CQZ.fastq.gz
// wget https://www.encodeproject.org/files/ENCFF736ZWY/@@download/ENCFF736ZWY.fastq.gz
// wget https://www.encodeproject.org/files/ENCFF411MVZ/@@download/ENCFF411MVZ.fastq.gz

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// USAGE INSTRUCTIONS

// To run the full pipeline:
// nextflow run atac_seq.nf \
//   --input "path/to/*.fastq" \
//   --outdir "path/to/output" \
//   --genome_index "path/to/bowtie2/index/prefix" \
//   --sample_info "path/to/sample_info.csv" \
//   --gtf_file "path/to/annotation.gtf"

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// INPUT FILE DETAILS AND HOW THEY ARE USED

// 1. FASTQ files (*.fastq or *.fq.gz):
//    - Used in: FASTQC_BEFORE → TRIM_GALORE → FASTQC_AFTER → ALIGN_BOWTIE2
//    - Description: Raw sequencing reads to be processed and aligned.

// 2. Bowtie2 Genome Index:
//    - Used in: ALIGN_BOWTIE2
//    - Description: Prebuilt index files (.bt2) used by Bowtie2 to align reads.

// 3. Sample Information CSV (samples.txt):
//    - Used in: REORDER_COUNTS → DESEQ2_ANALYSIS
//    - Description: Maps each BAM file to a sample group (e.g. Leukemia or Healthy).

// 4. GTF Annotation File (Gencode):
//    - Used in: SORT_GTF → EXTRACT_GENES_FROM_GTF → BEDTOOLS_CLOSEST_ANNOTATION
//    - Description: Contains gene annotations for identifying closest genes to peaks.

// 5. MACS2 narrowPeak files (*.narrowPeak):
//    - Used in: MERGE_NARROWPEAK_FILES
//    - Description: Represent peaks of accessible chromatin from each BAM.

// 6. Merged Peaks BED (merged_peaks.bed):
//    - Used in: CONVERT_BED_TO_SAF
//    - Description: Union of all narrowPeak files, sorted and merged using bedtools.

// 7. SAF file (generated from merged BED):
//    - Used in: FEATURECOUNTS
//    - Description: Tabular format describing genomic intervals for read counting.

// 8. BAM files (from alignment):
//    - Used in: FEATURECOUNTS
//    - Description: Aligned read files used to quantify accessibility over merged peaks.

// 9. read_counts.txt:
//    - Used in: REORDER_COUNTS
//    - Description: Raw count matrix produced by featureCounts; columns = samples, rows = peaks.

// 10. processed_counts.csv:
//     - Used in: DESEQ2_ANALYSIS → HEATMAP_VST_PLOT → PCA_PLOT
//     - Description: Reordered count matrix used for downstream statistical analysis and plots.

// 11. deseq2_results.csv:
//     - Used in: FILTER_SIGNIFICANT_PEAKS → FINAL_MERGE_SIGPEAKS_GENES → VOLCANO_PLOT
//     - Description: DESeq2 output of differential peak accessibility (log2FC, padj, etc).

// 12. significant_peaks.csv:
//     - Used in: CREATE_SIGNIFICANT_BED
//     - Description: Filtered peaks with padj < 0.001 and |log2FC| > 2.

// 13. significant_peaks.bed:
//     - Used in: SORT_PEAKS_BED → ADD_CHR_TO_PEAKS
//     - Description: BED-formatted version of significant peaks.

// 14. significant_peaks.chr.sorted.bed:
//     - Used in: BEDTOOLS_CLOSEST_ANNOTATION
//     - Description: Final BED-formatted peaks with chromosome prefix for annotation.

// 15. genes.sorted.gtf:
//     - Used in: BEDTOOLS_CLOSEST_ANNOTATION
//     - Description: GTF containing only gene entries, needed for peak-to-gene mapping.

// 16. peaks_with_nearest_genes.csv:
//     - Used in: FINAL_MERGE_SIGPEAKS_GENES
//     - Description: Output of bedtools closest; links each peak to its nearest gene.

// 17. final_cleaned_peaks.csv:
//     - Used in: VOLCANO_PLOT, ENRICHMENT_ANALYSIS
//     - Description: Merged table with peak statistics and gene annotations.

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// PROCESS DESCRIPTIONS

// FASTQC_BEFORE: Quality control check before trimming
// TRIM_GALORE: Adapter trimming and quality filtering
// FASTQC_AFTER: Quality control check after trimming
// ALIGN_BOWTIE2: Align reads to genome using Bowtie2
// MACS2_PEAKCALLING: Call peaks from aligned BAM files using MACS2
// MERGE_NARROWPEAK_FILES: Combine all narrowPeak files into one unified peak set
// CONVERT_BED_TO_SAF: Convert merged peaks to SAF format for featureCounts
// FEATURECOUNTS: Count reads overlapping peaks
// REORDER_COUNTS: Format and reorder the count matrix to match sample info
// DESEQ2_ANALYSIS: Perform differential accessibility analysis using DESeq2
// FILTER_SIGNIFICANT_PEAKS: Select significantly differentially accessible peaks
// CREATE_SIGNIFICANT_BED: Convert significant peaks to BED format for annotation
// SORT_PEAKS_BED: Sort significant peaks
// ADD_CHR_TO_PEAKS: Add 'chr' prefix to peaks if not present
// SORT_GTF: Sort GTF file for annotation
// EXTRACT_GENES_FROM_GTF: Filter GTF to retain only gene entries
// BEDTOOLS_CLOSEST_ANNOTATION: Annotate peaks with nearest genes using bedtools
// FINAL_MERGE_SIGPEAKS_GENES: Merge DESeq2 and gene annotation results
// VOLCANO_PLOT: Plot volcano plot of DESeq2 results with gene labels
// HEATMAP_VST_PLOT: Create heatmap of top variable peaks using VST transformation
// PCA_PLOT: Create PCA plot based on chromatin accessibility
// ENRICHMENT_ANALYSIS: Perform GO/KEGG/Reactome enrichment analysis on significant genes
// PLOT_ENRICHMENT_BAR_PLOTS: Generate barplots and dotplots of top enriched pathways
// IGV_SNAPSHOT: Generate IGV snapshot of peak regions using a batch script

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------





nextflow.enable.dsl=2
params.input = "/Users/ramya/Documents/academics/semester/bio-informatics/ATACsequencing/Input-files/*.fastq"
params.outdir = "/Users/ramya/Documents/academics/semester/bio-informatics/ATACsequencing/ATAC-seq-results"
params.genome_index = "/Users/ramya/Documents/academics/semester/bio-informatics/ATACsequencing/bowtie2_index/GRCh38_noalt_as"
params.genome_size = 'hs' 
params.sample_info = "/Users/ramya/Documents/academics/semester/bio-informatics/ATACsequencing/Bioinfo H/samples.txt"
params.gtf_file = "/Users/ramya/Documents/academics/semester/bio-informatics/ATACsequencing/ATAC-seq-results/manual/gencode.v43.annotation.gtf"



//  FASTQC (Before) 
process FASTQC_BEFORE {
    tag "$sample.baseName"

    input:
    path sample

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html
    path sample, emit: original

    script:
    """
    fastqc ${sample} -o .
    """
}

//  TRIM_GALORE 
process TRIM_GALORE {
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    tag "$sample.baseName"

    input:
    path sample

    output:
    path "*_trimmed.fq.gz", emit: trimmed

    script:
    """
    trim_galore ${sample} --gzip --output_dir .
    """
}

//  FASTQC (After) 
process FASTQC_AFTER {
    tag "$trimmed.baseName"

    input:
    path trimmed

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html
    path trimmed, emit: cleaned

    script:
    """
    fastqc ${trimmed} -o .
    """
}

//  GENOME ALIGNMENT 
process ALIGN_BOWTIE2 {
    tag "$sample.simpleName"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    memory '4 GB' 

    input:
        path sample

    output:
        path "*.bam", emit: bam
        path "*.log", emit: log

    script:
        def base = sample.simpleName.replaceFirst(/_trimmed$/, "")
        """
        echo "Running Bowtie2 on ${sample}" > ${base}_align.log

        bowtie2 -x ${params.genome_index} -U ${sample} --sensitive -p 4 2>> ${base}_align.log | \
        samtools sort -o ${base}.bam -
        samtools index ${base}.bam

        echo "Finished ${base}" >> ${base}_align.log
        """
}


////  PEAK CALLING 
process MACS2_PEAKCALLING {
    tag "$bam.baseName"
    publishDir "${params.outdir}/peaks", mode: 'copy'

    input:
    path bam

    output:
    path "${bam.baseName}_peaks_peaks.xls"
    path "${bam.baseName}_peaks_peaks.narrowPeak", emit: narrowPeak
    path "${bam.baseName}_peaks_summits.bed"  

    script:
    """
    macs2 callpeak -t ${bam} -f BAM -g ${params.genome_size} --nomodel -n ${bam.baseName}_peaks
    """
}

////  MERGING ALL NARROWPEAK FILES 
process MERGE_NARROWPEAK_FILES {
    publishDir "${params.outdir}/merged_narrowpeaks", mode: 'copy'

    input:
    path narrowPeak_files

    output:
    path "merged_peaks.bed", emit: merged_bed

    script:
    """
    cat ${narrowPeak_files.join(' ')} | sort -k1,1 -k2,2n | bedtools merge -i - > merged_peaks.bed
    """
}



////  BED TO SAF CONVERSION FOR FEATURE COUNTS 
process CONVERT_BED_TO_SAF {
    publishDir "${params.outdir}/merged_peaks", mode: 'copy'

    tag "$bed_file.baseName"

    input:
    path bed_file

    output:
    path "${bed_file.simpleName}.saf", emit: saf

    script:
    """
    awk 'BEGIN {OFS="\\t"; print "GeneID","Chr","Start","End","Strand"}
         {print "peak_"NR, \$1, \$2, \$3, "."}' ${bed_file} > ${bed_file.simpleName}.saf
    """
}


////  READ COUNTS 
process FEATURECOUNTS {
    publishDir "${params.outdir}/featureCounts", mode: 'copy'
    tag "featureCounts"

    input:
    path bam_files
    path saf_file

    output:
    path "read_counts.txt", emit: count_matrix

    script:
    """
    featureCounts -a ${saf_file} -F SAF -o read_counts.txt ${bam_files.join(' ')}
    """
}


////  REORDERING THE COLUMNS  
process REORDER_COUNTS {
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    input:
    path count_matrix

    output:
    path "processed_counts.csv", emit: reordered_counts

    script:
    """
    Rscript -e '
    counts <- read.delim("${count_matrix}", comment.char = "#", check.names = FALSE)
    counts <- counts[, c("Geneid", "ENCFF186CQZ_leukemia.bam", "ENCFF391BFJ_leukemia.bam", "ENCFF736ZWY.bam", "ENCFF411MVZ.bam")]
    rownames(counts) <- counts\$Geneid
    counts <- counts[, -1]
    write.csv(counts, file = "processed_counts.csv")
    '
    """
}


////  DIFFERENTIAL ANALYSIS 
process DESEQ2_ANALYSIS {
    publishDir "${params.outdir}/deseq2_results", mode: 'copy'

    conda 'bioconda::bioconductor-deseq2' 
    conda 'r-dplyr'

    input:
    path count_matrix
    path sample_info

    output:
    path "deseq2_results.csv", emit: deseq2_results

    script:
    """
  Rscript -e '
  library(DESeq2);
  library(dplyr);

  # Read counts
  counts <- read.csv("processed_counts.csv", row.names=1, check.names=FALSE)

  # Read sample info
  sample_info <- read.csv("samples.txt", row.names=1);
  sample_info <- sample_info[colnames(counts), , drop=FALSE];

  # DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = sample_info,
                                design = ~ Condition);
  dds <- DESeq(dds);

  res <- results(dds);

  write.csv(as.data.frame(res), file="deseq2_results.csv");
'

    """
}


////  FILTERING PEAKS padj<0.001 & log2FC>2 
process FILTER_SIGNIFICANT_PEAKS {
    publishDir "${params.outdir}/filtered_peaks", mode: 'copy'

    input:
    path deseq2_results

    output:
    path "significant_peaks.csv", emit: sig_peaks

    script:
    """
    Rscript -e '
      res <- read.csv("${deseq2_results}", row.names=1);

      # Filter significant peaks
      sig_res <- res[which(res\$padj < 0.001 & abs(res\$log2FoldChange) > 2), ];
      write.csv(sig_res, file="significant_peaks.csv", row.names=TRUE);
    '
    """
}


process CREATE_SIGNIFICANT_BED {
    publishDir "${params.outdir}/filtered_peaks", mode: 'copy'

    input:
    path sig_peaks
    path saf_file

    output:
    path "significant_peaks.bed", emit: sig_bed

    script:
    """
    cut -d',' -f1 ${sig_peaks} | tail -n +2 | sed 's/"//g' > sig_ids.txt

    #Grep matching entries from SAF
    grep -Ff sig_ids.txt ${saf_file} > tmp_grep_output.txt

    #Generating BED format
    awk 'BEGIN{OFS="\\t"} {print \$2, \$3, \$4, \$1}' tmp_grep_output.txt > significant_peaks.bed
    """

}


////  GTF SORTING 
process SORT_GTF {
    publishDir "${params.outdir}/genes", mode: 'copy'

    input:
    path gtf_file

    output:
    path "gencode.sorted.gtf", emit: sorted_gtf

    script:
    """
    sort -k1,1 -k4,4n ${gtf_file} > gencode.sorted.gtf
    """
}


process EXTRACT_GENES_FROM_GTF {
    publishDir "${params.outdir}/genes", mode: 'copy'

    input:
    path sorted_gtf

    output:
    path "genes.sorted.gtf", emit: genes_gtf

    script:
    """
    awk '\$3 == "gene"' ${sorted_gtf} > genes.sorted.gtf
    """
}


process ADD_CHR_TO_PEAKS {
    publishDir "${params.outdir}/filtered_peaks", mode: 'copy'

    input:
    path sorted_peaks

    output:
    path "significant_peaks.chr.sorted.bed", emit: chr_peaks

    script:
    """
    sed 's/^/chr/' ${sorted_peaks} > significant_peaks.chr.sorted.bed
    """
}

//////  SORTING SIGNIFICANT PEAKS 
process SORT_PEAKS_BED {
    publishDir "${params.outdir}/filtered_peaks", mode: 'copy'

    input:
    path peaks_bed

    output:
    path "significant_peaks.sorted.bed", emit: sorted_peaks

    script:
    """
    sort -k1,1 -k2,2n ${peaks_bed} > significant_peaks.sorted.bed
    """
}


////  GENE ANNOTATION 
process BEDTOOLS_CLOSEST_ANNOTATION {
    publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

    input:
    path chr_peaks
    path genes_gtf

    output:
    path "peaks_with_nearest_genes.csv", emit: annotated_peaks

    script:
    """
    bedtools closest -a ${chr_peaks} -b ${genes_gtf} -D a > peaks_with_nearest_genes.csv
    """
}




process FINAL_MERGE_SIGPEAKS_GENES {
    publishDir "${params.outdir}/final_table", mode: 'copy'

    input:
    path sig_peaks
    path peaks_with_genes

    output:
    path "final_cleaned_peaks.csv", emit: final_merged

    script:
    """
    Rscript -e '
    #Read significant peaks
    deseq <- read.csv("${sig_peaks}", row.names = 1)
    deseq\$PeakID <- rownames(deseq)

    #Read bedtools output
    genes <- read.table("${peaks_with_genes}", sep = "\\t", header = FALSE, stringsAsFactors = FALSE)
    #Extract necessary columns
    genes\$PeakID <- genes\$V4
    genes\$GeneInfo <- genes\$V13
    genes\$Distance <- genes\$V14

    #Extract gene name using string slicing
    genes\$Gene <- sapply(strsplit(genes\$GeneInfo, ";"), function(x) {
      gene_part <- x[grepl("gene_name", x)]
      if (length(gene_part) == 0) return(NA_character_)
      words <- strsplit(gene_part, " ")[[1]]
      gene <- gsub("\\\"", "", words[length(words)])
      return(gene)
    })

    gene_info <- genes[, c("PeakID", "Gene", "Distance")]
    #Merge with DESeq2 results
    merged <- merge(deseq, gene_info, by = "PeakID", all.x = TRUE)

    merged\$PeakNum <- as.numeric(sub("peak_", "", merged\$PeakID))
merged <- merged[order(merged\$PeakNum), ]
merged\$PeakNum <- NULL

    #Save result
    write.csv(merged, "final_cleaned_peaks.csv", row.names = FALSE)
    '
    """
}




////  VOLCANO PLOT 
process VOLCANO_PLOT {
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path final_output

    output:
    path "volcano_plot_with_gene_names_updated.pdf", emit: volcano_pdf
    path "volcano_plot.png", emit: volcano_png

    script:
    """
    
    Rscript -e '
    library(ggplot2)
    library(ggrepel)

    merged <- read.csv("${final_output}")
    merged <- merged[-log10(merged\$padj) <= 400, ]
    merged\$Significant <- ifelse(merged\$padj < 0.0001 & abs(merged\$log2FoldChange) > 2, "Yes", "No")

    significant_points <- merged[merged\$Significant == "Yes", ]
    top100 <- head(significant_points[order(significant_points\$padj), ], 100)

    # Volcano plot
    volcano_plot <- ggplot(merged, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
    theme_bw() +
    labs(title = "Volcano Plot: Differential Chromatin Accessibility",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-Value") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlim(c(-7, 7)) +
    ylim(c(0, 600)) +
    ggrepel::geom_text_repel(data = top100, aes(label = Gene), size = 2, max.overlaps = 100)

    ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
    ggsave("volcano_plot_with_gene_names_updated.pdf", plot = volcano_plot, width = 8, height = 6, dpi = 300)

    print(volcano_plot)
    '

    """
}

////  HEATMAP PLOT 
process HEATMAP_VST_PLOT {
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path processed_counts

    output:
    path "heatmap.png", emit: heatmap_png
    path "heatmap.pdf", emit: heatmap_pdf

    script:
    """
    Rscript -e '
    library(DESeq2)
    library(pheatmap)
    counts_raw <- read.csv("${processed_counts}", row.names = 1)
    sample_info <- data.frame(
      condition = c("Leukemia", "Leukemia", "Healthy", "Healthy"),
      row.names = colnames(counts_raw)
    )
    sample_info\$condition <- factor(sample_info\$condition)
    dds <- DESeqDataSetFromMatrix(countData = counts_raw,
                                  colData = sample_info,
                                  design = ~ condition)

    vsd <- vst(dds, blind = TRUE)
    vsd_mat <- assay(vsd)

    #top 100 variable peaks
    top_var <- head(order(rowVars(vsd_mat), decreasing = TRUE), 100)
    top_vsd <- vsd_mat[top_var, ]

    # Heatmap
    png("heatmap.png", width = 800, height = 600, res = 120)
    pheatmap(top_vsd,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             annotation_col = sample_info,
             main = "Top 100 Variable Peaks (VST)",
             scale = "row",
             angle_col = 0,
             fontsize_col = 4,
         fontsize_row = 4)
    dev.off()

    pdf("heatmap.pdf", width = 8, height = 6)
    pheatmap(top_vsd,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             annotation_col = sample_info,
             main = "Top 100 Variable Peaks (VST)",
             scale = "row",
             angle_col = 0,
             fontsize_col = 4,
         fontsize_row = 4)
    dev.off()
    '
    """
}


////  ENRICHMENT 

process ENRICHMENT_ANALYSIS {
    publishDir "${params.outdir}/pathway_analysis", mode: 'copy'
    conda 'r-base'

    input:
    path final_peaks_table

    output:
    path "enrichment_results.csv", emit: pathway_results
    path "enrichment_plot.pdf", optional: true, emit: pathway_plot_pdf
    path "enrichment_plot.png", optional: true, emit: pathway_plot_png
    path "enrichment_dotplot.pdf", optional: true, emit: dotplot_pdf
    path "enrichment_dotplot.png", optional: true, emit: dotplot_png
    path "kegg_dotplot.png", optional: true, emit: kegg_dotplot
    path "reactome_dotplot.png", optional: true, emit: reactome_dotplot
    path "log.txt", emit: enrichment_log

    script:
    """
    Rscript -e '
    if (!requireNamespace("gprofiler2", quietly = TRUE)) install.packages("gprofiler2", repos="https://cloud.r-project.org")
    if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos="https://cloud.r-project.org")

    library(gprofiler2)
    library(ggplot2)

    df <- read.csv("${final_peaks_table}")
    gene_list <- unique(na.omit(df\$Gene))
    gene_list <- gene_list[gene_list != "" & gene_list != "." & grepl("^[A-Za-z]", gene_list)]

    cat(paste("Number of genes submitted:", length(gene_list), "\\n"), file="log.txt")

    if (length(gene_list) >= 5) {
  
                       
    gene_list_named <- list(query_1 = gene_list)

    gost_res <- gost(query = gene_list_named,
                    organism = "hsapiens",
                    sources = c("GO:BP", "KEGG", "REAC"),
                    user_threshold = 0.05)

    res <- gost_res\$result

    if (!is.null(gost_res\$meta\$genes_metadata\$query_1\$intersections)) {
    intersections <- gost_res\$meta\$genes_metadata\$query_1\$intersections
    res\$intersection <- sapply(intersections, function(x) paste(x, collapse = ","))
    } else {
    res\$intersection <- NA  # Or you can skip this column entirely
    }


        if (!is.null(gost_res) && !is.null(gost_res\$result) && nrow(gost_res\$result) > 0) {
        res <- gost_res\$result

    # Convert all list-type columns to comma-separated strings
    res[] <- lapply(res, function(x) {
    if (is.list(x)) sapply(x, function(i) paste(i, collapse = ",")) else x
    })

    write.csv(res, "enrichment_results.csv", row.names = FALSE)

        top10 <- head(gost_res\$result, 10)
        top10\$log10_pval <- -log10(top10\$p_value)

        # General dot plot
        dot <- ggplot(top10, aes(x = log10_pval,
                                 y = reorder(term_name, log10_pval),
                                 size = intersection_size,
                                 color = source)) +
          geom_point(alpha = 0.7) +
          labs(title = "Top 10 Enriched Pathways",
               x = "-log10(p-value)",
               y = "Pathway",
               size = "Gene Count",
               color = "Source") +
          theme_minimal()

        ggsave("enrichment_dotplot.png", plot = dot, width = 10, height = 6)
        ggsave("enrichment_dotplot.pdf", plot = dot, width = 10, height = 6)

        # Bar plot
        p <- ggplot(top10, aes(x = reorder(term_name, -log10(p_value)),
                               y = -log10(p_value),
                               fill = source)) +
             geom_bar(stat = "identity") +
             coord_flip() +
             labs(title = "Top 10 Enriched Pathways",
                  x = "Pathway",
                  y = "-log10(p-value)") +
             theme_minimal()

        ggsave("enrichment_plot.pdf", plot = p, width = 10, height = 6)
        ggsave("enrichment_plot.png", plot = p, width = 10, height = 6)

        # KEGG-specific plot
        kegg <- subset(gost_res\$result, source == "KEGG")
        top_kegg <- head(kegg[order(kegg\$p_value), ], 10)
        if (nrow(top_kegg) > 0) {
          top_kegg\$log10_pval <- -log10(top_kegg\$p_value)
          kegg_plot <- ggplot(top_kegg, aes(x = log10_pval,
                                            y = reorder(term_name, log10_pval),
                                            size = intersection_size,
                                            color = source)) +
            geom_point(alpha = 0.7) +
            labs(title = "Top 10 Enriched KEGG Pathways",
                 x = "-log10(p-value)",
                 y = "Pathway",
                 size = "Gene Count",
                 color = "Source") +
            theme_minimal()
          ggsave("kegg_dotplot.png", plot = kegg_plot, width = 10, height = 6)
        }

        # Reactome-specific plot
        react <- subset(gost_res\$result, source == "REAC")
        top_react <- head(react[order(react\$p_value), ], 10)
        if (nrow(top_react) > 0) {
          top_react\$log10_pval <- -log10(top_react\$p_value)
          react_plot <- ggplot(top_react, aes(x = log10_pval,
                                              y = reorder(term_name, log10_pval),
                                              size = intersection_size,
                                              color = source)) +
            geom_point(alpha = 0.7) +
            labs(title = "Top 10 Enriched Reactome Pathways",
                 x = "-log10(p-value)",
                 y = "Pathway",
                 size = "Gene Count",
                 color = "Source") +
            theme_minimal()
          ggsave("reactome_dotplot.png", plot = react_plot, width = 10, height = 6)
        }

        cat("Enrichment successful.\\n", file="log.txt", append=TRUE)
      } else {
        write.csv(data.frame(Message = "No enrichment results found."), "enrichment_results.csv", row.names = FALSE)
        png("enrichment_plot.png")
        plot.new()
        text(0.5, 0.5, "No enrichment results", cex = 1.2)
        dev.off()
        cat("No enrichment results found.\\n", file="log.txt", append=TRUE)
      }
    } else {
      write.csv(data.frame(Message = "Too few genes for enrichment."), "enrichment_results.csv", row.names = FALSE)
      png("enrichment_plot.png")
      plot.new()
      text(0.5, 0.5, "Too few genes", cex = 1.2)
      dev.off()
      cat("Not enough valid gene names for enrichment.\\n", file="log.txt", append=TRUE)
    }
    '
    """
}

////  ENRICHMENT PLOTS 
process PLOT_ENRICHMENT_BAR_PLOTS {
    publishDir "${params.outdir}/pathway_analysis", mode: 'copy'
    conda 'r-base'

    input:
    path enrichment_results

    output:
    path "barplot_all.png", emit: barplot_all
    path "barplot_by_source.png", emit: barplot_by_source
    path "log_plotting.txt", emit: plot_log

    script:
    """
    Rscript -e '
    required_packages <- c("ggplot2", "tidyr", "dplyr")
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }

    library(ggplot2)
    library(dplyr)
    library(tidyr)

    df <- read.csv("${enrichment_results}", stringsAsFactors = FALSE)

    df\$p_value <- suppressWarnings(as.numeric(df\$p_value))
    df <- df[!is.na(df\$p_value), ]
    df\$log10p <- -log10(df\$p_value)
    df <- df[order(df\$log10p, decreasing = TRUE), ]

    top10 <- head(df, 10)
    bp1 <- ggplot(top10, aes(x = reorder(term_name, log10p), y = log10p, fill = source)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = "Top 10 Enriched Pathways", x = "Pathway", y = "-log10(p-value)") +
      theme_minimal(base_size = 10)
    ggsave("barplot_all.png", plot = bp1, width = 10, height = 6)
    top_by_source <- df %>%
      group_by(source) %>%
      slice_max(order_by = log10p, n = 10) %>%
      ungroup()

    bp2 <- ggplot(top_by_source, aes(x = reorder(term_name, log10p), y = log10p, fill = source)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      facet_wrap(~source, scales = "free_y") +
      labs(title = "Top 10 Enriched Pathways per Source", x = "Pathway", y = "-log10(p-value)") +
      theme_minimal(base_size = 8)
    ggsave("barplot_by_source.png", plot = bp2, width = 12, height = 6)

    cat("Plotting completed successfully.\\n", file = "log_plotting.txt")
    '
    """
}

////  NARROW PEAKS IGC SCREENSHOT 
process IGV_SNAPSHOT {
    publishDir "${params.outdir}/igv_snapshots", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path peak_files

    output:
    path "atac_snapshot.png"

    script:
    """
    mkdir -p snapshots

    cat <<EOF > igv_batch.bat
    new
    genome hg38
    snapshotDirectory snapshots
    EOF

    for peak in ${peak_files.join(' ')}; do
        echo "load \$peak" >> igv_batch.bat
        echo "expand" >> igv_batch.bat
    done

    echo "viewaspanel" >> igv_batch.bat
    echo "genomeview" >> igv_batch.bat
    echo "wait 10" >> igv_batch.bat
    echo "snapshot atac_snapshot.png" >> igv_batch.bat
    echo "exit" >> igv_batch.bat

    bash /Users/ramya/tools/igv/IGV_2.16.0/igv.sh -b igv_batch.bat

    cp snapshots/atac_snapshot.png .
    """
}


////  PCA PLOT 
process PCA_PLOT {
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path processed_counts

    output:
    path "pca_plot.pdf", emit: pca_pdf
    path "pca_plot.png", emit: pca_png

    script:
    """
    Rscript -e '
    library(DESeq2)
    library(ggplot2)

    # Load counts
    counts_raw <- read.csv("${processed_counts}", row.names = 1)

    sample_info <- data.frame(
      condition = c("Leukemia", "Leukemia", "Healthy", "Healthy"),
      row.names = colnames(counts_raw)
    )
    sample_info\$condition <- factor(sample_info\$condition)

    dds <- DESeqDataSetFromMatrix(countData = counts_raw, colData = sample_info, design = ~condition)
    vsd <- vst(dds, blind = TRUE)

    # PCA
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pca_data, "percentVar"))

    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
         geom_point(size = 4) +
         xlab(paste0("PC1: ", percentVar[1], "% variance")) +
         ylab(paste0("PC2: ", percentVar[2], "% variance")) +
         ggtitle("PCA of Chromatin Accessibility (VST)") +
         theme_minimal()

    ggsave("pca_plot.pdf", plot = p, width = 8, height = 6)
    ggsave("pca_plot.png", plot = p, width = 8, height = 6)
    '
    """
}




//  WORKFLOW 

workflow {
    samples = Channel.fromPath(params.input, checkIfExists: true)
   
    // Step 1: FastQC before trimming
    qc_raw = FASTQC_BEFORE(samples)

    // Step 2: Trim Galore
    trimmed = TRIM_GALORE(qc_raw.original)

    // Step 3: FastQC after trimming
    qc_clean = FASTQC_AFTER(trimmed.trimmed)

    // Step 4: Alignment
    aligned = ALIGN_BOWTIE2(qc_clean.cleaned)

    // Step 5: Peak calling
    peaks = MACS2_PEAKCALLING(aligned.bam)

    // Step 6: narrowPeak to saf Merge all bed files
    narrowPeak_files = Channel.fromPath("${params.outdir}/peaks/*.narrowPeak", checkIfExists: true)
    merged_bed = MERGE_NARROWPEAK_FILES(narrowPeak_files)
    saf_file = CONVERT_BED_TO_SAF(merged_bed)
    // Step 7: Read counts
    bam_files = Channel.fromPath("${params.outdir}/aligned/*.bam", checkIfExists: true)
    count_matrix = FEATURECOUNTS(bam_files.collect(), saf_file)
    reordered_counts = REORDER_COUNTS(count_matrix)

    // Step 8: Differential Analysis
    sample_info = file(params.sample_info)
    deseq2_results = DESEQ2_ANALYSIS(reordered_counts, sample_info)

    // Step 9 : Filtering Significant Peaks
    sig_peaks = FILTER_SIGNIFICANT_PEAKS(deseq2_results)

    // Step 10 : Significant peaks to BED Files
    sig_bed = CREATE_SIGNIFICANT_BED(sig_peaks, saf_file)
    sorted_peaks = SORT_PEAKS_BED(sig_bed)
    chr_peaks = ADD_CHR_TO_PEAKS(sorted_peaks)

    // Pre req steps for gene annotation
    gtf_file = file(params.gtf_file)
    sorted_gtf = SORT_GTF(gtf_file)
    genes_gtf = EXTRACT_GENES_FROM_GTF(sorted_gtf)
    // Step 11 : Gene Annotation  
    annotated_peaks = BEDTOOLS_CLOSEST_ANNOTATION(chr_peaks, genes_gtf)

    // Step 12 : Merge Peaks with nearest genes
    final_output = FINAL_MERGE_SIGPEAKS_GENES(sig_peaks, annotated_peaks)

    // PLOTS - [HEAT, VOLCANO, PCA, IGV Screenshot of narrow peak files]
    peak_inputs = Channel.fromPath("${params.outdir}/peaks/*_peaks_peaks.narrowPeak", checkIfExists: true).collect()
    IGV_SNAPSHOT(peak_inputs)    
    PCA_PLOT(reordered_counts)
    VOLCANO_PLOT(final_output)
    HEATMAP_VST_PLOT(reordered_counts)

    // Enrichment
    enrichment_outputs = ENRICHMENT_ANALYSIS(final_output)
    PLOT_ENRICHMENT_BAR_PLOTS(enrichment_outputs.pathway_results)

}
