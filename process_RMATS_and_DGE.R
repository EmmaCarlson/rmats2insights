############################################################
# Title: Pipeline for Downstream Analysis of rMATS Results
# Author: Emma Carlson
# Description:
#  Pipeline for rMATS downstream processing.
############################################################

#############################
# 1) Packages
#############################

# Unique, de-duplicated package list
libraries <- c(
  "argparse","patchwork","GenomicRanges","GenomicFeatures","ggplot2","ggprism",
  "dplyr","EnhancedVolcano","clusterProfiler","BSgenome","rtracklayer","VarCon",
  "IRanges","genomation","ggpubr","enrichplot","DESeq2","grDevices",
  "RColorBrewer","pheatmap","ggrepel","here","org.Hs.eg.db","tidyverse",
  "fgsea","reshape2","data.table"
)

# Install (if needed) and load packages
check_install_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (pkg %in% BiocManager::available()) BiocManager::install(pkg) else install.packages(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
invisible(lapply(unique(libraries), check_install_load))
message("Packages loaded.")

#############################
# 2) CLI Arguments
#############################

# Parse command line arguments 
parse_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument("-i","--input_dir", type="character", help="Path to rMATS input directory")
  parser$add_argument("-o","--output_dir", type="character", default="NONE", help="Path to output directory")
  parser$add_argument("-c","--controlCondition", type="character", default="Control", help="Control/reference condition label")
  parser$add_argument("-m","--caseCondition", type="character", default="Case", help="Case/mutant condition label")
  parser$add_argument("-f","--FDR_cutoff", type="numeric", default=0.05, help="FDR cutoff for DSE")
  parser$add_argument("-d","--ILD_cutoffs", type="character", default="0.05", help="ILD cutoff(s) comma-separated")
  parser$add_argument("-j","--IJC_cutoffs", type="character", default="3", help="IJC cutoff(s) comma-separated")
  parser$add_argument("-s","--SJC_cutoffs", type="character", default="3", help="SJC cutoff(s) comma-separated")
  parser$add_argument("-g","--genome_build", type="character", default="UNKNOWN", help="Genome build: hg38 or mm10")
  parser$add_argument("-id","--input_dir_DGE", type="character", default="NONE", help="Path for DGE count matrix + metadata")
  parser$add_argument("-fc","--fold_change_cutoff", type="numeric", default=1.5, help="Fold-change cutoff (absolute)")
  parser$add_argument("-mc","--min_read_count", type="numeric", default=20, help="Minimum total count for gene retention")
  parser$add_argument("-ap","--adj_pvalue_cutoff", type="numeric", default=0.05, help="Adjusted p-value cutoff for DGE")
  parser$parse_args()
}

# Read args + store as variables
opt <- parse_args()
input_dir <- opt$input_dir
output_dir <- ifelse(opt$output_dir == "NONE", file.path(input_dir, "output"), opt$output_dir)
controlCondition <- opt$controlCondition
caseCondition <- opt$caseCondition
FDR_cutoff <- opt$FDR_cutoff
ILD_list <- as.numeric(unlist(strsplit(opt$ILD_cutoffs, ",")))
IJC_list <- as.numeric(unlist(strsplit(opt$IJC_cutoffs, ",")))
SJC_list <- as.numeric(unlist(strsplit(opt$SJC_cutoffs, ",")))
genome_build <- opt$genome_build
input_dir_DGE <- opt$input_dir_DGE
fold_change_cutoff <- opt$fold_change_cutoff
min_read_count <- opt$min_read_count
adj_pvalue_cutoff <- opt$adj_pvalue_cutoff

# Prepare output directory and parameter log
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
message("Output directory: ", output_dir)

condsStr <- paste0(
  "This analysis was run with the following parameters:\n",
  "Conditions: ", controlCondition, ", ", caseCondition, "\n",
  "Reference Condition: ", controlCondition, "\n",
  "fold_change_cutoff: ", fold_change_cutoff, "\n",
  "Adjusted P-value: ", adj_pvalue_cutoff, "\n",
  "Minimum Read Cutoff: ", min_read_count, "\n",
  "FDR: ", FDR_cutoff, "\n",
  "ILD cutoffs: ", paste(ILD_list, collapse = ","), "\n",
  "IJC cutoffs: ", paste(IJC_list, collapse = ","), "\n",
  "SJC cutoffs: ", paste(SJC_list, collapse = ","), "\n",
  "Genome Build: ", genome_build
)
write(condsStr, file = file.path(output_dir, "analysis_parameters_info.txt"))

#############################
# 3) Utilities (I/O, plotting, helpers)
#############################

# Save ggplot in PNG + PDF (single definition; pass outdir when needed)
save_plot <- function(plot, filename, outdir = output_dir, width = 10, height = 8, dpi = 600) {
  ggsave(file.path(outdir, paste0(filename, ".png")), plot = plot, dpi = dpi, width = width, height = height)
  ggsave(file.path(outdir, paste0(filename, ".pdf")), plot = plot, dpi = dpi, width = width, height = height)
}

# Write table as TSV (Excel-friendly)
write_results <- function(data, filename) {
  write.table(data, file = here::here(output_dir, filename), sep = "\t", row.names = FALSE, quote = FALSE)
}

# Significance stars for p-values
get_significance <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

#############################
# 4) Enrichment helpers (GO/KEGG + plotting)
#############################

# Save bar/dot plots and CSV for enrichment results (GO/KEGG)
save_enrichment_plots <- function(results, direction, analysis_type, ont_type, output_path) {
  # Create enrichment subdir if needed
  enrich_dir <- file.path(output_path, "enrichment_plots")
  if (!dir.exists(enrich_dir)) dir.create(enrich_dir, recursive = TRUE)
  
  # For GO, simplify to reduce redundancy; for KEGG, plot directly
  if (ont_type != "KEGG") {
    SimResults <- clusterProfiler::simplify(
      results, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL
    )
    bar_plot <- barplot(SimResults, showCategory = 15)
    dot_plot <- dotplot(SimResults, showCategory = 15)
    write.csv(SimResults@result,
              file = file.path(enrich_dir, paste0(analysis_type,
                                                  ifelse(!is.null(direction), paste0("_", direction), ""),
                                                  ifelse(!is.null(ont_type), paste0("_", ont_type), ""),
                                                  "_results.csv")),
              row.names = FALSE)
  } else {
    bar_plot <- barplot(results, showCategory = 15)
    dot_plot <- dotplot(results, showCategory = 15)
    write.csv(results,
              file = file.path(enrich_dir, paste0(analysis_type,
                                                  ifelse(!is.null(direction), paste0("_", direction), ""),
                                                  ifelse(!is.null(ont_type), paste0("_", ont_type), ""),
                                                  "_results.csv")),
              row.names = FALSE)
  }
  
  # Save plots into enrichment dir
  save_plot(bar_plot, paste0("barPlot_", analysis_type,
                             ifelse(!is.null(direction), paste0("_", direction), ""),
                             ifelse(!is.null(ont_type), paste0("_", ont_type), "")),
            outdir = enrich_dir)
  save_plot(dot_plot, paste0("dotPlot_", analysis_type,
                             ifelse(!is.null(direction), paste0("_", direction), ""),
                             ifelse(!is.null(ont_type), paste0("_", ont_type), "")),
            outdir = enrich_dir)
}

# GO enrichment wrapper (takes ENTREZ IDs, org_db is OrgDb)
perform_GO_enrichment <- function(gene_list, direction, output_path, adj_pvalue_cutoff = 0.1, org_db) {
  GO_results <- enrichGO(
    gene            = gene_list,
    universe        = background,
    OrgDb           = org_db,
    keyType         = "ENTREZID",
    ont             = "ALL",
    pAdjustMethod   = "fdr",
    minGSSize       = 10,
    maxGSSize       = 2000,
    pvalueCutoff    = 1,
    qvalueCutoff    = 1,
    readable        = TRUE
  )
  
  if (!is.null(GO_results) && nrow(GO_results) > 0) {
    # Split by ontology and plot if significant at provided cutoff
    bpResults <- dplyr::filter(GO_results, ONTOLOGY == "BP", p.adjust < adj_pvalue_cutoff)
    if (nrow(bpResults) > 0) save_enrichment_plots(bpResults, direction, "GO", "BP", output_path)
    
    mfResults <- dplyr::filter(GO_results, ONTOLOGY == "MF", p.adjust < adj_pvalue_cutoff)
    if (nrow(mfResults) > 0) save_enrichment_plots(mfResults, direction, "GO", "MF", output_path)
    
    ccResults <- dplyr::filter(GO_results, ONTOLOGY == "CC", p.adjust < adj_pvalue_cutoff)
    if (nrow(ccResults) > 0) save_enrichment_plots(ccResults, direction, "GO", "CC", output_path)
    
    # Also write the full GO table for transparency
    write.csv(GO_results, file = file.path(output_path, paste0("GO_enrichment_results_", direction, ".csv")), row.names = FALSE)
  }
}

# KEGG enrichment wrapper (organism: "hsa" or "mmu")
perform_KEGG_enrichment <- function(gene_list, direction, output_path, adj_pvalue_cutoff = 0.1, organism_code) {
  KEGG_results <- enrichKEGG(
    gene              = gene_list,
    organism          = organism_code,
    keyType           = "kegg",
    pvalueCutoff      = adj_pvalue_cutoff,
    pAdjustMethod     = "BH",
    universe          = background,
    minGSSize         = 10,
    maxGSSize         = 2000,
    qvalueCutoff      = adj_pvalue_cutoff,
    use_internal_data = FALSE
  )
  if (!is.null(KEGG_results) && nrow(KEGG_results) > 0) {
    save_enrichment_plots(KEGG_results, direction, "KEGG", "KEGG", output_path)
  }
}

#############################
# 5) Genome annotations and resources
#############################

# Soft path to bundled resources (update as needed)
includedDataPath <- "/Users/ecarlson/Desktop/mastersProject_data/RMATS_analysis/raw_data/"

# Load build-specific packages and files
if (genome_build == "mm10") {
  add_libs <- c("TxDb.Mmusculus.UCSC.mm10.knownGene","BSgenome.Mmusculus.UCSC.mm10","org.Mm.eg.db")
  invisible(lapply(add_libs, check_install_load))
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  gtf_file <- rtracklayer::readGFF(list.files(file.path(includedDataPath, "mm10"), pattern = "gtf", full.names = TRUE))
  ce_df <- read.delim(list.files(file.path(includedDataPath, "mm10"), pattern = "ConstitutiveExons", full.names = TRUE))
  ci_df <- read.delim(list.files(file.path(includedDataPath, "mm10"), pattern = "ConstitutiveIntrons", full.names = TRUE))
  banchpoint_anno <- readBed(list.files(file.path(includedDataPath, "mm10"), pattern = "predictedBP_top", full.names = TRUE), zero.based = FALSE)
} else if (genome_build == "hg38") {
  add_libs <- c("TxDb.Hsapiens.UCSC.hg38.knownGene","BSgenome.Hsapiens.UCSC.hg38","org.Hs.eg.db")
  invisible(lapply(add_libs, check_install_load))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  gtf_file <- rtracklayer::readGFF(list.files(file.path(includedDataPath, "hg38"), pattern = "gtf", full.names = TRUE))
  ce_df <- read.delim(list.files(file.path(includedDataPath, "hg38"), pattern = "ConstitutiveExons", full.names = TRUE))
  ci_df <- read.delim(list.files(file.path(includedDataPath, "hg38"), pattern = "ConstitutiveIntrons", full.names = TRUE))
  banchpoint_anno <- readBed(list.files(file.path(includedDataPath, "hg38"), pattern = "predictedBP_top", full.names = TRUE), zero.based = FALSE)
} else {
  message("Genome build '", genome_build, "' not supported. Only filtering will be performed.")
}

# If build is supported, build GRanges resources + branchpoint table
if (genome_build %in% c("mm10","hg38")) {
  threeUTR <- threeUTRsByTranscript(txdb)
  fiveUTR  <- fiveUTRsByTranscript(txdb)
  
  ce_bed <- GRanges(seqnames = ce_df[,1], strand = ce_df[,7], ranges = IRanges(start = ce_df[,4], end = ce_df[,5]), names = ce_df[,9])
  ci_bed <- GRanges(seqnames = ci_df[,1], strand = ci_df[,7], ranges = IRanges(start = ci_df[,4], end = ci_df[,5]), names = ci_df[,9])
  
  # Normalize naming: keep original 'banchpoint_anno' for backward compatibility
  branchpoint_anno <- banchpoint_anno
  branchpoint_anno_df <- as.data.frame(branchpoint_anno)
  colnames(branchpoint_anno_df) <- c("Chrom","bp_start","bp_end","width","strand","bp_score","bp_pos")
}

#############################
# 6) Splice-site/PTC/branchpoint helpers
#############################

# Compute PTC/CE + splice site windows + MaxEntScan scores
PTC_CE_MaxEntScan_Calc <- function(modEventFile, event) {
  grange_eventFile <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                              ranges = IRanges(start = modEventFile[,6], end = modEventFile[,7]),
                              names  = modEventFile$GeneSymbol)
  
  grange_CE_overlap_equal <- subsetByOverlaps(grange_eventFile, ce_bed, type = "equal", maxgap = 1)
  
  # For SE/A3SS/A5SS/MXE, compute PTC status; skip for MXE explicitly
  if (event != "MXE") {
    grange_fiveUTR_overlap  <- subsetByOverlaps(grange_eventFile, fiveUTR)
    grange_threeUTR_overlap <- subsetByOverlaps(grange_eventFile, threeUTR)
    modEventFile$PTC_status <- ifelse(
      (modEventFile$GeneSymbol %in% mcols(grange_fiveUTR_overlap)$names |
         modEventFile$GeneSymbol %in% mcols(grange_threeUTR_overlap)$names),
      "AMBIGUOUS",
      ifelse(modEventFile[,8] %% 3 == 0, "nonPTC", "PTC")
    )
  }
  
  # CE status
  modEventFile$CE_status <- ifelse(modEventFile$GeneSymbol %in% mcols(grange_CE_overlap_equal)$names, "CE", "nonCE")
  
  # Define windows for MaxEntScan around 5' and 3' based on strand
  modEventFile$start_5prime <- ifelse(modEventFile$strand == "+", modEventFile[,6] - 6, modEventFile[,7] - 2)
  modEventFile$end_5prime   <- ifelse(modEventFile$strand == "+", modEventFile[,6] + 2, modEventFile[,7] + 6)
  modEventFile$start_3prime <- ifelse(modEventFile$strand == "+", modEventFile[,7] - 2, modEventFile[,6] - 20)
  modEventFile$end_3prime   <- ifelse(modEventFile$strand == "+", modEventFile[,7] + 20, modEventFile[,6] + 2)
  
  grange_5prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                           ranges  = IRanges(start = modEventFile$start_5prime, end = modEventFile$end_5prime),
                           names   = modEventFile$GeneID)
  grange_3prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                           ranges  = IRanges(start = modEventFile$start_3prime, end = modEventFile$end_3prime),
                           names   = modEventFile$GeneID)
  
  # Get sequence windows
  if (genome_build == "mm10") {
    modEventFile$seq_3prime <- getSeq(Mmusculus, grange_3prime, as.character = TRUE)
    modEventFile$seq_5prime <- getSeq(Mmusculus, grange_5prime, as.character = TRUE)
  } else {
    modEventFile$seq_3prime <- getSeq(Hsapiens, grange_3prime, as.character = TRUE)
    modEventFile$seq_5prime <- getSeq(Hsapiens, grange_5prime, as.character = TRUE)
  }
  
  # MaxEntScan scores
  modEventFile$maxEntScanScore_5prime <- unlist(lapply(modEventFile$seq_5prime, calculateMaxEntScanScore, ssType = 5))
  modEventFile$maxEntScanScore_3prime <- unlist(lapply(modEventFile$seq_3prime, calculateMaxEntScanScore, ssType = 3))
  
  modEventFile
}

# Branchpoint overlap + nearest BP per gene (by absolute distance)
branchPoint_calc <- function(modEventFile, event) {
  modEventFile <- modEventFile[!grepl("_", modEventFile$Chrom),]
  
  # Build event bed by event type
  if (event != "RI") {
    grange_bed <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                          ranges = IRanges(start = modEventFile$UpstreamIS, end = modEventFile$UpstreamIE),
                          names  = modEventFile$GeneSymbol)
  } else {
    grange_bed <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                          ranges = IRanges(start = modEventFile$IntronStart, end = modEventFile$IntronEnd),
                          names  = modEventFile$GeneSymbol)
  }
  
  event_bed_df <- as.data.frame(grange_bed)
  colnames(event_bed_df) <- c("Chrom","start","end","width","strand","GeneSymbol")
  
  # Overlap with predicted branchpoints
  hits <- findOverlaps(branchpoint_anno, grange_bed)
  hits <- as.data.frame(hits)
  
  bp_event_df <- data.frame()
  if (nrow(hits) > 0) {
    for (num in seq_len(nrow(hits))) {
      bp_pos_index    <- hits[num,]$queryHits
      event_bed_index <- hits[num,]$subjectHits
      bp_event_df <- rbind(bp_event_df,
                           merge(branchpoint_anno_df[bp_pos_index,],
                                 event_bed_df[event_bed_index,],
                                 by = c("Chrom","strand")))
    }
    bp_event_df$dist_to_bp <- bp_event_df$bp_pos - bp_event_df$bp_start
    # Keep nearest BP per gene
    bp_event_df <- bp_event_df %>%
      dplyr::group_by(GeneSymbol) %>%
      dplyr::slice_min(abs(dist_to_bp)) %>%
      dplyr::ungroup()
    modEventFile <- merge(modEventFile,
                          bp_event_df[,c("GeneSymbol","bp_start","bp_pos","dist_to_bp","bp_score")],
                          by = "GeneSymbol", all.x = TRUE)
  }
  
  modEventFile
}

#############################
# 7) Event filtering + per-event processors
#############################

# Filter rMATS event table by IJC/SJC/FDR/ILD
filterEventFile <- function(eventFile, event, ILD_cutoff, IJC_cutoff, SJC_cutoff) {
  
  # If counts are CSV-like strings, convert to per-group means
  if (is.character(eventFile$IJC_SAMPLE_1)) {
    eventFile$IJC_SAMPLE_1_mean <- unlist(lapply(type.convert(strsplit(eventFile$IJC_SAMPLE_1, ","), as.is = TRUE), mean))
    eventFile$SJC_SAMPLE_1_mean <- unlist(lapply(type.convert(strsplit(eventFile$SJC_SAMPLE_1, ","), as.is = TRUE), mean))
    eventFile$IJC_SAMPLE_2_mean <- unlist(lapply(type.convert(strsplit(eventFile$IJC_SAMPLE_2, ","), as.is = TRUE), mean))
    eventFile$SJC_SAMPLE_2_mean <- unlist(lapply(type.convert(strsplit(eventFile$SJC_SAMPLE_2, ","), as.is = TRUE), mean))
    eventFile$IncLevel1_mean    <- unlist(lapply(type.convert(strsplit(eventFile$IncLevel1, ","), as.is = TRUE), mean))
    eventFile$IncLevel2_mean    <- unlist(lapply(type.convert(strsplit(eventFile$IncLevel2, ","), as.is = TRUE), mean))
  }
  
  # Non-MXE events require either group to pass read cutoffs
  if (event != "MXE") {
    keep_idx <- (
      (eventFile$IJC_SAMPLE_1_mean > IJC_cutoff & eventFile$SJC_SAMPLE_1_mean > SJC_cutoff) |
        (eventFile$IJC_SAMPLE_2_mean > IJC_cutoff & eventFile$SJC_SAMPLE_2_mean > SJC_cutoff)
    ) & (eventFile$FDR < FDR_cutoff) & (abs(eventFile$IncLevelDifference) > ILD_cutoff)
  } else {
    # MXE special case: strand-aware filtering (kept your original logic)
    keep_idx <- ifelse(eventFile$strand == "+",
                       (eventFile$IJC_SAMPLE_1_mean > IJC_cutoff & eventFile$SJC_SAMPLE_1_mean > SJC_cutoff),
                       (eventFile$IJC_SAMPLE_2_mean > IJC_cutoff & eventFile$SJC_SAMPLE_2_mean > SJC_cutoff)) &
      (eventFile$FDR < FDR_cutoff) & (abs(eventFile$IncLevelDifference) > ILD_cutoff)
  }
  
  filtered <- eventFile[keep_idx, , drop = FALSE]
  
  # Assign Condition by sign of ILD based on event class
  if (event %in% c("SE","MXE")) {
    filtered$Condition <- ifelse(filtered$IncLevelDifference < -ILD_cutoff, controlCondition, caseCondition)
  } else {
    filtered$Condition <- ifelse(filtered$IncLevelDifference >  ILD_cutoff, controlCondition, caseCondition)
  }
  
  filtered
}

# === Per-event processors (kept your original math/columns) ===
process_se_mats <- function(eventFile, event) {
  ExonLength   <- eventFile[,7]  - eventFile[,6]
  UpstreamEL   <- (eventFile[,9] - eventFile[,8]) + 1
  UpstreamIL   <-  eventFile[,6] - eventFile[,9]
  DownstreamEL <- (eventFile[,11] - eventFile[,10]) + 1
  DownstreamIL <- (eventFile[,10] - eventFile[,7]) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4],
                             eventFile[,6], eventFile[,7], ExonLength,
                             eventFile[,8], eventFile[,9], UpstreamEL, eventFile[,9], eventFile[,6], UpstreamIL,
                             eventFile[,10], eventFile[,11], DownstreamEL, eventFile[,7], eventFile[,10], DownstreamIL,
                             eventFile[,20], eventFile[,23], eventFile$Condition)
  colnames(modEventFile) <- c("ID","GeneID","GeneSymbol","strand","Chrom","ExonStart","ExonEnd","ExonLength",
                              "UpstreamES","UpstreamEE","UpstreamEL","UpstreamIS","UpstreamIE","UpstreamIL",
                              "DownstreamES","DownstreamEE","DownstreamEL","DownstreamIS","DownstreamIE","DownstreamIL",
                              "FDR","ILD","Condition")
  
  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile, event)
  modEventFile <- branchPoint_calc(modEventFile, event)
  modEventFile
}

process_ri_mats <- function(eventFile, event) {
  IntronLength  <-  eventFile[,10] - eventFile[,9]
  UpstreamEL    <- (eventFile[,9]  - eventFile[,8])  + 1
  DownstreamEL  <- (eventFile[,11] - eventFile[,10]) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4],
                             eventFile[,9], eventFile[,10], IntronLength,
                             eventFile[,8], eventFile[,9], UpstreamEL, eventFile[,10], eventFile[,11], DownstreamEL,
                             eventFile[,20], eventFile[,23], eventFile$Condition)
  colnames(modEventFile) <- c("ID","GeneID","GeneSymbol","strand","Chrom","IntronStart","IntronEnd","IntronLength",
                              "UpstreamES","UpstreamEE","UpstreamEL","DownstreamES","DownstreamEE","DownstreamEL",
                              "FDR","ILD","Condition")
  
  modEventFile <- modEventFile[!grepl("_", modEventFile$Chrom),]
  
  # Splice windows for MaxEntScan
  modEventFile$start_5prime <- ifelse(modEventFile$strand == "+", modEventFile$DownstreamES - 6, modEventFile$UpstreamEE - 2)
  modEventFile$end_5prime   <- ifelse(modEventFile$strand == "+", modEventFile$DownstreamES + 2, modEventFile$UpstreamEE + 6)
  modEventFile$start_3prime <- ifelse(modEventFile$strand == "+", modEventFile$UpstreamEE - 2, modEventFile$DownstreamES - 20)
  modEventFile$end_3prime   <- ifelse(modEventFile$strand == "+", modEventFile$UpstreamEE + 20, modEventFile$DownstreamES + 2)
  
  grange_5prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                           ranges  = IRanges(start = modEventFile$start_5prime, end = modEventFile$end_5prime),
                           names   = modEventFile$GeneID)
  grange_3prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                           ranges  = IRanges(start = modEventFile$start_3prime, end = modEventFile$end_3prime),
                           names   = modEventFile$GeneID)
  
  if (genome_build == "mm10") {
    modEventFile$seq_3prime <- getSeq(Mmusculus, grange_3prime, as.character = TRUE)
    modEventFile$seq_5prime <- getSeq(Mmusculus, grange_5prime, as.character = TRUE)
  } else {
    modEventFile$seq_3prime <- getSeq(Hsapiens, grange_3prime, as.character = TRUE)
    modEventFile$seq_5prime <- getSeq(Hsapiens, grange_5prime, as.character = TRUE)
  }
  modEventFile$maxEntScanScore_5prime <- unlist(lapply(modEventFile$seq_5prime, calculateMaxEntScanScore, ssType = 5))
  modEventFile$maxEntScanScore_3prime <- unlist(lapply(modEventFile$seq_3prime, calculateMaxEntScanScore, ssType = 3))
  
  # CI overlap
  grange_eventFile <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                              ranges = IRanges(start = modEventFile[,6], end = modEventFile[,7]),
                              names  = modEventFile$GeneSymbol)
  grange_CI_overlap_equal <- subsetByOverlaps(grange_eventFile, ci_bed, type = "equal", maxgap = 1)
  modEventFile$CI_status <- ifelse(modEventFile$GeneSymbol %in% mcols(grange_CI_overlap_equal)$names, "CI", "nonCI")
  
  modEventFile <- branchPoint_calc(modEventFile, event)
  modEventFile
}

process_a3ss_mats <- function(eventFile, event) {
  ExonStart   <- ifelse((eventFile[,23] > 0), eventFile[,6], eventFile[,8])
  ExonEnd     <- ifelse((eventFile[,23] > 0), eventFile[,7], eventFile[,9])
  ExonLength  <- ExonEnd - ExonStart
  UpstreamIS  <- ifelse((ExonEnd > eventFile[,10]), eventFile[,11], ExonEnd)
  UpstreamIE  <- ifelse((ExonEnd > eventFile[,10]), ExonStart, eventFile[,10])
  UpstreamIL  <- UpstreamIE - UpstreamIS
  FlankingEL  <- (eventFile[,11] - eventFile[,10]) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4],
                             ExonStart, ExonEnd, ExonLength,
                             UpstreamIS, UpstreamIE, UpstreamIL, eventFile[,10], eventFile[,11], FlankingEL,
                             eventFile[,20], eventFile[,23], eventFile$Condition)
  
  colnames(modEventFile) <- c("ID","GeneID","GeneSymbol","strand","Chrom","ExonStart","ExonEnd","ExonLength",
                              "UpstreamIS","UpstreamIE","UpstreamIL","FlankingES","FlankingEE","FlankingEL",
                              "FDR","ILD","Condition")
  
  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile, event)
  modEventFile <- branchPoint_calc(modEventFile, event)
  modEventFile
}

process_a5ss_mats <- function(eventFile, event) {
  ExonStart   <- ifelse((eventFile[,23] > 0), eventFile[,6], eventFile[,8])
  ExonEnd     <- ifelse((eventFile[,23] > 0), eventFile[,7], eventFile[,9])
  ExonLength  <- ExonEnd - ExonStart
  UpstreamIS  <- ifelse((ExonEnd > eventFile[,10]), eventFile[,11], ExonEnd)
  UpstreamIE  <- ifelse((ExonEnd > eventFile[,10]), ExonStart, eventFile[,10])
  UpstreamIL  <- UpstreamIE - UpstreamIS
  FlankingEL  <- (eventFile[,11] - eventFile[,10]) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4],
                             ExonStart, ExonEnd, ExonLength,
                             UpstreamIS, UpstreamIE, UpstreamIL, eventFile[,10], eventFile[,11], FlankingEL,
                             eventFile[,20], eventFile[,23], eventFile$Condition)
  
  colnames(modEventFile) <- c("ID","GeneID","GeneSymbol","strand","Chrom","ExonStart","ExonEnd","ExonLength",
                              "UpstreamIS","UpstreamIE","UpstreamIL","FlankingES","FlankingEE","FlankingEL",
                              "FDR","ILD","Condition")
  
  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile, event)
  modEventFile <- branchPoint_calc(modEventFile, event)
  modEventFile
}

process_mxe_mats <- function(eventFile, event) {
  ExonInUseStart <- ifelse((eventFile[,5] == "+" & eventFile[,25] > 0), eventFile[,6],
                           ifelse((eventFile[,5] == "+" & eventFile[,25] < 0), eventFile[,8],
                                  ifelse((eventFile[,5] == "-" & eventFile[,25] > 0), eventFile[,8], eventFile[,6])))
  ExonInUseEnd <- ifelse((eventFile[,5] == "+" & eventFile[,25] > 0), eventFile[,7],
                         ifelse((eventFile[,5] == "+" & eventFile[,25] < 0), eventFile[,9],
                                ifelse((eventFile[,5] == "-" & eventFile[,25] > 0), eventFile[,9], eventFile[,7])))
  
  ExonInUseLength <- ExonInUseEnd - ExonInUseStart
  UpstreamEL <- (eventFile[,11] - eventFile[,10]) + 1
  UpstreamIS <- eventFile[,11]
  UpstreamIE <- ExonInUseStart
  UpstreamIL <- ExonInUseStart - eventFile[,11]
  DownstreamEL <- (eventFile[,13] - eventFile[,12]) + 1
  DownstreamIS <- ExonInUseEnd
  DownstreamIE <- eventFile[,12]
  DownstreamIL <- (eventFile[,12] - ExonInUseEnd) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4],
                             ExonInUseStart, ExonInUseEnd, ExonInUseLength,
                             eventFile[,10], eventFile[,11], UpstreamEL, UpstreamIS, UpstreamIE, UpstreamIL,
                             eventFile[,12], eventFile[,13], DownstreamEL, DownstreamIS, DownstreamIE, DownstreamIL,
                             eventFile[,22], eventFile[,25], eventFile$Condition)
  
  colnames(modEventFile) <- c("ID","GeneID","GeneSymbol","strand","Chrom","ExonInUseStart","ExonInUseEnd","ExonInUseLength",
                              "UpstreamES","UpstreamEE","UpstreamEL","UpstreamIS","UpstreamIE","UpstreamIL",
                              "DownstreamES","DownstreamEE","DownstreamEL","DownstreamIS","DownstreamIE","DownstreamIL",
                              "FDR","ILD","Condition")
  
  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile, event)
  modEventFile <- branchPoint_calc(modEventFile, event)
  modEventFile
}

#############################
# 8) Plotting helpers
#############################

# Condition colors used across plots
myColors <- c("cornflowerblue","red")
names(myColors) <- c(controlCondition, caseCondition)

# CE/PTC bar charts
plotEventBarChart <- function(modEventFile, output_dir, event) {
  barChart <- ggplot(as.data.frame(with(modEventFile, table(Condition, CE_status))),
                     aes(x = CE_status, y = Freq, fill = factor(Condition, names(myColors)))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = Freq), position = position_dodge(width = 0.9), vjust = -0.25, size = 20) +
    scale_fill_manual(name = "Condition", values = myColors) +
    labs(x = "CE Status", y = "Number of Events") +
    theme_prism() +
    theme(text = element_text(family = "Helvetica"),
          axis.line = element_line(colour = "black", linewidth = 5)) +
    theme(axis.text  = element_text(size = 70, color = "black"),
          axis.title = element_text(size = 70, face = "bold", color = "black"),
          legend.text  = element_text(size = 50, color = "black"),
          legend.title = element_text(size = 40, face = "bold", color = "black"),
          plot.title   = element_text(hjust = 0.5, size = 100, face = "bold", color = "black"))
  
  ggsave("CEstatus_bar.svg", device = "svg", plot = barChart, path = output_dir, width = 25, height = 30, dpi = 600)
  ggsave("CEstatus_bar.pdf", device = "pdf", plot = barChart, path = output_dir, width = 25, height = 30, dpi = 600)
  
  if (event != "MXE") {
    barChart2 <- ggplot(as.data.frame(with(modEventFile, table(Condition, PTC_status))),
                        aes(x = PTC_status, y = Freq, fill = factor(Condition, names(myColors)))) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_text(aes(label = Freq), position = position_dodge(width = 0.9), vjust = -0.25, size = 20) +
      scale_fill_manual(name = "Condition", values = myColors) +
      labs(x = "PTC Status", y = "Number of Events") +
      theme_prism() +
      theme(text = element_text(family = "Helvetica"),
            axis.line = element_line(colour = "black", linewidth = 5)) +
      theme(axis.text  = element_text(size = 70, color = "black"),
            axis.title = element_text(size = 70, face = "bold", color = "black"),
            legend.text  = element_text(size = 50, color = "black"),
            legend.title = element_text(size = 40, face = "bold", color = "black"),
            plot.title   = element_text(hjust = 0.5, size = 100, face = "bold", color = "black"))
    ggsave("PTCstatus_bar.svg", device = "svg", plot = barChart2, path = output_dir, width = 25, height = 30, dpi = 600)
    ggsave("PTCstatus_bar.pdf", device = "pdf", plot = barChart2, path = output_dir, width = 25, height = 30, dpi = 600)
  }
}

# Volcano for DSE (EnhancedVolcano handles -log10 transform internally; no manual ylim)
plotVolcano_events <- function(eventFile, eventFile_filtered, output_dir, event, FDR_cutoff, ILD_cutoff) {
  keyvals <- ifelse(
    (eventFile$ID %in% eventFile_filtered[eventFile_filtered$Condition == controlCondition, ]$ID), 'blue',
    ifelse((eventFile$ID %in% eventFile_filtered[eventFile_filtered$Condition == caseCondition, ]$ID), 'red', 'grey')
  )
  keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == 'red']  <- caseCondition
  names(keyvals)[keyvals == 'blue'] <- controlCondition
  names(keyvals)[keyvals == 'grey'] <- 'Not Significant'
  
  volcanoPlot <- EnhancedVolcano(eventFile,
                                 lab = "",
                                 x = 'IncLevelDifference',
                                 y = 'FDR',
                                 title = '',
                                 subtitle = paste(controlCondition, ':', sum(keyvals == 'blue'),
                                                  ' -- ', caseCondition, ':', sum(keyvals == 'red')),
                                 caption  = paste("ILD cutoff =", ILD_cutoff, "; FDR cutoff =", FDR_cutoff),
                                 xlab = "Inclusion Level Difference",
                                 xlim = c(-1, 1),
                                 pointSize = 2.0,
                                 colCustom = keyvals,
                                 legendPosition = 'bottom',
                                 legendLabSize = 14,
                                 legendIconSize = 4.0,
                                 gridlines.major = TRUE,
                                 gridlines.minor = FALSE,
                                 border = 'partial',
                                 borderWidth = 1.2,
                                 borderColour = 'black',
                                 max.overlaps = Inf)
  ggsave(paste0(event, "events_volcano.svg"), device = "svg", plot = volcanoPlot, path = output_dir, width = 25, height = 30, dpi = 600)
  ggsave(paste0(event, "events_volcano.pdf"), device = "pdf", plot = volcanoPlot, path = output_dir, width = 25, height = 30, dpi = 600)
}

# Violin + density side-by-side for a numeric column
plotViolin_and_Density <- function(modEventFile, col_name, y_label, event, output_dir) {
  clean_col <- paste0(col_name, "_num")
  
  # Clean, numeric, and mild outlier trimming
  filtered_file <- modEventFile %>%
    dplyr::filter(!!rlang::sym(col_name) != "-") %>%
    dplyr::mutate(
      !!clean_col := abs(as.numeric(!!rlang::sym(col_name))),
      Condition   = factor(Condition, levels = c(caseCondition, controlCondition))
    ) %>%
    dplyr::filter(dplyr::between(
      !!rlang::sym(clean_col),
      stats::quantile(!!rlang::sym(clean_col), 0.25, na.rm = TRUE) - 1.5 * IQR(!!rlang::sym(clean_col), na.rm = TRUE),
      stats::quantile(!!rlang::sym(clean_col), 0.75, na.rm = TRUE) + 1.5 * IQR(!!rlang::sym(clean_col), na.rm = TRUE)
    ))
  
  # Require at least two per group
  if (dplyr::n_distinct(filtered_file$Condition) > 1 &&
      sum(filtered_file$Condition == caseCondition) >= 2 &&
      sum(filtered_file$Condition == controlCondition) >= 2) {
    
    t_test  <- t.test(as.formula(paste(clean_col, "~ Condition")), data = filtered_file)
    p_value <- t_test$p.value
    sig     <- get_significance(p_value)
    max_val <- max(filtered_file[[clean_col]], na.rm = TRUE)
    y_pos   <- max_val * 1.25
    
    mean_df <- filtered_file %>%
      dplyr::group_by(Condition) %>%
      dplyr::summarize(mean_val = if (col_name %in% c("ILD","BPscore","Log2FC","maxEntScanScore_3prime","maxEntScanScore_5prime")) {
        round(mean(!!rlang::sym(clean_col), na.rm = TRUE), 4)
      } else round(mean(!!rlang::sym(clean_col), na.rm = TRUE)))
    
    violin <- ggplot(filtered_file, aes(x = Condition, y = !!rlang::sym(clean_col))) +
      geom_violin(aes(fill = Condition), trim = FALSE, width = 0.6) +
      geom_segment(data = mean_df,
                   aes(x = as.numeric(Condition) - 0.1, xend = as.numeric(Condition) + 0.1,
                       y = mean_val, yend = mean_val),
                   color = "darkblue", linewidth = 2) +
      geom_text(data = mean_df,
                aes(x = Condition, y = mean_val + (max_val * 0.05), label = mean_val),
                color = "darkblue", size = 14/.pt, fontface = "bold") +
      ggpubr::geom_bracket(xmin = controlCondition, xmax = caseCondition, y.position = y_pos,
                           label = sprintf("P-val = %.2g (%s)", p_value, sig),
                           label.size = 14/.pt, size = 1.3) +
      theme_bw() +
      scale_fill_manual(name = "Condition", values = myColors) +
      theme(axis.text = element_text(size = rel(3)),
            axis.title = element_text(size = rel(3.5)),
            legend.position = "none") +
      ylab(y_label) +
      labs(title = paste0(event, " Events in ", controlCondition, ": ", sum(modEventFile$Condition == controlCondition),
                          "\n", event, " Events in ", caseCondition, ": ", sum(modEventFile$Condition == caseCondition)))
    
    density <- ggplot(filtered_file, aes(x = !!rlang::sym(clean_col), fill = Condition)) +
      geom_density(alpha = 0.6, color = NA) +
      geom_vline(data = mean_df, aes(xintercept = mean_val), color = "darkblue", linewidth = 1, linetype = "dashed") +
      geom_text(data = mean_df, aes(x = mean_val, y = 0, label = mean_val), color = "darkblue", size = 14/.pt, fontface = "bold", vjust = -1) +
      theme_bw() +
      scale_fill_manual(name = "Condition", values = myColors) +
      theme(axis.text = element_text(size = rel(3.5)),
            axis.title = element_text(size = rel(3.5)),
            legend.position = "none") +
      xlab(y_label) + ylab("")
    
    combined <- violin + density + patchwork::plot_layout(ncol = 2)
    ggsave(paste0(event, "_", col_name, "_ViolinDensity_plots.svg"), device = "svg", plot = combined, path = output_dir, width = 20, height = 8)
    ggsave(paste0(event, "_", col_name, "_ViolinDensity_plots.pdf"), device = "pdf", plot = combined, path = output_dir, width = 20, height = 8)
  }
}

# Summary bar charts for event/gene counts
plotSummaryBarChart <- function(count, events, condition, gene_count, output_dir) {
  df <- data.frame(count, events, condition)
  bar1 <- ggplot(df, aes(x = events, y = count, fill = factor(condition, names(myColors)))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.25, size = 20) +
    scale_fill_manual(name = "Condition", values = myColors) +
    labs(x = "Event Type", y = "Event Count") +
    theme_prism() +
    theme(text = element_text(family = "Helvetica-Bold"),
          axis.line = element_line(colour = "black", linewidth = 5),
          axis.text = element_text(size = 70, color = "black"),
          axis.title = element_text(size = 70, face = "bold", color = "black"),
          legend.text = element_text(size = 50, color = "black"),
          legend.title = element_text(size = 40, face = "bold", color = "black"),
          plot.title = element_text(hjust = 0.5, size = 100, face = "bold", color = "black"))
  ggsave("eventCount_bar.png", plot = bar1, path = output_dir, width = 25, height = 30, dpi = 600)
  
  df2 <- df; df2$gene_count <- gene_count
  bar2 <- ggplot(df2, aes(x = events, y = gene_count, fill = factor(condition, names(myColors)))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = gene_count), position = position_dodge(width = 0.9), vjust = -0.25, size = 20) +
    scale_fill_manual(name = "Condition", values = myColors) +
    labs(x = "Event Type", y = "Gene Count") +
    theme_prism() +
    theme(text = element_text(family = "Helvetica-Bold"),
          axis.line = element_line(colour = "black", linewidth = 5),
          axis.text = element_text(size = 70, color = "black"),
          axis.title = element_text(size = 70, face = "bold", color = "black"),
          legend.text = element_text(size = 50, color = "black"),
          legend.title = element_text(size = 40, face = "bold", color = "black"),
          plot.title = element_text(hjust = 0.5, size = 100, face = "bold", color = "black"))
  ggsave("geneCount_bar.png", plot = bar2, path = output_dir, width = 25, height = 30, dpi = 600)
}

#############################
# 9) DGE (DESeq2) + Enrichment (genes)
#############################

output_dir_plots <- file.path(output_dir, "figures")
dir.create(output_dir_plots, recursive = TRUE, showWarnings = FALSE)

if (input_dir_DGE != "NONE") {
  countMatrixFile <- list.files(path = input_dir_DGE, pattern = "count|counts", ignore.case = TRUE, full.names = TRUE)
  sampleInfoFile  <- list.files(path = input_dir_DGE, pattern = "info|metadata|Info", ignore.case = TRUE, full.names = TRUE)
  
  if (length(countMatrixFile) != 1 || length(sampleInfoFile) != 1) {
    if (length(countMatrixFile) > 1) message("ERROR: Multiple count matrix files found (pattern '*count*').")
    if (length(sampleInfoFile)  > 1) message("ERROR: Multiple sample info files found (pattern '*info*|*metadata*').")
    if (length(countMatrixFile) == 0) message("ERROR: No count matrix file found.")
    if (length(sampleInfoFile)  == 0) message("ERROR: No sample info file found.")
    message("Skipping DGE analysis...")
  } else {
    message("Starting DGE analysis...")
    countMatrix <- read.delim(countMatrixFile, header = TRUE, row.names = 1, check.names = FALSE)
    sampleInfo  <- read.csv(sampleInfoFile, header = TRUE)
    colnames(sampleInfo) <- c("id","sample_names","condition","description")  # assumes four columns
    
    # Keep only the two conditions of interest
    sampleInfo <- sampleInfo[sampleInfo$condition %in% c(controlCondition, caseCondition), ]
    
    # Keep numeric columns only in countMatrix and match sample order
    countMatrix <- countMatrix %>% dplyr::select(where(is.numeric))
    countMatrix <- countMatrix[, colnames(countMatrix) %in% sampleInfo$sample_names, drop = FALSE]
    countMatrix <- countMatrix[, order(names(countMatrix)), drop = FALSE]
    sampleInfo  <- sampleInfo[order(sampleInfo$sample_names), ]
    
    # Run DESeq2 
    dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleInfo, design = ~ condition)
    keep <- rowSums(DESeq2::counts(dds)) >= min_read_count
    dds <- dds[keep,]
    dds$condition <- relevel(dds$condition, ref = controlCondition)
    dds <- DESeq(dds)
    
    res <- results(dds); res$geneName <- rownames(res)
    resSig <- res[!is.na(res$padj) & (res$padj < adj_pvalue_cutoff) & (abs(res$log2FoldChange) > log2(fold_change_cutoff)), ]
    
    # Up/Down lists
    upregulatedResSig   <- resSig[resSig$log2FoldChange >  log2(fold_change_cutoff), ]
    downregulatedResSig <- resSig[resSig$log2FoldChange < -log2(fold_change_cutoff), ]
    
    write_results(upregulatedResSig,   "upregulated_results.txt")
    write_results(downregulatedResSig, "downregulated_results.txt")
    write_results(resSig,              "significant_results.txt")
    write_results(res,                 "unfiltered_results.txt")
    
    # Correlation heatmap on normalized counts
    normalized_counts <- counts(dds, normalized = TRUE)
    sample_cor_matrix <- cor(normalized_counts)
    Samples <- pheatmap(sample_cor_matrix,
                        clustering_distance_rows = "correlation",
                        clustering_distance_cols = "correlation",
                        color = colorRampPalette(c("blue","white","red"))(20),
                        main = "Sample Correlation Heatmap")
    save_plot(Samples, "correlationPlot", outdir = output_dir_plots)
    
    # PCA
    vsd <- varianceStabilizingTransformation(dds)
    rv  <- matrixStats::rowVars(assay(vsd))
    select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
    pc <- prcomp(t(assay(vsd)[select, ]))
    condition <- sampleInfo$condition
    scores <- data.frame(pc$x, condition)
    samples <- sampleInfo$sample_names
    explained_variance <- summary(pc)$importance[2, ] * 100
    
    AllSamples <- ggplot(scores, aes(x = PC1, y = PC2, col = factor(condition))) +
      xlab(paste0("PC1 Variance Explained: ", round(explained_variance[1], 1), "%")) +
      ylab(paste0("PC2 Variance Explained: ", round(explained_variance[2], 1), "%")) +
      geom_point(size = 5) +
      ggrepel::geom_label_repel(aes(label = samples), show.legend = FALSE, max.overlaps = 50) +
      ggtitle("Principal Component Analysis") +
      scale_colour_brewer(name = " ", palette = "Set1") +
      theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
            legend.position = "right",
            legend.key = element_rect(fill = 'NA'),
            legend.text = element_text(size = 10, face = "bold"),
            axis.text  = element_text(colour = "black"),
            axis.title = element_text(face = "bold"),
            panel.background = element_rect(color = 'black', fill = NA),
            axis.line  = element_line(linewidth = 0.5, color = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major  = element_line(color = "grey", linetype = "dashed"))
    save_plot(AllSamples, "pcaPlot", outdir = output_dir_plots)
    
    # Volcano for DEGs
    keyvals <- ifelse(
      (res$log2FoldChange < -log2(fold_change_cutoff) & res$padj < adj_pvalue_cutoff), 'blue',
      ifelse((res$log2FoldChange >  log2(fold_change_cutoff) & res$padj < adj_pvalue_cutoff), 'red','grey')
    )
    keyvals[is.na(keyvals)] <- 'grey'
    names(keyvals)[keyvals == 'red']  <- 'Up-regulated'
    names(keyvals)[keyvals == 'blue'] <- 'Down-regulated'
    names(keyvals)[keyvals == 'grey'] <- 'NS'
    
    DEG <- EnhancedVolcano(res,
                           lab = "",
                           x = 'log2FoldChange',
                           y = 'padj',
                           title = 'Differentially Expressed Genes',
                           subtitle = paste('Down:', sum(keyvals == 'blue'),
                                            '  Up:',   sum(keyvals == 'red')),
                           caption = paste("FC cutoff =", fold_change_cutoff, "; p-adj cutoff =", adj_pvalue_cutoff),
                           xlab = bquote(~Log[2]~ 'Fold Change'),
                           xlim = c(-4, 4),
                           pCutoff = 0.1,  # leave defaults sensible; you filter earlier
                           FCcutoff = 0.0,
                           pointSize = 2.0,
                           colCustom = keyvals,
                           legendPosition = 'bottom',
                           legendLabSize = 14,
                           legendIconSize = 4.0,
                           gridlines.major = TRUE,
                           gridlines.minor = FALSE,
                           border = 'partial',
                           borderWidth = 1.2,
                           borderColour = 'black',
                           max.overlaps = Inf) +
      theme(plot.subtitle = element_text(hjust = 0, size = 12))
    save_plot(DEG, "volcanoPlot", outdir = output_dir_plots)
    
    # Enrichment run, currently only when genome build supported
    if (genome_build %in% c("mm10","hg38")) {
      background <- gsub("\\..*", "", rownames(res)); background <- background[!is.na(background)]
      listDN    <- gsub("\\..*", "", rownames(downregulatedResSig)); listDN <- listDN[!is.na(listDN)]
      listUP    <- gsub("\\..*", "", rownames(upregulatedResSig));   listUP <- listUP[!is.na(listUP)]
      
      if (genome_build == "mm10") {
        background_entrez <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")
        listDN_entrez     <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = listDN,  column = "ENTREZID", keytype = "ENSEMBL")
        listUP_entrez     <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = listUP,  column = "ENTREZID", keytype = "ENSEMBL")
        assign("background", background_entrez, inherits = TRUE)
        perform_GO_enrichment(listDN_entrez, "Down", output_dir, adj_pvalue_cutoff, org.Mm.eg.db)
        perform_KEGG_enrichment(listDN_entrez, "Down", output_dir, adj_pvalue_cutoff, "mmu")
        perform_GO_enrichment(listUP_entrez, "Up", output_dir, adj_pvalue_cutoff, org.Mm.eg.db)
        perform_KEGG_enrichment(listUP_entrez, "Up", output_dir, adj_pvalue_cutoff, "mmu")
      } else {
        background_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")
        listDN_entrez     <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = listDN,  column = "ENTREZID", keytype = "ENSEMBL")
        listUP_entrez     <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = listUP,  column = "ENTREZID", keytype = "ENSEMBL")
        assign("background", background_entrez, inherits = TRUE)
        perform_GO_enrichment(listDN_entrez, "Down", output_dir, adj_pvalue_cutoff, org.Hs.eg.db)
        perform_KEGG_enrichment(listDN_entrez, "Down", output_dir, adj_pvalue_cutoff, "hsa")
        perform_GO_enrichment(listUP_entrez, "Up", output_dir, adj_pvalue_cutoff, org.Hs.eg.db)
        perform_KEGG_enrichment(listUP_entrez, "Up", output_dir, adj_pvalue_cutoff, "hsa")
      }
    }
  }
}

#############################
# 10) Gene set enrichment using ILD (per-event)
#############################

# GO + GSEA using ILD ranking
GSEA_and_PE <- function(backgroundList, modEventFile, condition, output_dir, p_value_cutoff = 0.1) {
  # Define OrgDb from build
  if (genome_build == "hg38") {
    
    org_db <- org.Hs.eg.db
    
  } else if (genome_build == "mm10"){
    
    org_db <- org.Mm.eg.db
    
  }
  
  # GO over SYMBOL list
  geneList <- modEventFile[modEventFile$Condition == condition, ]$GeneSymbol
  GO_results <- enrichGO(
    gene = geneList,
    universe = backgroundList,
    OrgDb = org_db,
    keyType = "SYMBOL",
    ont = "ALL",
    pAdjustMethod = "fdr",
    minGSSize = 10,
    maxGSSize = 2000,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  # Write and plot GO if not empty results
  if (!is.null(GO_results) && nrow(GO_results) > 0) {
    bpResults <- dplyr::filter(GO_results, ONTOLOGY == "BP", p.adjust < p_value_cutoff)
    ccResults <- dplyr::filter(GO_results, ONTOLOGY == "CC", p.adjust < p_value_cutoff)
    mfResults <- dplyr::filter(GO_results, ONTOLOGY == "MF", p.adjust < p_value_cutoff)
    write.table(GO_results, file = file.path(output_dir, paste0("GOresultsUNFILTERED_", event, "_", condition, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
    GOresults <- dplyr::filter(GO_results, p.adjust < p_value_cutoff)
    write.table(GOresults, file = file.path(output_dir, paste0("GOresultsFILTERED_p_value_cutoff_", event, "_", condition, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
    if (nrow(bpResults) > 0) save_plot(dotplot(bpResults, showCategory = 15), paste0("BP_DotPlot_", event, "_", condition), outdir = output_dir)
    if (nrow(ccResults) > 0) save_plot(dotplot(ccResults, showCategory = 15), paste0("CC_DotPlot_", event, "_", condition), outdir = output_dir)
    if (nrow(mfResults) > 0) save_plot(dotplot(mfResults, showCategory = 15), paste0("MF_DotPlot_", event, "_", condition), outdir = output_dir)
  }
  
  # GSEA ranked by ILD (ENSEMBL IDs from GeneID)
  gene_list <- modEventFile[modEventFile$Condition == condition, ]$ILD
  names(gene_list) <- gsub("\\..*","", modEventFile[modEventFile$Condition == condition, ]$GeneID)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- stats::na.omit(gene_list)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  if (length(gene_list) > 1) {
    gse <- gseGO(geneList = gene_list,
                 ont = "ALL",
                 keyType = "ENSEMBL",
                 minGSSize = 3,
                 scoreType = "pos",
                 maxGSSize = 800,
                 pvalueCutoff = 1,
                 verbose = TRUE,
                 OrgDb = org_db,
                 pAdjustMethod = "fdr")
    if (!is.null(gse) && nrow(gse) > 0) {
      write.table(gse, file = file.path(output_dir, paste0("GSEAresultsUNFILTERED_", event, "_", condition, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
      sigGSE_res <- dplyr::filter(gse, p.adjust < p_value_cutoff)
      if (nrow(sigGSE_res) > 0) {
        write.table(sigGSE_res, file = file.path(output_dir, paste0("GSEAresultsFILTERED_", event, "_", condition, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
        save_plot(dotplot(sigGSE_res, showCategory = 15), paste0("GSEAresults_DotPlot_", event, "_", condition), outdir = output_dir)
      }
    }
  }
}

#############################
# 11) Main: process AS events for all cutoff combos
#############################

AS_events <- c("SE","RI")  # (keep MXE/A5SS/A3SS commented unless needed)
# AS_events <- c("SE","RI","MXE","A5SS","A3SS")

for (ILD_cutoff in ILD_list) {
  for (IJC_cutoff in IJC_list) {
    for (SJC_cutoff in SJC_list) {
      
      # Create a run label + init containers for summary bars
      run <- paste("run", "ILD", ILD_cutoff, "FDR", FDR_cutoff, "IJC", IJC_cutoff, "SJC", SJC_cutoff, sep = "_")
      count <- events <- condition <- gene_count <- c()
      
      for (event in AS_events) {
        inputFile <- list.files(paste0(input_dir, "/"), pattern = paste0(event, ".MATS.JC.txt"), full.names = TRUE)
        if (file.exists(inputFile)) {
          message("Processing: ", inputFile)
          eventFile <- read.delim(inputFile, header = TRUE, sep = "\t")
          
          # Filter by thresholds
          eventFile_filtered <- filterEventFile(eventFile, event, ILD_cutoff, IJC_cutoff, SJC_cutoff)
          message(event, " filtered; ", nrow(eventFile), " events originally; ",
                  nrow(eventFile_filtered), " events remaining.")
          
          # Update summary counters
          count <- append(count, c(sum(eventFile_filtered$Condition == controlCondition),
                                   sum(eventFile_filtered$Condition == caseCondition)))
          events <- append(events, rep(event, 2))
          condition <- append(condition, c(controlCondition, caseCondition))
          # rMATS JC files use 'geneSymbol'
          gene_count <- append(
            gene_count,
            c(length(unique(eventFile_filtered[eventFile_filtered$Condition == controlCondition, ]$geneSymbol)),
              length(unique(eventFile_filtered[eventFile_filtered$Condition == caseCondition, ]$geneSymbol)))
          )
          
          # Define output subdirs
          outputPath       <- file.path(output_dir, run, event)
          outputPath_plots <- file.path(output_dir, run, event, "figures")
          dir.create(outputPath, recursive = TRUE, showWarnings = FALSE)
          dir.create(outputPath_plots, recursive = TRUE, showWarnings = FALSE)
          
          # Write filtered rMATS event file
          write.table(eventFile_filtered,
                      file = file.path(outputPath, paste0(event, ".MATS.filtered.JC.txt")),
                      sep = "\t", row.names = FALSE, quote = FALSE)
          
          # If annotations loaded (hg38/mm10) and we have filtered events, process and visualize
          if (genome_build %in% c("mm10","hg38") && nrow(eventFile_filtered) > 0) {
            
            # Build annotated modEventFile per event type
            modEventFile <- switch(event,
                                   "SE"  = process_se_mats(eventFile_filtered, event),
                                   "RI"  = process_ri_mats(eventFile_filtered, event),
                                   "A3SS"= process_a3ss_mats(eventFile_filtered, event),
                                   "A5SS"= process_a5ss_mats(eventFile_filtered, event),
                                   "MXE" = process_mxe_mats(eventFile_filtered, event)
            )
            
            # Normalize Ensembl IDs and save annotated table
            modEventFile$GeneID <- gsub("\\..*","", modEventFile$GeneID)
            write.table(modEventFile,
                        file = file.path(outputPath, paste0(event, ".MATS.modified.JC.txt")),
                        quote = FALSE, row.names = FALSE, sep = "\t")
            
            # Background for GO over SYMBOLs in this event class
            gene_list_background <- eventFile$geneSymbol
            
            # GO + GSEA for both conditions
            GSEA_and_PE(gene_list_background, modEventFile, controlCondition, outputPath_plots)
            GSEA_and_PE(gene_list_background, modEventFile, caseCondition,    outputPath_plots)
            
            # Volcano (raw eventFile vs filtered)
            plotVolcano_events(eventFile, eventFile_filtered, outputPath_plots, event, FDR_cutoff, ILD_cutoff)
            
            # Common metrics
            plotViolin_and_Density(modEventFile, "ILD",                   "Inclusion Level Difference",  event, outputPath_plots)
            plotViolin_and_Density(modEventFile, "maxEntScanScore_5prime","5' Splice Site Strength",    event, outputPath_plots)
            plotViolin_and_Density(modEventFile, "maxEntScanScore_3prime","3' Splice Site Strength",    event, outputPath_plots)
            plotViolin_and_Density(modEventFile, "dist_to_bp",            "Distance to Branchpoint",    event, outputPath_plots)
            plotViolin_and_Density(modEventFile, "bp_score",              "Branchpoint Score",          event, outputPath_plots)
            
            # Event-specific lengths
            if (event == "RI") {
              plotViolin_and_Density(modEventFile, "IntronLength", "Intron Length", event, outputPath_plots)
            } else if (event == "MXE") {
              plotViolin_and_Density(modEventFile, "ExonInUseLength", "Exon Length",              event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "UpstreamIL",      "Upstream Intron Length",   event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "DownstreamIL",    "Downstream Intron Length", event, outputPath_plots)
            } else if (event == "SE") {
              plotViolin_and_Density(modEventFile, "ExonLength",      "Exon Length",              event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "UpstreamIL",      "Upstream Intron Length",   event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "DownstreamIL",    "Downstream Intron Length", event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "UpstreamEL",      "Upstream Exon Length",     event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "DownstreamEL",    "Downstream Exon Length",   event, outputPath_plots)
              plotEventBarChart(modEventFile, outputPath_plots, event)
            } else {
              plotViolin_and_Density(modEventFile, "UpstreamIL", "Upstream Intron Length", event, outputPath_plots)
            }
          }
        } else {
          message("Skipping (not found): ", paste0(event, ".MATS.JC.txt"))
        }
      } # end for event
      
      # Summary bars for this run
      plotSummaryBarChart(count, events, condition, gene_count, file.path(output_dir, run))
    }
  }
}

message("Complete.")
