
############# Import Libraries #############
# list of required packages
libraries <- c("argparse",
               "patchwork",
               "GenomicRanges",
               "GenomicFeatures",
               "ggplot2",
               "ggprism",
               "dplyr",
               "EnhancedVolcano",
               "clusterProfiler",
               "BSgenome",
               "rtracklayer",
               "VarCon",
               "IRanges",
               "genomation",
               "ggpubr",
               "clusterProfiler",
               "enrichplot",
               "ggplot2",
               "DESeq2",
               "EnhancedVolcano",
               "grDevices",
               "RColorBrewer",
               "pheatmap",
               "ggrepel",
               "dplyr",
               "argparse",
               "here",
               "org.Hs.eg.db",
               "tidyverse",
               "fgsea",
               "ggrepel",
               "reshape2",
               "here",
               "data.table",
               "grDevices", 
               "here")

# function to check, install, and load packages
check_install_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    # install the package if not installed
    if (package %in% BiocManager::available()) {
      BiocManager::install(package)
    } else {
      install.packages(package)
    }
  }
  # load the package
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

# apply the function to each package in the list
lapply(libraries, check_install_load)

save_plot <- function(plot, filename, output_dir) {
  png_filename <- paste0(filename, ".png")
  #  pdf_filename <- paste0(filename, ".pdf")
  
  # Save as PNG
  print(output_dir)
  ggsave(png_filename, plot = plot, path = output_dir, dpi = 600, device = 'png')
  
  # Save as PDF
  svg_filename <- paste0(filename, ".svg")
  ggsave(svg_filename, plot = plot, path = output_dir, dpi = 600, device = 'svg')
}

############# Function Definitions #############

# Function to load data from mm10 or hg38
load_required_data <- function(genome_build, includedDataPath) {
  
  if (genome_build == "mm10") {
    
    additonal_libraries <- c("TxDb.Mmusculus.UCSC.mm10.knownGene",
                             "BSgenome.Mmusculus.UCSC.mm10",
                             "org.Mm.eg.db")
    
    lapply(additonal_libraries, check_install_load)
    
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    gtf_file <- rtracklayer::readGFF(list.files(paste0(includedDataPath, "/mm10"), pattern = "gtf", full.names = TRUE))
    ce_df <- read.delim(list.files(paste0(includedDataPath, "/mm10"), pattern = "ConstitutiveExons", full.names = TRUE))
    ci_df <- read.delim(list.files(paste0(includedDataPath, "/mm10"), pattern = "ConstitutiveIntrons", full.names = TRUE))
    branchpoint_anno <- readBed(list.files(paste0(includedDataPath, "/mm10"), pattern = "predictedBP_top", full.names = TRUE), zero.based = FALSE)
    
  } else if (genome_build == "hg38") {
    
    additonal_libraries <- c("TxDb.Hsapiens.UCSC.hg38.knownGene",
                             'BSgenome.Hsapiens.UCSC.hg38',
                             "org.Hs.eg.db")
    
    lapply(additonal_libraries, check_install_load)
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    gtf_file <- rtracklayer::readGFF(list.files(paste0(includedDataPath, "/hg38"), pattern = "gtf", full.names = TRUE))
    ce_df <- read.delim(list.files(paste0(includedDataPath, "/hg38"), pattern = "ConstitutiveExons", full.names = TRUE))
    ci_df <- read.delim(list.files(paste0(includedDataPath, "/hg38"), pattern = "ConstitutiveIntrons", full.names = TRUE))
    branchpoint_anno <- readBed(list.files(paste0(includedDataPath, "/hg38"), pattern = "predictedBP_top", full.names = TRUE),  zero.based = FALSE)
    
    
  }
  
  threeUTR <- threeUTRsByTranscript(txdb)
  fiveUTR <- fiveUTRsByTranscript(txdb)
  
  ce_bed <- GRanges(seqnames = ce_df[,1], strand = ce_df[,7], ranges = IRanges(start = ce_df[,4], end = ce_df[,5]), names = ce_df[,9])
  ci_bed <- GRanges(seqnames = ci_df[,1], strand = ci_df[,7], ranges = IRanges(start = ci_df[,4], end = ci_df[,5]), names = ci_df[,9])
  
  branchpoint_anno_df <- as.data.frame(branchpoint_anno)
  colnames(branchpoint_anno_df) <- c("Chrom", "bp_start", "bp_end", "width", "strand", "bp_score", "bp_pos")
  
  return(list(branchpoint_anno_df, branchpoint_anno, threeUTR, fiveUTR, txdb, gtf_file, ce_bed, ci_bed))
  
}

# Function to filter the event file on ILD, FDR, IJC, and SJC
filterEventFile <- function(eventFile, eventType, ILD_cutoff, FDR_cutoff, IJC_cutoff, SJC_cutoff, controlCondition, caseCondition) {
  
  # checks if column is indeed a list a values
  if( class(eventFile$IJC_SAMPLE_1) == "character") {
    eventFile$IJC_SAMPLE_1_mean <- unlist(lapply(type.convert(strsplit(eventFile$IJC_SAMPLE_1, ","), as.is = TRUE), mean))
    eventFile$SJC_SAMPLE_1_mean <- unlist(lapply(type.convert(strsplit(eventFile$SJC_SAMPLE_1, ","), as.is = TRUE), mean))
    eventFile$IJC_SAMPLE_2_mean <- unlist(lapply(type.convert(strsplit(eventFile$IJC_SAMPLE_2, ","), as.is = TRUE), mean))
    eventFile$SJC_SAMPLE_2_mean <- unlist(lapply(type.convert(strsplit(eventFile$SJC_SAMPLE_2, ","), as.is = TRUE), mean))
    
    eventFile$IncLevel1_mean <- unlist(lapply(type.convert(strsplit(eventFile$IncLevel1, ","), as.is = TRUE), mean))
    eventFile$IncLevel2_mean <- unlist(lapply(type.convert(strsplit(eventFile$IncLevel2, ","), as.is = TRUE), mean))
  }
  
  if (eventType != "MXE") {
    eventFile_filtered <- eventFile[((eventFile$IJC_SAMPLE_1_mean > IJC_cutoff & eventFile$SJC_SAMPLE_1_mean > SJC_cutoff) | (eventFile$IJC_SAMPLE_2_mean > IJC_cutoff & eventFile$SJC_SAMPLE_2_mean > SJC_cutoff))
                                    & eventFile$FDR < FDR_cutoff & abs(eventFile$IncLevelDifference) > ILD_cutoff,]
  } else {
    # TO DO: Check if filtering for MXE is correct
    eventFile_filtered <- eventFile[ifelse(eventFile$strand=="+", (eventFile$IJC_SAMPLE_1_mean > IJC_cutoff & eventFile$SJC_SAMPLE_1_mean > SJC_cutoff),
                                           (eventFile$IJC_SAMPLE_2_mean > IJC_cutoff & eventFile$SJC_SAMPLE_2_mean > SJC_cutoff))
                                    & eventFile$FDR < FDR_cutoff & abs(eventFile$IncLevelDifference) > ILD_cutoff,]
  }
  
  if (eventType == "SE" | eventType == "MXE") {
    eventFile_filtered$Condition_Type <- ifelse(eventFile_filtered$IncLevelDifference < -ILD_cutoff, "control", "case")
    eventFile_filtered$Condition <- ifelse(eventFile_filtered$IncLevelDifference < -ILD_cutoff, controlCondition, caseCondition)
  } else {
    eventFile_filtered$Condition_Type <- ifelse(eventFile_filtered$IncLevelDifference > ILD_cutoff, "control", "case")
    eventFile_filtered$Condition <- ifelse(eventFile_filtered$IncLevelDifference > ILD_cutoff, controlCondition, caseCondition)
  }
  
  return(eventFile_filtered)
}

# Write results to Excel-compatible files with raw counts and Gene Symbol
write_results <- function(data, filename, output_dir) {
  write.table(data, file = here::here(output_dir, filename), sep = "\t", row.names = FALSE, quote = FALSE)
}

branchPoint_calc <- function(modEventFile, branchpoint_anno, branchpoint_anno_df, eventType) {
  
  # Remove nonstandard chromosomes
  modEventFile <- modEventFile[!grepl("_", modEventFile$Chrom), ]
  
  # Build GRanges for event regions depending on event type
  if (eventType != "RI") {
    grange_bed <- GRanges(
      seqnames = modEventFile$Chrom,
      strand   = modEventFile$strand,
      ranges   = IRanges(start = modEventFile$UpstreamIS, end = modEventFile$UpstreamIE),
      names    = modEventFile$ID       # <-- attach ID directly
    )
  } else {
    grange_bed <- GRanges(
      seqnames = modEventFile$Chrom,
      strand   = modEventFile$strand,
      ranges   = IRanges(start = modEventFile$IntronStart, end = modEventFile$IntronEnd),
      names    = modEventFile$ID
    )
  }
  
  # Convert event ranges to dataframe with ID included
  event_bed_df <- as.data.frame(grange_bed)
  colnames(event_bed_df) <- c("Chrom", "start", "end", "width", "strand", "ID")
  
  # Find overlaps
  hits <- as.data.frame(findOverlaps(branchpoint_anno, grange_bed))
  
  # If no overlaps at all, just return modEventFile with NA columns
  if (nrow(hits) == 0) {
    modEventFile$bp_start    <- NA
    modEventFile$bp_pos      <- NA
    modEventFile$dist_to_bp  <- NA
    modEventFile$bp_score    <- NA
    return(modEventFile)
  }
  
  # Build merged annotation table
  bp_event_df <- data.frame()
  for (i in 1:nrow(hits)) {
    
    bp_idx    <- hits$queryHits[i]
    event_idx <- hits$subjectHits[i]
    
    bp_row    <- branchpoint_anno_df[bp_idx, ]
    event_row <- event_bed_df[event_idx, c("Chrom", "strand", "ID")]
    
    merged <- merge(
      bp_row,
      event_row,
      by = c("Chrom", "strand"),
      all.x = TRUE
    )
    
    bp_event_df <- rbind(bp_event_df, merged)
  }
  
  # Compute distance to BP
  bp_event_df$dist_to_bp <- bp_event_df$bp_pos - bp_event_df$bp_start
  
  # Keep the closest BP hit for each event ID
  bp_event_df <- bp_event_df %>%
    group_by(ID) %>%
    slice_min(abs(dist_to_bp), with_ties = FALSE) %>%
    ungroup()
  
  # Merge back into modEventFile using ID (correct, stable mapping)
  modEventFile <- merge(
    modEventFile,
    bp_event_df[, c("ID", "bp_start", "bp_pos", "dist_to_bp", "bp_score")],
    by = "ID",
    all.x = TRUE
  )
  
  return(modEventFile)
}


PTC_CE_MaxEntScan_Calc <- function(modEventFile, ce_bed, fiveUTR, threeUTR, eventType, genome_build) {
  grange_eventFile <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                              ranges = IRanges(start = modEventFile[,6], end = modEventFile[,7]),
                              names = modEventFile$GeneSymbol)
  
  grange_CE_overlap_equal <- subsetByOverlaps(grange_eventFile, ce_bed, type = "equal", maxgap = 1)
  
  if (eventType != "MXE") {
    
    grange_fiveUTR_overlap <- subsetByOverlaps(grange_eventFile, fiveUTR)
    grange_threeUTR_overlap <- subsetByOverlaps(grange_eventFile, threeUTR)
    
    print("PTC status calculated...")
    modEventFile$PTC_status <- ifelse((modEventFile$GeneSymbol %in% mcols(grange_fiveUTR_overlap)$names |
                                         modEventFile$GeneSymbol %in% mcols(grange_threeUTR_overlap)$names),
                                      "AMBIGUOUS", ifelse(modEventFile[,8]%%3 == 0, "nonPTC", "PTC" ))
  }
  
  print("CE status calculated...")
  modEventFile$CE_status <- ifelse((modEventFile$GeneSymbol %in% mcols(grange_CE_overlap_equal)$names),
                                   "CE", "nonCE")
  
  modEventFile$start_5prime <- ifelse(modEventFile$strand == "+", modEventFile[,6] - 6, modEventFile[,7] - 2)
  modEventFile$end_5prime <- ifelse(modEventFile$strand == "+", modEventFile[,6] + 2, modEventFile[,7] + 6)
  modEventFile$start_3prime <- ifelse(modEventFile$strand == "+", modEventFile[,7] - 2, modEventFile[,6] - 20)
  modEventFile$end_3prime <- ifelse(modEventFile$strand == "+", modEventFile[,7] + 20, modEventFile[,6] + 2)
  
  grange_5prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand, ranges = IRanges(start = modEventFile$start_5prime, end = modEventFile$end_5prime), names = modEventFile$GeneID)
  grange_3prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand, ranges = IRanges(start = modEventFile$start_3prime, end = modEventFile$end_3prime), names = modEventFile$GeneID)
  
  if (genome_build == "mm10") {
    modEventFile$seq_3prime <- getSeq(Mmusculus, grange_3prime, as.character=TRUE)
    modEventFile$seq_5prime <- getSeq(Mmusculus, grange_5prime, as.character=TRUE)
  } else if  (genome_build == "hg38"){
    modEventFile$seq_3prime <- getSeq(Hsapiens, grange_3prime, as.character=TRUE)
    modEventFile$seq_5prime <- getSeq(Hsapiens, grange_5prime, as.character=TRUE)
  }
  
  print("5' maxEntScanScore calculated...")
  modEventFile$maxEntScanScore_5prime <- unlist(lapply(modEventFile$seq_5prime, calculateMaxEntScanScore, ssType = 5))
  print("3' maxEntScanScore calculated...")
  modEventFile$maxEntScanScore_3prime <- unlist(lapply(modEventFile$seq_3prime, calculateMaxEntScanScore, ssType = 3))
  
  
  
  return(modEventFile)
  
}
### Functions to process events ###
process_se_mats <- function(eventFile, full_annotation, genome_build, required_data_list) {
  
  ExonLength = eventFile[,7] - eventFile[,6]
  UpstreamEL = (eventFile[,9] - eventFile[,8]) + 1
  UpstreamIL = eventFile[,6] - eventFile[,9]
  DownstreamEL = (eventFile[,11] - eventFile[,10]) + 1
  DownstreamIL = (eventFile[,10] - eventFile[,7]) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4], eventFile[,6], eventFile[,7], ExonLength,
                             eventFile[,8], eventFile[,9], UpstreamEL, eventFile[,9], eventFile[,6], UpstreamIL,
                             eventFile[,10], eventFile[,11], DownstreamEL, eventFile[,7], eventFile[,10], DownstreamIL,
                             eventFile[,20], eventFile[,23], eventFile$Condition)
  print(dim(modEventFile))
  colnames(modEventFile) <- c("ID", "GeneID", "GeneSymbol", "strand", "Chrom", "ExonStart", "ExonEnd", "ExonLength",
                              "UpstreamES", "UpstreamEE", "UpstreamEL", "UpstreamIS", "UpstreamIE", "UpstreamIL",
                              "DownstreamES", "DownstreamEE", "DownstreamEL", "DownstreamIS", "DownstreamIE", "DownstreamIL",
                              "FDR", "ILD", "Condition")
  print("ModEvent - Before PTC")
  print(dim(modEventFile))
  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile, required_data_list[[7]], required_data_list[[3]], required_data_list[[4]], "SE", genome_build)
  print("ModEvent - After PTC")
  print(dim(modEventFile))
  modEventFile <- branchPoint_calc(modEventFile, required_data_list[[2]], required_data_list[[1]], "SE")
  print("ModEvent - After BP")
  print(dim(modEventFile))
  return(modEventFile)
  
}

process_ri_mats <- function(eventFile, full_annotation, genome_build, required_data_list) {
  
  IntronLength = eventFile[,10] - eventFile[,9]
  UpstreamEL = (eventFile[,9] - eventFile[,8]) + 1
  DownstreamEL = (eventFile[,11] - eventFile[,10]) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4], eventFile[,9], eventFile[,10], IntronLength,
                             eventFile[,8], eventFile[,9], UpstreamEL, eventFile[,10], eventFile[,11], DownstreamEL,
                             eventFile[,20], eventFile[,23], eventFile$Condition)
  
  colnames(modEventFile) <- c("ID", "GeneID", "GeneSymbol", "strand", "Chrom", "IntronStart", "IntronEnd", "IntronLength",
                              "UpstreamES", "UpstreamEE", "UpstreamEL", "DownstreamES", "DownstreamEE", "DownstreamEL",
                              "FDR", "ILD", "Condition")
  
  modEventFile <- modEventFile[!grepl("_", modEventFile$Chrom),] # removes random chromosomes
  
  if (isTRUE(full_annotation)) {
    
    ci_bed = required_data_list[[7]]
    
    modEventFile$start_5prime <- ifelse(modEventFile$strand == "+", modEventFile$DownstreamES - 6, modEventFile$UpstreamEE - 2)
    modEventFile$end_5prime <- ifelse(modEventFile$strand == "+", modEventFile$DownstreamES + 2, modEventFile$UpstreamEE + 6)
    modEventFile$start_3prime <- ifelse(modEventFile$strand == "+", modEventFile$UpstreamEE - 2, modEventFile$DownstreamES - 20)
    modEventFile$end_3prime <- ifelse(modEventFile$strand == "+", modEventFile$UpstreamEE + 20, modEventFile$DownstreamES + 2)
    
    grange_5prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand, ranges = IRanges(start = modEventFile$start_5prime, end = modEventFile$end_5prime), names = modEventFile$GeneID)
    grange_3prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand, ranges = IRanges(start = modEventFile$start_3prime, end = modEventFile$end_3prime), names = modEventFile$GeneID)
    
    if (genome_build == "mm10") {
      modEventFile$seq_3prime <- getSeq(Mmusculus, grange_3prime, as.character=TRUE)
      modEventFile$seq_5prime <- getSeq(Mmusculus, grange_5prime, as.character=TRUE)
    } else if (genome_build == "hg38") {
      modEventFile$seq_3prime <- getSeq(Hsapiens, grange_3prime, as.character=TRUE)
      modEventFile$seq_5prime <- getSeq(Hsapiens, grange_5prime, as.character=TRUE)
    }
    
    print("5' maxEntScanScore calculated...")
    modEventFile$maxEntScanScore_5prime <- unlist(lapply(modEventFile$seq_5prime, calculateMaxEntScanScore, ssType = 5))
    print("3' maxEntScanScore calculated...")
    modEventFile$maxEntScanScore_3prime <- unlist(lapply(modEventFile$seq_3prime, calculateMaxEntScanScore, ssType = 3))
    
    grange_eventFile <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                                ranges = IRanges(start = modEventFile[,6], end = modEventFile[,7]),
                                names = modEventFile$GeneSymbol)
    
    grange_CI_overlap_equal <- subsetByOverlaps(grange_eventFile, ci_bed, type = "equal", maxgap = 1)
    
    modEventFile$CI_status <- ifelse((modEventFile$GeneSymbol %in% mcols(grange_CI_overlap_equal)$names),
                                     "CI", "nonCI")
    print("HERE-BP_RI")
    
    modEventFile <- branchPoint_calc(modEventFile, required_data_list[[2]], required_data_list[[1]], "RI")
  }
  
  return(modEventFile)
  
}

process_a3ss_mats <- function(eventFile, full_annotation, genome_build, required_data_list) {
  
  ExonStart = ifelse((eventFile[,23] > 0),eventFile[,6], eventFile[,8])
  ExonEnd = ifelse((eventFile[,23] > 0), eventFile[,7], eventFile[,9])
  ExonLength = ExonEnd - ExonStart
  UpstreamIS = ifelse((ExonEnd > eventFile[,10]), eventFile[,11], ExonEnd)
  UpstreamIE = ifelse((ExonEnd > eventFile[,10]), ExonStart, eventFile[,10])
  UpstreamIL = UpstreamIE - UpstreamIS
  FlankingEL = (eventFile[,11] - eventFile[,10]) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4], ExonStart, ExonEnd, ExonLength,
                             UpstreamIS, UpstreamIE, UpstreamIL, eventFile[,10], eventFile[,11], FlankingEL,
                             eventFile[,20], eventFile[,23], eventFile$Condition)
  
  colnames(modEventFile) <- c("ID", "GeneID", "GeneSymbol", "strand", "Chrom", "ExonStart", "ExonEnd", "ExonLength",
                              "UpstreamIS", "UpstreamIE", "UpstreamIL", "FlankingES", "FlankingEE", "FlankingEL",
                              "FDR", "ILD", "Condition")
  
  if (isTRUE(full_annotation)) {
    modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile, required_data_list[[7]], required_data_list[[3]], required_data_list[[4]], "A3SS", genome_build)
    modEventFile <- branchPoint_calc(modEventFile, required_data_list[[2]], required_data_list[[1]], "A3SS")
  }
  
  
  return(modEventFile)
  
}

process_a5ss_mats <- function(eventFile, full_annotation, genome_build, required_data_list) {
  
  ExonStart = ifelse((eventFile[,23] > 0), eventFile[,6], eventFile[,8])
  ExonEnd = ifelse((eventFile[,23] > 0), eventFile[,7], eventFile[,9])
  ExonLength = ExonEnd - ExonStart
  UpstreamIS = ifelse((ExonEnd > eventFile[,10]), eventFile[,11], ExonEnd)
  UpstreamIE = ifelse((ExonEnd > eventFile[,10]), ExonStart, eventFile[,10])
  UpstreamIL = UpstreamIE - UpstreamIS
  FlankingEL = (eventFile[,11] - eventFile[,10]) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4], ExonStart, ExonEnd, ExonLength,
                             UpstreamIS, UpstreamIE, UpstreamIL, eventFile[,10], eventFile[,11], FlankingEL,
                             eventFile[,20], eventFile[,23], eventFile$Condition)
  
  colnames(modEventFile) <- c("ID", "GeneID", "GeneSymbol", "strand", "Chrom", "ExonStart", "ExonEnd", "ExonLength",
                              "UpstreamIS", "UpstreamIE", "UpstreamIL", "FlankingES", "FlankingEE", "FlankingEL",
                              "FDR", "ILD", "Condition")
  
  if (isTRUE(full_annotation)) {
    modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile, required_data_list[[7]], required_data_list[[3]], required_data_list[[4]], "A5SS", genome_build)
    modEventFile <- branchPoint_calc(modEventFile, required_data_list[[2]], required_data_list[[1]], "A5SS")
  }
  
  
  return(modEventFile)
  
}

process_mxe_mats <- function(eventFile, full_annotation, genome_build, required_data_list) {
  
  ExonInUseStart <- ifelse((eventFile[,5] == "+" & eventFile[,25] > 0), eventFile[,6],
                           ifelse((eventFile[,5] == "+" & eventFile[,25] < 0), eventFile[,8],
                                  ifelse((eventFile[,5] == "-" & eventFile[,25] > 0), eventFile[,8], eventFile[,6])))
  
  ExonInUseEnd <- ifelse((eventFile[,5] == "+" & eventFile[,25] > 0), eventFile[,7],
                         ifelse((eventFile[,5] == "+" & eventFile[,25] < 0), eventFile[,9],
                                ifelse((eventFile[,5] == "-" & eventFile[,25] > 0), eventFile[,9], eventFile[,7])))
  
  ExonInUseLength = ExonInUseEnd - ExonInUseStart
  UpstreamEL = (eventFile[,11] - eventFile[,10]) + 1
  UpstreamIS = eventFile[,11]
  UpstreamIE = ExonInUseStart
  UpstreamIL = ExonInUseStart - eventFile[,11]
  DownstreamEL = (eventFile[,13] - eventFile[,12]) + 1
  DownstreamIS = ExonInUseEnd
  DownstreamIE = eventFile[,12]
  DownstreamIL = (eventFile[,12] - ExonInUseEnd) + 1
  
  modEventFile <- data.frame(eventFile[,1], eventFile[,2], eventFile[,3], eventFile[,5], eventFile[,4], ExonInUseStart, ExonInUseEnd, ExonInUseLength,
                             eventFile[,10], eventFile[,11], UpstreamEL, UpstreamIS, UpstreamIE, UpstreamIL,
                             eventFile[,12], eventFile[,13], DownstreamEL, DownstreamIS, DownstreamIE, DownstreamIL,
                             eventFile[,22], eventFile[,25], eventFile$Condition)
  
  colnames(modEventFile) <- c("ID", "GeneID", "GeneSymbol", "strand", "Chrom", "ExonInUseStart", "ExonInUseEnd", "ExonInUseLength",
                              "UpstreamES", "UpstreamEE", "UpstreamEL", "UpstreamIS", "UpstreamIE", "UpstreamIL",
                              "DownstreamES", "DownstreamEE", "DownstreamEL", "DownstreamIS", "DownstreamIE", "DownstreamIL",
                              "FDR", "ILD", "Condition")
  
  if (isTRUE(full_annotation)) {
    modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile, required_data_list[[7]], required_data_list[[3]], required_data_list[[4]], "MXE", genome_build)
    modEventFile <- branchPoint_calc(modEventFile, required_data_list[[2]], required_data_list[[1]], "MXE")
  }
  
  return(modEventFile)
  
}

get_significance <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

saveOutputFiles <- function(eventFile_filtered, modEventFile_filtered, eventType, outputPath, outputPath_plots, controlCondition, caseCondition) {
  
  write.table(modEventFile_filtered, paste0(outputPath, "/", eventType, ".MATS.filtered.JC.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(modEventFile_filtered, paste0(outputPath, "/", eventType, ".MATS.modified.JC.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  write.table(modEventFile_filtered[modEventFile_filtered$Condition == "control",], paste0(outputPath, "/", eventType, ".MATS.CTL.JC.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(modEventFile_filtered[modEventFile_filtered$Condition == "case",], paste0(outputPath, "/", eventType, ".MATS.CASE.JC.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  print("HERE")
  print(eventType)
  plotViolin_and_Density(modEventFile_filtered, "ILD", "Inclusion Level Difference", eventType, outputPath_plots, controlCondition, caseCondition)
  plotViolin_and_Density(modEventFile_filtered, "maxEntScanScore_5prime", "5' Splice Site Strength Score", eventType, outputPath_plots, controlCondition, caseCondition)
  plotViolin_and_Density(modEventFile_filtered, "maxEntScanScore_3prime", "3' Splice Site Strength Score", eventType, outputPath_plots, controlCondition, caseCondition)
  plotViolin_and_Density(modEventFile_filtered, "dist_to_bp", "Distance to Branchpoint", eventType, outputPath_plots, controlCondition, caseCondition)
  plotViolin_and_Density(modEventFile_filtered, "bp_score", "Branchpoint Score", eventType, outputPath_plots, controlCondition, caseCondition)
  
  if ("log2FC" %in% colnames(modEventFile_filtered)) {
    
    plotViolin_and_Density(modEventFile_filtered, "log2FC", "Log2 Fold Change", eventType, outputPath_plots, controlCondition, caseCondition)
    
  }
  
  
  if (eventType != "RI") {
    plotEventBarChart(modEventFile_filtered, outputPath_plots, eventType, controlCondition, caseCondition)
  }
  
  if (eventType == "RI") {
    
    plotViolin_and_Density(modEventFile_filtered, "IntronLength", "Intron Length", eventType, outputPath_plots, controlCondition, caseCondition)
    
  } else if (eventType == "MXE") {
    
    plotViolin_and_Density(modEventFile_filtered, "ExonInUseLength", "Exon Length", eventType, outputPath_plots, controlCondition, caseCondition)
    plotViolin_and_Density(modEventFile_filtered, "UpstreamIL", " Upstream Intron Length", eventType, outputPath_plots, controlCondition, caseCondition)
    plotViolin_and_Density(modEventFile_filtered, "DownstreamIL", " Downstream Intron Length", eventType, outputPath_plots, controlCondition, caseCondition)
    
  } else if (eventType == "SE") {
    plotViolin_and_Density(modEventFile_filtered, "ExonLength", "Exon Length", eventType, outputPath_plots, controlCondition, caseCondition)
    plotViolin_and_Density(modEventFile_filtered, "UpstreamIL", " Upstream Intron Length", eventType, outputPath_plots, controlCondition, caseCondition)
    plotViolin_and_Density(modEventFile_filtered, "DownstreamIL", " Downstream Intron Length", eventType, outputPath_plots, controlCondition, caseCondition)
    plotViolin_and_Density(modEventFile_filtered, "UpstreamEL", " Upstream Exon Length", eventType, outputPath_plots, controlCondition, caseCondition)
    plotViolin_and_Density(modEventFile_filtered, "DownstreamEL", " Downstream Exon Length", eventType, outputPath_plots, controlCondition, caseCondition)
    
    # To be added:(BPscore, BPDist, Log2FC)
  } else {
    plotViolin_and_Density(modEventFile_filtered, "UpstreamIL", " Upstream Intron Length", eventType, outputPath_plots, controlCondition, caseCondition)
  }
}

plotViolin_and_Density <- function(modEventFile, col_name, y_label, event, output_dir, controlCondition, caseCondition) {
  
  myColors <- c("cornflowerblue", "red")
  names(myColors) <- c(controlCondition,caseCondition)
  
  clean_col <- paste0(col_name, "_num")
  
  filtered_file <- modEventFile %>%
    filter(!!sym(col_name) != "-") %>%
    mutate(
      !!clean_col := abs(as.numeric(!!sym(col_name))),
      Condition = factor(Condition, levels = c(caseCondition, controlCondition))
    ) %>%
    filter(between(
      !!sym(clean_col),
      quantile(!!sym(clean_col), 0.25, na.rm = TRUE) - 1.5 * IQR(!!sym(clean_col), na.rm = TRUE),
      quantile(!!sym(clean_col), 0.75, na.rm = TRUE) + 1.5 * IQR(!!sym(clean_col), na.rm = TRUE)
    ))
  if (n_distinct(filtered_file$Condition) > 1 & sum(filtered_file$Condition == caseCondition) >= 2 & sum(filtered_file$Condition == controlCondition) >=  2) {
    t_test <- t.test(as.formula(paste(clean_col, "~ Condition")), data = filtered_file)
    p_value <- t_test$p.value
    sig <- get_significance(p_value)
    
    # Calculate means and positions
    max_val <- max(filtered_file[[clean_col]], na.rm = TRUE)
    y_pos <- max_val * 1.25
    
    mean_df <- filtered_file %>%
      group_by(Condition) %>%
      summarize(mean_val = if (col_name %in% c("ILD", "BPscore", "log2FC", "maxEntScanScore_3prime", "maxEntScanScore_5prime")) {
        round(mean(!!sym(clean_col), na.rm = TRUE), 4)
      } else {
        round(mean(!!sym(clean_col), na.rm = TRUE))
      })
    
    
    violin <- ggplot(filtered_file, aes(x = Condition, y = !!sym(clean_col)), palette = c("cornflowerblue", "red")) +
      geom_violin(aes(fill = Condition), trim = FALSE, width = 0.6) +
      #scale_x_discrete(guide = guide_axis(n.dodge=3))+
      geom_segment(
        data = mean_df,
        aes(x = as.numeric(Condition) - 0.1, xend = as.numeric(Condition) + 0.1,
            y = mean_val, yend = mean_val),
        color = "darkblue", linewidth = 2) +
      geom_text(
        data = mean_df,
        aes(x = Condition, y = mean_val + (max_val * 0.05), label = mean_val),
        color = "darkblue", size = 14/.pt, fontface = "bold") +
      geom_bracket(
        xmin = controlCondition, xmax = caseCondition,
        y.position = y_pos,
        label = sprintf("P-val = %.2g (%s)", p_value, sig),
        label.size = 14/.pt, size = 1.3) +
      theme_bw() +
      scale_fill_manual(name="Condition", values = myColors) +
      theme(
        axis.text = element_text(size = rel(3)),
        axis.title = element_text(size = rel(3.5)),
        legend.position = "none") +
      ylab(y_label) + labs(title= paste0(event, " Events in ", controlCondition, ": ", sum(modEventFile$Condition == controlCondition),
                                         "\n", event, " Events in ", caseCondition, ": ", sum(modEventFile$Condition == caseCondition)))
    
    # Create density plot
    density <- ggplot(filtered_file, aes(x = !!sym(clean_col), fill = Condition), palette = c("cornflowerblue", "red")) +
      geom_density(alpha = 0.6, color = NA) +
      geom_vline(
        data = mean_df,
        aes(xintercept = mean_val),
        color = "darkblue", linewidth = 1, linetype = "dashed") +
      geom_text(
        data = mean_df,
        aes(x = mean_val, y = 0, label = mean_val),
        color = "darkblue", size = 14/.pt, fontface = "bold", vjust = -1) +
      theme_bw() +
      scale_fill_manual(name="Condition", values = myColors) +
      theme(
        axis.text = element_text(size = rel(3.5)),
        axis.title = element_text(size = rel(3.5)),
        legend.position = "none") +
      xlab(y_label) +
      ylab("")
    
    
    # Combine and save
    combined <- violin + density + plot_layout(ncol = 2)
    
    ggsave(paste0(event, "_", col_name, "_ViolinDensity_plots.svg"), device = "svg", plot = combined, path = output_dir, width = 20, height = 8)
    ggsave(paste0(event, "_", col_name, "_ViolinDensity_plots.pdf"), device = "pdf", plot = combined, path = output_dir, width = 20, height = 8)
  }
  
}

# Function to perform GO enrichment analysis
perform_GO_enrichment <- function(gene_list, direction, background, output_dir, adj_pvalue_cutoff = 0.1, org_type) {
  GO_results <- enrichGO(
    gene = gene_list,
    universe = background,
    OrgDb = org_type,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "fdr",
    minGSSize = 10,
    maxGSSize = 2000,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  #Change so GO_results is run once
  if (!is.null(GO_results)) {
    if (dim(GO_results)[1] > 0) {
      # Plot and save charts...
      write.csv(GO_results, file = paste0(output_dir, "/GO_enrichment_results_", direction, ".csv"), row.names = FALSE)
      bpResults <- filter(GO_results, ONTOLOGY == "BP", p.adjust < adj_pvalue_cutoff)
      if (dim(bpResults)[1] > 0) {
        save_enrichment_plots(bpResults, direction, "GO", "BP", output_dir)
      }
      mfResults <- filter(GO_results, ONTOLOGY == "MF", p.adjust < adj_pvalue_cutoff)
      if (dim(mfResults)[1] > 0) {
        save_enrichment_plots(mfResults, direction, "GO", "MF", output_dir)
      }
      ccResults <- filter(GO_results, ONTOLOGY == "CC", p.adjust < adj_pvalue_cutoff)
      if (dim(ccResults)[1] > 0) {
        save_enrichment_plots(ccResults, direction, "GO", "CC", output_dir)
      }
    }
  }
  
}

perform_KEGG_enrichment <- function(gene_list, direction, background, output_dir, adj_pvalue_cutoff = 0.1, org_type) {
  KEGG_results <- enrichKEGG(
    gene = gene_list,
    organism = org_type,
    keyType = "kegg",
    pvalueCutoff = adj_pvalue_cutoff,
    pAdjustMethod = "BH",
    universe = background,
    minGSSize = 10,
    maxGSSize = 2000,
    qvalueCutoff = adj_pvalue_cutoff,
    use_internal_data = FALSE
  )
  #print(KEGG_results)
  
  # Plot and save charts...
  if(!is.null(KEGG_results)) {
    if (dim(KEGG_results)[1] > 0) {
      save_enrichment_plots(KEGG_results, direction, "KEGG", "KEGG", output_dir) 
    }
  }
}

save_enrichment_plots <- function(results, direction, analysis_type, ont_type, output_dir) {
  #print(results)
  #print(direction)
  #print(analysis_type)
  #print(ont_type)
  if (ont_type != "KEGG") {
    SimResults <- clusterProfiler::simplify(results, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
    bar_plot <- barplot(SimResults, showCategory = 15)
    dot_plot <- dotplot(SimResults, showCategory = 15)
    write.csv(SimResults@result, file = paste0(output_dir, analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_results.csv"), row.names = FALSE)
  } else {
    bar_plot <- barplot(results, showCategory = 15)
    dot_plot <- dotplot(results, showCategory = 15)
    write.csv(results, file = paste0(output_dir, analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_results.csv"), row.names = FALSE)
  }
  
  save_plot(bar_plot, paste0("barPlot_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")), output_dir)
  save_plot(dot_plot, paste0("dotPlot_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")), output_dir)
}


### Figures

plotSummaryBarChart <- function(df, output_dir, controlCondition, caseCondition) {
  
  #df <- data.frame(count,events,condition)
  
  myColors <- c("cornflowerblue", "red")
  names(myColors) <- c(controlCondition,caseCondition)
  barChart <- ggplot(df, aes(x=events, y=count, fill=factor(condition, names(myColors)))) +
    geom_bar(stat="identity", position = position_dodge()) +
    geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.25, size = 20) +
    scale_fill_manual(name="Condition", values = myColors) +
    labs(x = "Event Type", y = "Event Count") +
    theme_prism() +
    theme(text = element_text(family = "Helvetica-Bold"), axis.line = element_line(colour = "black",
                                                                                   linewidth = 5, linetype = "solid")) +
    theme(axis.text=element_text(size=70, color = "black"),
          axis.title = element_text(size=70, face= "bold", color = "black"),
          legend.text = element_text(size = 50, color = "black"),
          legend.title = element_text(size = 40, face = "bold", color = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 100, face= "bold", color = "black"))
  
  ggsave(paste0("eventCount_bar.png"), plot = barChart, path = output_dir, width = 25, height = 30, dpi = 600)
  
  barChart <- ggplot(df, aes(x=events, y=gene_count, fill=factor(condition, names(myColors)))) +
    geom_bar(stat="identity", position = position_dodge()) +
    geom_text(aes(label=gene_count), position=position_dodge(width=0.9), vjust=-0.25, size = 20) +
    scale_fill_manual(name="Condition", values = myColors) +
    labs(x = "Event Type", y = "Gene Count") +
    theme_prism() +
    theme(text = element_text(family = "Helvetica-Bold"), axis.line = element_line(colour = "black",
                                                                                   linewidth = 5, linetype = "solid")) +
    theme(axis.text=element_text(size=70, color = "black"),
          axis.title = element_text(size=70, face= "bold", color = "black"),
          legend.text = element_text(size = 50, color = "black"),
          legend.title = element_text(size = 40, face = "bold", color = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 100, face= "bold", color = "black"))
  
  ggsave(paste0("geneCount_bar.png"), plot = barChart, path = output_dir, width = 25, height = 30, dpi = 600)
  
}

plotEventBarChart <- function(modEventFile, output_dir, event, controlCondition, caseCondition) {
  
  myColors <- c("cornflowerblue", "red")
  names(myColors) <- c(controlCondition,caseCondition)
  
  barChart <- ggplot(as.data.frame(with(modEventFile, table(Condition, CE_status))), aes(x=CE_status, y = Freq, fill=factor(Condition, names(myColors)))) +
    geom_bar(stat="identity", position = position_dodge()) +
    geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25, size = 20) +
    scale_fill_manual(name="Condition", values = myColors) +
    labs(x = "CE Status", y = "Number of Events") +
    theme_prism() +
    theme(text = element_text(family = "Helvetica"), axis.line = element_line(colour = "black",
                                                                              linewidth = 5, linetype = "solid")) +
    theme(axis.text=element_text(size=70, color = "black"),
          axis.title = element_text(size=70, face= "bold", color = "black"),
          legend.text = element_text(size = 50, color = "black"),
          legend.title = element_text(size = 40, face = "bold", color = "black")) +
    theme(plot.title = element_text(hjust = 0.5, size = 100, face= "bold", color = "black"))
  
  ggsave(paste0("CEstatus_bar.svg"), device = "svg", plot = barChart, path = output_dir, width = 25, height = 30, dpi = 600)
  ggsave(paste0("CEstatus_bar.pdf"), device = "pdf", plot = barChart, path = output_dir, width = 25, height = 30, dpi = 600)
  
  if (event != "MXE") {
    barChart <- ggplot(as.data.frame(with(modEventFile, table(Condition, PTC_status))), aes(x=PTC_status, y = Freq, fill=factor(Condition, names(myColors)))) +
      geom_bar(stat="identity", position = position_dodge()) +
      geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25, size = 20) +
      scale_fill_manual(name="Condition", values = myColors) +
      labs(x = "PTC Status", y = "Number of Events") +
      theme_prism() +
      theme(text = element_text(family = "Helvetica"), axis.line = element_line(colour = "black",
                                                                                linewidth = 5, linetype = "solid")) +
      theme(axis.text=element_text(size=70, color = "black"),
            axis.title = element_text(size=70, face= "bold", color = "black"),
            legend.text = element_text(size = 50, color = "black"),
            legend.title = element_text(size = 40, face = "bold", color = "black")) +
      theme(plot.title = element_text(hjust = 0.5, size = 100, face= "bold", color = "black"))
    
    ggsave(paste0("PTCstatus_bar.svg"), device = "svg", plot = barChart, path = output_dir, width = 25, height = 30, dpi = 600)
    ggsave(paste0("PTCstatus_bar.pdf"), device = "pdf", plot = barChart, path = output_dir, width = 25, height = 30, dpi = 600)
  }
  
}

###

run_DGEA <- function(input_dir_DGE, fc_cutoff, min_read_count, adjusted_p_value_cutoff, controlCondition, caseCondition, genome_build, output_dir, run_Enrichment) {
  
  countMatrixFile <- list.files(path = input_dir_DGE, pattern = "count|counts", ignore.case = TRUE, full.names = TRUE)
  sampleInfoFile <- list.files(path = input_dir_DGE, pattern = "info|metadata|Info", ignore.case = TRUE, full.names = TRUE)
  
  print(countMatrixFile)
  print(sampleInfoFile)
  
  if (length(countMatrixFile) != 1 | length(sampleInfoFile) != 1) {
    
    if (length(countMatrixFile) > 1) { 
      print("ERROR: Multiple files with pattern for count matrix found! (*count*) Please remove or rename files.")
    }
    if (length(sampleInfoFile) > 1) {
      print("ERROR: Multiple files with pattern for sample info found! (*info* or *metadata*) Please remove or rename files.")
    }
    if (length(countMatrixFile) == 0) {
      print("ERROR: NO files with pattern for count matrix found! (*count*) Please check that file exists.")
    }
    if (length(sampleInfoFile) == 0) {
      print("ERROR: NO files with pattern for sample info found! (*info* or *metadata*) Please check that file exists.")
    }
    
    print("Skipping DGE analysis ...")
    
  } else {
    
    print("Starting DGE analysis ...")
    
    countMatrix <- read.delim(countMatrixFile, header = TRUE, row.names = 1, check.names = FALSE)
    sampleInfo <- read.csv(sampleInfoFile, header = TRUE)
    
    colnames(sampleInfo) <- c("id", "sample_names", "condition", "description") # TO DO: Add check for file format
    
    sampleInfo <- sampleInfo[sampleInfo$condition == controlCondition | sampleInfo$condition == caseCondition,]
    print(dim(sampleInfo))
    
    print(colnames(countMatrix))
    countMatrix <- countMatrix %>% select(where(is.numeric)) # removes non-numerical columns
    
    countMatrix <- countMatrix[,colnames(countMatrix) %in% sampleInfo$sample_names]
    
    countMatrix <- countMatrix[ , order(names(countMatrix))]
    sampleInfo <- sampleInfo[order(sampleInfo$sample_names),]
    
    dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleInfo, design = ~ condition)
    
    keep <- rowSums(counts(dds)) >= min_read_count
    dds <- dds[keep,]
    
    dds$condition <- relevel(dds$condition, ref = controlCondition)
    
    dds <- DESeq(dds)
    
    res <- results(dds)
    res$geneName <- rownames(res)
    
    resSig <- res[!is.na(res$padj) & (res$padj < adjusted_p_value_cutoff) & (abs(res$log2FoldChange) > log2(fc_cutoff)), ]
    
    gene_status <- data.frame(res$geneName, res$log2FoldChange, res$padj)
    colnames(gene_status) <- c("geneName", "log2FC", "padj")
    
    gene_status$DGE_status <- ifelse((is.na(gene_status$padj) | gene_status$padj >= 0.05 | (abs(res$log2FoldChange) < log2(fc_cutoff))), "NS", ifelse(gene_status$log2FC > log2(fc_cutoff), "UP", "DOWN"))
    
    # Separate upregulated and downregulated genes
    upregulatedResSig <- resSig[resSig$log2FoldChange > log2(fc_cutoff), ]
    downregulatedResSig <- resSig[resSig$log2FoldChange < -log2(fc_cutoff), ]
    write_results(upregulatedResSig, "upregulated_results.txt", output_dir)
    print("WRITING RESULTS")
    write_results(downregulatedResSig, "downregulated_results.txt", output_dir)
    write_results(resSig, "significant_results.txt", output_dir)
    write_results(res, "unfiltered_results.txt", output_dir)
    
    # Access the normalized counts from DESeq2
    normalized_counts <- counts(dds, normalized = TRUE)
    
    # Calculate sample correlation matrix
    sample_cor_matrix <- cor(normalized_counts)
    
    # Create a heatmap using pheatmap
    Samples <- pheatmap(sample_cor_matrix,
                        clustering_distance_rows = "correlation",
                        clustering_distance_cols = "correlation",
                        color = colorRampPalette(c("blue", "white", "red"))(20),
                        main = "Sample Correlation Heatmap")
    
    save_plot(Samples, "correlationPlot", output_dir)
    
    # PCA
    vst_dds <- vst(dds)
    rv <- rowVars(assay(vst_dds))
    vsd <- varianceStabilizingTransformation(dds)
    select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
    pc <- prcomp(t(assay(vsd)[select, ]))
    condition <- sampleInfo$condition
    scores <- data.frame(pc$x, condition)
    samples <- sampleInfo$sample_names
    explained_variance <- summary(pc)$importance[2, ] * 100
    
    AllSamples <- ggplot(scores, aes(x = PC1, y = PC2, col = factor(condition))) +
      xlab(paste0("PC1 Variance Explained: ", explained_variance[1], "%")) +
      ylab(paste0("PC2 Variance Explained: ", explained_variance[2], "%")) +
      geom_point(size = 5) +
      geom_label_repel(aes(label = samples), show.legend = FALSE, max.overlaps = 50) +
      ggtitle("Principal Component Analysis") +
      scale_colour_brewer(name = " ", palette = "Set1") +
      theme(
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        legend.position = "right",
        legend.key = element_rect(fill = 'NA'),
        legend.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = 'bold'),
        panel.background = element_rect(color = 'black', fill = NA),
        axis.line.y = element_line(linewidth = 0.5, color = "black"),
        axis.line.x = element_line(linewidth = 0.5, color = "black"),
        axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),
        panel.grid.major.x = element_line(color = "grey", linetype = "dashed")
      )
    
    save_plot(AllSamples, "pcaPlot", output_dir)
    
    # Volcano plot of up and down regulated genes
    keyvals <- ifelse(
      (res$log2FoldChange < -log2(fc_cutoff) & res$padj < adjusted_p_value_cutoff), 'blue',
      ifelse((res$log2FoldChange > log2(fc_cutoff) & res$padj < adjusted_p_value_cutoff), 'red',
             'grey')
    )
    keyvals[is.na(keyvals)] <- 'grey'
    names(keyvals)[keyvals == 'red'] <- 'Up-regulated'
    names(keyvals)[keyvals == 'grey'] <- 'NS'
    names(keyvals)[keyvals == 'blue'] <- 'Down-regulated'
    
    # Set xlim and ylim dynamically
    DEG <- EnhancedVolcano(res,
                           lab = "",
                           x = 'log2FoldChange',
                           y = 'padj',
                           title = 'Differentially Expressed Genes',
                           subtitle = paste('Down-regulated Genes:', sum(keyvals == 'blue'),
                                            ' -- Up-regulated Genes:', sum(keyvals == 'red')),
                           caption = paste("FC cutoff =", fc_cutoff, "; p-adj cutoff =", adjusted_p_value_cutoff),
                           xlab = bquote(~Log[2]~ 'Fold Change'),
                           xlim = c(-4,4),
                           ylim = c(0,max(res$padj)),
                           pCutoff = 0.1,
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
                           borderColour = 'black', max.overlaps = Inf
    )
    
    DEG <- DEG + theme(plot.subtitle = element_text(hjust = 0, size = 12))
    
    save_plot(DEG, "volcanoPlot", output_dir)
    if (run_Enrichment)
    if (genome_build == "mm10" | genome_build == "hg38") {
      
      background <- gsub("\\..*","", rownames(res))
      
      background <- background[!is.na(background)]
      
      listDN <- gsub("\\..*","", rownames(downregulatedResSig))
      listDN <- listDN[!is.na(listDN)]
      
      
      listUP <- gsub("\\..*","", rownames(upregulatedResSig))
      listUP <- listUP[!is.na(listUP)]
      
      if (run_Enrichment) {
        if (genome_build == "mm10") {
          
          #background <- mapIds(org.Mm.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")
          listDN <- mapIds(org.Mm.eg.db, keys = listDN, column = "ENTREZID", keytype = "ENSEMBL")
          listUP <- mapIds(org.Mm.eg.db, keys = listUP, column = "ENTREZID", keytype = "ENSEMBL")
          
          # Main analysis for downregulated genes
          perform_GO_enrichment(listDN, "Down", background, output_dir, adjusted_p_value_cutoff, "org.Mm.eg.db")
          perform_KEGG_enrichment(listDN, "Down", background, output_dir, adjusted_p_value_cutoff, "mmu")
          
          # Define other input parameters as needed
          
          # Main analysis for upregulated genes
          perform_GO_enrichment(listUP, "Up", background, output_dir, adjusted_p_value_cutoff, "org.Mm.eg.db")
          perform_KEGG_enrichment(listUP, "Up", background, output_dir,adjusted_p_value_cutoff, "mmu")
          
        }
        if (genome_build == "hg38") {
          background <- mapIds(org.Hs.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")
          listDN <- mapIds(org.Hs.eg.db, keys = listDN, column = "ENTREZID", keytype = "ENSEMBL")
          listUP <- mapIds(org.Hs.eg.db, keys = listUP, column = "ENTREZID", keytype = "ENSEMBL")
          
          # Main analysis for downregulated genes
          perform_GO_enrichment(listDN, "Down", background, output_dir, adjusted_p_value_cutoff, org.Hs.eg.db)
          perform_KEGG_enrichment(listDN, "Down", background, output_dir, adjusted_p_value_cutoff, "hsa")
          
          # Define other input parameters as needed
          
          # Main analysis for upregulated genes
          perform_GO_enrichment(listUP, "Up", background, output_dir, adjusted_p_value_cutoff, org.Hs.eg.db)
          perform_KEGG_enrichment(listUP, "Up", background, output_dir,adjusted_p_value_cutoff, "hsa")
          
        }
      }
    }
    
  }
  
  return(gene_status)
  
}

####


# rmats2insights_run: function that executes main workflow of package
rmats2insights_run <- function(ILD_list = 0.05, IJC_list = 3, SJC_list = 3, FDR_list = 0.05, input_dir, 
                               output_dir = paste0(input_dir, "/output"), genome_build = "UNKNOWN",
                               controlCondition = "CTL", caseCondition = "MUT",
                               input_dir_DGE = "NONE", fc_cutoff = 1.5, min_read_count = 20, adjusted_p_value_cutoff = 0.05,
                               gtf_file = "NONE", CE_file = "NONE", CI_file = "NONE", branchpoint_anno = "NONE", which_txdb = "NONE", which_BSgenome = "NONE", run_Enrichment = FALSE) {
  
  # retrieves minimum value for each filtering variable list, to perform filtering and annotation
  # after filtering and annotation, filtered event files, modified event files, and data visualizations
  # are generated for each combination of filtering criteria
  
  ILD_min_cutoff <- min(ILD_list)
  FDR_min_cutoff <- min(FDR_list)
  IJC_min_cutoff <- min(IJC_list)
  SJC_min_cutoff <- min(SJC_list)
  
  input_files <- list.files(input_dir, pattern = ".MATS.JC.txt", full.names = TRUE)
  
  dir.create(output_dir, recursive = TRUE)
  
  
  if (input_dir_DGE != "NONE") {
    
    DGE_file <- run_DGEA(input_dir_DGE, fc_cutoff, min_read_count, adjusted_p_value_cutoff, controlCondition, caseCondition, genome_build, output_dir, run_Enrichment)
    print(DGE_file)
  }
  
  if (length(input_files) == 0) {
    
    print("NO EVENT FILES FOUND. INPUT DIRECTORY SHOULD CONTAIN FILES WITH NAMING PATTERN OF *.MATS.JC.txt\nEXITING FUNCTION")
    return(NULL)
    
  }
  print(paste0("FOUND EVENT FILES: ", input_files))
  
  # checks if genome build is specified. If it isn't, will only perform partial annotation. If it isn't mm10 or hg38, user must provide
  # additional data (TO BE IMPLEMENTED)
  if (genome_build == "UNKNOWN") {
    
    # if genome build is unspecified, will only modify file to include exon/intron length, upstream and downstream lengths
    full_annotation = FALSE
    required_data_list <- "NONE"
    
  } else if (genome_build == "hg38" | genome_build == "mm10") {
    full_annotation = TRUE
    #required_data_list <- load_required_data(genome_build, here("rawData")) #TO DO, MAKE SOFT PATH OR LOAD FROM PACKAGE
    required_data_list <- load_required_data(genome_build, "/Users/ecarlson/Desktop/mastersProject_data/RMATS_analysis/raw_data/")
    print("FULL ANNO")

    print("PASSED")
    
  } else {
    # TO DO: Run partially if TxDB and GTF provided, handle loading failures
    if (gtf_file != "NONE" & CE_file != "NONE" & CI_file != "NONE" & branchpoint_anno != "NONE" 
        & which_txdb != "NONE" & which_BSgenome != "NONE") {
      
      additonal_libraries <- c(which_txdb,
                               which_BSgenome)
      
      lapply(additonal_libraries, check_install_load)
      
      txdb <- get(which_txdb)
      gtf_file <- rtracklayer::readGFF(gtf_file)
      branchpoint_anno <- readBed(branchpoint_anno, zero.based = FALSE)
      ce_df <- read.delim(CE_file)
      ci_df <- read.delim(CI_file)
      
      threeUTR <- threeUTRsByTranscript(txdb)
      fiveUTR <- fiveUTRsByTranscript(txdb)
      
      ce_bed <- GRanges(seqnames = ce_df[,1], strand = ce_df[,7], ranges = IRanges(start = ce_df[,4], end = ce_df[,5]), names = ce_df[,9])
      ci_bed <- GRanges(seqnames = ci_df[,1], strand = ci_df[,7], ranges = IRanges(start = ci_df[,4], end = ci_df[,5]), names = ci_df[,9])
      
      branchpoint_anno_df <- as.data.frame(banchpoint_anno)
      colnames(branchpoint_anno_df) <- c("Chrom", "bp_start", "bp_end", "width", "strand", "bp_score", "bp_pos")
      
      full_annotation = list(branchpoint_anno_df, branchpoint_anno, threeUTR, fiveUTR, txdb, gtf_file, ce_bed, ci_bed)
      required_data_list <- "NONE"
      
    } else {
      print("Files missing - only basic annotation provided")
      full_annotation = FALSE
      required_data_list <- "NONE"
    }
    # USER NEEDS TO SPECIFY TXDB, PROVIDE GTF, CE/CI FILES, AND BRANCHPOINT ANNO
    # IF TXDB AND CE/CI NOT PROVIDED - SKIP CE/CI ANNOTATION
    # IF TXDB, BRANCHPOINT ANNO NOT PROVIDED - SKIP BP ANNOTATION
    # IF TXDB NOT PROVIDED, SKIP MAXENTSCAN
    # IF GTF NOT PROVIDED, SKIP GTF ANNOTATIONS
  }
  
  count <- c()
  events <- c()
  condition <- c()
  gene_count <- c()
  run_list <- c()
  
  for (file_n in input_files) {
    
    print(paste0("Starting analysis for ", file_n))
    
    eventType <- sub("\\..*$", "", basename(file_n))
    
    print(eventType)
    
    eventFile <- read.delim(file_n, header = TRUE, sep = "\t") # reads in event file
    
    eventFile_filtered <- filterEventFile(eventFile, eventType, ILD_min_cutoff, 
                                          FDR_min_cutoff, IJC_min_cutoff, SJC_min_cutoff,
                                          controlCondition, caseCondition)
    
    if (nrow(eventFile_filtered) > 0) {
      if (eventType == "SE") {
        modEventFile <- process_se_mats(eventFile_filtered, full_annotation, genome_build, required_data_list)
      } else if (eventType == "RI") {
        modEventFile <- process_ri_mats(eventFile_filtered, full_annotation, genome_build, required_data_list)
      } else if (eventType == "A3SS") {
        modEventFile <- process_a3ss_mats(eventFile_filtered, full_annotation, genome_build, required_data_list)
      } else if (eventType == "A5SS") {
        modEventFile <- process_a5ss_mats(eventFile_filtered, full_annotation, genome_build, required_data_list)
      } else if (eventType == "MXE") {
        modEventFile <- process_mxe_mats(eventFile_filtered, full_annotation, genome_build, required_data_list)
      }
      
      modEventFile$GeneID <- gsub("\\..*","", modEventFile$GeneID) # removes decimal at end of GeneID
      
      if (input_dir_DGE != "NONE") {
        print(dim(modEventFile))
        modEventFile <- merge(modEventFile, DGE_file, by.x = "GeneID", by.y = "geneName", all.x = TRUE)
        print(modEventFile)
        print(dim(modEventFile))
        
      }
      
      for (FDR_cutoff in FDR_list) {
        for (ILD_cutoff in ILD_list) {
          for (IJC_cutoff in IJC_list) {
            for (SJC_cutoff in SJC_list) {
              
              run <- paste("run", "ILD", ILD_cutoff, "FDR", FDR_cutoff, "IJC", IJC_cutoff, "SJC", SJC_cutoff, sep = "_")
              outputPath <- file.path(output_dir, run, eventType)
              dir.create(outputPath, recursive = TRUE)
              
              outputPath_plots <- file.path(output_dir, run, eventType, "figures") # creates output path to write output to
              dir.create(outputPath_plots, recursive = TRUE)
              
              eventFile_filtered_SUB <- filterEventFile(eventFile_filtered, eventType, ILD_cutoff, FDR_cutoff, IJC_cutoff, SJC_cutoff, controlCondition, caseCondition)
              print("Event File Filtered")
              print(dim(eventFile_filtered_SUB))
              modEventFile_SUB <- modEventFile[modEventFile$ID %in% eventFile_filtered_SUB$ID,]
              print("Mod Event File Filtered")
              print(dim(modEventFile_SUB))
              saveOutputFiles(eventFile_filtered_SUB, modEventFile_SUB, eventType, outputPath, outputPath_plots, controlCondition, caseCondition)
              
              count <- append(count, c(sum(eventFile_filtered_SUB$Condition_Type == "control"), sum(eventFile_filtered_SUB$Condition_Type == "case")))
              events <- append(events, rep(eventType, 2))
              condition <- append(condition, c(controlCondition, caseCondition))
              gene_count <- append(gene_count, c(length(unique(eventFile_filtered_SUB[eventFile_filtered_SUB$Condition_Type == "control",]$geneSymbol)),
                                                 length(unique(eventFile_filtered_SUB[eventFile_filtered_SUB$Condition_Type == "case",]$geneSymbol))))
              run_list <- append(run_list, rep(run, 2))
              
            }
          }
          
        }
      }
      
    }
    
  }
  summary_df <- data.frame(count, gene_count, events, condition, run_list)
  for (unique_run in unique(run_list))
  {
    summary_df_sub <- summary_df[summary_df$run_list == unique_run,]
    plotSummaryBarChart(summary_df_sub, file.path(output_dir, unique_run), controlCondition, caseCondition) 
  }
}
