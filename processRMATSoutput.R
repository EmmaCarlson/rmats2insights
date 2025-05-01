
##################################################### import libraries

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
               "ggpubr")

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

##################################################### parse arguments

parse_args <- function() {

  # only input path (directory to RMATS output) is required
  parser <- ArgumentParser()

  parser$add_argument("-i", "--input_dir", type = "character", help = "Path to the input directory")
  parser$add_argument("-o", "--output_dir", type = "character", help = "Path to the output directory", default = "NONE")

  parser$add_argument("-c", "--controlCondition", type = "character", help = "Name for control condition", default = "Control")
  parser$add_argument("-m", "--caseCondition", type = "character", help = "Name for case/mutant condition", default = "Case")

  parser$add_argument("-f", "--FDR_cutoff", type = "numeric", help = "FDR cutoff", default = 0.05)

  parser$add_argument("-d", "--ILD_cutoffs", type = "character", default = "0.05",
                      help = "Values for filtering on inclusion level difference, can be a single value or a comma-separated list of values")
  parser$add_argument("-j", "--IJC_cutoffs", type = "character", default = "3",
                      help = "Values for filtering on included juction count, can be a single value or a comma-separated list of values")
  parser$add_argument("-s", "--SJC_cutoffs", type = "character", default = "3",
                      help = "Values for filtering on skipped juction count, can be a single value or a comma-separated list of values")

  parser$add_argument("-g", "--genome_build", type = "character", help = "genome build - either hg38 or mm10", default = "mm10")

  args <- parser$parse_args()
  return(args)

}

opt <- parse_args()

input_dir <- opt$input_dir
output_dir <- opt$output_dir

controlCondition <- opt$controlCondition
caseCondition <- opt$caseCondition

FDR_cutoff <- opt$FDR_cutoff
ILD_list <- as.numeric(unlist(strsplit(opt$ILD_cutoffs, ",")))
IJC_list <- as.numeric(unlist(strsplit(opt$IJC_cutoffs, ",")))
SJC_list <- as.numeric(unlist(strsplit(opt$SJC_cutoffs, ",")))

genome_build <- opt$genome_build

# if no output directory is provided, create new output directory in input directory
if (output_dir == "NONE") {
  output_dir <- file.path(input_dir, "output")
}

##################################################### load required data
# TO DO: Allow for input of GTF, CE/CI, and genome for unsupported builds

includedDataPath <- "/Users/ecarlson/Desktop/mastersProject/RMATS_analysis/raw_data/" # path to data used for annotation for hg38 and mm10; TO DO: make soft path

if (genome_build == "mm10") {

  additonal_libraries <- c("TxDb.Mmusculus.UCSC.mm10.knownGene",
                           "BSgenome.Mmusculus.UCSC.mm10",
                            "org.Mm.eg.db")

  lapply(additonal_libraries, check_install_load)

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  gtf_file <- rtracklayer::readGFF(list.files(paste0(includedDataPath, "/mm10"), pattern = "gtf", full.names = TRUE))
  ce_df <- read.delim(list.files(paste0(includedDataPath, "/mm10"), pattern = "ConstitutiveExons", full.names = TRUE))
  ci_df <- read.delim(list.files(paste0(includedDataPath, "/mm10"), pattern = "ConstitutiveIntrons", full.names = TRUE))
  banchpoint_anno <- readBed(list.files(paste0(includedDataPath, "/mm10"), pattern = "predictedBP_top", full.names = TRUE), zero.based = FALSE)

} else if (genome_build == "hg38") {

  additonal_libraries <- c("TxDb.Hsapiens.UCSC.hg38.knownGene",
                           'BSgenome.Hsapiens.UCSC.hg38',
                           "org.Hs.eg.db")

  lapply(additonal_libraries, check_install_load)

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  gtf_file <- rtracklayer::readGFF(list.files(paste0(includedDataPath, "/hg38"), pattern = "gtf", full.names = TRUE))
  ce_df <- read.delim(list.files(paste0(includedDataPath, "/hg38"), pattern = "ConstitutiveExons", full.names = TRUE))
  ci_df <- read.delim(list.files(paste0(includedDataPath, "/hg38"), pattern = "ConstitutiveIntrons", full.names = TRUE))
  banchpoint_anno <- readBed(list.files(paste0(includedDataPath, "/hg38"), pattern = "predictedBP_top", full.names = TRUE),  zero.based = FALSE)


} else {
  print(paste0("Genome build ", genome_build, " is currently not supported. Only filtering will be performed."))
}

if (genome_build == "mm10" | genome_build == "hg38") {
  threeUTR <- threeUTRsByTranscript(txdb)
  fiveUTR <- fiveUTRsByTranscript(txdb)

  ce_bed <- GRanges(seqnames = ce_df[,1], strand = ce_df[,7], ranges = IRanges(start = ce_df[,4], end = ce_df[,5]), names = ce_df[,9])
  ci_bed <- GRanges(seqnames = ci_df[,1], strand = ci_df[,7], ranges = IRanges(start = ci_df[,4], end = ci_df[,5]), names = ci_df[,9])

  branchpoint_anno_df <- as.data.frame(banchpoint_anno)
  colnames(branchpoint_anno_df) <- c("Chrom", "bp_start", "bp_end", "width", "strand", "bp_score", "bp_pos")

}

PTC_CE_MaxEntScan_Calc <- function(modEventFile) {
  grange_eventFile <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                              ranges = IRanges(start = modEventFile[,6], end = modEventFile[,7]),
                              names = modEventFile$GeneSymbol)

  grange_CE_overlap_equal <- subsetByOverlaps(grange_eventFile, ce_bed, type = "equal", maxgap = 1)

  if (event != "MXE") {

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
  } else {
    modEventFile$seq_3prime <- getSeq(Hsapiens, grange_3prime, as.character=TRUE)
    modEventFile$seq_5prime <- getSeq(Hsapiens, grange_5prime, as.character=TRUE)
  }

  print("5' maxEntScanScore calculated...")
  modEventFile$maxEntScanScore_5prime <- unlist(lapply(modEventFile$seq_5prime, calculateMaxEntScanScore, ssType = 5))
  print("3' maxEntScanScore calculated...")
  modEventFile$maxEntScanScore_3prime <- unlist(lapply(modEventFile$seq_3prime, calculateMaxEntScanScore, ssType = 3))



  return(modEventFile)

}

plotEventBarChart <- function(modEventFile, outputPath, event) {

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

  ggsave(paste0("CEstatus_bar.svg"), device = "svg", plot = barChart, path = outputPath, width = 25, height = 30, dpi = 600)
  print("ici")
  ggsave(paste0("CEstatus_bar.pdf"), device = "pdf", plot = barChart, path = outputPath, width = 25, height = 30, dpi = 600)

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

    ggsave(paste0("PTCstatus_bar.svg"), device = "svg", plot = barChart, path = outputPath, width = 25, height = 30, dpi = 600)
    ggsave(paste0("PTCstatus_bar.pdf"), device = "pdf", plot = barChart, path = outputPath, width = 25, height = 30, dpi = 600)
  }

}

plotVolcano_events <- function(eventFile, eventFile_filtered, outputPath, event, FDR_cutoff, ILD_cutoff) {

  keyvals <- ifelse(
    (eventFile$ID %in% eventFile_filtered[eventFile_filtered$Condition == controlCondition,]$ID), 'blue',
    ifelse((eventFile$ID %in% eventFile_filtered[eventFile_filtered$Condition == caseCondition,]$ID), 'red',
           'grey'))
  keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == 'red'] <- caseCondition
  names(keyvals)[keyvals == 'grey'] <- 'Not Significant'
  names(keyvals)[keyvals == 'blue'] <- controlCondition

  print("ici 2")
  print(table(keyvals))
  print(colnames(eventFile))
  volcanoPlot <- EnhancedVolcano(eventFile,
                                 lab = "",
                                 x = 'IncLevelDifference',
                                 y = 'FDR',
                                 title = '',
                                 subtitle = paste(controlCondition, ':', sum(keyvals == 'blue'),
                                                  ' -- ', caseCondition, ':', sum(keyvals == 'red')),
                                 caption = paste("ILD cutoff =", ILD_cutoff, "; FDR cutoff =", FDR_cutoff),
                                 xlab = "Inclusion Level Difference",
                                 xlim = c(-1,1),
                                 ylim = c(0, -log10(min(eventFile[eventFile$FDR != 0,]$FDR) * 10^-1)),
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
  print("ici 3")
  ggsave(paste0(event, "events_volcano.svg"), device = "svg", plot = volcanoPlot, path = outputPath, width = 25, height = 30, dpi = 600)
  ggsave(paste0(event, "events_volcano.pdf"), device = "pdf", plot = volcanoPlot, path = outputPath, width = 25, height = 30, dpi = 600)
}

save_plot <- function(plot, filename) {
  png_filename <- paste0(filename, ".png")

  enrichPlot_dir <- paste0(outputPath_plots, "/enrichment_plots")
  if (!dir.exists(enrichPlot_dir)) {
    dir.create(enrichPlot_dir, recursive = TRUE)
  } else {
    cat("Output directory already exists. Skipping creation.\n")
  }

  ggsave(png_filename, plot = plot, path =  enrichPlot_dir, dpi = 600, device = 'png')
}

GSEA_and_PE <- function(backgroundList, modEventFile, condition, outputPath) {

  geneList <- modEventFile[modEventFile$Condition == condition,]$GeneSymbol

  p_value_cutoff <- FDR_cutoff

  if (genome_build == "hg38") {

    org_db <- org.Hs.eg.db

  } else {

    org_db <- org.Mm.eg.db

  }

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
    readable = TRUE )

  if (!is.null(GO_results)) {
    if (dim(GO_results)[1] > 0) {
      bpResults <- filter(GO_results, ONTOLOGY == "BP", p.adjust < p_value_cutoff)
      ccResults <- filter(GO_results, ONTOLOGY == "CC", p.adjust < p_value_cutoff)
      mfResults <- filter(GO_results, ONTOLOGY == "MF", p.adjust < p_value_cutoff)
      write.table(GO_results, file = paste0(outputPath, "/GOresultsUNFILTERED_", event, "_", condition, ".txt"), sep = "\t", quote = FALSE)
      GOresults <- filter(GO_results, p.adjust < p_value_cutoff)
      write.table(GOresults, file = paste0(outputPath, "/GOresultsFILTERED_p_value_cutoff_", event, "_", condition, ".txt"), sep = "\t", quote = FALSE)
      if (dim(bpResults)[1] > 0) {
        save_plot(dotplot(bpResults, showCategory = 15), paste0("BP_DotPlot_", event, "_", condition))
      }
      if (dim(ccResults)[1] > 0) {
        save_plot(dotplot(ccResults, showCategory = 15), paste0("CC_DotPlot_", event, "_", condition))
      }
      if (dim(mfResults)[1] > 0) {
        save_plot(dotplot(mfResults, showCategory = 15), paste0("MF_DotPlot_", event, "_", condition))
      }
    }
  }

  gene_list <- modEventFile[modEventFile$Condition == condition,]$ILD
  names(gene_list) <- gsub("\\..*","", modEventFile[modEventFile$Condition == condition,]$GeneID)
  print(gene_list)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list<-na.omit(gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)

  if (length(gene_list) > 1) {
    print("Here - gene-list")
    print(gene_list)
    gse <- gseGO(geneList=gene_list,
                 ont ="ALL",
                 keyType = "ENSEMBL",
                 minGSSize = 3,
                 scoreType = "pos",
                 maxGSSize = 800,
                 pvalueCutoff = 1,
                 verbose = TRUE,
                 OrgDb = org_db,
                 pAdjustMethod = "fdr")
    print("Here = passed")
    if (!is.null(gse)) {
      if(dim(gse)[1] > 0) {
        write.table(gse, file = paste0(outputPath, "/GSEAresultsUNFILTERED_", event, "_", condition, ".txt"), sep = "\t", quote = FALSE)
        sigGSE_res <- filter(gse, p.adjust < p_value_cutoff)

        if (dim(sigGSE_res)[1] > 0) {
          write.table(sigGSE_res, file = paste0(outputPath, "/GSEAresultsFILTERED_", event, "_", condition, ".txt"), sep = "\t", quote = FALSE)
          save_plot(dotplot(sigGSE_res, showCategory = 15), paste0("GSEAresults_DotPlot_", event, "_", condition))
        }
      }
    }
  }
}

# Function to get significance symbol
get_significance <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

plotViolin_and_Density <- function(modEventFile, col_name, y_label, event, outputPath) {

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

  print(filtered_file)
  if (n_distinct(filtered_file$Condition) > 1 & sum(filtered_file$Condition == caseCondition) >= 2 & sum(filtered_file$Condition == controlCondition) >=  2) {
    print("Got through")
    t_test <- t.test(as.formula(paste(clean_col, "~ Condition")), data = filtered_file)
    p_value <- t_test$p.value
    sig <- get_significance(p_value)

    # Calculate means and positions
    max_val <- max(filtered_file[[clean_col]], na.rm = TRUE)
    y_pos <- max_val * 1.25

    mean_df <- filtered_file %>%
      group_by(Condition) %>%
      summarize(mean_val = if (col_name %in% c("ILD", "BPscore", "Log2FC", "maxEntScanScore_3prime", "maxEntScanScore_5prime")) {
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

    ggsave(paste0(event, "_", col_name, "_ViolinDensity_plots.svg"), device = "svg", plot = combined, path = outputPath, width = 20, height = 8)
    ggsave(paste0(event, "_", col_name, "_ViolinDensity_plots.pdf"), device = "pdf", plot = combined, path = outputPath, width = 20, height = 8)
  }

}

branchPoint_calc <- function(modEventFile) {
  modEventFile <- modEventFile[!grepl("_", modEventFile$Chrom),]

  if (event != "RI") {
    print("here")
    grange_bed <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                          ranges = IRanges(start = modEventFile$UpstreamIS, end = modEventFile$UpstreamIE),
                          names = modEventFile$GeneSymbol)
  } else {
    grange_bed <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand,
                          ranges = IRanges(start = modEventFile$IntronStart, end = modEventFile$IntronEnd),
                          names = modEventFile$GeneSymbol)
  }
  event_bed_df <- as.data.frame(grange_bed)
  colnames(event_bed_df) <- c("Chrom", "start", "end", "width", "strand", "GeneSymbol")

  hits <- findOverlaps(banchpoint_anno, grange_bed)
  hits <- as.data.frame(hits)
  print(hits)
  bp_event_df <- data.frame()

  for (num in 1:nrow(hits)) {
    bp_pos_index <- hits[num,]$queryHits
    event_bed_index <- hits[num,]$subjectHits
    print(event_bed_df[event_bed_index,])
    print("_______________")
    print(branchpoint_anno_df[bp_pos_index,])

    bp_event_df <- rbind(bp_event_df,merge(branchpoint_anno_df[bp_pos_index,], event_bed_df[event_bed_index,], by = c("Chrom", "strand")))
  }

  print("HERE")
  bp_event_df$dist_to_bp <- bp_event_df$bp_pos - bp_event_df$bp_start
  print(bp_event_df)
  bp_event_df <- bp_event_df %>%
    group_by(GeneSymbol) %>%
    slice_min(abs(dist_to_bp)) %>%
    ungroup

  modEventFile <- merge(modEventFile, bp_event_df[,c("GeneSymbol", "bp_start", "bp_pos", "dist_to_bp", "bp_score")], by = "GeneSymbol")
  print("HERE_2")
  return(modEventFile)
}
#####################################################

filterEventFile <- function(eventFile, ILD_cutoff, IJC_cutoff, SJC_cutoff) {

  # checks if column is indeed a list a values
  if( class(eventFile$IJC_SAMPLE_1) == "character") {
    eventFile$IJC_SAMPLE_1_mean <- unlist(lapply(type.convert(strsplit(eventFile$IJC_SAMPLE_1, ","), as.is = TRUE), mean))
    eventFile$SJC_SAMPLE_1_mean <- unlist(lapply(type.convert(strsplit(eventFile$SJC_SAMPLE_1, ","), as.is = TRUE), mean))
    eventFile$IJC_SAMPLE_2_mean <- unlist(lapply(type.convert(strsplit(eventFile$IJC_SAMPLE_2, ","), as.is = TRUE), mean))
    eventFile$SJC_SAMPLE_2_mean <- unlist(lapply(type.convert(strsplit(eventFile$SJC_SAMPLE_2, ","), as.is = TRUE), mean))

    eventFile$IncLevel1_mean <- unlist(lapply(type.convert(strsplit(eventFile$IncLevel1, ","), as.is = TRUE), mean))
    eventFile$IncLevel2_mean <- unlist(lapply(type.convert(strsplit(eventFile$IncLevel2, ","), as.is = TRUE), mean))
  }

  if (event != "MXE") {
    eventFile_filtered <- eventFile[((eventFile$IJC_SAMPLE_1_mean > IJC_cutoff & eventFile$SJC_SAMPLE_1_mean > SJC_cutoff) | (eventFile$IJC_SAMPLE_2_mean > IJC_cutoff & eventFile$SJC_SAMPLE_2_mean > SJC_cutoff))
                                    & eventFile$FDR < FDR_cutoff & abs(eventFile$IncLevelDifference) > ILD_cutoff,]
  } else {
    # TO DO: Check if filtering for MXE is correct
    eventFile_filtered <- eventFile[ifelse(eventFile$strand=="+", (eventFile$IJC_SAMPLE_1_mean > IJC_cutoff & eventFile$SJC_SAMPLE_1_mean > SJC_cutoff),
                                           (eventFile$IJC_SAMPLE_2_mean > IJC_cutoff & eventFile$SJC_SAMPLE_2_mean > SJC_cutoff))
                                    & eventFile$FDR < FDR_cutoff & abs(eventFile$IncLevelDifference) > ILD_cutoff,]
  }

  if (event == "SE" | event == "MXE") {
    eventFile_filtered$Condition <- ifelse(eventFile_filtered$IncLevelDifference < -ILD_cutoff, controlCondition, caseCondition)
  } else {
    eventFile_filtered$Condition <- ifelse(eventFile_filtered$IncLevelDifference > ILD_cutoff, controlCondition, caseCondition)
  }

  return(eventFile_filtered)
}

### Functions to process events ###
process_se_mats <- function(eventFile) {

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

  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile)
  modEventFile <- branchPoint_calc(modEventFile)

  return(modEventFile)

}

process_ri_mats <- function(eventFile) {

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

  modEventFile$start_5prime <- ifelse(modEventFile$strand == "+", modEventFile$DownstreamES - 6, modEventFile$UpstreamEE - 2)
  modEventFile$end_5prime <- ifelse(modEventFile$strand == "+", modEventFile$DownstreamES + 2, modEventFile$UpstreamEE + 6)
  modEventFile$start_3prime <- ifelse(modEventFile$strand == "+", modEventFile$UpstreamEE - 2, modEventFile$DownstreamES - 20)
  modEventFile$end_3prime <- ifelse(modEventFile$strand == "+", modEventFile$UpstreamEE + 20, modEventFile$DownstreamES + 2)

  grange_5prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand, ranges = IRanges(start = modEventFile$start_5prime, end = modEventFile$end_5prime), names = modEventFile$GeneID)
  grange_3prime <- GRanges(seqnames = modEventFile$Chrom, strand = modEventFile$strand, ranges = IRanges(start = modEventFile$start_3prime, end = modEventFile$end_3prime), names = modEventFile$GeneID)

  if (genome_build == "mm10") {
    modEventFile$seq_3prime <- getSeq(Mmusculus, grange_3prime, as.character=TRUE)
    modEventFile$seq_5prime <- getSeq(Mmusculus, grange_5prime, as.character=TRUE)
  } else {
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

  modEventFile <- branchPoint_calc(modEventFile)

  return(modEventFile)

}

process_a3ss_mats <- function(eventFile) {

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

  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile)
  modEventFile <- branchPoint_calc(modEventFile)


  return(modEventFile)

}

process_a5ss_mats <- function(eventFile) {

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

  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile)
  modEventFile <- branchPoint_calc(modEventFile)


  return(modEventFile)

}

process_mxe_mats <- function(eventFile) {

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

  modEventFile <- PTC_CE_MaxEntScan_Calc(modEventFile)
  modEventFile <- branchPoint_calc(modEventFile)

  return(modEventFile)

}

### Functions to plot summary visualizations ###

plot_gene_type_pie <- function() {
  gtf_file <- as.data.frame(gtf_file)
  gtf_file <- gtf_file[,c("gene_id", "gene_biotype")]
  gtf_file <- gtf_file[!duplicated(gtf_file),]

  mm10_SE_events <- merge(mm10_SE_events, gtf_file, by.x = "GeneID", by.y = "gene_id")
  mm10_SE_events$gene_type_cond <- paste0(mm10_SE_events$gene_biotype, "_", mm10_SE_events$Condition)

  mm10_sub <- data.frame(table((mm10_SE_events$gene_type_cond)))
  mm10_sub <- merge(mm10_sub, mm10_SE_events[,c("Condition", "gene_type_cond", "gene_biotype")], by.x = "Var1", by.y = "gene_type_cond")
  mm10_sub <- mm10_sub[!duplicated(mm10_sub),]

  ggplot(mm10_sub, aes(x = factor(Condition), y = Freq, fill = factor(gene_biotype))) +
    geom_col() +
    scale_x_discrete(limits = unique(mm10_SE_events$Condition)) +
    coord_polar("y")
}
plotSummaryBarChart <- function(count, events, condition, gene_count, outputPath) {

  df <- data.frame(count,events,condition)

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

  ggsave(paste0("eventCount_bar.png"), plot = barChart, path = outputPath, width = 25, height = 30, dpi = 600)

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

  ggsave(paste0("geneCount_bar.png"), plot = barChart, path = outputPath, width = 25, height = 30, dpi = 600)

}

##################################################### run script

AS_events <- c("SE", "RI", "MXE", "A5SS", "A3SS")

# will run for all combinations of ILD cutoffs, IJC cutoffs, and SJC cutoffs

for (ILD_cutoff in ILD_list) {

  for (IJC_cutoff in IJC_list) {

    for (SJC_cutoff in SJC_list) {

      # name of run, specifies which cutoffs are used, will be used to create subdirectory of output directory
      run <- paste("run", "ILD", ILD_cutoff, "FDR", FDR_cutoff, "IJC", IJC_cutoff, "SJC", SJC_cutoff, sep = "_")

      # initialize empty lists to store data to create summary barcharts for each run
      count <- c()
      events <- c()
      condition <- c()
      gene_count <- c()

      for (event in AS_events) {

        inputFile <- list.files(paste0(input_dir, "/"), pattern = paste0(event, ".MATS.JC.txt"), full.names = TRUE) # TO DO: Provide option for JCEC

        print(inputFile)
        # checks if input file exists, if not, proceeds onto the next event
        if (file.exists(inputFile)) {

          print(paste0("Processing ",  inputFile, "..."))
          eventFile <- read.delim(inputFile, header = TRUE, sep = "\t") # reads in input file

          eventFile_filtered <- filterEventFile(eventFile, ILD_cutoff, IJC_cutoff, SJC_cutoff) # filters events on current cutoffs
          cat(event, " filtered; ", nrow(eventFile), "events originally, ", nrow(eventFile_filtered), " events remaining\n")

          # appends data to lists for creation of summary barcharts
          # count is number of all events that pass filtering, gene_count is number of unique genes with events that pass filtering
          count <- append(count, c(sum(eventFile_filtered$Condition == controlCondition), sum(eventFile_filtered$Condition == caseCondition)))
          events <- append(events, rep(event, 2))
          condition <- append(condition, c(controlCondition, caseCondition))
          gene_count <- append(gene_count, c(length(unique(eventFile_filtered[eventFile_filtered$Condition == controlCondition,]$geneSymbol)),
                                             length(unique(eventFile_filtered[eventFile_filtered$Condition == caseCondition,]$geneSymbol))))


          outputPath <- file.path(output_dir, run, event) # creates output path to write output to
          dir.create(outputPath, recursive = TRUE)

          outputPath_plots <- file.path(output_dir, run, event, "figures") # creates output path to write output to
          dir.create(outputPath_plots, recursive = TRUE)

          write.table(eventFile_filtered, paste0(outputPath, "/", event, ".MATS.filtered.JC.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

          if (genome_build == "mm10" | genome_build == "hg38") {
            if (nrow(eventFile_filtered) > 0) {
              # call event-specific processing function

              if (event == "SE") {
                modEventFile <- process_se_mats(eventFile_filtered)
              } else if (event == "RI") {
                modEventFile <- process_ri_mats(eventFile_filtered)
              } else if (event == "A3SS") {
                modEventFile <- process_a3ss_mats(eventFile_filtered)
              } else if (event == "A5SS") {
                modEventFile <- process_a5ss_mats(eventFile_filtered)
              } else if (event == "MXE") {
                modEventFile <- process_mxe_mats(eventFile_filtered)
              }

              modEventFile$GeneID <- gsub("\\..*","", modEventFile$GeneID)

              #modEventFile <- merge(modEventFile, gtf_file[c("gene_id", "gene_biotype"),], by.x = "GeneID", by.y = "gene_id")

              write.table(modEventFile, file = paste0(outputPath, "/", event, ".MATS.modified.JC.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

              # for each event - plot chart of ILD, 5' and 3' maxEntScanScore, intron/exon length, BPscore, BPdist

              gene_list_background <- eventFile$geneSymbol # TO DO: Change background
              GSEA_and_PE(gene_list_background, modEventFile, controlCondition, outputPath_plots)
              GSEA_and_PE(gene_list_background, modEventFile, caseCondition, outputPath_plots)

              plotViolin_and_Density(modEventFile, "ILD", "Inclusion Level Difference", event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "maxEntScanScore_5prime", "5' Splice Site Strength Score", event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "maxEntScanScore_3prime", "3' Splice Site Strength Score", event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "dist_to_bp", "Distance to Branchpoint", event, outputPath_plots)
              plotViolin_and_Density(modEventFile, "bp_score", "Branchpoint Score", event, outputPath_plots)

              print(head(modEventFile))
              if (event != "RI") {
                print(outputPath_plots)
                plotEventBarChart(modEventFile, outputPath_plots, event)
              }

              if (event == "RI") {

                plotViolin_and_Density(modEventFile, "IntronLength", "Intron Length", event, outputPath_plots)

              } else if (event == "MXE") {

                plotViolin_and_Density(modEventFile, "ExonInUseLength", "Exon Length", event, outputPath_plots)
                plotViolin_and_Density(modEventFile, "UpstreamIL", " Upstream Intron Length", event, outputPath_plots)
                plotViolin_and_Density(modEventFile, "DownstreamIL", " Downstream Intron Length", event, outputPath_plots)

              } else if (event == "SE") {
                plotViolin_and_Density(modEventFile, "ExonLength", "Exon Length", event, outputPath_plots)
                plotViolin_and_Density(modEventFile, "UpstreamIL", " Upstream Intron Length", event, outputPath_plots)
                plotViolin_and_Density(modEventFile, "DownstreamIL", " Downstream Intron Length", event, outputPath_plots)
                plotViolin_and_Density(modEventFile, "UpstreamEL", " Upstream Exon Length", event, outputPath_plots)
                plotViolin_and_Density(modEventFile, "DownstreamEL", " Downstream Exon Length", event, outputPath_plots)

                # To be added:(BPscore, BPDist, Log2FC)
              } else {
                plotViolin_and_Density(modEventFile, "UpstreamIL", " Upstream Intron Length", event, outputPath_plots)
              }

              }
          }

        } else {

          print(paste0("Skipping ", inputFile, " (not found)."))

        }

      }
      plotSummaryBarChart(count, events, condition, gene_count, file.path(output_dir, run)) # creates summary barcharts

    }
  }
}

