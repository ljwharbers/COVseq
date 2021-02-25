## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30

samples = c("MS126_S1", "MS127_S2", "MS128_S3", "MS129_S4", "MS130_S5", "MS131_S6", "MS132_S7", "MS133_S8")

for(sample in samples) {
  libraries = list.files(paste0("/mnt/AchTeraD/data/BICRO264/", sample), pattern = ".bam$", recursive = T, full.names = T)
  annot = fread("/mnt/AchTeraD/data/BICRO264/barcodes-annot.txt", header = F)
  ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)
  
  # annot[, V2 := c(paste0("sample", 1:55), "N1", "P1", "N2", "P2", "N3")]
  #libraries = libraries[grepl("merged", libraries)]
  # libraries = libraries[c(1:2, 6)]
  #libraries = libraries[grepl("sample|pos|neg", libraries)]
  
  # params
  param = ScanBamParam(what=c("qname","flag"))
  s_region = GRanges(seqnames = "NC_045512.2", IRanges(start = 21563, end = 25384))
  s_length = 25384 - 21562
  
  # go through list of bam files and calculate coverage in the full genome and in the S region of SARS-CoV-2
  res = pblapply(libraries, function(library) {
    # Load bamfile as granges
    gr = granges(readGAlignments(library, param=param))
    
    # Get coverage full genome
    cov = coverage(gr)
    full_1 = sum(cov$NC_045512.2@lengths[cov$NC_045512.2@values > 0]) / sum(cov$NC_045512.2@lengths)
    full_5 = sum(cov$NC_045512.2@lengths[cov$NC_045512.2@values > 4]) / sum(cov$NC_045512.2@lengths)
    full_10 = sum(cov$NC_045512.2@lengths[cov$NC_045512.2@values > 9]) / sum(cov$NC_045512.2@lengths)
    
    # Subset for S-region  
    gr_sregion = subsetByOverlaps(gr, s_region, type = "within")
    
    # Get coverage S-region
    cov_s = coverage(gr_sregion)
    
    s_1 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 0]) / s_length
    s_5 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 4]) / s_length
    s_10 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 9]) / s_length 
    
    return(data.table(sample = gsub(".*/|\\..*", "", library),
                      aligned_reads = length(gr),
                      full_1x = full_1,
                      full_5x = full_5,
                      full_10x = full_10,
                      Sregion_1x = s_1,
                      Sregion_5x = s_5,
                      Sregion_10x = s_10,
                      average_coverage = mean(cov)))
  }, cl = num_threads)
  
  total = rbindlist(res, use.names = T)
  total = merge(total, annot, by.x = "sample", by.y = "V1")
  # Set new names
  total[, sample := factor(V2, levels = annot$V2)]
  # total[, sample := factor(paste0("Sample ", 1:30), levels = paste0("Sample ", 1:30))]
  #total[, sample := factor(sample, levels = c(paste0("sample", 1:30), "pos", "neg"))]
  total_m = melt(total, id.vars = c("sample", "V2"))
  # Order by ct
  total_m[, sample := factor(sample, levels = c(ct$sample[order(ct$ct)], "Sample Neg", "Sample Pos"))]
  
  
  #total_m[, sample := factor(sample, levels = c(paste0("sample", 1:30), "pos", "neg"))]
  
  plt1 = ggplot(total_m[variable == "aligned_reads"], aes(x = sample, y = value, fill = sample)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_viridis_d() +
    labs(title = "Total reads") +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none")
  
  plt2 = ggplot(total_m[variable == "average_coverage"], aes(x = sample, y = value, fill = sample)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_viridis_d() +
    labs(title = "Average coverage") +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none")
  
  plt3 = ggplot(total_m[!variable %in% c("aligned_reads", "average_coverage"), ], aes(x = sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_viridis_d() +
    labs(title = "Percentage of genome covered") +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom")
  
  plt = grid.arrange(plt1, plt2, plt3)
  save_and_plot(plot_grid(plt), paste0("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/reads-cov-boc/", sample),
                height = 14, width = 20)
}

