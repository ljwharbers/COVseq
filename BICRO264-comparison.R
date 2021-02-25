## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=32

libraries = c("MS126_S1", "MS127_S2", "MS128_S3", "MS129_S4", "MS130_S5", "MS131_S6", "MS132_S7", "MS133_S8")
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")
annot = fread("/mnt/AchTeraD/data/BICRO264/barcodes-annot.txt", header = F)

results = lapply(libraries, function(lib) {
  files = list.files(paste0("/mnt/AchTeraD/data/BICRO264/", lib, "/trimmed/"), pattern = ".bam$", full.names = T)
  
  samplenames = gsub(".trimmed.sorted.bam", "", basename(files))
  # params
  param = ScanBamParam(what=c("qname","flag"))
  s_region = GRanges(seqnames = "NC_045512.2", IRanges(start = 21563, end = 25384))
  s_length = 25384 - 21562
  
  # go through list of bam files and calculate coverage in the full genome and in the S region of SARS-CoV-2
  res = pblapply(files, function(library) {
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
    
    return(data.table(aligned_reads = length(gr),
                      full_1x = full_1,
                      full_5x = full_5,
                      full_10x = full_10,
                      Sregion_1x = s_1,
                      Sregion_5x = s_5,
                      Sregion_10x = s_10,
                      average_coverage = mean(cov)))
  }, cl = num_threads)
  
  total = rbindlist(res)
  total[, sample := samplenames]
  
  # Merge with annotation
  total = merge(total, annot, by.x = "sample", by.y = "V1")
  total[, sample := V2]
  total[, library := lib]
  return(total)
})

total = rbindlist(results)

# Melt with only 10x columns
total_m = melt(total[, c(1, 5, 8, 11)], id.vars = c("sample", "library"))
total_m[, depth := gsub(".*_", "", variable)]
total_m[, region := gsub("_.*", "", variable)]

# Order by ct
total_m[, sample := factor(sample, levels = c(ct$sample[order(ct$ct)], "Sample Neg", "Sample Pos"))]


plt1 = ggplot(total_m[region == "full"], aes(x = sample, y = value, group = library, color = library)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_hue() +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage (10x)", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/lineplots-boc/library-comparisons_10xBoC-full",
              height = 8, width = 14)

plt2 = ggplot(total_m[region == "Sregion"], aes(x = sample, y = value, group = library, color = library)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_hue() +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[5]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage (10x)", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/lineplots-boc/library-comparisons_10xBoC-sregion",
              height = 8, width = 14)

# Plot mean per condition
conditions = data.table(library = libraries, condition = c("EDTA", "EDTA",
                                                           "95c", "95c",
                                                           "65c", "65c",
                                                           "std-ligase", "std-ligase"))

total_m = merge(total_m, conditions, by = "library")
total_m = total_m[, .(mean_val = mean(value)), by = .(sample, variable, region, condition)]

plt1 = ggplot(total_m[region == "full"], aes(x = sample, y = mean_val, group = condition, color = condition)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_hue() +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage (10x)", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/lineplots-boc/condition_mean-comparisons_10xBoC-full",
              height = 8, width = 14)

plt2 = ggplot(total_m[region == "Sregion"], aes(x = sample, y = mean_val, group = condition, color = condition)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_hue() +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[5]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage (10x)", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/lineplots-boc/condition_mean-comparisons_10xBoC-sregion",
              height = 8, width = 14)
