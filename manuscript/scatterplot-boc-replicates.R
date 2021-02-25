## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2", "scales")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/55-patients_CT-clean.tsv")
annot = fread("/mnt/AchTeraD/data/BICRO251/TN21-annot_barcodes.txt", header = F)
annot[, V1 := paste0("Sample ", V1)]

files = list.files("/mnt/AchTeraD/data/BICRO251/TN21/trimmed/", pattern = ".bam$", full.names = T)
files = c(files, list.files("/mnt/AchTeraD/data/BICRO254/", pattern = ".bam$", full.names = T, recursive = T))

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

samplenames = paste0(gsub(".*BICRO254\\/\\/|.*BICRO251\\/|\\/trimmed.*", "", files), "-", gsub(".trimmed.*", "", basename(files)))

total[, barcode := gsub(".*-", "", samplenames)]
total[, library := gsub("-.*", "", samplenames)]

total = merge(annot, total, by.x = "V2", by.y = "barcode")

# Remove P/N
total = total[!grepl("N|P", V1)]

# Correlate scatterplots
total_subs = total[, .(V1, library, full_10x)]
total_wide = dcast(total_subs, V1~library)

# Make plots

plt1 = ggplot(total_wide, aes(x = TN21, y = MS82_S3)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2) +
  scale_y_continuous(labels = number_format()) + 
  scale_x_continuous() +
  stat_cor()

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/scatterplots-boc/MS82-ref-scatterplot",
              width = 7, height = 7)

plt2 = ggplot(total_wide, aes(x = TN21, y = MS83_S4)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2) +
  scale_y_continuous(labels = number_format()) + 
  scale_x_continuous() +
  stat_cor()

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/scatterplots-boc/MS83-ref-scatterplot",
              width = 7, height = 7)

plt3 = ggplot(total_wide, aes(x = TN21, y = MS84_S5)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2) +
  scale_y_continuous(labels = number_format()) + 
  scale_x_continuous() +
  stat_cor()

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/scatterplots-boc/MS84-ref-scatterplot",
              width = 7, height = 7)

plt4 = ggplot(total_wide, aes(x = TN21, y = NZ190_S1)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2) +
  scale_y_continuous(labels = number_format()) + 
  scale_x_continuous() +
  stat_cor()

save_and_plot(plt4, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/scatterplots-boc/NZ190-ref-scatterplot",
              width = 7, height = 7)

plt5 = ggplot(total_wide, aes(x = TN21, y = NZ191_S2)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2) +
  scale_y_continuous(labels = number_format()) + 
  scale_x_continuous() +
  stat_cor()

save_and_plot(plt5, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/scatterplots-boc/NZ191-ref-scatterplot",
              width = 7, height = 7)
