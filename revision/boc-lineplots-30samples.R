## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")
#annot = fread("/mnt/AchTeraD/data/BICRO268/MS147+148-barcodes-annot.txt", header = F)

files = list.files("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/MS147/variants/bam/", pattern = ".bam$", full.names = T)

samplenames = gsub(".trim.sorted.bam", "", basename(files))
samplenames = gsub("_", " ", samplenames)
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
#total = merge(total, annot, by.x = "sample", by.y = "V1")

# Melt and prepare for plotting
total_m = melt(total[, c(9, 4, 7)], id.vars = "sample")
total_m[, depth := gsub(".*_", "", variable)]
total_m[, region := gsub("_.*", "", variable)]

# Order by ct
total_m[, sample := factor(sample, levels = c(ct$sample[order(ct$ct)], "Sample Neg", "Sample Pos"))]

plt1 = ggplot(total_m[region == "full"], aes(x = sample, y = value, group = depth, color = depth)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_viridis_d() +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/lineplots-boc/MS147_fullregion-30+ctrls",
              height = 7, width = 8)

plt2 = ggplot(total_m[region == "Sregion"], aes(x = sample, y = value, group = depth, color = depth)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_viridis_d() +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/lineplots-boc/MS147_sregion-30+ctrls",
              height = 7, width = 8)