## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/CCI-96-ct.tsv", header = T)
ct[, sample := gsub("_", " ", sample)]
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")
#annot = fread("/mnt/AchTeraD/data/BICRO268/MS147+148-barcodes-annot.txt", header = F)

files = list.files("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/MS169/variants/bam/", pattern = ".bam$", full.names = T)

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
total = merge(total, ct, by = "sample")

# Melt and prepare for plotting
#total_m = melt(total[, c(9, 4, 7)], id.vars = "sample")
total_m = melt(total[, c(1, 5, 8, 11)], id.vars = c("sample", "ct"))
total_m[, depth := gsub(".*_", "", variable)]
total_m[, region := gsub("_.*", "", variable)]

# Order by ct
setorder(total_m, ct, -value)
total_m = total_m[!grepl("Neg", sample)]
total_m[, sample := factor(sample, levels = total_m[region == "full"]$sample)]

plt1 = ggplot(total_m[region == "full"], aes(x = sample, y = value, group = depth, fill = ct)) +
  geom_col() +
  scale_fill_viridis_c(direction = -1) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limit = c(0, 1)) +
  #geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/barplots-boc/MS169-96samples-full",
              height = 6, width = 18)

total_m[, sample := factor(sample, levels = total_m[region == "Sregion"]$sample)]
plt2 = ggplot(total_m[region == "Sregion"], aes(x = sample, y = value, group = depth, fill = value)) +
  geom_col() +
  scale_fill_viridis_c(direction = -1) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limit = c(0, 1)) +
  #geom_hline(yintercept = c(theoret$`NlaIII + MseI`[5]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/barplots-boc/MS169-96samples-Sregion",
              height = 6, width = 18)

# Dotplot
total_m[, sample := factor(sample, levels = total_m[region == "full"]$sample)]

plt3 = ggplot(total_m[region == "full"], aes(x = sample, y = value)) +
  geom_point(size = 3) +
  labs(y = "Breadth of Coverage (10x) ", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/barplots-boc/MS169-96samples-scatter",
              height = 7, width = 8)
