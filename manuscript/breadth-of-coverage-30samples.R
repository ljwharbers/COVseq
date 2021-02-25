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

files = list.files("/mnt/AchTeraD/data/BICRO251/TN11/trimmed/", pattern = ".bam$", full.names = T)
files = c(files, list.files("/mnt/AchTeraD/data/BICRO248+249/TN11_S2/trimmed/", pattern = ".bam$", full.names = T))
files = files[grepl("sample", files)]

samplenames = gsub(".bam", "", basename(files))
samplenames[1:30] = gsub("sample" ,"SE150 - Sample ", samplenames[1:30])
samplenames[31:60] = gsub("sample", "SE75 - Sample ", samplenames[31:60])

names(files) = samplenames

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
total[, sample := names(files)]

total_m = melt(total[, c(1, 3:6)], id.vars = "sample")
total_m[, readlength := ifelse(grepl("SE75", sample), "SE75", "SE150")]
total_m[, sample := gsub(".*- ", "", sample)]
total_m[, sample := factor(sample, levels = ct$sample[order(ct$ct)])]
total_m[, depth := gsub(".*_", "", variable)]
total_m[, region := gsub("_.*", "", variable)]

plt1 = ggplot(total_m[readlength == "SE75" & region == "full"], aes(x = sample, y = value, group = depth, color = depth)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_viridis_d() +
  scale_y_continuous(labels = percent_format()) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[1]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/lineplots-boc/SE75-COvseq-30samples-color-depth-genome",
              height = 7, width = 8)

plt2 = ggplot(total_m[readlength == "SE150" & region == "full"], aes(x = sample, y = value, group = depth, color = depth)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_viridis_d() +
  scale_y_continuous(labels = percent_format()) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/lineplots-boc/SE150-COvseq-30samples-color-depth-genome",
              height = 7, width = 8)

plt3 = ggplot(total_m[readlength == "SE75" & region == "full"], aes(x = sample, y = value, group = depth, color = depth)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_viridis_d() +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[1]), color = "red", linetype = "dashed", size = 1) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/lineplots-boc/SE75-COvseq-30samples-color-depth-genome_fixedY",
              height = 7, width = 8)

plt4 = ggplot(total_m[readlength == "SE150" & region == "full"], aes(x = sample, y = value, group = depth, color = depth)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  scale_color_viridis_d() +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt4, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/lineplots-boc/SE150-COvseq-30samples-color-depth-genome_fixedY",
              height = 7, width = 8)
