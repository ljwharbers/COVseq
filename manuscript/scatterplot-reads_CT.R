## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for correlating number of reads vs CT

## Load/install packages
packages = c("data.table", "pbapply", "scales", "Rsamtools", "GenomicRanges", "GenomicAlignments")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load files
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)

# files = list.files("/mnt/AchTeraD/data/BICRO251/TN11/trimmed/", pattern = ".bam$", full.names = T)
# files = files[grepl("sample", files)]
files = list.files("/mnt/AchTeraD/data/BICRO237/covseq/", pattern = ".bam$", full.names = T, recursive = T)
files = c(files, list.files("/mnt/AchTeraD/data/BICRO240/covseq/", pattern = ".bam$", full.names = T, recursive = T))
files = files[!grepl("dedup", files)]
# samplenames = gsub(".bam", "", basename(files))
# samplenames = gsub("sample", "Sample ", samplenames)
samplenames = paste0("Sample ", 1:30)

# Get read counts
param = ScanBamParam(what=c("qname"))

reads = pblapply(1:length(files), function(i) {
  gr = granges(readGAlignments(files[i], param=param))
  return(data.table(sample = samplenames[i], reads = length(gr)))
}, cl = 30)

tot_reads = rbindlist(reads)
tot_reads = merge(tot_reads, ct, by = "sample")

# Plot scatterplot
plt = ggplot(tot_reads, aes(x = ct, y = reads)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2) +
  scale_y_continuous(labels = number_format()) + 
  scale_x_continuous() +
  stat_cor()

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/scatterplot-reads/SE150-reads_CT-scatterplot-NEBNext",
              height = 6, width = 7)
