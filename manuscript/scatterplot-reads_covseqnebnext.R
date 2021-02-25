## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for correlating number of reads vs CT

## Load/install packages
packages = c("data.table", "pbapply", "scales")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load files
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)

files = list.files("/mnt/AchTeraD/data/BICRO251/TN11/trimmed/", pattern = ".bam$", full.names = T)
files = files[grepl("sample", files)]
files = c(files, list.files("/mnt/AchTeraD/data/BICRO237/covseq/", recursive = T, pattern = "trimmed.sorted.bam$", full.names = T),
          list.files("/mnt/AchTeraD/data/BICRO240/covseq/", recursive = T, pattern = "trimmed.sorted.bam$", full.names = T))

samplenames = gsub(".bam", "", basename(files[1:30]))
samplenames = c(gsub("sample", "COVseq - Sample ", samplenames), paste0("NEBNext - Sample ", 1:30))

# Get read counts
param = ScanBamParam(what=c("qname"))

reads = pblapply(1:length(files), function(i) {
  gr = granges(readGAlignments(files[i], param=param))
  return(data.table(sample = samplenames[i], reads = length(gr)))
}, cl = 30)

tot_reads = rbindlist(reads)
tot_reads[, method := gsub(" -.*", "", sample)]
tot_reads[, sample := gsub(".* -", "", sample)]
tot_reads = reshape(tot_reads, idvar = "sample", timevar = "method", direction = "wide")
setnames(tot_reads, c("sample", "COVseq", "NEBNext"))

# Plot scatterplot
plt = ggplot(tot_reads, aes(x = NEBNext, y = COVseq)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 1) +
  scale_y_continuous(labels = number_format()) + 
  scale_x_continuous(labels = number_format()) +
  stat_cor()

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/scatterplot-reads/SE150-reads_covseq_nebnext-scatterplot",
              height = 6, width = 7)
