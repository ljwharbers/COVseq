## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for correlating number of reads vs CT

## Load/install packages
packages = c("data.table", "pbapply", "scales")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load files
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)

files = list.files("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/MS147_miseq/variants/bam/", pattern = ".bam$", full.names = T)
files = c(files, list.files("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/NEBNext/variants/bam/", pattern = ".bam$", full.names = T))
files = files[!grepl("Pos|Neg|Sample_12", files)]


samplenames = gsub(".trim.sorted.bam", "", basename(files))
samplenames = paste0(c(rep("COVseq - ", 29), rep("NEBNext - ", 29)), samplenames)

# Get read counts
param = ScanBamParam(what=c("qname"))

reads = pblapply(1:length(files), function(i) {
  gr = granges(readGAlignments(files[i], param=param))
  return(data.table(sample = samplenames[i], reads = length(gr)))
}, cl = 30)

tot_reads = rbindlist(reads)
tot_reads[, method := gsub(" -.*", "", sample)]
tot_reads[, sample := gsub(".* - ", "", sample)]
tot_reads = reshape(tot_reads, idvar = "sample", timevar = "method", direction = "wide")
setnames(tot_reads, c("sample", "COVseq", "NEBNext"))

tot_reads[, sample := gsub("_", " ", sample)]

# Add Ct value
tot_reads = merge(tot_reads, ct, by = "sample", all.x = T)

# Plot scatterplot
plt = ggplot(tot_reads, aes(x = NEBNext, y = COVseq, color = ct)) +
  geom_point(size = 4) +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2, size = 2) +
  scale_color_viridis_c(option = "E", direction = -1) +
  scale_y_continuous(labels = scales::number_format()) + 
  scale_x_continuous(labels = scales::number_format()) +
  labs(y = "# Reads (COVseq)", x = "# Reads (NEBNext)") +
  stat_cor()

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-reads/SE300-reads_covseq_nebnext-scatterplot",
              height = 9, width = 9)
