## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ggplot2", "pbapply", "Rsamtools", "GenomicAlignments")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

files = list.files(paste0("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/MS147/variants/bam/"), 
                   pattern = ".bam$", full.names = T)

ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv")

param = ScanBamParam(what=c("qname","flag"))

reads = pblapply(files, function(file) {
  gr = granges(readGAlignments(file, param=param))
  return(data.table(sample = gsub(".*\\/|.trim.sorted.bam", "", file),
                    reads = length(gr)))
}, cl = 32)

reads = rbindlist(reads)
reads[, sample := gsub("_", " ", sample)]

# Remove controls and sample 12
reads = reads[!grepl("Neg|Pos|Sample 12", sample)]
ct = ct[!grepl("Sample 12", sample)]

total = merge(reads, ct)

plt = ggplot(total, aes(x = ct, y = reads)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = F, color = "red", linetype = 2, size = 2, fullrange = T) +
  scale_y_continuous(breaks = c(0, 2.5e5, 5e5, 7.5e5, 1e6), labels = scales::label_comma()) +
  stat_cor(method = "pearson") +
  labs(y = "Reads",
       x = "Ct value")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-reads/MS147_PE150-reads-ct",
              height = 9, width = 9)
