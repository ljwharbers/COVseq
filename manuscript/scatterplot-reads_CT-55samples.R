## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for correlating number of reads vs CT

## Load/install packages
packages = c("data.table", "pbapply", "scales", "Rsamtools", "GenomicRanges")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load files
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/55-patients_CT-clean.tsv", header = T)
annot = fread("/mnt/AchTeraD/data/BICRO251/TN21-annot_barcodes.txt")
annot[, V1 := paste0("Sample ", V1)]
ct = merge(annot, ct, by.x = "V1", by.y = "sample", all.x = T)

dirs = list.dirs("/mnt/AchTeraD/data/BICRO258/", recursive = F, full.names = F)
#dirs = dirs[!grepl("fastq", dirs)]
dirs = dirs[grepl("TN22|TN23", dirs)]

res = lapply(dirs, function(dir){
  files = list.files(paste0("/mnt/AchTeraD/data/BICRO258/", dir, "/trimmed"), pattern = ".bam$", full.names = T)
  
  filenames = gsub(".trimmed.*|.*\\/", "", files)
  
  # Get read counts
  param = ScanBamParam(what=c("qname"))
  
  reads = pblapply(1:length(files), function(i) {
    gr = granges(readGAlignments(files[i], param=param))
    return(data.table(sample = filenames[i], reads = length(gr), run = dir))
  }, cl = 30)
  
  tot_reads = rbindlist(reads)
})
# Merge all data
tot_reads = rbindlist(res)
tot_reads = merge(tot_reads, ct, by.x = "sample", by.y = "V2")

# Plot scatterplot
plt = ggplot(tot_reads[run == "TN22"], aes(x = N, y = reads)) +
  geom_point() +
  # geom_text_repel(aes(label = V1), max.overlaps = 15) +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 1) +
  scale_color_viridis_d() + 
  scale_y_continuous(labels = number_format()) + 
  scale_x_continuous() +
  stat_cor() +
  labs(y = "Reads", x = "Ct")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/scatterplot-reads/TN22-reads_CT-scatterplot-55samples-unpurified_10nl",
              height = 6, width = 7)
