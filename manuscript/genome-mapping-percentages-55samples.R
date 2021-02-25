## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# files = c(list.files("/mnt/AchTeraD/data/BICRO237/covseq/fastq/", pattern = "screen.txt", full.names = T),
#           list.files("/mnt/AchTeraD/data/BICRO240/covseq/fastq/", pattern = "screen.txt", full.names = T))
files = list.files("/mnt/AchTeraD/data/BICRO254/NZ191_S2/demultiplexed/", pattern = "screen.txt", full.names = T)
#files = files[grepl("sample", files)]
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/55-patients_CT-clean.tsv")
annot = fread("/mnt/AchTeraD/data/BICRO251/TN21-annot_barcodes.txt")

res = lapply(1:length(files), function(x) {
  # samplename = x
  samplename = gsub(".*/|_screen.*", "", files[x])
  samplename = gsub("sample" ,"Sample ", samplename)
  dt = fread(files[x], select = c(1, 2, 4, 6, 8, 10, 12))
  data.table(sample = samplename,
             genome = dt$Genome,
             hits = dt$`%One_hit_one_genome` + dt$`%Multiple_hits_one_genome`)
})

total = rbindlist(res)
total = merge(total, annot, by.x = "sample", by.y = "V2")
total = total[!grepl("N|P", V1)]
total[, sample := paste0("Sample ", V1)]
total$V1 = NULL

total[genome == "Human", colors := "human"]
total[genome == "SARS-CoV-2", colors := "SARS-CoV-2"]
total[is.na(colors), colors := "other"]

# add unmapped
unmapped = total[, list(genome = "unmapped", hits=100-sum(hits), colors = "unmapped"), by = sample]

total = rbind(total, unmapped)
total_m = melt(total, id.vars = c("sample", "genome", "colors"))
total_m[, colors := factor(colors, levels = c("other", "unmapped", "human", "SARS-CoV-2"))]

# Merge with sample names
total_m = total_m[sample != "unassigned"]
#total_m[, sample := factor(paste0("Sample ", V1))]

# Order both on ct value
total_m[, sample := factor(sample, levels = ct$sample[order(ct$N)])]

ct[, sample := factor(sample, levels = levels(total_m$sample))]
ct[, ct := N*2]
plt = ggplot(total_m, aes(x = sample, y = value)) +
  geom_col(aes(fill = colors)) +
  geom_line(data = ct, aes(x = sample, y = ct, group = 1), color = "red", size = 1.5) +
  scale_fill_viridis_d(option = "D", direction = -1) +
  labs(y = "Percentage of reads", x = "", fill = "") +
  scale_y_continuous(expand = c(0, 0), sec.axis = sec_axis(~.*0.5, name = "CT value")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mapping/SE150-COVseq-55-samples-genome-percentages-CT_NZ191",
              width = 12, height = 6)
