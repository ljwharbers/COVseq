## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# files = c(list.files("/mnt/AchTeraD/data/BICRO237/covseq/fastq/", pattern = "screen.txt", full.names = T),
#           list.files("/mnt/AchTeraD/data/BICRO240/covseq/fastq/", pattern = "screen.txt", full.names = T))
files = list.files("/mnt/AchTeraD/data/BICRO268/MS147/demultiplexed/", pattern = "screen.txt", full.names = T)
files = files[!grepl("unassigned", files)]
annot = fread("/mnt/AchTeraD/data/BICRO268/MS146+147-barcodes-annot.txt", header = F)

ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv")

res = lapply(1:length(files), function(x) {
  # samplename = x
  samplename = gsub(".*/|_screen.*", "", files[x])
  dt = fread(files[x], select = c(1, 2, 4, 6, 8, 10, 12))
  # Remove adapters
  dt = dt[Genome != "Adapters"]
  
  data.table(sample = samplename,
             genome = dt$Genome,
             hits = dt$`%One_hit_one_genome` + dt$`%Multiple_hits_one_genome`)
})

total = rbindlist(res)

total[genome == "Human", colors := "human"]
total[genome == "SARS-CoV-2", colors := "SARS-CoV-2"]
total[is.na(colors), colors := "other"]

# Merge with samplenames
total[, sample := gsub("_R.", "", sample)]
total = total[, .(hits = mean(hits)), by = .(sample, genome, colors)]

# add unmapped
unmapped = total[, list(genome = "unmapped", hits=100-sum(hits), colors = "unmapped"), by = sample]

total = rbind(total, unmapped)


# Get samplenames
total = merge(total, annot, by.x = "sample", by.y = "V1")
total[, sample := V2]
total = total[!grepl("Neg|Pos|Sample 12", sample)]

# Order both on ct value
total[, sample := factor(sample, levels = ct$sample[order(ct$ct)])]
total[, colors := factor(colors, levels = c("other", "unmapped", "human", "SARS-CoV-2"))]

# Prepare Ct line
ct[, ct := ct*2]
ct = ct[!grepl("Sample 12", sample)]
ct[, sample := factor(sample, levels = levels(total$sample))]


plt = ggplot(total, aes(x = sample, y = hits)) +
  geom_col(aes(fill = colors)) +
  geom_line(data = ct, aes(x = sample, y = ct, group = 1), color = "red", size = 1.5) +
  scale_fill_viridis_d(option = "D", direction = -1) +
  labs(y = "Percentage of reads", x = "", fill = "") +
  scale_y_continuous(expand = c(0, 0), sec.axis = sec_axis(~.*0.5, name = "CT value")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/mapping/COVseq-PE150-genome-percentages-CT_MS147",
              width = 9, height = 6)
