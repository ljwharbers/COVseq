## Author: Luuk Harbers
## Date: 2021-05-17
## Script for plotting coverage at VoC locations

## Load/install packages
packages = c("data.table", "pbapply", "GenomicAlignments")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

muts = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/mutation-sourceinfo-updated3march.csv")

# Prepare mutations for overlapping
setkey(muts, chr, start, end)

files = list.files(paste0("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/MS147_miseq/variants/bam/"), 
                   pattern = ".bam$", full.names = T)


param = ScanBamParam(what=c("qname"))

res = pblapply(files, function(file) {
  gr = granges(readGAlignments(file, param=param))
  dt = as.data.table(gr)
  setkey(dt, seqnames, start, end)
  # Get overlaps
  overlaps = foverlaps(muts, dt)
  
  # Get reads per overlap and normalize by library size
  counts = overlaps[, .N, by = Event]
  counts[, normalized := (N / nrow(dt)) * 1e6]
  return(data.table(sample = gsub(".*\\/|.trim.sorted.bam", "", file),
                    mutation = counts$Event,
                    coverage = counts$N,
                    normalized_coverage = counts$normalized))
}, cl = 32)
total = rbindlist(res)

# Remove files from negative samples
total = total[!grepl("Sample_12|Sample_Pos|Sample_Neg|Sample_21|Sample_1$|Sample_30|Sample_3$|Sample_17|Sample_29|Sample_2$|Sample_27|Sample_18", sample),]
total = merge(total, muts, by.x = "mutation", by.y = "Event")
total[, mutation_ann := paste0(mutation, " (", ann, ", ", gene, ")")]

# Get mutation order
mut_order = total[, mean(normalized_coverage), by = mutation_ann]
setorder(mut_order, V1)
total[, mutation_ann := factor(mutation_ann, levels = mut_order$mutation_ann)]
total[, sample := gsub("_", " ", sample)]

muts[, mutation_ann := paste0(Event, " (", ann, ", ", gene, ")")]
muts[, mutation_ann := factor(mutation_ann, levels = mut_order$mutation)]

# Remove Sample prefix and set coverage < 10 to NA
total[, sample := gsub("Sample ", "", sample)]
total[, sample := factor(sample, levels = as.character(sort(unique(as.numeric(total$sample)))))]
total[coverage < 15, coverage := NA]

annot_plt = ggplot(muts, aes(x = 1, y = mutation_ann, fill = source)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  labs(fill = "Variant") +
  theme_void() +
  theme(legend.position = "left",
        axis.text.y = element_text(hjust = .95))

# plot
heatmap = ggplot(total, aes(x = sample, y = mutation_ann, fill = coverage)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c(limits = c(0, 5e2), oob = scales::squish, labels = c(0, 100, 200, 300, 400, "> 500")) +
  labs(y = "Mutation", x = "OAS-29 samples", fill = "Variant coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank())

plt = plot_grid(annot_plt, heatmap, ncol = 2, rel_widths = c(.225, 1), align = "h")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/mutations/heatmaps/PE300-non-normalized_Coverage-VoC-muts_15+",
              width = 20, height = 12)

means = total[, mean(normalized_coverage), by = mutation]
sds = total[, sd(normalized_coverage), by = mutation]
