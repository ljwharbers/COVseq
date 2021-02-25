## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "pbapply", "GenomicAlignments")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

#libraries = c("TN11", "TN21", "NZ190_S1", "NZ191_S2", "MS82_S3", "MS83_S4", "MS84_S5")
libraries = c("TN11", "TN21")
#runs = c(rep("BICRO251", 2), rep("BICRO254", 5))
runs = c("BICRO251", "BICRO251")
annot = fread("/mnt/AchTeraD/data/BICRO251/TN21-annot_barcodes.txt")
muts = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/mutation-sourceinfo.tsv")

# Prepare mutations for overlapping
setkey(muts, chr, start, end)

files = lapply(1:length(libraries), function(lib) {
  files = list.files(paste0("/mnt/AchTeraD/data/", runs[lib], "/", libraries[lib], "/trimmed"), pattern = "bam$", full.names = T)
  data.table(files = files, library = libraries[lib])
})
files = rbindlist(files)
files = files[library != "TN11" | grepl("sample", files)]
files = files[!grepl("fixed", files)]
files[, name := gsub(".trimmed.*|.bam", "", basename(files))]

# Remove files from negative samples
discard_rest = c("Sample 22", "Sample 24", "Sample N1", "Sample N2", "Sample N3", "Sample P1", "Sample P2")
discard_tn11 = c("sample17", "sample21", "sample29", "sample3", "sample30", "sample1", "sample2", "sample27", "sample18")

files = merge(files, annot, by.x = "name", by.y = "V2", all.x = T, sort = F)
files[, V1 := paste0("Sample ", V1)]
files = files[!name %in% discard_tn11 & !V1 %in% discard_rest]
files[, name := gsub("sample", "Sample ", name)]
files[library != "TN11", name := V1]

param = ScanBamParam(what=c("qname"))

res = pblapply(1:nrow(files), function(i) {
  gr = granges(readGAlignments(files$files[i], param=param))
  dt = as.data.table(gr)
  setkey(dt, seqnames, start, end)
  # Get overlaps
  overlaps = foverlaps(muts, dt)
  
  # Get reads per overlap and normalize by library size
  counts = overlaps[, .N, by = Event]
  counts[, normalized := (N / nrow(dt)) * 1e6]
  return(data.table(sample = files$name[i], library = files$library[i],
                    mutation = counts$Event, normalized_coverage = counts$normalized))
}, cl = 30)
total = rbindlist(res)
total[, id := paste(library, sample, sep = "_")]

# Prepare total dt for annotation
total[, cohort := ifelse(library == "TN11", "30 samples", "55 samples")]
# total[library == "TN21" | library == "TN11", replicate := "Replicate 1"]
# total[library == "NZ190_S1", replicate := "Replicate 2"]
# total[library == "NZ191_S2", replicate := "Replicate 3"]
# total[library == "MS82_S3", replicate := "Replicate 4"]
# total[library == "MS83_S4", replicate := "Replicate 5"]
# total[library == "MS84_S5", replicate := "Replicate 6"]

# Reorder
setorder(total, cohort)
total[, id := factor(id, levels = unique(id))]

# Get mean for replicates
total_mean = total[, .(normalized_coverage = mean(normalized_coverage)), by = .(sample, mutation, cohort)]
total_mean[, id := paste(sample, cohort)]
setorder(total, cohort)
total_mean[, id := factor(id, levels = unique(id))]

# Merge with mutation annot
total_mean = merge(total_mean, muts, by.x = "mutation", by.y = "Event")
#total_mean = merge(total_mean, muts, by.x = "mutation", by.y = "Event", allow.cartesian = T)

total_mean[, mutation_ann := paste0(mutation, " (", ann, ", ", gene, ")")]

# Get mutation order
mut_order = total_mean[, mean(normalized_coverage), by = mutation_ann]
setorder(mut_order, V1)
total_mean[, mutation_ann := factor(mutation_ann, levels = mut_order$mutation_ann)]
muts[, Event := factor(Event, levels = mut_order$mutation)]



annot_plt1 = ggplot(total_mean, aes(x = id, y = 1, fill = cohort)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d(end = 0.4) +
  theme_void() 

# annot_plt2 = ggplot(total, aes(x = id, y = 1, fill = replicate)) +
#   geom_bar(stat = "identity",
#            width = 1) +
#   scale_fill_viridis_d() +
#   theme_void() 
# 
annot_plt3 = ggplot(muts, aes(x = 1, y = Event, fill = source)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d(begin = 0.4, end = 0.8) +
  theme_void() +
  theme(legend.position = "left")

# plot
heatmap = ggplot(total_mean, aes(x = id, y = mutation_ann, fill = normalized_coverage)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1e3), oob = scales::squish, labels = c("0", "250", "500", "750", "> 1000")) +
  labs(y = "Mutation", fill = "Normalized Coverage") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

plt = plot_grid(heatmap, annot_plt1, ncol = 1, rel_widths = c(0.5, 1), 
                rel_heights = c(1, 0.05), align = "hv", axis = "lrtb")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/potential-coverage-topmuts-noreps",
              width = 20, height = 12)

means = total_mean[, mean(normalized_coverage), by = mutation]
sds = total_mean[, sd(normalized_coverage), by = mutation]
