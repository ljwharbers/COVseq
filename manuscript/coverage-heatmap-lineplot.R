## Author: Luuk Harbers
## Date: 2020-11-24
## Script for plotting Coverage of SARS-CoV-2 genomes

## Load/install packages
packages = c("data.table", "ggplot2", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 32
binsize = 100
# 
files = list.files("/mnt/AchTeraD/data/BICRO248+249/TN11_S2/trimmed/", pattern = ".bam$", full.names = T)
files = files[grepl("sample", files)]
files = c(files, list.files("/mnt/AchTeraD/data/BICRO237/covseq/", recursive = T, pattern = "trimmed.sorted.bam$", full.names = T),
          list.files("/mnt/AchTeraD/data/BICRO240/covseq/", recursive = T, pattern = "trimmed.sorted.bam$", full.names = T))

ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)

samplenames = gsub(".*/|.bam", "", files)
samplenames = gsub("sample" ,"COVseq - Sample ", samplenames)
samplenames[31:60] = paste0("NEBNext - Sample ", 1:30)

# Set params
param = ScanBamParam(what=c("qname","flag"))

res = pblapply(1:length(files), function(x) {
  samplename = samplenames[x]
  gr = granges(readGAlignments(files[x], param=param))
  nreads = length(gr)
  cov = unlist(coverage(gr))
  # Normalize by nreads
  cov@values = cov@values / nreads
  dt = data.table(cov@values, cov@lengths)
  values = unlist(apply(dt, 1, function(x) rep(x[1], x[2])))
  
  res = data.table(coverage = values)
  setnames(res, samplename)
  
  return(res)
}, cl = nthreads)


total = do.call(cbind, res)
# binned = total[, as.list(colSums(.SD)), by = gl(ceiling(nrow(total)/binsize), binsize, nrow(total))]
# 
# binned_m = melt(binned, id.vars = "gl")
# binned_m[, variable := factor(variable, levels = rev(ct$sample[order(ct$ct)]))]
# 
# plt1 = ggplot(binned_m, aes(x = gl, y = variable, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient2(high = "red") +
#   labs(y = "", x = "SARS-CoV-2", fill = "Coverage per bin") +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# 
# save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/coverage-heatmap/NEBNext-heatmap-100bp-normalized",
#               width = 13, height = 7)

# calculate ct for ordering
# Select samples
samples = c("COVseq - Sample 5" , "NEBNext - Sample 5", "COVseq - Sample 17", "NEBNext - Sample 17")
total = total[, ..samples]

total[, base := factor(1:nrow(total), levels = 1:nrow(total))]

total_m = melt(total, id.vars = "base")
total_m[, method := ifelse(grepl("COVseq", variable), "COVseq", "NEBNext")]
total_m[, quality := ifelse(grepl("Sample 5", variable), "Likely negative", "Likely positive")]

#total_m[, variable := factor(variable, levels = rev(ct$sample[order(ct$ct)]))]

ggplot(total_m, aes(x = base, y = value, color = method, linetype = quality, group = variable)) +
  geom_path(size = 1.5) +
  scale_color_viridis_d() +
  labs(y = "Coverage", x = "SARS-CoV-2", color = "Sample") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/coverage-lineplots/COVseq-coverage-non-normalized-log2",
              width = 15, height = 7)
