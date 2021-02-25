## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "VennDiagram")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

muts = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/top20AAmutations-frontimicrob.tsv", header = T)

files = list.files("/mnt/AchTeraD/Documents/Projects/COVseq/data/calls", pattern = ".bed", recursive = T, full.names = T)
samples = gsub(".*calls\\/|.bed", "", files)
samples = gsub("\\/", "-", samples)

discard = c("sample17", "sample21", "sample29", "sample3", "sample30", "sample1", "sample2", "sample27", "sample18")
discard_mut = c("T16649A", "C16650A")
keep_ind = which(!gsub(".*-", "", samples) %in% discard)

samples = samples[keep_ind]
files = files[keep_ind]

# Get mutations info
res = lapply(1:length(files), function(x) {
  #dt = fread(files[x], select = c(1:3, 5:7), col.names = c("chr", "start", "end", "qual", "ref", "alt"))
  dt = fread(files[x])
  dt[, method := gsub("-.*", "", samples[x])]
  dt[, sample := gsub(".*-", "", samples[x])]
  dt[, depth := as.numeric(gsub(".*DP=|;.*", "", dt$V9))]
  
  # Get HQ alt/ref allele counts
  counts = dt[, gsub(".*DP4=|;.*", "", V9)]
  counts = tstrsplit(counts, ",")
  ref = as.numeric(counts[[1]]) + as.numeric(counts[[2]])
  alt = as.numeric(counts[[3]]) + as.numeric(counts[[4]])
  
  dt[, nalt := alt]
  dt[, nref := ref]
  dt = dt[, .(V1, V2, V3, V5, V6, V7, depth, method, sample, nalt, nref)]
  setnames(dt, c("chr", "start", "end", "qual", "ref", "alt", "depth", "method", "sample", "nalt", "nref"))
  return(dt)
})
total = rbindlist(res)

# Match actual starts/end to 1 indexed start site
total[, start := start + 1]
total[, end := end + 1]
total[, vaf := nalt / (nalt + nref)]

total[, id := factor(paste0(ref, start, alt))]

# Remove muts from sequencing artifacts
total = total[!id %in% discard_mut, ]

# Filter out indels (only keep snps) and low quality snps
total = total[nchar(ref) == 1 & nchar(alt) == 1 & qual >= 30 & depth >= 10 & nalt > nref]


total_mat = dcast(total, id ~ sample, fun.aggregate = function(x) paste(x, collapse = ", "), value.var = "method")
total_mat[total_mat == ""] = NA
total_m = melt(total_mat, id.vars = "id")
total_m = total_m[complete.cases(total_m)]

counts = total_m[, .N, by = .(id, value)]
order = counts[, sum(N), by = .(id)]
setorder(order, -V1)

counts[, id := factor(id, levels = order$id)]


# plt = ggplot(counts, aes(x = id, y = N, fill = value)) +
#   geom_bar(stat = "identity") +
#   scale_fill_viridis_d() +
#   labs(x = "Mutation", y = "Number of times mutation is called", fill = "Called by:") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/barplot-mutationscalled-posSamples-new",
#               width = 24, height = 12)

counts_muts = total[, .N, by = .(sample, method)]
counts_muts = dcast(counts_muts, sample~method)

plt1 = ggplot(counts_muts, aes(x = nebnext, y = covseq)) +
  geom_jitter() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2) +
  stat_cor() +
  coord_fixed(xlim = c(4, 15), ylim = c(4, 15))

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/scatterplot-numSNPs-COV-NEB-jitter-fixed",
              height = 7, width = 9)

# # Get VAFs
# total_vaf = merge(total_m, total[, .(nalt, nref, sample, id, method, vaf)], by.x = c("variable", "id"),  by.y = c("sample", "id"), all.x = T)
# 
# plt2 = ggplot(total_vaf, aes(x = value, y = vaf, color = method)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
#   labs(y = "VAF", x = "Method", color = "Method") +
#   scale_color_viridis_d()
# 
# save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/boxplot-mutationscalled-VAF",
#               width = 9, height = 7)
# 
# plt3 = ggplot(total_vaf, aes(x = value, y = depth, color = method)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
#   labs(y = "Depth", x = "Method", color = "Method") +
#   scale_color_viridis_d()
# 
# save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/boxplot-mutationscalled-depth",
#               width = 9, height = 7)
# log_ind = duplicated(total, by = c("chr", "start", "end", "ref", "alt", "sample"), fromLast = T) | duplicated(total, by = c("chr", "start", "end", "ref", "alt", "sample"))
# 
# common = total[log_ind, c("chr", "start", "end", "ref", "alt", "sample", "nalt", "nref", "method")]
# unique = total[!log_ind, c("chr", "start", "end", "ref", "alt", "method", "qual", "sample")]
# 
#
data = total_m[, .N, by = .(value)]
venn_data = data.table(covseq_nebnext = data$N[1],
                       covseq = data$N[1] + data$N[3],
                       nebnext = data$N[1] + data$N[2])

# Load VennDiagram library
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/VennDiagram_cairo-posSamples-new_ps.eps",
         onefile = TRUE, height=7, width=7, family="Helvetica",
         pointsize=8, antialias="none")
# Render a Venn Diagram
draw.pairwise.venn(
  venn_data$covseq,                            # Size of the left circle.
  venn_data$nebnext,                            # Size of the right circle.
  venn_data$covseq_nebnext,                        # Size of the overlapping area.
  category = c("COVseq", "NEBNext"),  # Label text.
  cat.pos = c(0, 0),                   # Position of labels.
  scaled=TRUE,                         # Scale the circle size or not.
)
dev.off()
# 
# common[, vaf := nalt / (nalt + nref)]
# common_m = melt(common[, c(1:6, 9:10)], id.vars = c("chr", "start", "end", "ref", "alt", "sample", "method"))
# common_m[, id := factor(paste0(ref, start, alt))]
# 
# freq = table(common_m$id)
# neworder = order(freq, decreasing = T)
# common_m[, id := factor(id, levels = names(freq)[neworder])]
# 
# plt = ggplot(common_m, aes(x = id, y = value, color = method)) +
#   stat_summary(fun = mean,
#                fun.min = function(x) mean(x) - sd(x), 
#                fun.max = function(x) mean(x) + sd(x), 
#                geom = "pointrange",
#                position = position_dodge(width = 0.5)) + 
#   scale_color_viridis_d() +
#   labs(y = "Fraction alt allele", x = "Mutation", color = "Method") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/dot-plot-VAF-commonmuts",
#               height = 7, width = 9)
