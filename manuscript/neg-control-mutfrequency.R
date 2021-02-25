## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "VennDiagram", "pbapply", "UpSetR")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 32

muts = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/top20AAmutations-frontimicrob.tsv", header = T)
annot = fread("/mnt/AchTeraD/data/BICRO251/TN21-annot_barcodes.txt")

dirs = list.dirs("/mnt/AchTeraD/data/BICRO254/", recursive = F)
dirs = dirs[!grepl("fastq", dirs)]


res = lapply(dirs, function(x){
  samplename = gsub(".*\\/", "", x)
  files = list.files(x, pattern = ".bed", recursive = T, full.names = T)
  samples = gsub(".trimmed.*|.*\\/", "", files)
  
  subres = pblapply(1:length(files), function(x) {
    #dt = fread(files[x], select = c(1:3, 5:7), col.names = c("chr", "start", "end", "qual", "ref", "alt"))
    dt = fread(files[x])
    dt[, run := samplename]
    dt[, sample := samples[x]]
    dt[, depth := as.numeric(gsub(".*DP=|;.*", "", dt$V9))]
    
    # Get HQ alt/ref allele counts
    counts = dt[, gsub(".*DP4=|;.*", "", V9)]
    counts = tstrsplit(counts, ",")
    ref = as.numeric(counts[[1]]) + as.numeric(counts[[2]])
    alt = as.numeric(counts[[3]]) + as.numeric(counts[[4]])
    
    dt[, nalt := alt]
    dt[, nref := ref]
    dt = dt[, .(V1, V2, V3, V5, V6, V7, depth, run, sample, nalt, nref)]
    setnames(dt, c("chr", "start", "end", "qual", "ref", "alt", "depth", "run", "sample", "nalt", "nref"))
    return(dt)
  }, cl = nthreads)
  total = rbindlist(subres)
})

total = rbindlist(res)

# Merge with sample names
total = merge(total, annot, by.x = "sample", by.y = "V2")
total[, sample := paste0("Sample ", V1)]
total$V1 = NULL

# Match actual starts/end to 1 indexed start site
total[, start := start + 1]
total[, end := end + 1]

# Filter out indels (only keep snps) and low quality snps
#total = total[nchar(ref) == 1 & nchar(alt) == 1 & qual >= 30 & depth >= 10 & nalt > nref]
total = total[nchar(ref) == 1 & nchar(alt) == 1]

total[, id := factor(paste0(ref, start, alt))]
total[, vaf := nalt / (nalt + nref)]
total[, count := 1]

# Filter pos/neg and get only neg 
total_samples = total[!grepl("N|P", sample)]
total_neg = total[grepl("N", sample)]

total_samples = total_samples[nchar(ref) == 1 & nchar(alt) == 1 & qual >= 30 & depth >= 10 & nalt > nref]
counts = total_samples[, .N, by = .(id, run)]

total_neg_merged = merge(total_neg, counts, by = c("id", "run"))

ggplot(total_neg_merged, aes(x = vaf, y = N)) +
  facet_wrap(~run + sample, ncol = 3) +
  geom_point(size = 2) +
  scale_color_viridis_d() +
  stat_cor() +
  labs(y = "Frequency of mutation in COVseq 55 samples", x = "VAF in negative control")




order = counts[, sum(N), by = .(id)]
# Get matrix
total_mat = dcast(total, id ~ run, fun.aggregate = function(x) sum(x), value.var = "count")
upset(as.data.frame(total_mat), nsets = 5)


counts = total[, .N, by = .(id, run)]
setorder(counts, -N)
counts[, id := factor(id, levels = id)]


plt = ggplot(counts, aes(x = id, y = N, fill = value)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(x = "Mutation", y = "Number of times mutation is called", fill = "Called by:") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/barplot-mutationscalled",
              width = 24, height = 12)



# Get VAFs
total_vaf = merge(total_m, total[, .(nalt, nref, sample, id, method)], by.x = c("variable", "id"),  by.y = c("sample", "id"), all.x = T)
total_vaf[, vaf := nalt / (nalt + nref)]
total_vaf[, depth := nalt+nref]

plt2 = ggplot(total_vaf, aes(x = value, y = vaf, color = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
  labs(y = "VAF", x = "Method", color = "Method") +
  scale_color_viridis_d()

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/boxplot-mutationscalled-VAF",
              width = 9, height = 7)

plt3 = ggplot(total_vaf, aes(x = value, y = depth, color = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
  labs(y = "Depth", x = "Method", color = "Method") +
  scale_color_viridis_d()

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/boxplot-mutationscalled-depth",
              width = 9, height = 7)
# log_ind = duplicated(total, by = c("chr", "start", "end", "ref", "alt", "sample"), fromLast = T) | duplicated(total, by = c("chr", "start", "end", "ref", "alt", "sample"))
# 
# common = total[log_ind, c("chr", "start", "end", "ref", "alt", "sample", "nalt", "nref", "method")]
# unique = total[!log_ind, c("chr", "start", "end", "ref", "alt", "method", "qual", "sample")]
# 
#
data = total_m[, .N, by = .(value)]
venn_data = data.table(covseq_nebnext = data$N[2],
                       covseq = data$N[2] + data$N[1],
                       nebnext = data$N[2] + data$N[3])

# Load VennDiagram library
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/VennDiagram_cairo_ps.eps",
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
