## Author: Luuk Harbers
## Date: 2020-02-16
## Script for analysis of variants

## Load/install packages
packages = c("data.table", "UpSetR", "stringr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
libraries = libraries[grepl("NZ265|NZ266|NZ267", libraries)]

lin = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/NZ266/lineage_report.csv")
name_relabel = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/relabeling-table-AOS95.tsv")

# Load in variants
variants = lapply(libraries, function(lib) {
  files = list.files(paste0(basedir, "/", lib, "/variants/ivar/"), full.names = T, pattern = ".tsv")
  variants = lapply(files, function(file) {
    dt = fread(file)
    dt[, sample := gsub(".tsv", "", basename(file))]
    return(dt)
  })
  variants = rbindlist(variants)
  variants[, library := lib]
  return(variants)
})

variants = rbindlist(variants)

# Get unique variants
variants = unique(variants, by = c("POS", "REF", "ALT", "sample", "library"))
variants = variants[PASS == TRUE & ALT_FREQ > 0.75 & nchar(ALT) == 1 & TOTAL_DP > 15,]
variants = variants[!grepl("Sample_Neg", sample)]
variants[, id := factor(paste0(REF, POS, ALT))]
variants[, sample_id := paste0(sample, "_", id)]
variants[, sample := gsub("_", " ", sample)]
lin = lin[!grepl("Sample_Neg", taxon)]
lin[, taxon := gsub("_", " ", taxon)]


# RELABEL
name_relabel[, Oldname := paste0("Sample ", Oldname)]
lin = merge(lin, name_relabel, by.x = "taxon", by.y = "Oldname")
variants = merge(variants, name_relabel, by.x = "sample", by.y = "Oldname")
variants[, sample := Newname]
lin[, taxon := Newname]

variants_collapse = variants[, .(library = paste(library, collapse = ", "),
                                 sample = sample,
                                 id = id), by = sample_id]

# Reorder factors
variants_collapse[, library := factor(library, levels = unique(library))]
mut_order = variants_collapse[, .N, by = id]
setorder(mut_order, -N)
variants_collapse[, id := factor(id, levels = mut_order$id)]

# Set NAs
variants_collapse = tidyr::complete(variants_collapse, sample, id)
setDT(variants_collapse)

# Count number of libraries
variants_collapse[, nlibraries := sapply(variants_collapse$library, function(x) str_count(x, ",")) + 1]
variants_collapse[is.na(nlibraries), nlibraries := 0]

mat = dcast(variants_collapse, sample ~ id, value.var = "nlibraries", fun.aggregate = function(x) sum(as.numeric(as.character(x))))
col_clust = hclust(dist(mat[, -1]))$order
row_clust = hclust(dist(t(mat[, -1])))$order

# Set order
variants_collapse[, sample := factor(sample, mat$sample[col_clust])]
variants_collapse[, id := factor(id, colnames(mat[, -1])[row_clust])]

# Make annotation plot
lin[, taxon := factor(taxon, levels = levels(variants_collapse$sample))]



lin_plt = ggplot(lin, aes(x = taxon, y = 1, fill = lineage)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  labs(fill = "Pangolin Lineage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

variants_collapse[nlibraries == 0, nlibraries := NA]
variants_collapse[, nlibraries := factor(nlibraries, levels = c("3", "2", "1"))]

plt = ggplot(variants_collapse, aes(x = sample, y = id, fill = nlibraries)) + 
  geom_tile(color = "black", size = .5) +
  scale_fill_viridis_d(option = "E", end = 1, begin = 0.2, na.value = "white", drop = F) +
  #scale_fill_hue(na.value = "white") +
  labs(y = "Mutation ID", x = "Sample", fill = "Detected in:") +
  theme(axis.line = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

combined = plot_grid(plt, lin_plt, align = "v", axis = "lr", ncol = 1, rel_heights = c(1, 0.16))

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/mutations/heatmaps/NZ265-NZ267-96samples_triple-replicate-nlibraries",
              width = 24, height = 16)

# Plot scatterplots and correlation heatmap
variants_num = variants[, .N, by = .(library, sample)]
variants_num = dcast(variants_num, sample ~ library, value.var = "N")

plt1 = ggplot(variants_num, aes(x = NZ265, y = NZ266)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(method = "lm", linetype = 2, color = "red", se = F) +
  labs(x = "Number of SNVs (NZ265)", y = "Number of SNVs (NZ266)") +
  stat_cor()

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-muts/96-triplicate-NZ265-NZ266",
              height = 7, width = 8)

plt2 = ggplot(variants_num, aes(x = NZ265, y = NZ267)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(method = "lm", linetype = 2, color = "red", se = F) +
  labs(x = "Number of SNVs (NZ265)", y = "Number of SNVs (NZ267)") +
  stat_cor()

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-muts/96-triplicate-NZ265-NZ267",
              height = 7, width = 8)

plt3 = ggplot(variants_num, aes(x = NZ266, y = NZ267)) +
  geom_point(size = 3, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_smooth(method = "lm", linetype = 2, color = "red", se = F) +
  labs(x = "Number of SNVs (NZ266)", y = "Number of SNVs (NZ267)") +
  stat_cor()

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-muts/96-triplicate-NZ266-NZ267",
              height = 7, width = 8)

# Plot correlation heatmap
cors = cor(variants_num[, 2:4], method = "pearson")
cors = melt(cors)

plt4 = ggplot(cors, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_x_discrete("", position = "top") +
  scale_fill_viridis("Pearson's\ncorrelation", option="inferno", begin = .8, direction = -1) +
  labs(y = "", x = "") +
  theme(axis.line = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.text.x.top = element_text(angle = 45, hjust = 0))

save_and_plot(plt4, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/corr-snvs/96-triplicate-NZ266-NZ267-numsnvs",
              height = 7, width = 8)
