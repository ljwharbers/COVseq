## Author: Luuk Harbers
## Date: 2020-02-16
## Script for analysis of variants

## Load/install packages
packages = c("data.table", "UpSetR", "stringr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
libraries = libraries[grepl("MS166|NZ248", libraries)]

# Load annot
old = fread("/mnt/AchTeraD/data/BICRO271/MS166-annot_old.tsv", header = F, col.names = c("barcode", "old"))
new = fread("/mnt/AchTeraD/data/BICRO271/MS166-annot.tsv", header = F, col.names = c("barcode", "new", "location", "person"))
annot = merge(old, new)

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

# Give new names
variants = merge(variants, annot, by.x = "sample", by.y = "old")
variants[, sample := new]

# Get unique variants
variants = unique(variants, by = c("POS", "REF", "ALT", "sample", "library"))
variants = variants[PASS == TRUE & ALT_FREQ > 0.75 & nchar(ALT) == 1 & TOTAL_DP > 15,]
variants = variants[!grepl("Sample_2$|Sample_Neg", sample)]



variants[, id := factor(paste0(REF, POS, ALT))]
variants[, sample_id := paste0(sample, "_", id)]

variants_collapse = variants[, .(library = paste(library, collapse = ", "),
                                 sample = sample,
                                 id = id), by = sample_id]

# Reorder factors
variants_collapse[, library := factor(library, levels = unique(library))]
mut_order = variants_collapse[, .N, by = id]
setorder(mut_order, -N)
variants_collapse[, id := factor(id, levels = mut_order$id)]
variants_collapse[, sample := gsub("_", " ", sample)]

# Set NAs
variants_collapse = tidyr::complete(variants_collapse, sample, id)
setDT(variants_collapse)

# Count number of libraries
variants_collapse[, nlibraries := sapply(variants_collapse$library, function(x) str_count(x, ",")) + 1]
variants_collapse[is.na(nlibraries), nlibraries := 0]
variants_collapse[, nlibraries := factor(nlibraries, levels = c("4", "3", "2", "1", "0"))]

mat = dcast(variants_collapse, sample ~ id, value.var = "nlibraries", fun.aggregate = function(x) sum(as.numeric(as.character(x))))
col_clust = hclust(dist(mat[, -1]))$order
row_clust = hclust(dist(t(mat[, -1])))$order

# Set order
variants_collapse[, sample := factor(sample, mat$sample[col_clust])]
variants_collapse[, id := factor(id, colnames(mat[, -1])[row_clust])]

# Make annotation plot
annot = annot[!grepl("Sample_2$|Sample_Neg", new)]
annot[, sample := gsub("_", " ", new)]
annot[, sample := factor(sample, levels = mat$sample[col_clust])]
annot_plt_loc = ggplot(annot, aes(x = sample, y = 1, fill = location)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

annot_plt_aff = ggplot(annot, aes(x = sample, y = 1, fill = person)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")

plt = ggplot(variants_collapse, aes(x = sample, y = id, fill = library)) + 
  geom_tile(color = "black", size = .5) +
  scale_fill_viridis_d(option = "E", end = 1, begin = 0.2, na.value = "white", drop = F) +
  labs(y = "Mutation ID", x = "Sample", fill = "Detected in:") +
  theme(axis.line = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

combined = plot_grid(plt, annot_plt_loc, annot_plt_aff, align = "v", axis = "lr", ncol = 1, rel_heights = c(1, 0.065, 0.15))

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/mutations/heatmaps/MS166-NZ248_overlap-newnames-location-person",
              width = 24, height = 16)
