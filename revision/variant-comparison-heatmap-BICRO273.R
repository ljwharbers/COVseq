## Author: Luuk Harbers
## Date: 2020-02-16
## Script for analysis of variants

## Load/install packages
packages = c("data.table", "UpSetR", "stringr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
libraries = libraries[grepl("NZ244|NZ245", libraries)]

# Lineage info
lin = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/pangolin/NZ244/lineage_report.csv")

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

variants[, id := factor(paste0(REF, POS, ALT))]
variants[, sample_id := paste0(sample, "_", id)]
variants = variants[!grepl("Sample_N|Sample_P|Sample_52|Sample_54", sample)]
lin = lin[!grepl("Sample_N|Sample_P|Sample_52|Sample_54", taxon)]

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

# Plot annotation
lin[lineage == "None", lineage := NA]
lin[, taxon := gsub("_", " ", taxon)]
lin[, taxon := factor(taxon, levels = mat$sample[col_clust])]
annot_plt = ggplot(lin, aes(x = taxon, y = 1, fill = lineage)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  theme_void()

plt = ggplot(variants_collapse, aes(x = sample, y = id, fill = library)) + 
  geom_tile(color = "black", size = .5) +
  scale_fill_viridis_d(option = "E", end = 1, begin = 0.2, na.value = "white", drop = F) +
  labs(y = "Mutation ID", x = "Sample", fill = "Detected in:") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_blank())

# Combine plots and save
combined_noleg = plot_grid(plt, annot_plt + theme(legend.position = ""),
                           align = "v", axis = "lr", rel_heights = c(1, 0.1), ncol = 1)

legend = cowplot::get_legend(
  annot_plt + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

combined = plot_grid(combined_noleg, legend, ncol = 1, rel_heights = c(1, 0.05))

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/mutations/heatmaps/BICRO273-NZ244+NZ245-55samples-lineage",
              width = 18, height = 14)

plt1 = upset(variants_mat, sets = libraries, nintersect = 20, keep.order = T, order.by = "freq", point.size = 3,
             line.size = 1)
save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/upsetPlots/BICRO273-all_libs",
              width = 12, height = 8)
