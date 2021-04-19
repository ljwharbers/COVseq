## Author: Luuk Harbers
## Date: 2020-02-16
## Script for analysis of variants

## Load/install packages
packages = c("data.table", "UpSetR", "VennDiagram")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
libraries = libraries[grepl("MS147$|NEBNext", libraries)]

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
variants[, library := gsub("MS147", "COVseq", library)]
# Get unique variants
variants = unique(variants, by = c("POS", "REF", "ALT", "sample", "library"))
variants = variants[PASS == TRUE & ALT_FREQ > 0.75 & nchar(ALT) == 1 & TOTAL_DP > 15,]
#variants = variants[PASS == TRUE & nchar(ALT) == 1 & TOTAL_DP > 20,]

# Remove positive and negative samples (and perhaps samples <= 35ct)
variants = variants[!grepl("Sample_12|Sample_Pos|Sample_Neg|Sample_21|Sample_1$|Sample_30|Sample_3$|Sample_17|Sample_29|Sample_2$|Sample_27|Sample_18", sample),]
#variants = variants[!grepl("Sample_Pos|Sample_Neg", sample),]

# Get id and dcast
variants[, id := factor(paste0(REF, POS, ALT))]
variants[, sample_id := paste0(sample, "_", id)]

variants_collapse = variants[, .(library = paste(library, collapse = ", "),
                                 sample = sample,
                                 id = id), by = sample_id]

# Reorder factors
variants_collapse[, library := factor(library, levels = c("COVseq, NEBNext", "NEBNext", "COVseq"))]
mut_order = variants_collapse[, .N, by = id]
setorder(mut_order, -N)
variants_collapse[, id := factor(id, levels = mut_order$id)]
variants_collapse[, sample := gsub("_", " ", sample)]

# Set NAs
variants_collapse = tidyr::complete(variants_collapse, sample, id)
setDT(variants_collapse)

# Set order
order = variants_collapse[!is.na(library), .N, by = sample]
setorder(order, -N)
variants_collapse[, sample := factor(sample, levels = order$sample)]

plt = ggplot(variants_collapse, aes(x = sample, y = id, fill = library)) + 
  geom_tile(color = "black", size = .5) +
  scale_fill_viridis_d(option = "E", end = 0.8, begin = 0.2, na.value = "white", drop = F) +
  labs(y = "Mutation ID", x = "Sample", fill = "Detected in:") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_blank())

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/mutations/heatmaps/MS147-NEBNext_overlap",
              width = 12, height = 10)
