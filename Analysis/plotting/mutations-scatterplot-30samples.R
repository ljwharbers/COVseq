## Author: Luuk Harbers
## Date: 2020-02-16
## Script for analysis of variants

## Load/install packages
packages = c("data.table", "UpSetR", "VennDiagram")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
libraries = libraries[grepl("MS147_miseq$|NEBNext", libraries)]

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
variants[, library := gsub("MS147_miseq", "COVseq", library)]
# Get unique variants
variants = unique(variants, by = c("POS", "REF", "ALT", "sample", "library"))
variants = variants[PASS == TRUE & ALT_FREQ > 0.75 & nchar(ALT) == 1 & TOTAL_DP > 15,]
#variants = variants[PASS == TRUE & nchar(ALT) == 1 & TOTAL_DP > 20,]

# Remove positive and negative samples and perhaps samples <= 35ct
variants = variants[!grepl("Sample_12|Sample_Pos|Sample_Neg|Sample_21|Sample_1$|Sample_30|Sample_3$|Sample_17|Sample_29|Sample_2$|Sample_27|Sample_18", sample),]

# Get counts
num_snv = variants[, .N, by = .(sample, library)]
num_snv = dcast(num_snv, sample ~ library)

plt = ggplot(num_snv, aes(x = NEBNext, y = COVseq)) + 
  geom_jitter(size = 4, width = .15, height = .15) +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2, size = 2, fullrange = TRUE) +
  stat_cor() +
  labs(y = "Number of SNVs detected (COVseq PE300)", x = "Number of SNVs detected (NEBNext)")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-muts/COVseq300PE-NEBNext-numMutations",
              height = 7, width = 7)
