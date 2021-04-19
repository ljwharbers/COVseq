
## Load/install packages
packages = c("data.table", "UpSetR", "VennDiagram")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = c("NZ265", "NZ266", "NZ267")

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
variants[, tosum := 1]
variants[, sample_id := paste0(sample, "_", id)]

# Select <=35 Ct samples
variants = variants[!grepl("Neg", sample)]

variants_mat = dcast(variants, sample_id ~ library, fun.aggregate = function(x) sum(x), value.var = "tosum")


plt1 = upset(variants_mat, sets = libraries, nintersect = 20, keep.order = T, order.by = "freq", point.size = 3,
             line.size = 1)

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/upsetPlots/NZ265-NZ267-upsetR",
              width = 8, height = 5)
