## Author: Luuk Harbers
## Date: 2020-02-16
## Script for analysis of variants

## Load/install packages
packages = c("data.table", "UpSetR", "VennDiagram")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
#libraries = libraries[!grepl("nextflow|work", libraries)]
libraries = c("NEBNext", "MS146", "MS147")

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
#variants = variants[PASS == TRUE & nchar(ALT) == 1 & TOTAL_DP > 20,]

# Remove positive and negative samples (and perhaps samples <= 35ct)
#variants = variants[!grepl("Sample_Pos|Sample_Neg|Sample_10", sample),]
variants = variants[!grepl("Sample_Pos|Sample_Neg|Sample_21|Sample_1|Sample_30|Sample_3|Sample_17|Sample_29|Sample_2|Sample_27|Sample_18", sample),]

# Get id and dcast
variants[, id := factor(paste0(REF, POS, ALT))]
variants[, tosum := 1]
variants[, sample_id := paste0(sample, "_", id)]
variants_mat = dcast(variants, sample_id ~ library, fun.aggregate = function(x) sum(x), value.var = "tosum")

# check exclusive variants
nebonly = variants_mat[NEBNext == 1 & MS146 == 0 & MS147 == 0,]
replicates_total = variants_mat[MS146 == 1 | MS147 == 1 ,]

plt1 = upset(variants_mat, sets = libraries, nintersect = 10, keep.order = T, order.by = "freq", point.size = 3,
             line.size = 1)

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/upsetPlots/MS146-MS147-replicates-overlap_ct<=35",
              height = 5, width = 8)

# plt2 = upset(variants_mat[, c(1, 3:4)], sets = libraries, nintersect = 10, keep.order = T, order.by = "freq", point.size = 3,
#              line.size = 1)
# 
# save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/upsetPlots/MS147-replicates-overlap_ct<=35",
#               height = 5, width = 8)

# Make Venn diagrams
MS146 = variants_mat[, .(sample_id, NEBNext, MS146)]
MS147 = variants_mat[, .(sample_id, NEBNext, MS147)]


MS146_dt = data.table(NEBNext = nrow(MS146[NEBNext == 1 & MS146 == 0, ]),
                      MS146 = nrow(MS146[NEBNext == 0 & MS146 == 1, ]),
                      shared = nrow(MS146[NEBNext == 1 & MS146 == 1, ]))

MS147_dt = data.table(NEBNext = nrow(MS147[NEBNext == 1 & MS147 == 0, ]),
                      MS147 = nrow(MS147[NEBNext == 0 & MS147 == 1, ]),
                      shared = nrow(MS147[NEBNext == 1 & MS147 == 1, ]))


dts = list(MS146_dt, MS147_dt)
invisible(lapply(dts, function(dt) {
  lib = names(dt)[2]
  venn_data = data.table(ref = dt$NEBNext + dt$shared,
                         replicate = dt[, 2] + dt$shared,
                         ref_replicate = dt$shared)
  setnames(venn_data, c("ref", "replicate", "ref_replicate"))
  
  setEPS()
  cairo_ps(paste0("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/mutations/venndiagrams/venndiagramsVennDiagram_cairo-ref-",lib, "_ct<=35_ps.eps"),
           onefile = TRUE, height=7, width=7, family="Helvetica",
           pointsize=8, antialias="none")
  draw.pairwise.venn(
    venn_data$ref,                            # Size of the left circle.
    venn_data$replicate,                            # Size of the right circle.
    venn_data$ref_replicate,                        # Size of the overlapping area.
    category = c("NEBNext", lib),  # Label text.
    cat.pos = c(0, 0),                   # Position of labels.
    scaled=TRUE,                         # Scale the circle size or not.
  )
  dev.off()
}))
