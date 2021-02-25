## Author: Luuk Harbers
## Date: 2020-02-16
## Script for analysis of variants

## Load/install packages
packages = c("data.table", "UpSetR", "VennDiagram")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
libraries = libraries[!grepl("nextflow|work|MS146|MS147|NZ", libraries)]

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
variants = variants[!grepl("Sample_Pos|Sample_Neg|Sample_10|Sample_21|Sample_1|Sample_30|Sample_3|Sample_17|Sample_29|Sample_2|Sample_27|Sample_18", sample),]

# Get id and dcast
variants[, id := factor(paste0(REF, POS, ALT))]
variants[, tosum := 1]
variants[, sample_id := paste0(sample, "_", id)]
variants_mat = dcast(variants, sample_id ~ library, fun.aggregate = function(x) sum(x), value.var = "tosum")

# check exclusive variants
nebonly = variants_mat[NEBNext == 1 & MS126 == 0 & MS127 == 0 & MS128 == 0 & MS129 == 0 & MS130 == 0 & MS131 == 0,]
replicates_total = variants_mat[MS126 == 1 | MS127 == 1 | MS128 == 1 | MS129 == 1 | MS130 == 1 | MS131 == 1,]

plt1 = upset(variants_mat, sets = libraries[1:6], nintersect = 10, keep.order = T, order.by = "freq", point.size = 3,
            line.size = 1)

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/upsetPlots/MS126-MS131-replicates-overlap_ct<=35",
              height = 5, width = 8)

plt2 = upset(variants_mat, sets = libraries[1:7], nintersect = 10, keep.order = T, order.by = "freq", point.size = 3,
            line.size = 1)

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/upsetPlots/MS126-MS131-replicates-NEBNext-overlap_ct<=35",
              height = 5, width = 8)

# Make Venn diagrams
MS126 = variants_mat[, .(sample_id, NEBNext, MS126)]
MS127 = variants_mat[, .(sample_id, NEBNext, MS127)]
MS128 = variants_mat[, .(sample_id, NEBNext, MS128)]
MS129 = variants_mat[, .(sample_id, NEBNext, MS129)]
MS130 = variants_mat[, .(sample_id, NEBNext, MS130)]
MS131 = variants_mat[, .(sample_id, NEBNext, MS131)]

MS126_dt = data.table(NEBNext = nrow(MS126[NEBNext == 1 & MS126 == 0, ]),
                     MS126 = nrow(MS126[NEBNext == 0 & MS126 == 1, ]),
                     shared = nrow(MS126[NEBNext == 1 & MS126 == 1, ]))

MS127_dt = data.table(NEBNext = nrow(MS127[NEBNext == 1 & MS127 == 0, ]),
                      MS127 = nrow(MS127[NEBNext == 0 & MS127 == 1, ]),
                     shared = nrow(MS127[NEBNext == 1 & MS127 == 1, ]))

MS128_dt = data.table(NEBNext = nrow(MS128[NEBNext == 1 & MS128 == 0, ]),
                      MS128 = nrow(MS128[NEBNext == 0 & MS128 == 1, ]),
                     shared = nrow(MS128[NEBNext == 1 & MS128 == 1, ]))

MS129_dt = data.table(NEBNext = nrow(MS129[NEBNext == 1 & MS129 == 0, ]),
                      MS129 = nrow(MS129[NEBNext == 0 & MS129 == 1, ]),
                      shared = nrow(MS129[NEBNext == 1 & MS129 == 1, ]))

MS130_dt = data.table(NEBNext = nrow(MS130[NEBNext == 1 & MS130 == 0, ]),
                      MS130 = nrow(MS130[NEBNext == 0 & MS130 == 1, ]),
                      shared = nrow(MS130[NEBNext == 1 & MS130 == 1, ]))

MS131_dt = data.table(NEBNext = nrow(MS131[NEBNext == 1 & MS131 == 0, ]),
                      MS131 = nrow(MS131[NEBNext == 0 & MS131 == 1, ]),
                      shared = nrow(MS131[NEBNext == 1 & MS131 == 1, ]))

dts = list(MS126_dt, MS127_dt, MS128_dt, MS129_dt, MS130_dt, MS131_dt)
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
