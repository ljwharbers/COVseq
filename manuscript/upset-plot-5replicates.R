## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "UpSetR", "VennDiagram")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

muts = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/top20AAmutations-frontimicrob.tsv", header = T)

# Load files
annot = fread("/mnt/AchTeraD/data/BICRO251/TN21-annot_barcodes.txt")
annot[, V1 := paste0("Sample ", V1)]

dirs = list.dirs("/mnt/AchTeraD/data/BICRO254/", recursive = F, full.names = T)
dirs = dirs[!grepl("fastq", dirs)]
dirs = c(dirs, "/mnt/AchTeraD/data/BICRO251/TN21/")


res = lapply(dirs, function(dir){
  files = list.files(paste0(dir, "/called"), full.names = T)
  
  filenames = gsub(".trimmed.*", "", basename(files))
  
  res = lapply(1:length(files), function(x) {
    #dt = fread(files[x], select = c(1:3, 5:7), col.names = c("chr", "start", "end", "qual", "ref", "alt"))
    dt = fread(files[x])
    dt[, sample := filenames[x]]
    dt[, run := basename(dir)]
    dt[, depth := as.numeric(gsub(".*DP=|;.*", "", dt$V9))]
    
    # Get HQ alt/ref allele counts
    counts = dt[, gsub(".*DP4=|;.*", "", V9)]
    counts = tstrsplit(counts, ",")
    ref = as.numeric(counts[[1]]) + as.numeric(counts[[2]])
    alt = as.numeric(counts[[3]]) + as.numeric(counts[[4]])
    
    dt[, nalt := alt]
    dt[, nref := ref]
    dt = dt[, .(V1, V2, V3, V5, V6, V7, depth, sample, run, nalt, nref)]
    setnames(dt, c("chr", "start", "end", "qual", "ref", "alt", "depth", "sample", "run", "nalt", "nref"))
    return(dt)
  })
  total = rbindlist(res)
})
total = rbindlist(res)

# Merge with annotation file
total = merge(total, annot, by.x = "sample", by.y = "V2")
total[, sample := V1]
total$V1 = NULL

# discard =  c("Sample 22", "Sample 24", "Sample N1", "Sample N2", "Sample N3", "Sample P1", "Sample P2")
discard =  c("Sample 22", "Sample 24", "Sample N1", "Sample N2", "Sample N3", "Sample P1", "Sample P2",
             "Sample 30", "Sample 23", "Sample 48", "Sample 44", "Sample 32", "Sample 42", "Sample 37",
             "Sample 46", "Sample 35", "Sample 11", "Sample 28", "Sample 21")

total = total[!sample %in% discard]
# discard_mut = c("T16649A", "C16650A")

# Match actual starts/end to 1 indexed start site
total[, start := start + 1]
total[, end := end + 1]
total[, vaf := nalt / (nalt + nref)]

total[, id := factor(paste0(ref, start, alt))]

# Remove muts from sequencing artifacts
# total = total[!id %in% discard_mut, ]

# Filter out indels (only keep snps) and low quality snps
total = total[nchar(ref) == 1 & nchar(alt) == 1 & qual >= 30 & depth >= 10 & nalt > nref]
total[, present := 1]
total[, sample_id := paste0(sample, "_", id)]
total_mat = dcast(total, sample_id ~ run, fun.aggregate = function(x) sum(x), value.var = "present")

plt = upset(total_mat, sets = basename(dirs)[1:3], nintersect = 30, keep.order = T, order.by = "freq", point.size = 3,
            line.size = 1)

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/upsetR_55samples-replicates(MS)<=30",
              height = 8, width = 12)

total_mat[total_mat == ""] = NA
total_m = melt(total_mat, id.vars = "sample_id")
total_m = total_m[complete.cases(total_m)]
# 
# counts = total_m[, .N, by = .(id, value)]
# order = counts[, sum(N), by = .(id)]
# setorder(order, -V1)
# 
# counts[, id := factor(id, levels = order$id)]
# 
# 
# plt = ggplot(counts, aes(x = id, y = N, fill = value)) +
#   geom_bar(stat = "identity") +
#   scale_fill_viridis_d() +
#   labs(x = "Mutation", y = "Number of times mutation is called", fill = "Called by:") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/barplot-mutationscalled-posSamples",
#               width = 24, height = 12)
# 
# Get overlap for each combination
ms82 = total_mat[, .(sample_id, TN21, MS82_S3)]
ms83 = total_mat[, .(sample_id, TN21, MS83_S4)]
ms84 = total_mat[, .(sample_id, TN21, MS84_S5)]
nz190 = total_mat[, .(sample_id, TN21, NZ190_S1)]
nz191 = total_mat[, .(sample_id, TN21, NZ191_S2)]

ms82_dt = data.table(TN21 = nrow(ms82[TN21 == 1 & MS82_S3 == 0, ]),
                     MS82 = nrow(ms82[TN21 == 0 & MS82_S3 == 1, ]),
                     shared = nrow(ms82[TN21 == 1 & MS82_S3 == 1, ]))

ms83_dt = data.table(TN21 = nrow(ms83[TN21 == 1 & MS83_S4 == 0, ]),
                     MS83 = nrow(ms83[TN21 == 0 & MS83_S4 == 1, ]),
                     shared = nrow(ms83[TN21 == 1 & MS83_S4 == 1, ]))

ms84_dt = data.table(TN21 = nrow(ms84[TN21 == 1 & MS84_S5 == 0, ]),
                     MS84 = nrow(ms84[TN21 == 0 & MS84_S5 == 1, ]),
                     shared = nrow(ms84[TN21 == 1 & MS84_S5 == 1, ]))

nz190_dt = data.table(TN21 = nrow(nz190[TN21 == 1 & NZ190_S1 == 0, ]),
                      NZ190 = nrow(nz190[TN21 == 0 & NZ190_S1 == 1, ]),
                      shared = nrow(nz190[TN21 == 1 & NZ190_S1 == 1, ]))

nz191_dt = data.table(TN21 = nrow(nz191[TN21 == 1 & NZ191_S2 == 0, ]),
                      NZ191 = nrow(nz191[TN21 == 0 & NZ191_S2 == 1, ]),
                      shared = nrow(nz191[TN21 == 1 & NZ191_S2 == 1, ]))

dts = list(ms82_dt, ms83_dt, ms84_dt, nz190_dt, nz191_dt)
invisible(lapply(dts, function(dt) {
  lib = names(dt)[2]
  venn_data = data.table(ref = dt$TN21 + dt$shared,
                         replicate = dt[, 2] + dt$shared,
                         ref_replicate = dt$shared)
  setnames(venn_data, c("ref", "replicate", "ref_replicate"))

  setEPS()
  cairo_ps(paste0("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/VennDiagram_cairo-ref-",lib, "_ps.eps"),
           onefile = TRUE, height=7, width=7, family="Helvetica",
           pointsize=8, antialias="none")
  draw.pairwise.venn(
    venn_data$ref,                            # Size of the left circle.
    venn_data$replicate,                            # Size of the right circle.
    venn_data$ref_replicate,                        # Size of the overlapping area.
    category = c("Ref", lib),  # Label text.
    cat.pos = c(0, 0),                   # Position of labels.
    scaled=TRUE,                         # Scale the circle size or not.
  )
  dev.off()
}))



# Load VennDiagram library
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/VennDiagram_cairo-posSamples_ps.eps",
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