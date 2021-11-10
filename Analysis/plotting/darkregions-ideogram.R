## Date: xxxx-xx-xx
## Script for calculation of theoretical coverage

## Load/install packages
packages = c("data.table", "karyoploteR", "GenomicAlignments")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nla = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_NlaIII-cutsites.bed",
            col.names = c("chr", "start", "end"))
mse = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_MseI-cutsites.bed",
            col.names = c("chr",  "start", "end"))
muts = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/mutation-sourceinfo-updated3march.csv")


# Prepare ref
ref = data.table(region = c("SARS-CoV-2", "S-region"), start = c(0, 21563), end = c(29903, 25384))

nla[, enzyme := "NlaIII"]
mse[, enzyme := "MseI"]

# Check for each amplicon fragment the coverage
total = rbind(nla, mse)

total[, start_cover_300 := start - 300]
total[, end_cover_300 := end + 300]
total[start_cover_300 < 1, start_cover_300 := 1]
total[end_cover_300 > 29903, end_cover_300 := 29903]

# Reduce
ranges = reduce(GRanges(total$chr, IRanges(total$start_cover_300, total$end_cover_300)))

# Get non covered region
ranges = as.data.table(ranges)
ranges[, data.table::shift(start, type = "lead") - end]

uncovered = data.table(chr = ranges$seqnames, start = ranges$end,
                       end = ranges[, end + data.table::shift(start, type = "lead") - end])
uncovered = uncovered[1:nrow(uncovered)-1]

# SARS-COV-2 karyoploteR annot
covid_gr = GRanges(seqnames = "NC_045512.2", ranges = IRanges(start = 1, end = 29903))
annot = data.table(chr = "NC_045512.2",
                   start = c(1, 266, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558, 29675),
                   end = c(265, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674, 29903),
                   annot = c("5` UTR", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3` UTR"),
                   gieStain = c("gneg", "gpos50", "acen", "gpos75", "gneg", "gvar", "gpos25", "stalk", "gpos50", "gneg", "acen", "gpos50", "gpos75"))

# Plot ideogram
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/ideograms/darkregions_ideogram_cairo_ps.eps",
         onefile = TRUE, height=4, width=12, family="Helvetica",
         pointsize=8)
kp = plotKaryotype(genome = covid_gr, plot.type = 2, cytoband = annot, labels.plotter = NULL)
kp = kpAddChromosomeNames(kp, c("SARS-CoV-2"))
kp = kpPlotRegions(kp, data = makeGRangesFromDataFrame(uncovered), data.panel = 2, r0 = 0, r1 = 0.15)
kp = kpAddLabels(kp, labels = "Dark regions", data.panel = 2, r0 = 0, r1 = 0.15)
kp = kpPlotMarkers(kp, data = makeGRangesFromDataFrame(muts[, .(chr, start, end)]), data.panel = 2, r0 = 0.20, r1 = 0.35,
                   labels = "", marker.parts = c(1, 0, 0))
kp = kpAddLabels(kp, labels = "Variants of Concern", data.panel = 2, r0 = 0.2, r1 = 0.35)
dev.off()
