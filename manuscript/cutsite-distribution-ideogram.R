## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for creating plot showing cutsite locations of NlaIII and MseI in the SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "karyoploteR", "GenomicRanges")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nla = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_NlaIII-cutsites.bed",
            col.names = c("chr", "start", "end"))
mse = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_MseI-cutsites.bed",
            col.names = c("chr", "start", "end"))

# Transform to granges
nla[, chr := "SARS-CoV-2"]
mse[, chr := "SARS-CoV-2"]
nla = makeGRangesFromDataFrame(nla)
mse = makeGRangesFromDataFrame(mse)

# Create SARS-CoV-2 genome plot
# SARS-COV-2 karyoploteR annot
covid_gr = GRanges(seqnames = "SARS-CoV-2", ranges = IRanges(start = 1, end = 29903))
annot = data.table(chr = "SARS-CoV-2",
                   start = c(1, 266, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558, 29675),
                   end = c(265, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674, 29903),
                   annot = c("5` UTR", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3` UTR"),
                   gieStain = c("gneg", "gpos50", "acen", "gpos75", "gneg", "gvar", "gpos25", "stalk", "gpos50", "gneg", "acen", "gpos50", "gpos75"))

# Get uncovered regions (300bp)
covered = c(nla, mse)
covered = as.data.table(covered)
covered[, start := start - 281]
covered[, end := end + 281]
covered[, width := end-start]

# Fix before and after genome
covered[start < 1, start := 1]
covered[end > 29903, end := 29903]

# Transform back
covered = reduce(makeGRangesFromDataFrame(covered))

noncovered = as.data.table(covered)
noncovered = noncovered[, .(seqnames = seqnames, start = data.table::shift(end), end = start)]
noncovered = makeGRangesFromDataFrame(noncovered[complete.cases(noncovered)])

# Plot ideogram
pp = getDefaultPlotParams(plot.type=2)
pp$data1inmargin = 0
pp$data2inmargin = 0

# Plot
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/cutsite-distribution/NlaIII-MseI-cutsite_ideogram_cairo_ps.eps",
         onefile = TRUE, height=6, width=8, family="Helvetica",
         pointsize=8)
kp = plotKaryotype(genome = covid_gr, plot.type = 2, cytoband = annot, plot.params = pp)
kpPlotMarkers(kp, data=nla, labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0, r1 = 0.2)
kpPlotMarkers(kp, data=mse, labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.2, r1 = 0.4)
kpAddLabels(kp, labels = "NlaIII - Cutsites", data.panel = 2, r0 = 0, r1 = 0.15)
kpAddLabels(kp, labels = "MseI - Cutsites", data.panel = 2, r0 = 0.2, r1 = 0.35)

dev.off()


# Plot 'dark' regions
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/uncovered_ideogram/300bp-uncovered_ideogram_cairo_ps.eps",
         onefile = TRUE, height=6, width=8, family="Helvetica",
         pointsize=8)
kp = plotKaryotype(genome = covid_gr, plot.type = 2, cytoband = annot, plot.params = pp)
kpPlotRegions(kp, data=noncovered, data.panel = 2, r0 = 0.1, r1 = 0.2)
kpAddLabels(kp, labels = "Not covered (300bp)", data.panel = 2, r0 = 0, r1 = 0.15)

dev.off()
