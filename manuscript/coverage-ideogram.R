## Author: Luuk Harbers
## Date: 2020-07-29
## Script for plotting genome coverage on custom genome ideograms

## Load/install packages
packages = c("data.table", "karyoploteR", "GenomicAlignments", "Rsamtools")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

bam_1 = "/mnt/AchTeraD/data/BICRO231/MS47_S10/trimmed/ACGACCGC.trimmed.sorted.dedup.bam"
bam_2 = "/mnt/AchTeraD/data/BICRO236/MS48_S11/bamfiles/TGATGCGC.trimmed.sorted.dedup.bam"

param = ScanBamParam(what=c("qname","flag"))
gr_1 = granges(readGAlignments(bam_1, param=param))
gr_2 = granges(readGAlignments(bam_2, param=param))

# SARS-COV-2 karyoploteR annot
covid_gr = GRanges(seqnames = "NC_045512.2", ranges = IRanges(start = 1, end = 29903))
annot = data.table(chr = "NC_045512.2",
                   start = c(1, 266, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558, 29675),
                   end = c(265, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674, 29903),
                   annot = c("5` UTR", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3` UTR"),
                   gieStain = c("gneg", "gpos50", "acen", "gpos75", "gneg", "gvar", "gpos25", "stalk", "gpos50", "gneg", "acen", "gpos50", "gpos75"))

# Get coverage max for y-axis
max_cov1 = max(coverage(gr_1))
max_cov2 = max(coverage(gr_2))
max_cov = max(max_cov1, max_cov2)

# Plot
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/coverage/MS47_MS48_zoom-S-regioncairo_ps.eps",
         onefile = TRUE, height=6, width=8, family="Helvetica",
         pointsize=8)
kp = plotKaryotype(genome = covid_gr, plot.type = 2, cytoband = annot, zoom = GRanges(seqnames = "NC_045512.2", ranges = IRanges(start = 21563, end = 25384)))
kp = kpAddBaseNumbers(kp, tick.dist = 5000, add.units = TRUE, cex = 0.75)

# Plot coverage bam1
kp = kpPlotBAMCoverage(kp, data=bam_1, col="#FFD700", data.panel = 1)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1)
kpAddLabels(kp, "MS47", label.margin = 0.08, srt = 90, cex = 1.75, data.panel = 1)

# Plot coverage bam2
kp = kpPlotBAMCoverage(kp, data=bam_2, col="#1E90FF", data.panel = 2)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 2)
kpAddLabels(kp, "MS48", label.margin = 0.08, srt = 90, cex = 1.75, data.panel = 2)

dev.off()

