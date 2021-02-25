## Author: Luuk Harbers
## Date: 2020-07-29
## Script for plotting genome coverage on custom genome ideograms

## Load/install packages
packages = c("data.table", "karyoploteR", "GenomicAlignments", "Rsamtools")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

files = list.files("/mnt/AchTeraD/data/BICRO241/covseq/TN7_S2/trimmed/", pattern = "bam$", full.names = T)
plot_dir = "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/coverage/TN7/"


# SARS-COV-2 karyoploteR annot
covid_gr = GRanges(seqnames = "NC_045512.2", ranges = IRanges(start = 1, end = 29903))
annot = data.table(chr = "NC_045512.2",
                   start = c(1, 266, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558, 29675),
                   end = c(265, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674, 29903),
                   annot = c("5` UTR", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3` UTR"),
                   gieStain = c("gneg", "gpos50", "acen", "gpos75", "gneg", "gvar", "gpos25", "stalk", "gpos50", "gneg", "acen", "gpos50", "gpos75"))

for(bam in files){
  sample = gsub(".*/|\\..*", "", bam)
  param = ScanBamParam(what=c("qname","flag"))
  gr = granges(readGAlignments(bam, param=param))
  
  # Plot
  dir.create(plot_dir, showWarnings = F)
  setEPS()
  cairo_ps(paste0(plot_dir, sample, "_coverage_cairo.eps"),
           onefile = TRUE, height=6, width=8, family="Helvetica",
           pointsize=8)
  kp = plotKaryotype(genome = covid_gr, plot.type = 2, cytoband = annot)
  kp = kpAddBaseNumbers(kp, tick.dist = 5000, add.units = TRUE, cex = 0.75)
  
  # Plot coverage bam1
  kp = kpPlotBAMCoverage(kp, data=bam, col="#FFD700", data.panel = 1)
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage, data.panel = 1)
  kpAddLabels(kp, sample, label.margin = 0.08, srt = 90, cex = 1.75, data.panel = 1)
  
  dev.off()
}



