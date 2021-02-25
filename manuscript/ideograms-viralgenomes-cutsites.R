## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "karyoploteR")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Make annotations
sarscov2 = GRanges(seqnames = "NC_045512.2", ranges = IRanges(start = 1, end = 29903))
annot = data.table(chr = "NC_045512.2",
                   start = c(1, 266, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558, 29675),
                   end = c(265, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674, 29903),
                   annot = c("5` UTR", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3` UTR"),
                   gieStain = c("gneg", "gpos50", "acen", "gpos75", "gneg", "gvar", "gpos25", "stalk", "gpos50", "gneg", "acen", "gpos50", "gpos75"))

inf_a = GRanges(seqnames = c("NC_026438.1", "NC_026435.1", 
                             "NC_026437.1", "NC_026433.1", 
                             "NC_026436.1", "NC_026434.1", 
                             "NC_026431.1", "NC_026432.1"),
                ranges = c(IRanges(start = 1, end = 2368),
                           IRanges(start = 1, end = 2313),
                           IRanges(start = 1, end = 2204),
                           IRanges(start = 1, end = 1882),
                           IRanges(start = 1, end = 1841),
                           IRanges(start = 1, end = 1557),
                           IRanges(start = 1, end = 1191),
                           IRanges(start = 1, end = 1096)))

inf_b = GRanges(seqnames = c("NC_002204.1", "NC_002205.1", 
                             "NC_002206.1", "NC_002207.1", 
                             "NC_002208.1", "NC_002209.1", 
                             "NC_002210.1", "NC_002211.1"),
                ranges = c(IRanges(start = 1, end = 2280),
                           IRanges(start = 1, end = 2274),
                           IRanges(start = 1, end = 2151),
                           IRanges(start = 1, end = 1701),
                           IRanges(start = 1, end = 1497),
                           IRanges(start = 1, end = 1410),
                           IRanges(start = 1, end = 982),
                           IRanges(start = 1, end = 863)))

dengue = GRanges(seqnames = "NC_001474.2", ranges = IRanges(start = 1, end = 10723))

dirs = list.dirs("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/all-enzymes/", full.names = F, recursive = F)

# Load in all cutsite information
cutsites = lapply(dirs, function(dir) {
  files = list.files(paste0("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/all-enzymes/", dir), full.names = T)
  dt_list = lapply(files, function(file){
    dt = fread(file)
    dt[, V4 := gsub(".*_|-cut.*", "", basename(file))]
    setnames(dt, c("chr", "start", "end", "enzyme"))
    return(dt)
  })
  return(rbindlist(dt_list))
})
names(cutsites) = dirs

# Set params
pp = getDefaultPlotParams(plot.type=2)
pp$data1inmargin = 0
pp$data2inmargin = 0


# Plot ideograms
# Plot sars-cov-2
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/cutsite-distribution/sarscov2-all_4cutters-cutsite_ideogram_cairo_ps.eps",
         onefile = TRUE, height=6, width=16, family="Helvetica",
         pointsize=8)
kp = plotKaryotype(genome = sarscov2, plot.type = 2, cytoband = annot, plot.params = pp)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$`sars-cov-2`[enzyme == "BfaI"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.1, r1 = 0.3)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$`sars-cov-2`[enzyme == "NlaIII"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.3, r1 = 0.5)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$`sars-cov-2`[enzyme == "MseI"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.5, r1 = 0.7)

kpAddLabels(kp, labels = "BfaI - Cutsites", data.panel = 2, r0 = 0.1, r1 = 0.25)
kpAddLabels(kp, labels = "CviAII / FatI / NlaIII - Cutsites", data.panel = 2, r0 = 0.3, r1 = 0.45)
kpAddLabels(kp, labels = "MseI - Cutsites", data.panel = 2, r0 = 0.5, r1 = 0.65)

kp = kpAddBaseNumbers(kp, tick.dist = 1000, add.units = TRUE, cex = 0.75)

dev.off()

# Plot influenza a
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/cutsite-distribution/influenza_a-all_4cutters-cutsite_ideogram_cairo_ps.eps",
         onefile = TRUE, height=10, width=16, family="Helvetica",
         pointsize=8)
kp = plotKaryotype(genome = inf_a, plot.type = 2, plot.params = pp)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$influenza_a[enzyme == "BfaI"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.1, r1 = 0.3)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$influenza_a[enzyme == "NlaIII"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.4, r1 = 0.6)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$influenza_a[enzyme == "MseI"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.7, r1 = 0.9)

kpAddLabels(kp, labels = "BfaI - Cutsites", data.panel = 2, r0 = 0.1, r1 = 0.25)
kpAddLabels(kp, labels = "CviAII / FatI / NlaIII - Cutsites", data.panel = 2, r0 = 0.4, r1 = 0.55)
kpAddLabels(kp, labels = "MseI - Cutsites", data.panel = 2, r0 = 0.5, r1 = 0.85)

kp = kpAddBaseNumbers(kp, tick.dist = 1000, add.units = TRUE, cex = 0.75)

dev.off()

# Plot influenza b
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/cutsite-distribution/influenza_b-all_4cutters-cutsite_ideogram_cairo_ps.eps",
         onefile = TRUE, height=10, width=16, family="Helvetica",
         pointsize=8)
kp = plotKaryotype(genome = inf_b, plot.type = 2, plot.params = pp)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$influenza_b[enzyme == "BfaI"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.1, r1 = 0.3)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$influenza_b[enzyme == "NlaIII"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.4, r1 = 0.6)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$influenza_b[enzyme == "MseI"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.7, r1 = 0.9)

kpAddLabels(kp, labels = "BfaI - Cutsites", data.panel = 2, r0 = 0.1, r1 = 0.25)
kpAddLabels(kp, labels = "CviAII / FatI / NlaIII - Cutsites", data.panel = 2, r0 = 0.4, r1 = 0.55)
kpAddLabels(kp, labels = "MseI - Cutsites", data.panel = 2, r0 = 0.7, r1 = 0.85)

kp = kpAddBaseNumbers(kp, tick.dist = 1000, add.units = TRUE, cex = 0.75)

dev.off()

# Plot dengue
setEPS()
cairo_ps("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/cutsite-distribution/dengue-all_4cutters-cutsite_ideogram_cairo_ps.eps",
         onefile = TRUE, height=6, width=16, family="Helvetica",
         pointsize=8)
kp = plotKaryotype(genome = dengue, plot.type = 2, plot.params = pp)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$dengue[enzyme == "BfaI"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.1, r1 = 0.3)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$dengue[enzyme == "NlaIII"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.3, r1 = 0.5)
kpPlotMarkers(kp, data=makeGRangesFromDataFrame(cutsites$dengue[enzyme == "MseI"]), 
              labels = "", marker.parts = c(1, 0, 0), data.panel = 2, r0 = 0.5, r1 = 0.7)

kpAddLabels(kp, labels = "BfaI - Cutsites", data.panel = 2, r0 = 0.1, r1 = 0.25)
kpAddLabels(kp, labels = "CviAII / FatI / NlaIII - Cutsites", data.panel = 2, r0 = 0.3, r1 = 0.45)
kpAddLabels(kp, labels = "MseI - Cutsites", data.panel = 2, r0 = 0.5, r1 = 0.65)

kp = kpAddBaseNumbers(kp, tick.dist = 1000, add.units = TRUE, cex = 0.75)

dev.off()