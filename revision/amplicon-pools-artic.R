# custom COV-ID19 genome creation and primer allocation

require(karyoploteR)
require(data.table)
require(GenomicRanges)
require(dplyr)
require(RColorBrewer)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

covid_gr = GRanges(seqnames = "SARS-CoV-2", ranges = IRanges(start = 1, end = 29903))

annot = data.table(chr = rep("SARS-CoV-2", 12),
                   start = c(1, 266, 21563, 25393, 26245, 26523, 27202, 27394, 27894, 28274, 29558, 29675),
                   end = c(265, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 28259, 29533, 29674, 29903),
                   annot = c("5` UTR", "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N", "ORF10", "3` UTR"),
                   gieStain = c("gneg", "gpos50", "acen", "gpos75", "gneg", "gvar", "gpos25", "stalk", "gneg", "acen", "gpos50", "gpos75"))


#load primers
primers = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cov2-primers/CDC-SARS-CoV-2-Primers.csv")
pool1 = unlist(lapply(seq(1, nrow(primers), by = 2), function(x) paste0(primers$Pool1[x], "-", primers$Pool1[x+1])))
pool2 = unlist(lapply(seq(1, nrow(primers), by = 2), function(x) paste0(primers$Pool2[x], "-", primers$Pool2[x+1])))
pool3 = unlist(lapply(seq(1, nrow(primers), by = 2), function(x) paste0(primers$Pool3[x], "-", primers$Pool3[x+1])))
pool4 = unlist(lapply(seq(1, nrow(primers), by = 2), function(x) paste0(primers$Pool4[x], "-", primers$Pool4[x+1])))
pool5 = unlist(lapply(seq(1, nrow(primers), by = 2), function(x) paste0(primers$Pool5[x], "-", primers$Pool5[x+1])))
pool6 = unlist(lapply(seq(1, nrow(primers), by = 2), function(x) paste0(primers$Pool6[x], "-", primers$Pool6[x+1])))
pools = as.data.table(cbind(pool1, pool2, pool3, pool4, pool5, pool6))
pools[pools == "NA-NA"] <- NA
pool1 = GRanges(seqnames = "SARS-CoV-2", IRanges(pools$pool1))
pool2 = GRanges(seqnames = "SARS-CoV-2", IRanges(pools$pool2[!is.na(pools$pool2)]))
pool3 = GRanges(seqnames = "SARS-CoV-2", IRanges(pools$pool3[!is.na(pools$pool3)]))
pool4 = GRanges(seqnames = "SARS-CoV-2", IRanges(pools$pool4[!is.na(pools$pool4)]))
pool5 = GRanges(seqnames = "SARS-CoV-2", IRanges(pools$pool5[!is.na(pools$pool5)]))
pool6 = GRanges(seqnames = "SARS-CoV-2", IRanges(pools$pool6[!is.na(pools$pool6)]))

# Load artic primers
artic = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cov2-primers/artic_v3-inserts.bed")
artic_1 = GRanges(seqnames = "SARS-CoV-2", IRanges(start = artic[V5 == 1]$V2,
                                                   end = artic[V5 == 1]$V3))
artic_2 = GRanges(seqnames = "SARS-CoV-2", IRanges(start = artic[V5 == 2]$V2,
                                                   end = artic[V5 == 2]$V3))
# plot genome
cairo_ps(file = paste0("/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/primer-pools/artic+cdc-pools.eps"), height = 6, width = 12)
kp = plotKaryotype(genome = covid_gr, cytoband = annot)
kpPlotRegions(kp, data=pool1, col = brewer.pal(6, "Set1")[1], layer.margin = 0.01, border=NA, r0=0, r1=0.1)
kpPlotRegions(kp, data=pool2, col = brewer.pal(6, "Set1")[2], layer.margin = 0.01, border=NA, r0=0.15, r1=0.25)
kpPlotRegions(kp, data=pool3, col = brewer.pal(6, "Set1")[3], layer.margin = 0.01, border=NA, r0=0.30, r1=0.40)
kpPlotRegions(kp, data=pool4, col = brewer.pal(6, "Set1")[4], layer.margin = 0.01, border=NA, r0=0.45, r1=0.55)
kpPlotRegions(kp, data=pool5, col = brewer.pal(6, "Set1")[5], layer.margin = 0.01, border=NA, r0=0.6, r1=0.7)
kpPlotRegions(kp, data=pool6, col = brewer.pal(7, "Set1")[7], layer.margin = 0.01, border=NA, r0=0.75, r1=0.85)
kpPlotRegions(kp, data=artic_1, col = brewer.pal(9, "Set1")[8], layer.margin = 0.01, border=NA, r0=0.9, r1=1.05)
kpPlotRegions(kp, data=artic_2, col = brewer.pal(9, "Set1")[9], layer.margin = 0.01, border=NA, r0=1.1, r1=1.25)
kpAddLabels(kp, labels="CDC 1", r0=0, r1=0.1)
kpAddLabels(kp, labels="CDC 2", r0=0.15, r1=0.25)
kpAddLabels(kp, labels="CDC 3", r0=0.30, r1=0.40)
kpAddLabels(kp, labels="CDC 4", r0=0.45, r1=0.55)
kpAddLabels(kp, labels="CDC 5", r0=0.6, r1=0.7)
kpAddLabels(kp, labels="CDC 6", r0=0.75, r1=0.85)
kpAddLabels(kp, labels="Artic 1", r0=0.9, r1=1.05)
kpAddLabels(kp, labels="Artic 2", r0=1.1, r1=1.25)
dev.off()

png(file = paste0("/mnt/AchTeraD/Documents/Projects/COVID19/Plots/primer-pools.png"), height = 6, width = 12, units = "in", res = 300)
kp = plotKaryotype(genome = covid_gr, cytoband = annot)
kpPlotRegions(kp, data=pool1, col = brewer.pal(6, "Set1")[1], layer.margin = 0.01, border=NA, r0=0, r1=0.1)
kpPlotRegions(kp, data=pool2, col = brewer.pal(6, "Set1")[2], layer.margin = 0.01, border=NA, r0=0.15, r1=0.25)
kpPlotRegions(kp, data=pool3, col = brewer.pal(6, "Set1")[3], layer.margin = 0.01, border=NA, r0=0.30, r1=0.40)
kpPlotRegions(kp, data=pool4, col = brewer.pal(6, "Set1")[4], layer.margin = 0.01, border=NA, r0=0.45, r1=0.55)
kpPlotRegions(kp, data=pool5, col = brewer.pal(6, "Set1")[5], layer.margin = 0.01, border=NA, r0=0.6, r1=0.7)
kpPlotRegions(kp, data=pool6, col = brewer.pal(7, "Set1")[7], layer.margin = 0.01, border=NA, r0=0.75, r1=0.85)
kpAddLabels(kp, labels="Pool 1", r0=0, r1=0.1)
kpAddLabels(kp, labels="Pool 2", r0=0.15, r1=0.25)
kpAddLabels(kp, labels="Pool 3", r0=0.30, r1=0.40)
kpAddLabels(kp, labels="Pool 4", r0=0.45, r1=0.55)
kpAddLabels(kp, labels="Pool 5", r0=0.6, r1=0.7)
kpAddLabels(kp, labels="Pool 6", r0=0.75, r1=0.85)
dev.off()


# Size distribution
pool1_size = data.table(pool = "pool1", left = gsub("-.*", "", pools$pool1), right = gsub(".*-", "", pools$pool1))
pool2_size = data.table(pool = "pool2", left = gsub("-.*", "", pools$pool2), right = gsub(".*-", "", pools$pool2))
pool3_size = data.table(pool = "pool3", left = gsub("-.*", "", pools$pool3), right = gsub(".*-", "", pools$pool3))
pool4_size = data.table(pool = "pool4", left = gsub("-.*", "", pools$pool4), right = gsub(".*-", "", pools$pool4))
pool5_size = data.table(pool = "pool5", left = gsub("-.*", "", pools$pool5), right = gsub(".*-", "", pools$pool5))
pool6_size = data.table(pool = "pool6", left = gsub("-.*", "", pools$pool6), right = gsub(".*-", "", pools$pool6))

pools_size = rbind(pool1_size, pool2_size, pool3_size, pool4_size, pool5_size, pool6_size)
pools_size[, left := as.numeric(left)]
pools_size[, right := as.numeric(right)]
pools_size[, length := right - left]

plt = ggplot(pools_size, aes(x = length, color = pool)) + 
  geom_density(size = .8, adjust = 1) +
  scale_color_brewer(palette = "Set1")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVID19/Plots/length-distribution", height = 7, width = 7)
