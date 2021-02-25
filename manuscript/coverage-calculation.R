## Author: Luuk Harbers
## Date: 2020-07-29
## Script for calculating coverage and intersecting with known mutations

## Load/install packages
packages = c("data.table", "karyoploteR", "GenomicAlignments", "Rsamtools", "IRanges")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

bam_1 = "/mnt/AchTeraD/data/BICRO231/COVseq/trimmed/MS45_S4.trimmed.sorted.bam"
bam_2 = "/mnt/AchTeraD/data/BICRO231/COVseq/trimmed/MS46_S5.trimmed.sorted.bam"


param = ScanBamParam(what=c("qname","flag"))
gr_1 = granges(readGAlignments(bam_1, param=param))
gr_2 = granges(readGAlignments(bam_2, param=param))

cov1 = coverage(gr_1)
cov2 = coverage(gr_2)

# Get sum of bases covered by >n reads
sum(cov1$NC_045512.2@lengths[cov1$NC_045512.2@values > 0]) / sum(cov1$NC_045512.2@lengths)
sum(cov2$NC_045512.2@lengths[cov2$NC_045512.2@values > 0]) / sum(cov2$NC_045512.2@lengths)

# Subset
subsetrange = GRanges(seqnames = "NC_045512.2", IRanges(start = 21563, end = 25384))
gr_1sub = subsetByOverlaps(gr_1, subsetrange)
gr_2sub = subsetByOverlaps(gr_2, subsetrange)

cov1_sub = coverage(gr_1sub)
cov2_sub = coverage(gr_2sub)

#sublength
sublength = 25384 - 21563

# Get sum of bases covered by >n reads
sum(cov1_sub$NC_045512.2@lengths[cov1_sub$NC_045512.2@values > 1]) / sublength
sum(cov2_sub$NC_045512.2@lengths[cov2_sub$NC_045512.2@values > 1]) / sublength
