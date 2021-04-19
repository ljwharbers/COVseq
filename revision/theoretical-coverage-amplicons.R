## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculation of theoretical coverage

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nla = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_NlaIII-cutsites.bed",
            col.names = c("chr", "start", "end"))
mse = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_MseI-cutsites.bed",
            col.names = c("chr",  "start", "end"))


# Prepare ref
ref = data.table(region = c("SARS-CoV-2", "S-region"), start = c(0, 21563), end = c(29903, 25384))

nla[, enzyme := "NlaIII"]
mse[, enzyme := "MseI"]

# Check for each amplicon fragment the coverage
total = rbind(nla, mse)

total[, start_cover_150 := start - 150]
total[, end_cover_150 := end + 150]
total[, start_cover_300 := start - 280]
total[, end_cover_300 := end + 280]

total[start_cover_150 < 1, start_cover_150 := 1]
total[start_cover_300 < 1, start_cover_300 := 1]

total[end_cover_150 > 29903, end_cover_150 := 29903]
total[end_cover_300 > 29903, end_cover_300 := 29903]

# order
setorder(total, start, end)

dt = total[grepl("NlaIII", enzyme), list(chr = chr, start = start_cover_150, stop = end_cover_150)]
nla_150 = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / 29903

dt = total[grepl("NlaIII", enzyme), list(chr = chr, start = start_cover_300, stop = end_cover_300)]
nla_300 = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / 29903

dt = total[grepl("MseI", enzyme), list(chr = chr, start = start_cover_150, stop = end_cover_150)]
mse_150 = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / 29903

dt = total[grepl("MseI", enzyme), list(chr = chr, start = start_cover_300, stop = end_cover_300)]
mse_300 = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / 29903

dt = total[grepl("MseI|NlaIII", enzyme), list(chr = chr, start = start_cover_150, stop = end_cover_150)]
nla_mse_150 = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / 29903

dt = total[grepl("MseI|NlaIII", enzyme), list(chr = chr, start = start_cover_300, stop = end_cover_300)]
nla_mse_300 = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / 29903

# Same for S region
dt = total[grepl("NlaIII", enzyme) & start >= ref[region == "S-region"]$start & end <= ref[region == "S-region"]$end, 
           list(chr = chr, start = start_cover_150, stop = end_cover_150)]
nla_150_s = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / (ref[region == "S-region"]$end - ref[region == "S-region"]$start)

dt = total[grepl("NlaIII", enzyme) & start >= ref[region == "S-region"]$start & end <= ref[region == "S-region"]$end, 
           list(chr = chr, start = start_cover_300, stop = end_cover_300)]
nla_300_s = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / (ref[region == "S-region"]$end - ref[region == "S-region"]$start)

dt = total[grepl("MseI", enzyme) & start >= ref[region == "S-region"]$start & end <= ref[region == "S-region"]$end, 
           list(chr = chr, start = start_cover_150, stop = end_cover_150)]
mse_150_s = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / (ref[region == "S-region"]$end - ref[region == "S-region"]$start)

dt = total[grepl("MseI", enzyme) & start >= ref[region == "S-region"]$start & end <= ref[region == "S-region"]$end, 
           list(chr = chr, start = start_cover_300, stop = end_cover_300)]
mse_300_s = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / (ref[region == "S-region"]$end - ref[region == "S-region"]$start)

dt = total[start >= ref[region == "S-region"]$start & end <= ref[region == "S-region"]$end, 
           list(chr = chr, start = start_cover_150, stop = end_cover_150)]
nla_mse_150_s = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / (ref[region == "S-region"]$end - ref[region == "S-region"]$start)

dt = total[start >= ref[region == "S-region"]$start & end <= ref[region == "S-region"]$end, 
           list(chr = chr, start = start_cover_300, stop = end_cover_300)]
nla_mse_300_s = sum(width(reduce(GRanges(dt$chr, IRanges(dt$start, dt$stop))))) / (ref[region == "S-region"]$end - ref[region == "S-region"]$start)


theoretical = data.table(readlength = c(150, 300, 150, 300), 
                         region = c("SARS-CoV-2", "SARS-CoV-2", "S", "S"),
                         NlaIII = c(nla_150, nla_300, nla_150_s, nla_300_s), 
                         MseI = c(mse_150, mse_300, mse_150_s, mse_300_s),
                         "NlaIII + MseI" = c(nla_mse_150, nla_mse_300, nla_mse_150_s, nla_mse_300_s))
write.table(theoretical, "/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage-updated.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
