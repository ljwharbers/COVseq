## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")
annot = fread("/mnt/AchTeraD/data/BICRO251/TN21-annot_barcodes.txt")
annot[, V1 := paste0("Sample ", V1)]

files = list.files("/mnt/AchTeraD/data/BICRO251/TN11/trimmed/", pattern = ".bam$", full.names = T)
files = files[grepl("sample", files)]

files = c(files, list.files("/mnt/AchTeraD/data/BICRO251/TN21/trimmed/", pattern = ".bam$", full.names = T))
files = files[!grepl("fixed", files)]

samplenames = gsub(".bam|.trimmed.sorted", "", basename(files))

names(files) = samplenames

# params
param = ScanBamParam(what=c("qname","flag"))
s_region = GRanges(seqnames = "NC_045512.2", IRanges(start = 21563, end = 25384))
s_length = 25384 - 21562

# go through list of bam files and calculate coverage in the full genome and in the S region of SARS-CoV-2
res = pblapply(files, function(library) {
  # Load bamfile as granges
  gr = granges(readGAlignments(library, param=param))
  
  # Subset for S-region  
  gr_sregion = subsetByOverlaps(gr, s_region, type = "within")
  
  # Get coverage S-region
  cov_s = coverage(gr_sregion)
  return(cov_s)
}, cl = num_threads)


dt = pblapply(1:length(res), function(x) {
  
  data.table(sample = names(res)[x], 
             `10x` = sum(res[[x]]$NC_045512.2@lengths[res[[x]]$NC_045512.2@values > 9]) / s_length,
             `100x` = sum(res[[x]]$NC_045512.2@lengths[res[[x]]$NC_045512.2@values > 99]) / s_length
  )
})

total = rbindlist(dt)
total = merge(total, annot, by.x = "sample", by.y = "V2", all.x = T)
total[!is.na(V1), run := "TN21"]
total[is.na(V1), run := "TN11"]
total[!is.na(V1), sample := V1]

total$V1 = NULL

# Remove negative samples
discard = c("Sample 22", "Sample 24", "Sample N1", "Sample N2", "Sample N3", "Sample P1", "Sample P2",
            "sample17", "sample21", "sample29", "sample3", "sample30", "sample1", "sample2", "sample27", "sample18")

total = total[!sample %in% discard, ]
