## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply")
sapply(packages, require, character.only = T)

num_threads=30
basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/"
libraries = c("NZ265", "NZ266", "NZ267")

# params
param = ScanBamParam(what=c("qname","flag"))
s_region = GRanges(seqnames = "NC_045512.2", IRanges(start = 21563, end = 25384))
s_length = 25384 - 21562

res = lapply(libraries, function(x) {
  files = list.files(paste0(basedir, x, "/variants/bam/"), 
                     pattern = ".bam$", full.names = T)
  
  samplenames = gsub(".trim.sorted.bam", "", basename(files))
  samplenames = gsub("_", " ", samplenames)
  
  res = pblapply(files, function(library) {
    # Load bamfile as granges
    gr = granges(readGAlignments(library, param=param))
    
    # Get coverage full genome
    cov = coverage(gr)
    full_1 = sum(cov$NC_045512.2@lengths[cov$NC_045512.2@values > 0]) / sum(cov$NC_045512.2@lengths)
    full_5 = sum(cov$NC_045512.2@lengths[cov$NC_045512.2@values > 4]) / sum(cov$NC_045512.2@lengths)
    full_10 = sum(cov$NC_045512.2@lengths[cov$NC_045512.2@values > 9]) / sum(cov$NC_045512.2@lengths)
    
    # Subset for S-region  
    gr_sregion = subsetByOverlaps(gr, s_region, type = "within")
    
    # Get coverage S-region
    cov_s = coverage(gr_sregion)
    
    s_1 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 0]) / s_length
    s_5 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 4]) / s_length
    s_10 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 9]) / s_length 
    
    # Return all BoCs and #reads
    return(data.table(aligned_reads = length(gr),
                      full_1x = full_1,
                      full_5x = full_5,
                      full_10x = full_10,
                      Sregion_1x = s_1,
                      Sregion_5x = s_5,
                      Sregion_10x = s_10,
                      average_coverage = mean(cov)))
  }, cl = num_threads)
  
  total = rbindlist(res)
  total[, sample := samplenames]
  total[, library := x]
})

total = rbindlist(res)