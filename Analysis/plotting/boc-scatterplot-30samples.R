
## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values-cleaned.tsv", header = T)
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")
#annot = fread("/mnt/AchTeraD/data/BICRO268/MS147+148-barcodes-annot.txt", header = F)

libraries = c("MS147", "NEBNext")

res = lapply(libraries, function(lib) {
  files = list.files(paste0("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/", lib, "/variants/bam/"), pattern = ".bam$", full.names = T)
  
  samplenames = gsub(".trim.sorted.bam", "", basename(files))
  samplenames = gsub("_", " ", samplenames)
  # params
  param = ScanBamParam(what=c("qname","flag"))
  s_region = GRanges(seqnames = "NC_045512.2", IRanges(start = 21563, end = 25384))
  s_length = 25384 - 21562
  
  # go through list of bam files and calculate coverage in the full genome and in the S region of SARS-CoV-2
  res = pblapply(files, function(library) {
    # Load bamfile as granges
    gr = granges(readGAlignments(library, param=param))
    
    # Get coverage full genome
    cov = coverage(gr)
    full_10 = sum(cov$NC_045512.2@lengths[cov$NC_045512.2@values > 9]) / sum(cov$NC_045512.2@lengths)
    
    # Subset for S-region  
    gr_sregion = subsetByOverlaps(gr, s_region, type = "within")
    
    # Get coverage S-region
    cov_s = coverage(gr_sregion)
    
    s_10 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 9]) / s_length 
    
    return(data.table(aligned_reads = length(gr),
                      full_10x = full_10,
                      Sregion_10x = s_10,
                      average_coverage = mean(cov)))
  }, cl = num_threads)
  
  total = rbindlist(res)
  total[, sample := samplenames]
  total[, library := lib]
  return(total)
})

total = rbindlist(res)

# Remove controls + sample 12
total = total[!grepl("Sample 12|Neg|Pos", sample)]

# Make dt for scatterplot
dt = data.table("NEBNext" = total[library == "NEBNext"]$full_10x,
                "COVseq" = total[library == "MS147"]$full_10x,
                sample = total[library == "NEBNext"]$sample)

# Merge with Ct
dt = merge(dt, ct)

plt = ggplot(dt, aes(x = NEBNext, y = COVseq, color = ct)) + 
  geom_point(size = 4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_viridis_c(option = "E", direction = -1) +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 2, size = 2, fullrange = TRUE) +
  stat_cor() +
  labs(y = "COVseq - Breadth of Coverage (10x)", x = "NEBNext - Breadth of Coverage (10x)")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-boc/MS147-NEBNext-BoC(10x)-scatter",
              height = 7, width = 7)
