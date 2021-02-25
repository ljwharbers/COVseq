## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30

libraries = list.files("/mnt/AchTeraD/data/BICRO258/TN31/trimmed/", full.names = T)
libraries = libraries[!grepl("sorted|bai", libraries)]
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")

# params
param = ScanBamParam(what=c("qname","flag"))
s_region = GRanges(seqnames = "NC_045512.2", IRanges(start = 21563, end = 25384))
s_length = 25384 - 21562

# go through list of bam files and calculate coverage in the full genome and in the S region of SARS-CoV-2
res = pblapply(libraries, function(library) {
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
  
  return(data.table(sample = gsub(".*/|\\..*", "", library),
                    aligned_reads = length(gr),
                    full_1x = full_1,
                    full_5x = full_5,
                    full_10x = full_10,
                    Sregion_1x = s_1,
                    Sregion_5x = s_5,
                    Sregion_10x = s_10,
                    average_coverage = mean(cov)))
}, cl = num_threads)

total = rbindlist(res)
total[, sample := gsub(".bam", "", basename(libraries))]

total_m = melt(total[, c(1, 3:8)], id.vars = "sample")

total_m[, region := ifelse(grepl("Sregion", variable), "S-region", "SARS-CoV-2")]
total_m[grepl("1x", variable), depth := "1x"]
total_m[grepl("5x", variable), depth := "5x"]
total_m[grepl("10x", variable), depth := "10x"]
total_m[, depth := factor(depth, levels = c("1x", "5x", "10x"))]
total_m[, region := factor(region, levels = c("SARS-CoV-2", "S-region"))]
total_m[, sample := factor(sample, levels = c("50k", "10k", "2k", "400", "neg"))]

#total_m = total_m[grepl("COVseq", readlength)]

plt1 = ggplot(total_m[region == "SARS-CoV-2"], aes(x = depth, y = value, fill = sample)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_d(end = 0.5) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sequencing depth")

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/barplots-boc/TN31-BreadtOfCoverage-fullregion",
              height = 6, width = 8)

plt2 = ggplot(total_m[region == "S-region"], aes(x = depth, y = value, fill = sample)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_d(end = 0.5) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[5]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sequencing depth")

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/barplots-boc/TN31-BreadtOfCoverage-Sregion",
              height = 6, width = 8)

total[, sample := factor(sample, levels = c("50k", "10k", "2k", "400", "neg"))]
plt3 = ggplot(total, aes(x = sample, y = aligned_reads, fill = sample)) +
  geom_col() +
  scale_fill_viridis_d(end = 0.5) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(y = "Reads", x = "Sample")

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/barplots-reads/TN31-reads",
              height = 6, width = 8)
