## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=15

libraries = list.files("/mnt/AchTeraD/data/BICRO241/covseq/", pattern = ".bam$", recursive = T, full.names = T)
libraries = libraries[grepl("MseI|nlaIII|merged", libraries)]
#libraries = c(libraries, "/mnt/AchTeraD/data/BICRO237/covseq/MS52_S4/bamfiles/MS52_S4.trimmed.sorted.bam")

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
  full_10 = sum(cov$NC_045512.2@lengths[cov$NC_045512.2@values > 9]) / sum(cov$NC_045512.2@lengths)
  
  # Subset for S-region  
  gr_sregion = subsetByOverlaps(gr, s_region, type = "within")
  
  # Get coverage S-region
  cov_s = coverage(gr_sregion)
  
  s_1 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 0]) / s_length
  s_10 = sum(cov_s$NC_045512.2@lengths[cov_s$NC_045512.2@values > 9]) / s_length 
  
  return(data.table("Digestion" = gsub(".*/|\\..*", "", library),
                    "Aligned reads" = length(gr),
                    "SARS-CoV-2" = full_1,
                    "S-region" = s_1,
                    "Average coverage" = mean(cov)))
}, cl = num_threads)

total = rbindlist(res, use.names = T)
total = total[c(1, 4, 5, 6), ]
total[, Digestion := c("MseI", "NlaIII", "NlaIII + MseI", "NEBNext")]

total_m = data.table::melt(total, id.vars = "Digestion")
total_m[, Digestion := factor(Digestion, levels = c("NlaIII", "MseI", "NlaIII + MseI", "NEBNext"))]

plt1 = ggplot(total_m[variable == "Aligned reads",], aes(x = variable, y = value, fill = Digestion)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_viridis_d() +
  labs(y = "Reads") +
  theme(axis.title.x = element_blank())

plt2 = ggplot(total_m[variable == "Average coverage"], aes(x = variable, y = value, fill = Digestion)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_viridis_d() +
  labs(y = "Average coverage") +
  theme(axis.title.x = element_blank())

plt3 = ggplot(total_m[!variable %in% c("Aligned reads", "Average coverage"), ], aes(x = variable, y = value, fill = Digestion)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_viridis_d() +
  labs(y = "Breadth of coverage", x = "Genomic region", fill = "Library type")

grid.arrange(plt1, plt2, plt3)
save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/Coverage-digestion_vs_nebnext",
              height = 5, width = 7)
