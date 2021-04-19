## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/96-ct.tsv", header = T)
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")
#annot = fread("/mnt/AchTeraD/data/BICRO268/MS147+148-barcodes-annot.txt", header = F)

libraries = c("NZ265", "NZ266", "NZ267")
# params
param = ScanBamParam(what=c("qname","flag"))
s_region = GRanges(seqnames = "NC_045512.2", IRanges(start = 21563, end = 25384))
s_length = 25384 - 21562

res = lapply(libraries, function(x) {
  files = list.files(paste0("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/", x, "/variants/bam/"), 
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
total_ct = merge(total, ct, by.x = "sample", by.y = "id")

plt = ggplot(total_ct, aes(x = ct, y = full_10x, color = library)) +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "D") +
  labs(y = "Breadth of Coverage (10x)", x = "Ct", color = "Replicate") +
  scale_x_reverse() +
  scale_y_continuous(limits = c(0, 1))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-boc/NZ265-NZ267-96samplestriplicate",
              height = 7, width = 8)

# Plot corr heatmap
total_boc = total[, .(full_10x, sample, library)]
total_boc = dcast(total_boc, sample ~ library, value.var = "full_10x")
cors = cor(total_boc[, 2:4], method = "pearson")
cors = melt(cors)

plt2 = ggplot(cors, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_x_discrete("", position = "top") +
  scale_fill_viridis("Pearson's\ncorrelation", option="inferno", begin = .8, direction = -1) +
  labs(y = "", x = "") +
  theme(axis.line = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.text.x.top = element_text(angle = 45, hjust = 0))
  
save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/corr-boc/Correlation-BoC-96samples_NZ265-NZ267",
              height = 7, width = 8)

# Plot separate scatterplots
total_mat = dcast(total, sample~library, value.var = "full_10x")
total_mat = total_mat[sample != "Sample Neg"]

plt3 = ggplot(total_mat, aes(x = NZ265, y = NZ266)) +
  geom_point(size = 3) +
  labs(y = "Breadth of Coverage (10x) - NZ265", x = "Breadth of Coverage (10x) - NZ266") +
  geom_smooth(method = "lm", linetype = 2, color = "red", se = F) +
  stat_cor()

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-boc/NZ265+NZ266-96samplestriplicate",
              height = 7, width = 8)

plt4 = ggplot(total_mat, aes(x = NZ265, y = NZ267)) +
  geom_point(size = 3) +
  labs(y = "Breadth of Coverage (10x) - NZ265", x = "Breadth of Coverage (10x) - NZ267") +
  geom_smooth(method = "lm", linetype = 2, color = "red", se = F) +
  stat_cor()

save_and_plot(plt4, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-boc/NZ265+NZ267-96samplestriplicate",
              height = 7, width = 8)

plt5 = ggplot(total_mat, aes(x = NZ266, y = NZ267)) +
  geom_point(size = 3) +
  labs(y = "Breadth of Coverage (10x) - NZ266", x = "Breadth of Coverage (10x) - NZ267") +
  geom_smooth(method = "lm", linetype = 2, color = "red", se = F) +
  stat_cor()

save_and_plot(plt5, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-boc/NZ266+NZ267-96samplestriplicate",
              height = 7, width = 8)