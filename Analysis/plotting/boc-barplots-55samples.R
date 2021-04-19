## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating and plotting coverage of SARS-CoV-2 genome

## Load/install packages
packages = c("data.table", "GenomicAlignments", "Rsamtools", "IRanges", "pbapply", "ggplot2", "GGally")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

num_threads=30
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/55-patients_CT-clean-newnames.tsv", header = T)
theoret = fread("/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv")
#annot = fread("/mnt/AchTeraD/data/BICRO268/MS147+148-barcodes-annot.txt", header = F)

libraries = c("NZ242", "NZ243", "NZ244", "NZ245")
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

total_ct = merge(total, ct)

plt = ggplot(total_ct, aes(x = N, y = full_10x, color = library)) +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "D") +
  labs(y = "Breadth of Coverage (10x)", x = "Ct", color = "Replicate") +
  scale_x_reverse() +
  scale_y_continuous(limits = c(0, 1))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-boc/NZ242-NZ245_boc-ct_55samples",
              height = 7, width = 8)


# Melt
total_d = dcast(total[, .(full_10x, sample, library)], sample ~ library, value.var = "full_10x")
total_d = total_d[!grepl("Neg|Pos", sample)]


plt = ggpairs(total_d[, 2:ncol(total_d)]) +
  theme_base() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-boc/NZ242-NZ245_fullregion(10x)-55+ctrls",
              height = 11, width = 14)

# Cor
cors = melt(cor(total_d[, 2:5]))

plt3 = ggplot(cors, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_x_discrete("", position = "top") +
  scale_fill_viridis("Pearson's\ncorrelation", option="inferno", begin = .8, direction = -1) +
  labs(y = "", x = "") +
  theme(axis.line = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.text.x.top = element_text(angle = 45, hjust = 0))

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/corr-boc/Correlation-BoC-55samples_NZ242-NZ245",
              height = 7, width = 8)
# Merge with annotation
#total = merge(total, annot, by.x = "sample", by.y = "V1")

# Melt and prepare for plotting
total_m = melt(total[, c(9, 4, 7)], id.vars = "sample")
total_m[, depth := gsub(".*_", "", variable)]
total_m[, region := gsub("_.*", "", variable)]

# Order by ct
total_m[, sample := factor(sample, levels = c(ct$sample[order(ct$N)], 
                                              "Sample Neg1", "Sample Neg2", "Sample Neg3", 
                                              "Sample Pos1", "Sample Pos2"))]

plt1 = ggplot(total_m[region == "full"], aes(x = sample, y = value, group = depth, fill = value)) +
  geom_col() +
  scale_fill_viridis_c(direction = -1) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limit = c(0, 1)) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[2]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/barplots-boc/NZ225_fullregion-55+ctrls",
              height = 7, width = 14)

plt2 = ggplot(total_m[region == "Sregion"], aes(x = sample, y = value, group = depth, fill = value)) +
  geom_col() +
  scale_fill_viridis_c(direction = -1) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limit = c(0, 1)) +
  geom_hline(yintercept = c(theoret$`NlaIII + MseI`[5]), color = "red", linetype = "dashed", size = 1) +
  labs(y = "Breadth of Coverage", x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/barplots-boc/NZ225_sregion-55+ctrls",
              height = 7, width = 14)
