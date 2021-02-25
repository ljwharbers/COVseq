## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "UpSetR")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

muts = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/top20AAmutations-frontimicrob.tsv", header = T)

# Load files
annot = fread("/mnt/AchTeraD/data/BICRO251/TN21-annot_barcodes.txt")
annot[, V1 := paste0("Sample ", V1)]

dirs = list.dirs("/mnt/AchTeraD/data/BICRO254/", recursive = F, full.names = T)
dirs = dirs[!grepl("fastq", dirs)]
dirs = c(dirs, "/mnt/AchTeraD/data/BICRO251/TN21/")


res = lapply(dirs, function(dir){
  files = list.files(paste0(dir, "/called"), full.names = T)
  
  filenames = gsub(".trimmed.*", "", basename(files))
  
  res = lapply(1:length(files), function(x) {
    #dt = fread(files[x], select = c(1:3, 5:7), col.names = c("chr", "start", "end", "qual", "ref", "alt"))
    dt = fread(files[x])
    dt[, sample := filenames[x]]
    dt[, run := basename(dir)]
    dt[, depth := as.numeric(gsub(".*DP=|;.*", "", dt$V9))]
    
    # Get HQ alt/ref allele counts
    counts = dt[, gsub(".*DP4=|;.*", "", V9)]
    counts = tstrsplit(counts, ",")
    ref = as.numeric(counts[[1]]) + as.numeric(counts[[2]])
    alt = as.numeric(counts[[3]]) + as.numeric(counts[[4]])
    
    dt[, nalt := alt]
    dt[, nref := ref]
    dt = dt[, .(V1, V2, V3, V5, V6, V7, depth, sample, run, nalt, nref)]
    setnames(dt, c("chr", "start", "end", "qual", "ref", "alt", "depth", "sample", "run", "nalt", "nref"))
    return(dt)
  })
  total = rbindlist(res)
})
total = rbindlist(res)

# Merge with annotation file
total = merge(total, annot, by.x = "sample", by.y = "V2")
total[, sample := V1]
total$V1 = NULL

# discard = c("Sample 22", "Sample 42", "Sample 18", "Sample 23", "Sample 37", "Sample 44", 
#             "Sample 32", "Sample 24", "Sample 48", "Sample 35", "Sample 1")

discard = c("Sample 22", "Sample 24")
total = total[!sample %in% discard]
# discard_mut = c("T16649A", "C16650A")

# Match actual starts/end to 1 indexed start site
total[, start := start + 1]
total[, end := end + 1]
total[, vaf := nalt / (nalt + nref)]

total[, id := factor(paste0(ref, start, alt))]

# Remove muts from sequencing artifacts
# total = total[!id %in% discard_mut, ]

# Filter out indels (only keep snps) and low quality snps
total = total[nchar(ref) == 1 & nchar(alt) == 1 & qual >= 30 & depth >= 10 & nalt > nref]
total[, sample_id := paste0(sample, "_", id)]

ref = total[run == "TN21"]
reps = total[run != "TN21"]
keepmuts = reps$sample_id %in% ref$sample_id
reps = reps[keepmuts, ]

counts = reps[, .N, by = .(sample_id)]
sum_overlaps = counts[, .N, by = .(N)]
setnames(sum_overlaps, c("noverlaps", "N"))
sum_overlaps = rbind(sum_overlaps, data.table(noverlaps = c("TN21", 0), N = c(nrow(ref), nrow(ref) - sum(sum_overlaps$N))))
sum_overlaps[, fraction := N / nrow(ref)]

# Set factor
sum_overlaps[, noverlaps := factor(noverlaps, levels = c("TN21", as.character(5:0)))]


plt = ggplot(sum_overlaps, aes(x = noverlaps, y = fraction, fill = noverlaps)) +
  geom_col() +
  geom_text(aes(label = paste0("n = ", N)), vjust = -.5) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), labels = scales::label_percent()) +
  scale_fill_viridis_d() +
  labs(y = "Percentage of overlap", x = "Number of replicates with mutation") +
  theme(legend.position = "none")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/barplot-replicates-overlap-TN21",
              height = 7, width = 9)

# plot scatterplots
mut_counts = total[, .N, by = .(run, sample)]
mut_counts = dcast(mut_counts, sample~run)
mut_counts = mut_counts[!grepl("N|P", sample)]

plt1 = ggplot(mut_counts, aes(x = TN21, y = MS82_S3)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 1) +
  stat_cor()

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/scatterplot-numSNPs-TN21-MS82",
              height = 7, width = 9)

plt2 = ggplot(mut_counts, aes(x = TN21, y = MS83_S4)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 1) +
  stat_cor()

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/scatterplot-numSNPs-TN21-MS83",
              height = 7, width = 9)

plt3 = ggplot(mut_counts, aes(x = TN21, y = MS84_S5)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 1) +
  stat_cor()

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/scatterplot-numSNPs-TN21-MS84",
              height = 7, width = 9)

plt4 = ggplot(mut_counts, aes(x = TN21, y = NZ190_S1)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 1) +
  stat_cor()

save_and_plot(plt4, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/scatterplot-numSNPs-TN21-NZ190",
              height = 7, width = 9)

plt5 = ggplot(mut_counts, aes(x = TN21, y = NZ191_S2)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm", se = F, color = "red", linetype = 1) +
  stat_cor()

save_and_plot(plt5, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/scatterplot-numSNPs-TN21-NZ191",
              height = 7, width = 9)
