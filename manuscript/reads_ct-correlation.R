source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

options(scipen=999)
read_files = list.files("/mnt/AchTeraD/data/BICRO237/covseq/", pattern = ".bam$", recursive = T, full.names = T)
read_files = c(read_files,
               list.files("/mnt/AchTeraD/data/BICRO240/covseq/", pattern = "sorted.bam$", recursive = T, full.names = T))

ct_values = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/ct_values.tsv")

param = ScanBamParam(what=c("qname","flag"))

reads = unlist(pblapply(read_files, function(file) {
  gr = granges(readGAlignments(file, param=param))
  return(length(gr))
}, cl = 30))

total = data.table(ct = as.numeric(ct_values$V8),
                   reads = reads / 1e6)


plt = ggplot(total, aes(x = ct, y = reads)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = F, color = "red", linetype = 2) +
  stat_cor(method = "pearson") +
  labs(y = "Reads (millions)",
       x = "CT value")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/reads-ctvalues",
              height = 7, width = 7)
