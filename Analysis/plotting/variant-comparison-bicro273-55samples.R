
## Load/install packages
packages = c("data.table", "UpSetR", "GGally")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Ct values
ct = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/55-patients_CT-clean.tsv")
ct[, sample := gsub(" ", "_", )]

basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
#libraries = libraries[!grepl("nextflow|work", libraries)]
libraries = c("NZ242", "NZ243", "NZ244", "NZ245")

# Load in variants
variants = lapply(libraries, function(lib) {
  files = list.files(paste0(basedir, "/", lib, "/variants/ivar/"), full.names = T, pattern = ".tsv")
  variants = lapply(files, function(file) {
    dt = fread(file)
    dt[, sample := gsub(".tsv", "", basename(file))]
    return(dt)
  })
  variants = rbindlist(variants)
  variants[, library := lib]
  return(variants)
})

variants = rbindlist(variants)

# Get unique variants
variants = unique(variants, by = c("POS", "REF", "ALT", "sample", "library"))
variants = variants[PASS == TRUE & ALT_FREQ > 0.75 & nchar(ALT) == 1 & TOTAL_DP > 15,]

variants[, id := factor(paste0(REF, POS, ALT))]
variants[, tosum := 1]
variants[, sample_id := paste0(sample, "_", id)]

# Select <=35 Ct samples
variants = variants[!grepl("Sample_N|Sample_P|Sample_52|Sample_54", sample)]

variants_mat = dcast(variants, sample_id ~ library, fun.aggregate = function(x) sum(x), value.var = "tosum")


plt1 = upset(variants_mat, sets = libraries, nintersect = 20, keep.order = T, order.by = "freq", point.size = 3,
             line.size = 1)

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/upsetPlots/NZ242-NZ245_Ct<=35",
              width = 8, height = 5)

num_mat = variants[, .N, by = .(library, sample)]
num_mat = dcast(num_mat, sample ~ library, value.var = "N")

plt2 = ggpairs(num_mat[, 2:ncol(num_mat)]) +
  theme_base()

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/scatter-muts/55-replicates-NZ242-NZ245",
              width = 11, height = 14)

# cor
cors = melt(cor(num_mat[,2:5]))

plt3 = ggplot(cors, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3))) +
  scale_x_discrete("", position = "top") +
  scale_fill_viridis("Pearson's\ncorrelation", option="inferno", begin = .8, direction = -1) +
  labs(y = "", x = "") +
  theme(axis.line = element_blank(),
        text = element_text(family = "Helvetica"),
        axis.text.x.top = element_text(angle = 45, hjust = 0))

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/corr-snvs/Correlation-numsnvs-55samples_NZ242-NZ245",
              height = 7, width = 8)
