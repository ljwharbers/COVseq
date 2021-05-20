## Author: Luuk Harbers
## Date: 2020-11-24
## Script for plotting Coverage of SARS-CoV-2 genomes

## Load/install packages
packages = c("data.table", "ggplot2", "pbapply", "GenomicAlignments", "Rsamtools", "IRanges")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 32

# Set files
files = list.files("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/NEBNext/variants/bam", pattern = ".bam$", full.names = T)

# Set params
param = ScanBamParam(what=c("qname","flag"))

neb = pblapply(1:length(files), function(x) {
  samplename = gsub(".trim.sorted.bam", "", basename(files[x]))
  gr = granges(readGAlignments(files[x], param=param))
  nreads = length(gr)
  cov = unlist(coverage(gr))
  # Normalize by nreads
  # cov@values = cov@values / nreads
  dt = data.table(cov@values, cov@lengths)
  values = unlist(apply(dt, 1, function(x) rep(x[1], x[2])))
  
  res = data.table(coverage = values)
  setnames(res, samplename)
  
  return(res)
}, cl = nthreads)

# Bind
neb = do.call(cbind, neb)
neb[, pos := rownames(neb)]
neb = melt(neb, id.vars = "pos")
setnames(neb, c("pos", "sample", "NEBNext"))

# same for COVseq
files = list.files("/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon/MS147_miseq/variants/bam", pattern = ".bam$", full.names = T)

cov = pblapply(1:length(files), function(x) {
  samplename = gsub(".trim.sorted.bam", "", basename(files[x]))
  gr = granges(readGAlignments(files[x], param=param))
  nreads = length(gr)
  cov = unlist(coverage(gr))
  # Normalize by nreads
  # cov@values = cov@values / nreads
  dt = data.table(cov@values, cov@lengths)
  values = unlist(apply(dt, 1, function(x) rep(x[1], x[2])))
  
  res = data.table(coverage = values)
  setnames(res, samplename)
  
  return(res)
}, cl = nthreads)

# Bind
cov = do.call(cbind, cov)
cov[, pos := rownames(cov)]
cov = melt(cov, id.vars = "pos")
setnames(cov, c("pos", "sample", "COVseq"))

# Load variants
basedir = "/mnt/AchTeraD/Documents/Projects/COVseq/data/viralrecon"
libraries = list.dirs(basedir, recursive = F, full.names = F)
libraries = libraries[grepl("MS147_miseq$|NEBNext", libraries)]

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
variants[, library := gsub("MS147_miseq", "COVseq", library)]
# Get unique variants
variants = unique(variants, by = c("POS", "REF", "ALT", "sample", "library"))
variants = variants[PASS == TRUE & ALT_FREQ > 0.75 & nchar(ALT) == 1 & TOTAL_DP > 15,]
#variants = variants[PASS == TRUE & nchar(ALT) == 1 & TOTAL_DP > 20,]

# Remove positive and negative samples (and perhaps samples <= 35ct)
variants = variants[!grepl("Sample_12|Sample_Pos|Sample_Neg|Sample_21|Sample_1$|Sample_30|Sample_3$|Sample_17|Sample_29|Sample_2$|Sample_27|Sample_18", sample),]
#variants = variants[!grepl("Sample_Pos|Sample_Neg", sample),]

# Get id and dcast
variants[, id := factor(paste0(REF, POS, ALT))]
variants[, sample_id := paste0(sample, "_", id)]

variants_collapse = variants[, .(library = paste(library, collapse = ", "),
                                 sample = sample,
                                 id = id,
                                 pos = as.character(POS)), by = sample_id]

# Reorder factors
variants_collapse[, library := factor(library, levels = c("COVseq, NEBNext", "NEBNext", "COVseq"))]

total = merge(variants_collapse, neb, by = c("pos", "sample"))
total = merge(total, cov, by = c("pos", "sample"))
total = melt(total, id.vars = c("pos", "sample", "sample_id", "library", "id"))

# Plot
plt = ggplot(total, aes(x = library, y = value, color = variable)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .2)) +
  scale_color_npg() +
  labs(y = "Depth at SNV", x = "Library in which SNV is detected", color = "Method")
  
save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/reads-snv/reads_at_snv",
              width = 9, height = 7)
