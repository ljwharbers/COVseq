## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "gridExtra", "grid", "ggplot2", "RColorBrewer")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

#top folder
run_name = "BICRO236"
top_folder = "/mnt/AchTeraD/data/BICRO236/"
save_folder = "/mnt/AchTeraD/Documents/Projects/COVseq/QC/"

#list libraries
libraries = list.dirs(top_folder, recursive = F, full.names = F)
libraries = libraries[grepl("NZ|MS", libraries)] #can just change the grepl expression in most cases to match library directories


#get input and output reads after barcode/cutsite extraction
total_reads = data.table(library = character(), variable = numeric(), value = numeric())
sample_reads = data.table(V1 = character(), V2 = numeric(), sample = character(), library = character())

# Loop through libraries
for(library in libraries) {
  
  # Get total reads and demultiplexed reads
  count_file = list.files(paste0(top_folder, library, "/logs/"), pattern = "demultiplex.log", full.names = T)
  counts = fread(count_file, sep = ":", header = F)
  
  # Get number of mapped reads and reads mapped within range of cutsite
  cutsitefiles = list.files(paste0(top_folder, library, "/logs/"), pattern = "filter", full.names = T)
  cutsitefilter = lapply(cutsitefiles, function(file){
    data = fread(file, sep = ":")
    data[, sample := gsub("-.*|.*\\/", "", file)]
  })
  cutsite_dt = rbindlist(cutsitefilter)
  cutsite_dt[, V1 := ifelse(grepl("pre", V1), "mapped", "mapped in cutsite range")]
  
  # Get deduplicated reads
  dedup_file = list.files(paste0(top_folder, library, "/logs/"), pattern = "deduplication.log", full.names = T)
  dedup = fread(dedup_file, skip = "INFO Number of reads out", header = F)
  mapped_dedup = as.numeric(gsub(".*: ", "", dedup$V2[1]))
  
  
  # Make dt for total reads
  total = data.table(library = library, "total reads" = counts$V2[1], "with barcode" = counts$V2[2], 
                     "in cutsite range" = sum(cutsite_dt[V1 == "mapped in cutsite range"]$V2), 
                     deduplicated = mapped_dedup)
  total = melt(total, id.vars = "library")
  
  # reads to millions and set factors
  total[, value := value / 1e6]
  
  # rbind with total DTs
  total_reads = rbind(total_reads, total)
}

# Set factor levels for total_reads
total_reads[, variable := factor(variable, levels = c("total reads", "with barcode", "in cutsite range", "deduplicated"))]

# Plot
plt = ggplot(total_reads, aes(x = library, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Reads (millions)",
       x = "") +
  theme(legend.title = element_blank(),
        legend.position = "top")


# Save plots
if(!dir.exists(save_folder)) dir.create(save_folder, recursive = T)
save_and_plot(grid.draw(plt), paste0(save_folder, "sequence_reads_QC"), width = 12, height = 8)


