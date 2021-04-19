## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ggtree", "tidytree", "tidyr", "rjson", "colorblindr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

tree = read.tree("/home/luukharbers/ncov/results/italy_covseq/tree.nwk")
annot = fread("~/Downloads/Campioni_OAS_feb2021.csv", col.names = c("sample", "person", "date", "location", "S", "Orf", "N"))
annot[, sample := paste0("OAS-Sample_", sample)]

# Prepare clade information
clades = unlist(fromJSON(file = "/home/luukharbers/ncov/results/italy_covseq/clades.json")[2])
clades = clades[3:length(clades)]
clades = data.table(clade = clades,
                    node = names(clades))
clades[, node := gsub("nodes.|.clade.*", "", node)]

# Make tibble and merge with clade information and country/group
tree_data = tree %>% as.treedata %>% as_tibble
tree_data = full_join(tree_data, clades[!grepl("NODE", node),], by = c("label" = "node"))

list_groups = list(Italy = grep("^Italy", tree_data$label),
                   "OAS-29" = grep("^Sample_[1-2][0-9]", tree_data$label),
                   "OAS-95" = grep("AOC", tree_data$label),
                   "CCI-55" = grep("^Sample_[3-8][0-9]", tree_data$label))
list_groups$rest = setdiff(1:1982, unlist(list_groups))

clades = unique(tree_data$clade)
clades = gsub(" .*", "", clades)
clades = clades[!is.na(clades)]
list_clades = 
  lapply(clades, function(clade){
    grep(clade, tree_data$clade)
  })
names(list_clades) = clades

# Add sample information
tree_data = full_join(tree_data, annot[, .(sample, location)], by = c("label" = "sample"))
tree_data = tree_data %>% replace_na(list(location = "rest"))
list_location = 
  lapply(unique(tree_data$location), function(loc){
    grep(loc, tree_data$location)
  })
names(list_location) = unique(tree_data$location)


tree = groupOTU(tree, list_groups, group_name = "group_otu")
tree = groupOTU(tree, list_clades, group_name = "clade_otu")
#tree = groupOTU(tree, list_location, group_name = "location_otu")

tree = as.treedata(tree)

plt1 = ggtree(tree, aes(color = group_otu)) + 
  scale_color_viridis_d(option = "E", direction = -1) +
  labs(color = "Cohort")

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/phylogeny/phylo-tree_samples",
              height = 16, width = 14)

plt2 = ggtree(tree, aes(color = clade_otu)) + 
  scale_color_hue(direction = -1) +
  labs(color = "Clade")

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/phylogeny/phylo-tree_clades",
              height = 16, width = 14)


plt3 = ggtree(tree) + 
  geom_tiplab(align = T, size = 1, aes(subset = grepl("AOC", label))) +
  scale_color_npg() 

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/revision/phylogeny/phylo-tree_AOC--labels",
              height = 30, width = 14)
