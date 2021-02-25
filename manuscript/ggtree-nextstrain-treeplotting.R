## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ggtree", "tidytree", "dplyr", "rjson")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

tree = read.tree("/home/luukharbers/ncov/results/italy_covseq/tree_raw.nwk")

# Prepare clade information
clades = unlist(fromJSON(file = "/home/luukharbers/ncov/results/italy_covseq/clades.json")[2])
clades = clades[3:length(clades)]
clades = data.table(clade = clades,
                    node = names(clades))
clades[, node := gsub("nodes.|.clade.*", "", node)]

# Make tibble and merge with clade information and country/group
tree_data = tree %>% as.treedata %>% as_tibble
tree_data = full_join(tree_data, clades[!grepl("NODE", node),], by = c("label" = "node"))

#tree_data = unique(tree_data)
tree_data = tree_data %>% mutate(
  group = case_when(grepl("^Italy", label) ~ "Italy",
                    grepl("TN11", label) ~ "COVseq - 30",
                    grepl("TN21", label) ~ "COVseq - 50",
                    grepl("Wuhan", label) ~ "Wuhan",
                    grepl("NEB", label) ~ "NEBNext"),
  group = ifelse(is.na(group), "Rest", group),
  clade = ifelse(is.na(clade), "", clade)
)

list_groups = list(italy = grep("^Italy", tree_data$label),
  "COVseq - 30" = grep("TN11", tree_data$label),
  "COVseq - 50" = grep("TN21", tree_data$label),
  "Wuhan" = grep("Wuhan", tree_data$label),
  "NEBNext" = grep("NEB", tree_data$label))
list_groups$rest = setdiff(1:1086, unlist(list_groups))

list_clades = list("19A" = grep("19A", tree_data$clade),
                   "19B" = grep("19B", tree_data$clade),
                   "20A" = grep("20A", tree_data$clade),
                   "20B" = grep("20B", tree_data$clade),
                   "20C" = grep("20C", tree_data$clade))


tree_data = groupOTU(tree_data, list_groups, group_name = "group_otu")
tree_data = groupOTU(tree_data, list_clades, group_name = "clade_otu")

tree = as.treedata(tree_data)

plt1 = ggtree(tree, aes(color = group_otu)) + 
  scale_color_brewer(palette = "Dark2") +
  # geom_highlight(node = 1087, fill = brewer.pal(5, "Dark2")[1], alpha = 0.3) +
  # geom_highlight(node = 2081, fill = brewer.pal(5, "Dark2")[2], alpha = 0.3) +
  # geom_highlight(node = 1117, fill = brewer.pal(5, "Dark2")[3], alpha = 0.3) +
  # geom_highlight(node = 1955, fill = brewer.pal(5, "Dark2")[4], alpha = 0.3) +
  # geom_highlight(node = 1265, fill = brewer.pal(5, "Dark2")[5], alpha = 0.3) +
  geom_cladelabel(node = 2118, label = "19A") +
  geom_cladelabel(node = 2081, label = "19B") +
  geom_cladelabel(node = 2076, label = "19A") +
  geom_cladelabel(node = 1265, label = "20B") +
  geom_cladelabel(node = 1955, label = "20C") +
  geom_cladelabel(node = 1117, label = "20A")

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/phylogeny/IQTree-TN11-TN21-NEB",
              height = 16, width = 14)


# Plot without NEB
drop = tree@phylo$tip.label[grepl("NEB", tree@phylo$tip.label)]
tree_noneb = treeio::drop.tip(tree, drop)
tree_noneb = tree_noneb %>% as_tibble

list_groups = list(italy = grep("^Italy", tree_data$label),
                   "COVseq - 30" = grep("TN11", tree_data$label),
                   "COVseq - 50" = grep("TN21", tree_data$label),
                   "Wuhan" = grep("Wuhan", tree_data$label))
list_groups$rest = setdiff(1:1065, unlist(list_groups))

list_clades = list("19A" = grep("19A", tree_data$clade),
                   "19B" = grep("19B", tree_data$clade),
                   "20A" = grep("20A", tree_data$clade),
                   "20B" = grep("20B", tree_data$clade),
                   "20C" = grep("20C", tree_data$clade))


tree_noneb = groupOTU(tree_noneb, list_groups, group_name = "group_otu")
tree_noneb = groupOTU(tree_noneb, list_clades, group_name = "clade_otu")
tree_noneb = as.treedata(tree_noneb)

plt2 = ggtree(tree_noneb, aes(color = group_otu)) + 
  scale_color_brewer(palette = "Dark2", direction = -1) +
  # geom_highlight(node = 1087, fill = brewer.pal(5, "Dark2")[1], alpha = 0.3) +
  # geom_highlight(node = 2081, fill = brewer.pal(5, "Dark2")[2], alpha = 0.3) +
  # geom_highlight(node = 1117, fill = brewer.pal(5, "Dark2")[3], alpha = 0.3) +
  # geom_highlight(node = 1955, fill = brewer.pal(5, "Dark2")[4], alpha = 0.3) +
  # geom_highlight(node = 1265, fill = brewer.pal(5, "Dark2")[5], alpha = 0.3) +
  geom_cladelabel(node = 2076, label = "19A") +
  geom_cladelabel(node = 1086, label = "20A") +
  geom_cladelabel(node = 2037, label = "19B") +
  geom_cladelabel(node = 1946, label = "20C") +
  geom_cladelabel(node = 1243, label = "20B")
  # geom_cladelabel(node = 1117, label = "20A")

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/phylogeny/IQTree-TN11-TN21",
              height = 16, width = 14)
