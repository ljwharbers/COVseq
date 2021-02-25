## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for plotting distance density between cutsites

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nla = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_NlaIII-cutsites.bed")
mse = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_MseI-cutsites.bed")
nlamse = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_NlaIII-MseI-cutsites.bed")


nla = nla[, list(enzyme = "NlaIII", distance = data.table::shift(V2, type = "lead") - V2)]
mse = mse[, list(enzyme = "MseI", distance = data.table::shift(V2, type = "lead") - V2)]
nlamse = nlamse[, list(enzyme = "NlaIII + MseI", distance = data.table::shift(V2, type = "lead") - V2)]

cutsites = rbindlist(list(nla, mse, nlamse))
cutsites[, enzyme := factor(enzyme, levels = c("NlaIII" ,"MseI", "NlaIII + MseI"))]

plt = ggplot(cutsites, aes(x = distance, color = enzyme)) +
  geom_density(size = 1.7, adjust = 1) +
  scale_color_npg() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1050), breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)) +
  labs(y = "Density", x = "Distance between consecutive restriction sites (bp)", color = "Enzyme") +
  theme(legend.position = "top")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/cutsite-distribution/cutsite-density",
              height = 7, width = 7)
  
