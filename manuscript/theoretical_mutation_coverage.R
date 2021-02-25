## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculation of theoretical coverage

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nla = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_NlaIII-cutsites.bed",
            col.names = c("chr", "start", "end"))
mse = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_MseI-cutsites.bed",
            col.names = c("chr",  "start", "end"))

mutations = fread("/mnt/AchTeraD/Documents/Projects/COVseq/data/top20AAmutations-frontimicrob.tsv")

nla[, enzyme := "NlaIII"]
mse[, enzyme := "MseI"]

read_lengths = c("150" = 131, "75" = 56, "300" = 281)

total = rbind(nla, mse)
total[, start_cover_75 := start - read_lengths["75"]]
total[, end_cover_75 := end + read_lengths["75"]]
total[, start_cover_150 := start - read_lengths["150"]]
total[, end_cover_150 := end + read_lengths["150"]]
total[, start_cover_300 := start - read_lengths["300"]]
total[, end_cover_300 := end + read_lengths["300"]]

total[start_cover_75 < 0, start_cover_75 := 0]
total[start_cover_150 < 0, start_cover_150 := 0]
total[end_cover_75 > 29903, end_cover_75 := 29903]
total[end_cover_150 > 29903, end_cover_150 := 29903]

setorder(total, start, end)

# Get mutation location
mutations[, location := as.numeric(gsub("[^0-9]", "", Event))]
both_75 = unlist(lapply(mutations$location, function(x) {
  return(nrow(total[start_cover_75 < x & end_cover_75 > x,]))
}))

both_150 = unlist(lapply(mutations$location, function(x) {
  return(nrow(total[start_cover_150 < x & end_cover_150 > x,]))
}))

both_300 = unlist(lapply(mutations$location, function(x) {
  return(nrow(total[start_cover_300 < x & end_cover_300 > x,]))
}))
# nla_75 = unlist(lapply(mutations$location, function(x) {
#   return(nrow(total[start_cover_75 < x & end_cover_75 > x & enzyme == "NlaIII",]))
# }))
# 
# nla_150 = unlist(lapply(mutations$location, function(x) {
#   return(nrow(total[start_cover_150 < x & end_cover_150 > x & enzyme == "NlaIII",]))
# }))
# 
# mse_75 = unlist(lapply(mutations$location, function(x) {
#   return(nrow(total[start_cover_75 < x & end_cover_75 > x & enzyme == "MseI",]))
# }))
# 
# mse_150 = unlist(lapply(mutations$location, function(x) {
#   return(nrow(total[start_cover_150 < x & end_cover_150 > x & enzyme == "MseI",]))
# }))

res = data.table(mutations = mutations$Event, "75" = both_75, "150" = both_150, "300" = both_300, source = mutations$source)

res_m = melt(res, id.vars = c("mutations", "source"))          
setorder(res_m, -variable, -value)
res_m[, mutations := factor(mutations, levels = unique(mutations))]
res_m[, variable := factor(variable, levels = c("300", "150", "75"))]

plt1 = ggplot(res_m, aes(x = mutations, y = value, fill = variable)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_d(end = 0.6) +
  labs(y = "Number of cutsites covering mutation", x = "Mutations", fill = "Read length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

res[, mutations := factor(mutations, levels = levels(res_m$mutations))]
annot_plt = ggplot(res, aes(y = 1, x = mutations, fill = source)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d(begin = 0.4) +
  theme_void() 

plt = plot_grid(plt1, annot_plt, ncol = 1, rel_heights = c(1, 0.03), align = "v")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/COVseq/Plots/mutations/theoretical-mutation-coverage-inclUK-300-150-75",
              height = 7, width = 9)
