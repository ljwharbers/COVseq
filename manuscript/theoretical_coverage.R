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
bfa = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_BfaI-cutsites.bed",
            col.names = c("chr",  "start", "end"))

ref = data.table(region = c("SARS-CoV-2", "S-region"), start = c(0, 21563), end = c(29903, 25384))

nla[, enzyme := "NlaIII"]
mse[, enzyme := "MseI"]
bfa[, enzyme := "BfaI"]

read_lengths = c("300" = 280, "150" = 130, "75" = 56)


total = rbind(nla, mse, bfa)
total[, start_cover_75 := start - read_lengths["75"]]
total[, end_cover_75 := end + read_lengths["75"]]
total[, start_cover_150 := start - read_lengths["150"]]
total[, end_cover_150 := end + read_lengths["150"]]
total[, start_cover_300 := start - read_lengths["300"]]
total[, end_cover_300 := end + read_lengths["300"]]

total[start_cover_75 < 0, start_cover_75 := 0]
total[start_cover_150 < 0, start_cover_150 := 0]
total[start_cover_300 < 0, start_cover_300 := 0]

total[end_cover_75 > 29903, end_cover_75 := 29903]
total[end_cover_150 > 29903, end_cover_150 := 29903]
total[end_cover_300 > 29903, end_cover_300 := 29903]

setorder(total, start, end)

getCoverage = function(dt, total_length){
  x = dt[, .(start=min(start), stop=max(stop)),
         by=.(group=cumsum(c(1, tail(start, -1) > head(stop, -1))))]
  return(x[, sum(stop-start)/total_length])
}




# Take subsets from data with proper col names (start/stop mandatory) for full genome
nla_150 = getCoverage(dt = total[enzyme == "NlaIII", list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
mse_150 = getCoverage(dt = total[enzyme == "MseI", list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
bfa_150 = getCoverage(dt = total[enzyme == "BfaI", list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
nla_mse_150 = getCoverage(dt = total[enzyme == "NlaIII" | enzyme == "MseI", list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
nla_bfa_150 = getCoverage(dt = total[enzyme == "NlaIII" | enzyme == "BfaI", list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
mse_bfa_150 = getCoverage(dt = total[enzyme == "BfaI" | enzyme == "MseI", list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
triple_150 = getCoverage(dt = total[, list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])


s_length = ref$end[2] - ref$start[2]
nla_150_s = getCoverage(dt = total[enzyme == "NlaIII" & start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
mse_150_s = getCoverage(dt = total[enzyme == "MseI" & start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
bfa_150_s = getCoverage(dt = total[enzyme == "BfaI" & start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
nla_mse_150_s = getCoverage(dt = total[(enzyme == "NlaIII" | enzyme == "MseI") & start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
nla_bfa_150_s = getCoverage(dt = total[(enzyme == "NlaIII" | enzyme == "BfaI") & start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
mse_bfa_150_s = getCoverage(dt = total[(enzyme == "BfaI" | enzyme == "MseI") & start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
triple_150_s = getCoverage(dt = total[start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)

theoretical = data.table(readlength = c(150, 150), 
                         region = c("SARS-CoV-2", "S-region"), 
                         NlaIII = c(nla_150, nla_150_s), 
                         MseI = c(mse_150, mse_150_s),
                         BfaI = c(bfa_150, bfa_150_s),
                         "NlaIII + MseI" = c(nla_mse_150, nla_mse_150_s),
                         "NlaIII + BfaI" = c(nla_bfa_150, nla_bfa_150_s),
                         "MseI + BfaI" = c(mse_bfa_150, mse_bfa_150_s),
                         triple = c(triple_150, triple_150_s))
theoretical[triple > 1, triple := 1]

# both_75 = getCoverage(dt = total[, list(start = start_cover_75, stop = end_cover_75)], total_length = ref$end[1])
# both_150 = getCoverage(dt = total[, list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
# both_300 = getCoverage(dt = total[, list(start = start_cover_300, stop = end_cover_300)], total_length = ref$end[1])
# 
# nla_75 = getCoverage(dt = total[enzyme == "NlaIII", list(start = start_cover_75, stop = end_cover_75)], total_length = ref$end[1])
# nla_150 = getCoverage(dt = total[enzyme == "NlaIII", list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
# nla_300 = getCoverage(dt = total[enzyme == "NlaIII", list(start = start_cover_300, stop = end_cover_300)], total_length = ref$end[1])
# 
# mse_75 = getCoverage(dt = total[enzyme == "MseI", list(start = start_cover_75, stop = end_cover_75)], total_length = ref$end[1])
# mse_150 = getCoverage(dt = total[enzyme == "MseI", list(start = start_cover_150, stop = end_cover_150)], total_length = ref$end[1])
# mse_300 = getCoverage(dt = total[enzyme == "MseI", list(start = start_cover_300, stop = end_cover_300)], total_length = ref$end[1])
# 
# # Take subsets from data with proper col names (start/stop mandatory) for S-region
# s_length = ref$end[2] - ref$start[2]
# both_75_s = getCoverage(dt = total[start > 21563 & end < 25384, list(start = start_cover_75, stop = end_cover_75)], total_length = s_length)
# both_150_s = getCoverage(dt = total[start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
# both_300_s = getCoverage(dt = total[start > 21563 & end < 25384, list(start = start_cover_300, stop = end_cover_300)], total_length = s_length)
# 
# nla_75_s = getCoverage(dt = total[enzyme == "NlaIII" & start > 21563 & end < 25384, list(start = start_cover_75, stop = end_cover_75)], total_length = s_length)
# nla_150_s = getCoverage(dt = total[enzyme == "NlaIII" & start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
# nla_300_s = getCoverage(dt = total[enzyme == "NlaIII" & start > 21563 & end < 25384, list(start = start_cover_300, stop = end_cover_300)], total_length = s_length)
# 
# mse_75_s = getCoverage(dt = total[enzyme == "MseI" & start > 21563 & end < 25384, list(start = start_cover_75, stop = end_cover_75)], total_length = s_length)
# mse_150_s = getCoverage(dt = total[enzyme == "MseI" & start > 21563 & end < 25384, list(start = start_cover_150, stop = end_cover_150)], total_length = s_length)
# mse_300_s = getCoverage(dt = total[enzyme == "MseI" & start > 21563 & end < 25384, list(start = start_cover_300, stop = end_cover_300)], total_length = s_length)


# Combine and write output
theoretical = data.table(readlength = c(75, 150, 300, 75, 150, 300), 
                         region = c("SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "S-region", "S-region", "S-region"), 
                         NlaIII = c(nla_75, nla_150, nla_300, nla_75_s, nla_150_s, nla_300_s), 
                         MseI = c(mse_75, mse_150, mse_300, mse_75_s, mse_150_s, mse_300_s), 
                         "NlaIII + MseI" = c(both_75, both_150, both_300, both_75_s, both_150_s, both_300_s))
theoretical[`NlaIII + MseI` > 1, `NlaIII + MseI` := 1]
theoretical[MseI > 1, MseI := 1]

write.table(theoretical, "/mnt/AchTeraD/Documents/Projects/COVseq/theoretical_coverage.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")
