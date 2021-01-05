# Author: Luuk Harbers
# Date: 2020-07-18
# Email: l.j.w.harbers@gmail.com

from pysam import FastxFile
import pandas as pd

# Functions
# Hamming
def cHamming(s0, s1):
	assert len(s0) == len(s1), "read-barcode and barcode unequal length"
	N = len(s0)
	i = 0
	count = 0
	for i in range(N):
		count += (s0[i] != s1[i])
	return count

def hammingDistanceLoop(barcode, barcodelist, mismatches):
	for bc in barcodelist:
		if cHamming(bc, barcode) <= mismatches:
			return bc
	return "unassigned"

# Do record processing
def process_fastq(entry_r1, entry_r2, barcodes, mismatches,
	barcode_length, umi_length):
	umi = entry_r1.sequence[0:umi_length]
	barcode = entry_r1.sequence[umi_length:umi_length + barcode_length]
	match = hammingDistanceLoop(barcode, barcodes.iloc[:,0], mismatches)
	return (match,
		f"@{entry_r1.name}_{match}_{umi} {entry_r1.comment}\n{entry_r1.sequence[20:]}\n+\n{entry_r1.quality[20:]}\n",
		f"@{entry_r2.name}_{match}_{umi} {entry_r2.comment}\n{entry_r2.sequence[20:]}\n+\n{entry_r2.quality[20:]}\n")

def iterateFastq(in_r1, in_r2, barcodes, mismatches,
	update, outdir, log, barcode_length, umi_length):
	# Initialize counts
	match_count = 0
	unassigned_count = 0

	with FastxFile(in_r1) as r1, FastxFile(in_r2) as r2:
		for entry_r1, entry_r2 in zip(r1, r2):
			result = process_fastq(entry_r1, entry_r2, barcodes, mismatches,
				barcode_length, umi_length)

			# Update counts
			if result[0] == "unassigned":
				unassigned_count += 1
			else:
				match_count += 1

			# Progress updates
			if (unassigned_count + match_count) % update == 0:
				print(f"Reads processed: {unassigned_count + match_count}")

			# Write output to file
			with open(outdir + result[0] + "_R1.fastq", "a+") as out_r1:
				with open(outdir + result[0] + "_R2.fastq", "a+") as out_r2:
					out_r1.write(result[1])
					out_r2.write(result[2])

	# Write counts
	with open(log, "w") as logfile:
		logfile.write(
			f"total reads: {match_count + unassigned_count}\n"
			f"reads with barcode: {match_count}\n"
			f"read unassigned: {unassigned_count}")

# Read in barcode file
barcodes = pd.read_csv("/mnt/AchTeraD/data/BICRO205+206/MS18_barcodes.txt",
	delimiter = ",", header = None)

iterateFastq(
	"/media/bs2-seq/BICRO201/fastq/N__Crosetto_19_07/Sample_P14955_1001/P14955_1001_S1_R1_001.fastq.gz",
	"/media/bs2-seq/BICRO201/fastq/N__Crosetto_19_07/Sample_P14955_1001/P14955_1001_S1_R2_001.fastq.gz",
	barcodes, 1, 100000,
	"/mnt/AchTeraD/data/paired_test/", "/mnt/AchTeraD/data/paired_test/log",
	8, 8)