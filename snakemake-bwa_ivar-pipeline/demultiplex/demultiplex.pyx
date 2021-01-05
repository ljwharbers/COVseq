#!/usr/bin/python3

# Author: Luuk Harbers
# Date: 2020-07-18
# Email: l.j.w.harbers@gmail.com

cimport cython
from pysam import FastxFile

# Functions
# Hamming
def cHamming(str s0, str s1):
	assert len(s0) == len(s1), "read-barcode and barcode unequal length"
	cdef:
		int N = len(s0)
		int i, count = 0
	for i in range(N):
		count += (s0[i] != s1[i])
	return count

def hammingDistanceLoop(str barcode, barcodelist, int mismatches):
	for bc in barcodelist:
		if cHamming(bc, barcode) <= mismatches:
			return bc
	return "unassigned"

# Do record processing
def process_fastq(entry, barcodes, mismatches,
	int barcode_length, int umi_length):
	umi = entry.sequence[0:umi_length]
	barcode = entry.sequence[umi_length:umi_length + barcode_length]
	match = hammingDistanceLoop(barcode, barcodes.iloc[:,0], mismatches)
	return (match, f"@{entry.name}_{match}_{umi} {entry.comment}\n{entry.sequence[20:]}\n+\n{entry.quality[20:]}\n")

def iterateFastq(str fastq, barcodes, int mismatches,
	int update, str outdir, str log, int barcode_length, int umi_length):
	# Initialize counts
	cdef:
		int match_count = 0
		int unassigned_count = 0

	with FastxFile(fastq) as input_handle:
		for entry in input_handle:
			result = process_fastq(entry, barcodes, mismatches,
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
			with open(outdir + result[0] + ".fastq", "a+") as outfile:
				outfile.write(result[1])

	# Write counts
	with open(log, "w") as logfile:
		logfile.write(
			f"total reads: {match_count + unassigned_count}\n"
			f"reads with barcode: {match_count}\n"
			f"read unassigned: {unassigned_count}")

# Functions for paired processing
# Do record processing
def process_fastq_paired(entry_r1, entry_r2, barcodes, mismatches,
	barcode_length, umi_length):
	umi = entry_r1.sequence[0:umi_length]
	barcode = entry_r1.sequence[umi_length:umi_length + barcode_length]
	match = hammingDistanceLoop(barcode, barcodes.iloc[:,0], mismatches)
	return (match,
		f"@{entry_r1.name}_{match}_{umi} {entry_r1.comment}\n{entry_r1.sequence[20:]}\n+\n{entry_r1.quality[20:]}\n",
		f"@{entry_r2.name}_{match}_{umi} {entry_r2.comment}\n{entry_r2.sequence[20:]}\n+\n{entry_r2.quality[20:]}\n")

def iterateFastq_paired(in_r1, in_r2, barcodes, mismatches,
	update, outdir, log, barcode_length, umi_length):
	# Initialize counts
	match_count = 0
	unassigned_count = 0

	with FastxFile(in_r1) as r1, FastxFile(in_r2) as r2:
		for entry_r1, entry_r2 in zip(r1, r2):
			result = process_fastq_paired(entry_r1, entry_r2, barcodes, mismatches,
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
