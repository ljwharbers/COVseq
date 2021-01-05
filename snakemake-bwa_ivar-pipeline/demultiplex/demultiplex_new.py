#!/usr/bin/python3

# Author: Luuk Harbers
# Date: 2020-07-18
# Email: l.j.w.harbers@gmail.com

import pandas as pd
import gzip, argparse
from pysam import FastxFile
import multiprocessing as mp

# ARGPARSER time
parser = argparse.ArgumentParser(
	description="""Moves umi and barcode to header
	and demultiplexes fastq file into separate fastq files""")
parser.add_argument("-f", "--fastq", type=str, help="Path to fastq file.")
parser.add_argument(
	"-o", "--outdir", type=str,
	help="Path to output directory.")
parser.add_argument("-l", "--log", type=str, help="Path to logfile.")
parser.add_argument("-b", "--barcodes", type=str, help="Path to barcode file.")
parser.add_argument("-m", "--mismatches", type=int, default=1,
	help="Number of mismatches allowed for barcode matching, default = 1")
parser.add_argument(
	"-u", "--update", type=int, default=1e6,
	help="Update to stdout every n lines, default = 1.000.000")

args = parser.parse_args()


# Functions

# Open file
def process_lines(lines=None):
	ks = ['name', 'sequence', 'optional', 'quality']
	return {k: v for k, v in zip(ks, lines)}

def process_wrapper(
	lineByte, file, barcodes, mismatches,
	update, outdir, log):
	n = 4
	match_count = 0
	unassigned_count = 0

	# Process fastq
	with open(file) as f:
		lines = []
		f.seek(lineByte)
		for line in f:
			lines.append()
			if len(lines) == 4:
				record = process_lines(lines)
				result = process_fastq(record, barcodes, mismatches)

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
			lines = []


# Hamming
def cHamming(s0, s1):
	N = len(s0)
	i, count = 0
	for i in range(N):
		count += (s0[i] != s1[i])
	return count

def hammingDistanceLoop(barcode, barcodelist, mismatches):
	for bc in barcodelist:
		if cHamming(bc, barcode) <= mismatches:
			return bc
	return "unassigned"

# Do record processing
def process_fastq(record, barcodes, mismatches):
	umi = record.sequence[0:8]
	barcode = record.sequence[8:16]
	match = hammingDistanceLoop(barcode, barcodes.iloc[:,0], mismatches)

	# Get new string
	space_index = record.name.index(" ")
	new_name = f"{record.name[:space_index]}_{match}_{umi} {record.name[space_index:]}"
	return (match, f"@{new_name}\n{record.sequence[20:]}\n+\n{entry.quality[20:]}\n")
'''
def iterateFastq(str fastq, barcodes, int mismatches,
	int update, str outdir, str log):
	# Initialize counts

	int match_count = 0
	int unassigned_count = 0

	with FastxFile(fastq) as input_handle:
		for entry in input_handle:
			result = process_fastq(entry, barcodes, mismatches)

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
'''
# Read in barcode file
barcodes = pd.read_csv(args.barcodes, delimiter = ",")

pool = mp.Pool(4)
jobs = []

with open(args.fastq, "rb") as f:
	nextLineByte = f.tell()
	for line in f:
		jobs.append(pool.apply_async(process_fastq, (nextLineByte, args.fastq, barcodes, args.mismatches,
		args.update, args.outdir, args.log)))
		nextLineByte = f.tell()

for job in jobs:
	job.get()

pool.close()

