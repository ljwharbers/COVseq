# Author: Luuk Harbers
# Date: 2020-07-18
# Email: l.j.w.harbers@gmail.com

import pandas as pd
import argparse
import demultiplex

# ARGPARSER
parser = argparse.ArgumentParser(
	description="""Moves umi and barcode to header
	and demultiplexes fastq file into separate fastq files""")
parser.add_argument("-f", "--fastq", type=str, help="Path to fastq file.")
parser.add_argument("-f2", "--fastq2", type=str, help="Path to R2 fastq file")
parser.add_argument("-p", "--paired", action="store_true")
parser.add_argument(
	"-o", "--outdir", type=str,
	help="Path to output directory.")
parser.add_argument(
	"--barcode_length", type=int,
	help="length of barcode", default=8)
parser.add_argument("--umi_length", type=int, help="length of UMI", default=8)
parser.add_argument("-l", "--log", type=str, help="Path to logfile.")
parser.add_argument("-b", "--barcodes", type=str, help="Path to barcode file.")
parser.add_argument("-m", "--mismatches", type=int, default=1,
	help="Number of mismatches allowed for barcode matching, default = 1")
parser.add_argument(
	"-u", "--update", type=int, default=1e6,
	help="Update to stdout every n lines, default = 1.000.000")

args = parser.parse_args()


# Read in barcode file
barcodes = pd.read_csv(args.barcodes, delimiter = ",", header = None)

if args.paired:
	demultiplex.iterateFastq_paired(
		args.fastq, args.fastq2, barcodes, args.mismatches, args.update,
	args.outdir, args.log, args.barcode_length, args.umi_length)
else:
	demultiplex.iterateFastq(
		args.fastq, barcodes, args.mismatches, args.update,
		args.outdir, args.log, args.barcode_length, args.umi_length)
