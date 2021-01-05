#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: 2020-07-28
## Script Parse bam file and select read names from reads within n bases from a cutsite 

## Load/install packages
packages = c("data.table", "Rsamtools", "argparser")
invisible(lapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

## Parse arguments
parser = arg_parser("Parse bam file and select read names from reads within n bases from a cutsite")
parser = add_argument(parser, "--bam", help = "Path to bamfile")
parser = add_argument(parser, "--cutsites", help = "Path to bed file with locations of cutsites")
parser = add_argument(parser, "--distance", default = 20, help = "Distance treshold used (Default = 20)")
parser = add_argument(parser, "--readlength", default = 56, help = "Effective read length (taking into account removal of adapters etc.)")
parser = add_argument(parser, "--removecontigs", default = "MT|GL", 
                      help = "Contigs to remove, specified as simple regex (Default = MT|GL")
parser = add_argument(parser, "--outfile", help = "Path to outputfile")

# Parse arguments
argv = parse_args(parser)

# Check if arguments given
if(is.na(argv$bam)) stop("Filepath of --bam not specified")
if(is.na(argv$cutsites)) stop("Filepath of --cutsites not specified")
if(is.na(argv$outfile)) stop("Filepath of --outfile not specified")

# Load cutsites
sites = fread(argv$cutsites)
sites = sites[!grepl(argv$removecontigs, V1),]

# Load bam
bam = BamFile(argv$bam)

# Make bam data.table
aln = scanBam(bam)
bam_dt = data.table(chr = aln[[1]]$rname, 
                    pos = aln[[1]]$pos, 
                    strand = aln[[1]]$strand, 
                    readname = aln[[1]]$qname)

# Set keys
setkey(bam_dt, chr, pos)
setkey(sites, V1, V2)

# Overlap using 'roll = nearest' to get closest custite
bam_dt_merged = sites[bam_dt, roll = "nearest"]
bam_dt_merged[, distance := abs(V2 - V3)]

# Shift reverse reads
bam_dt_merged[strand == "-", distance := abs(distance - argv$readlength)]

# Select reads
bam_dt_keep = bam_dt_merged[distance <= argv$distance, ]

# Write readnames for gatk filtering
write.table(bam_dt_keep$readname, argv$outfile, quote = F, row.names = F, col.names = F, sep = "\t")
