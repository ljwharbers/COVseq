# COVseq
Scripts used for the analysis of COVseq

## Preprocessing
Use the snakemake pipeline 'snakemake-bwa_ivar-pipeline' to process COVseq libraries.
Use the snakemake pipeline 'snakemake-nebnext-bwa_ivar-pipeline' to process NEBNext libraries.

Requirements in pipeline:
bwa (0.7.17-r1188)
ivar (1.3)
samtools (1.10)
bcftools (1.10.2)
gatk (4.1.4.1)

For COVseq libraries make sure all the paths to the required scripts in the config are correct. 
Build the python/cython library using `$ python setup.py build_ext --inplace`, dependencies for demultiplexing are: pandas, argparse and pysam.

## Preparation for preprocessing COVseq libraries
For the demultiplexing of fastq files a custom python script is used. Main input required for this script is the fastq file, a list of barcodes used (no column names, just one barcode per row), the length of the barcode (default 8) and the number of mismatches allowed.
To filter out reads that map too far from cutsites you also need a bed file with all the cut site locations in the genome and the read length.
For the ivar pipeline the only extra file that is required is the file with the primer locations used for the amplicon sequencing.
