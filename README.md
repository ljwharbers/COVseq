# COVseq

Code used for "COVseq is a cost-effective workflow for mass-scale SARS-CoV-2 genomic surveillance"

Simonetti, Zhang, Harbers et al. (2021)

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

## Demultiplexing COVseq libraries
Since COVseq libraries are multiplexed libraries (multiple samples in one library). These libraries need to be demultiplexed. To demultiplex COVseq fastq files you can run the script in the `Demultiplex` folder.

Build the python/cython library using `$ python setup.py build_ext --inplace` and make sure the following dependencies are met: `pandas`, `argparse` and `pysam`.

Following this you can run the demultiplexing. An example command would be: `$ demultiplex_withcython.py -f {fastq1} -f2 {fastq2} --paired -o {output} -l {logfile-output} -b {list-of-barcodes} -m {mismatches}`
With the {list-of-barcodes} being a text file with 1 COVseq barcode per line (no headers). For more information regarding the different commands you can run `demultiplex_withcython.py --help`.

## Running FastQ-Screen
To get information regarding to the mapping percentages of sample specific fastq files (output of demultiplexing step) we used FastQ-Screen (version 0.14.1) https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/.

Simply add another entry in the database with the SARS-CoV-2 reference genome and run with default settings. 

Due to sensitive patient information we do not share any patient specific fastq files.

## Running nf-co.re/viralrecon
For further processing of fastq files we used the nextflow based pipeline from nf-core called viralrecon (version 1.1.0) https://nf-co.re/viralrecon/1.1.0. For any extra information or troubleshooting please check out their website and/or join the slack channel for the specific pipeline. 

Commands we used to run this pipeline are as follows: `$ nextflow run nf-core/viralrecon --input {samplesheet.csv} --genome 'NC_045512.2' --fasta {sarscov2-fastafile} --save_reference --protocol amplicon --amplicon_bed {ampliconbedfile} --skip_assembly --skip_markduplicates --skip_mosdepth --callers ivar --outdir {outdir} -profile docker --max-cpus 40 -r 1.1.0`

This generates all output used for further plotting and visualization in the manuscript.

## Pangolin lineage
To determine pangolin lineage in all the samples. Per library the consensus sequence was concatenated. These sequences are available in the viralrecon output folder `output/variants/ivar/consensus/`. Pangolin (version 2.3.2) can then be run with default settings.

For more information regarding pangolin please go to: https://github.com/cov-lineages/pangolin and https://cov-lineages.org/

## Further analysis
For further analysis of viralrecon output there are 2 basic scripts in the `Analysis` folder. These scripts are used to filter SNVs and combine SNVs from multiple samples/libraries and to retrieve bam-file metrics (Breadth of Coverage, number of reads). The output of this can be used to generate all the plots in the manuscript. These scripts are in `Analysis/plotting` but are not ready to be run due to personal paths. 

If you have any questions or if you are missing files feel free to email me at: luuk.harbers@scilifelab.se
