# COVseq
Code used for "COVseq is a cost-effective workflow for mass-scale SARS-CoV-2 genomic surveillance"

Simonetti, Zhang, Harbers et al. (2021)

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
