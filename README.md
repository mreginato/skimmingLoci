# skimmingLoci

## A wrapper for mapping, SNP calling and alignment generation   
   
#### Dependencies:						     

* bwa - https://github.com/lh3/bwa				     
* samtools - https://github.com/samtools/samtools		     
* bcftools - https://github.com/samtools/bcftools		     
* vcftools - https://github.com/vcftools/vcftools		     
* mafft - https://mafft.cbrc.jp/alignment/software/			     
* seqtk - https://github.com/lh3/seqtk			  


**********
####  Installation
**********


Download and install miniconda3:

- https://conda.io/miniconda.html

Create a conda environment:

`conda create --name skimmingLoci`

Activate the environment:

`conda activate skimmingLoci`

Install dependencies:

Download the script [install_deps_skimmingLoci.sh](https://github.com/mreginato/skimmingLoci/blob/main/install_deps_skimmingLoci.sh) (here) or:
  
  `wget https://raw.githubusercontent.com/mreginato/skimmingLoci/main/install_deps_skimmingLoci.sh`
 
Run the install script 

`bash install_deps_skimmingLoci.sh`
  
Download the pipeline script:

- [skimmingLoci.sh](https://raw.githubusercontent.com/mreginato/skimmingLoci/main/skimmingLoci.sh) (here) or:  
`wget https://raw.githubusercontent.com/mreginato/skimmingLoci/main/skimmingLoci.sh`

 Test if the pipeline is ok (the help should be printed if it is all good):
 
 `bash skimmingLoci.sh`
 
 * you can place the skimmingLoci.sh file in your path to make things easier.

**********
#### skimmingLoci
**********

The pipeline has 4 steps, which should be executed in the folder where 
the reference and fastq.gz files are located, in the following order:

* **Step 1. Map reads to reference**


This step uses the reference (*.fasta) and 
fastq.gz files, paired (*.R1.fastq.gz, *.R2.fastq.gz) or single 
(*.R1.fastq.gz):	

`skimmingLoci.sh -m`
	

* **Step 2. Call snps**


This step uses the reference (*.fasta) and the *.bam files generated 
in step 1:		
			
`skimmingLoci.sh -s`
	

* **Step 3. Generate consensus sequences**


This step uses the reference (*.fasta) and the *.vcf files generated 
in step 2:
	
`skimmingLoci.sh -c`
  

* **Step 4. Generate alignments**


This step uses the reference (*.fasta) and the *.fas files generated
in step 3:
  
`skimmingLoci.sh -a`
  

**********
#### Key parameters
**********

* **Step 3. Generate consensus sequences**

*--min-depth*: [numeric] mininum depth to export a base call to consensus


*--max-depth*: [numeric] maximum depth to export a base call to consensus


*Step 2 generates histograms of depth distribution for each and all 
samples together to help pick those numbers (*.depth.jpg).

In the following example bases with depth greater than 50 and smaller
than 4 would be excluded from the alignment (i.e., exported as "N"):

`skimmingLoci.sh -c --min-depth 4 --max-depth 50`


* **Step 4. Generate alignments**



*--min-cov*: [numeric] mininum coverage (% of the reference) of a sequence to be kept in the alignment


*--min-seq-number*: [numeric] mininum number of sequences to keep an alignment in the final data set



*Step 3 generates histograms of coverage distribution for each and all 
samples together to help pick those numbers (*.coverage.jpg).

In the following example consensus sequences with less than 20% of coverage from the reference will be excluded, and loci with less than 8 sequences will be excluded:

`skimmingLoci -c --min-cov 0.2 --min-seq.number 8`


**********

* **Statistics**

The pipeline generates basic statistics throughout the process. However, for big data sets this might slow the process. Also, if there is any problem with the R packages the pipeline might crash. To avoid this use the flag "--no-stats".

**********
#### Tutorial
**********

A detailed tutorial on how to use the pipeline (including an example data set) and how to filter outlier loci is available [here](http://htmlpreview.github.io/?https://github.com/mreginato/skimmingLoci/blob/main/skimmingLoci_Tutorial.html).


**********
#### Citation
**********

If you use this pipeline and or the companion R package please cite:

Reginato, M. (submitted). A pipeline for assembling low copy nuclear markers from plant genome skimming data for phylogenetic use. 

