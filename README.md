# skimmingLoci

## A pipeline for mapping, SNP calling and alignment generation   
   
#### Dependencies:						     

* bwa - https://github.com/lh3/bwa				     
* samtools - https://github.com/samtools/samtools		     
* bcftools - https://github.com/samtools/bcftools		     
* vcftools - https://github.com/vcftools/vcftools		     
* mafft - https://mafft.cbrc.jp/alignment/software/			     
* seqtk - https://github.com/lh3/seqtk	

**********
#### Tutorial
**********

A detailed tutorial on how to use the pipeline (including an example data set) and how to filter outlier loci is available [here](http://htmlpreview.github.io/?https://github.com/mreginato/skimmingLoci/blob/main/skimmingLoci_Tutorial.html).


**********
####  Installation
**********


Download and install miniconda3:

- https://conda.io/miniconda.html

Create a conda environment:

It is a good idea to create a conda enviroment for skimmingLoci (to install the dependencies and run the pipeline). To create the environment:

`conda create --name skimmingLoci`

Activate the environment:

`conda activate skimmingLoci`

Install dependencies:

Download the script install_deps_skimmingLoci.sh [(here)](https://github.com/mreginato/skimmingLoci/blob/main/install_deps_skimmingLoci.sh) or:
  
  `wget https://raw.githubusercontent.com/mreginato/skimmingLoci/main/install_deps_skimmingLoci.sh`
   
Run the install script 

`bash install_deps_skimmingLoci.sh`
  
Download the pipeline script:

- skimmingLoci.sh [(here)](https://raw.githubusercontent.com/mreginato/skimmingLoci/main/skimmingLoci.sh) or:  
`wget https://raw.githubusercontent.com/mreginato/skimmingLoci/main/skimmingLoci.sh`  
`chmod +x skimmingLoci.sh`

 Test if the pipeline is ok (the help should be printed if it is all good):
 
 `bash skimmingLoci.sh` or 
 `./skimmingLoci.sh`
 
 * You can place the skimmingLoci.sh file in your path to make things easier.

Deactivate skimmingLoci

To deactivate the skimmingLoci environment just type:

`conda deactivate`

* Remember to activate the environment before running the pipeline (otherwise you'll likely get an error).

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
#### Citation
**********

If you use this pipeline and or the companion R package please cite:

Reginato M. 2022. A pipeline for assembling low copy nuclear markers from plant genome skimming data for phylogenetic use. PeerJ 10:e14525 https://doi.org/10.7717/peerj.14525 

