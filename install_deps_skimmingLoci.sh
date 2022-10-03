### First install miniconda: https://conda.io/miniconda.html

### Create an environment

# conda create --name skimmingLoci
# activate
# conda activate skimmingLoci

## install dependencies

# wget https://raw.githubusercontent.com/mreginato/skimmingLoci/main/install_deps_skimmingLoci.sh
# bash install_deps_skimmingLoci.sh

# deactivate
# conda deactivate 

# to remove the environment
# conda remove --name skimmingLoci --all

### Install packages

conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -y -c bioconda bwa samtools bcftools vcftools mafft seqtk
conda install -y -c r r-base

### Install R packages

echo "install.packages('Rcpp', repos = 'https://cloud.r-project.org', dependencies=T); install.packages('ape', repos = 'https://cloud.r-project.org', dependencies=T); install.packages('https://github.com/mreginato/skimmingLoci/raw/main/skimmingLociR_1.0.tar.gz', repos = NULL, type='source')" > install_r_packs.R
R CMD BATCH install_r_packs.R
rm install_r_packs.R
