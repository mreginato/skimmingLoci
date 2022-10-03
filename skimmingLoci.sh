#!/bin/bash

# ------------------------------------------------------------------ #
#    skimmingLoci 			        Marcelo Reginato     #
# 								     #	
#    A wrapper for mapping, SNP calling and alignment generation     #
#     								     #
#    Dependencies:						     #
# 								     #	
#    bwa - https://github.com/lh3/bwa				     #
#    samtools - https://github.com/samtools/samtools		     #
#    bcftools - https://github.com/samtools/bcftools		     #
#    vcftools - https://github.com/vcftools/vcftools		     #
#    mafft - https://mafft.cbrc.jp/alignment/software/		     #
#    seqtk - https://github.com/lh3/seqtk			     #
# ------------------------------------------------------------------ #

VERSION=1.0.0

usage()
{
    echo -e "\n usage: skimmingLoci -r fasta.file -m -s -c -a | [parameters]\n"
    echo "where:"
    echo -e "general parameters:"
    echo -e "\t-r [--reference] 	A fasta file [reference.fasta]"
    echo -e "\t-i [--input-dir] 	Input folder [wd]. Directory where the reference and fastq.gz files are located"
    echo -e "\t-o [--output-dir] 	Output folder [wd]. Directory where the output will be exported to"
    echo -e "\t-x [--cores]		Number of cores [2] \n"
    echo -e "pipeline steps:"
    echo -e "\t-m [--mapping]   	Map reads - Generates *.bam files. Requires the reference and fastq.gz files"
    echo -e "\t-s [--snps]		Call snps - Generates *.vcf files. Requires the reference and bam files"
    echo -e "\t-c [--consensus] 	Consensus sequences - Generates *.fas files for each sample. Requires the reference and vcf files"
    echo -e "\t-a [--aligning] 	Loci alignments - Generates aligned *.fas files for each locus. Requires the samples *.fas files "
    echo -e "\t-N [--concatenate] 	Concatenate final alignments - Generates *.phy file and a map file"
    echo -e "\t-Z [--no-stats] 	Do not generate statistics"
    echo -e "reads:"
    echo -e "\t-p [--single] 		Paired or single reads [paired]\n"
    echo -e "mapping parameters:"
    echo -e "\t-A [--bwa-a] 		Score for a sequence match [1]"
    echo -e "\t-B [--bwa-b] 		Penalty for a mismatch [3]"
    echo -e "\t-O [--bwa-o ]		Gap open penalty for deletions and insertions [5]\n"
    echo -e "consensus parameters:"
    echo -e "\t-d [--min-depth]	Mininum depth to keep a base call in the consensus [3]"
    echo -e "\t-D [--max-depth]	Maximum depth to keep a base call in the consensus [100]"
    echo -e "\t-Q [--min-q] 		Minimum Q to keep a base call in the consensus [10]\n"
    echo -e "alignment parameters:"
    echo -e "\t-S [--min-seq-num]	Minimum number of sequences to keep an alignment [4]"
    echo -e "\t-C [--min-cov] 		Minimum coverage to keep a sequence in the alignment [0.1]
				values of --min-cov < 1 are taken as percent of the reference
				values --min-cov > 1 as number of bases"
    echo -e "\t-T [--trim]		Remove alignment positions with less than --min-seq-num\n"	
    echo -e "examples: \n"
    echo -e "\t# Mapping only "
    echo -e "\tskimmingLoci -r reference.fasta -m \n"
    echo -e "\t# SNP calling only "
    echo -e "\tskimmingLoci --reference reference.fasta --snps \n"
    echo -e "\t# Mapping, calling snps, generating consensus and loci alignments, with mininum depth of 6 for consensus"
    echo -e "\tskimmingLoci --reference reference.fasta -m -s -c -a --min-depth 6\n\n"
    echo -e "IMPORTANT: fastq.gz files should be named as:\n"
    echo -e "\tSomething.R1.fastq.gz , Something.R2.fastq.gz \n"
    echo -e "examples: \n"
    echo -e "\tLeandra_1.R1.fastq.gz , Leandra_1.R2.fastq.gz "
    echo -e "\tLeandra_melastomoides.R1.fastq.gz , Leandra_melastomoides.R2.fastq.gz "
    echo -e "\tLeandra_melastomoides_5.R1.fastq.gz , Leandra_melastomoides_5.R2.fastq.gz \n"
    echo -e "IMPORTANT: The fasta reference should be named as:\n"
    echo -e "\tsomething.fasta \n\n"
}

# ------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------ #

progress_bar() {
    local w=36 p=$1;  shift
    printf -v dashes "%*s" "$(( $p*$w/100 ))" ""; dashes=${dashes// /#};
    printf "\r\e[K[%-*s] %3d %% %s" "$w" "$dashes" "$p" "$*"; 
}

# ------------------------------------------------------------------ #
# Variables
# ------------------------------------------------------------------ #

REFERENCE=$(printf "$WD/reference.fasta\n")
MAPPING=0
CALLING=0
CONSENSUS=0
ALIGNING=0
CONCATENATE=0
STATS=1
READS=1
WD=$(pwd)
WD_IN=$(pwd)
CORES=2

### Bwa parameters 
BWA_A=1
BWA_B=3
BWA_O=5

### vcf2fq parameters
DEPTH_MIN=3
DEPTH_MAX=100
Q_MIN=10 

### alignment parameters
MIN_SEQ_NUM=4
MIN_COV=0.1
TRIM=0

# ------------------------------------------------------------------ #
# Options processing
# ------------------------------------------------------------------ #

### Get options

while [ "$1" != "" ]; do
    case $1 in
        -r | --reference )      shift
                                REFERENCE=$1
                                ;;
        -m | --mapping )        MAPPING=1
                                ;;
	-s | --snps )		CALLING=1
				;;
	-c | --consensus )	CONSENSUS=1
				;;
	-a | --aligning )	ALIGNING=1
				;;
	-i | --input-dir )	shift
				WD_IN=$1
				;;
	-o | --output-dir )	shift
				WD=$1
				;;
	-x | --cores )		shift
				CORES=$1
				;;
	-p | --single )		READS=0
				;;
	-A | --bwa-a )		shift
				BWA_A=$1
				;;
	-B | --bwa-b )		shift
				BWA_B=$1
				;;
	-O | --bwa-o )		shift
				BWA_O=$1
				;;
	-d | --min-depth )	shift
				DEPTH_MIN=$1
				;;
	-D | --max-depth )	shift
				DEPTH_MAX=$1
				;;
	-Q | --min-q )		shift
				Q_MIN=$1
				;;
	-S | --min-seq-num )	shift
				MIN_SEQ_NUM=$1
				;;
	-C | --min-cov )	shift
				MIN_COV=$1
				;;
	-T | --trim )	        shift
				TRIM=1
                                ;;
	-N | --concatenate )	CONCATENATE=1
                                ;;
	-Z | --no-stats )	STATS=0
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

# ------------------------------------------------------------------ #
# Dirs 
# ------------------------------------------------------------------ #

BAM_DIR=$(printf "$WD/1_mapping\n")
VCF_DIR=$(printf "$WD/2_calling\n")
CONSENSUS_DIR=$(printf "$WD/3_consensus\n")
ALIGNS_DIR=$(printf "$WD/4_aligns_all\n")
ALIGNS_F_DIR=$(printf "$WD/5_aligns_final\n")
CONC_DIR=$(printf "$WD/6_concatenated\n")

cd $WD

# ------------------------------------------------------------------ #
# Check dependencies
# ------------------------------------------------------------------ #


DEPS=(bwa samtools bcftools vcftools mafft seqtk R)
NUMDEP=0
for i in "${DEPS[@]}"
do
	if which $i &> /dev/null; then
		x=0
	else
    		echo "The dependency" $i "is not installed"
    		NUMDEP=$((NUMDEP + 1))
	fi
done

if [ $NUMDEP -gt 0 ]; then
	echo -e "\nPlease install all required dependencies"
	exit 1
fi


# ------------------------------------------------------------------ #
# Messages 
# ------------------------------------------------------------------ #

if [ "$MAPPING" = 0 ] && [ "$CALLING" = 0 ] && [ "$ALIGNING" = 0 ] && [ "$CONSENSUS" = 0 ] && [ "$CONCATENATE" = 0 ]; then
	usage
	exit 1
fi

echo -e "\n skimmingLoci" $VERSION "\n"
echo -e "Analysis started at: $(date) \n"
echo -e "Using $CORES cores \n"
echo -e "Input from: $WD_IN "
echo -e "Output to: $WD \n"

if [ -f $REFERENCE ]; then
	echo "# -------- General options --------- #"
	echo "Reference file is: $REFERENCE"
else
	echo "Reference file not found. Please provide a reference."
	usage
	exit 1
fi

if [ "$MAPPING" = 1 ]; then
	echo "Mapping reads: yes"
	if [ "$READS" = 1 ]; then
		echo "Reads: paired"
	else 
		echo "Reads: single"
	fi
	echo -e "Mapping output : $BAM_DIR"
else 
	echo "Mapping reads: no"
fi

if [ "$CALLING" = 1 ]; then
	echo "Calling snps: yes"
	echo -e "SNPs output : $VCF_DIR"
else 
	echo "Calling snps: no"
fi

if [ "$CONSENSUS" = 1 ]; then
	echo "Generating consensus: yes"
	echo -e "Consensus output : $CONSENSUS_DIR"
else 
	echo "Generating consensus: no"
fi

if [ "$ALIGNING" = 1 ]; then 
	echo -e "Generating alignments: yes"
	echo -e "Final alignments will be exported to: \n$ALIGNS_F_DIR"
else 
	echo -e "Generating alignments: no"
fi

if [ "$CONCATENATE" = 1 ]; then 
	echo -e "Concatenating: yes"
	echo -e "Concatenated alignment will be exported to: \n$CONC_DIR"
else 
	echo -e "Concatenating: no"
fi

if [ "$STATS" = 1 ]; then 
	echo -e "Generating stats: yes"
else 
	echo -e "Generating stats: no"
fi

echo -e "# ---------------------------------- # \n"

if [ "$MAPPING" = 1 ]; then
	echo -e "# ------- Mapping parameters ------- #"
	echo -e "--bwa-a : $BWA_A"
	echo -e "--bwa-b : $BWA_B"
	echo -e "--bwa-o : $BWA_O "
	echo -e "# ---------------------------------- # \n"
fi

if [ "$CONSENSUS" = 1 ]; then
	echo -e "# ------ Consensus parameters ------ #"
	echo -e "--min-depth : $DEPTH_MIN"
	echo -e "--max-depth : $DEPTH_MAX"
	echo -e "--min-q : $Q_MIN "
	echo -e "# ---------------------------------- # \n"
fi

if [ "$ALIGNING" = 1 ]; then 
	echo -e "# ------ Alignment parameters ------ #"
	echo -e "--min-seq-num : $MIN_SEQ_NUM"
	echo -e "--min-cov : $MIN_COV"
	if [ "$TRIM" = 1 ]; then
		echo -e "--trim : yes"
	else 
		echo -e "--trim : no"
	fi
	echo -e "# ---------------------------------- # \n"
fi


# ------------------------------------------------------------------ #
# Mapping
# ------------------------------------------------------------------ #

if [ "$MAPPING" = 1 ]; then
	
	### Listing files

	if [ "$(ls $WD_IN/*.fastq.gz 2>/dev/null | wc -l)" != "0" ]; then
		ls $WD_IN/*.R1.* | awk -F/ '{print $NF}' > fastq.files
	else 
		echo "Could not find the *.fastq.gz files."
		echo "See skimmingLoci --help"
		exit 1
	fi

	LABS=$(sed -e "s/.R1.fastq.gz//g" fastq.files)
	TOTAL_PP=$(ls $WD_IN/*.R1.fastq.gz | wc -l)
	COUNTER=1
	
	### Create index
	samtools faidx $REFERENCE 2>/dev/null
	bwa index $REFERENCE 2>/dev/null

	### Map 	
	echo -e "\n Mapping reads...\n"

	mkdir -p $BAM_DIR
	progress_bar 0
		
	for i in ${LABS}
		do
		if [ "$READS" = 1 ]; then
		bwa mem $REFERENCE $WD_IN/$i.R1.fastq.gz $WD_IN/$i.R2.fastq.gz "-t$CORES" -A $BWA_A -B $BWA_B -O $BWA_O -a -T 10 -M 2> $BAM_DIR/$i.bwa.log | samtools view -@$CORES -q 1 -b -F 4 -S -T $REFERENCE 2>/dev/null | samtools sort -@$CORES -o $BAM_DIR/$i.bam 2>/dev/null
		else
		bwa mem $REFERENCE $WD_IN/$i.R1.fastq.gz "-t$CORES" -A $BWA_A -B $BWA_B -O $BWA_O -a -T 10 -M 2> $BAM_DIR/$i.bwa.log | samtools view -@$CORES -q 1 -b -F 4 -S -T $REFERENCE 2>/dev/null | samtools sort -@$CORES -o $BAM_DIR/$i.bam 2>/dev/null
		fi 
		(( PP = 100* $COUNTER/$TOTAL_PP ))
		(( COUNTER = $COUNTER+1 ))
		progress_bar "$PP"
    		done
	rm fastq.files
	echo -e "\n"
fi


# ------------------------------------------------------------------ #
# SNP calling
# ------------------------------------------------------------------ #

if [ "$CALLING" = 1 ]; then
	
	### Listing files

	if [ "$(ls $BAM_DIR/*.bam 2>/dev/null | wc -l)" != "0" ]; then
		ls $BAM_DIR/*.bam | awk -F/ '{print $NF}' > bam.files
	else 
		echo "Could not find the *.bam files."
		echo "First run the mapping step."
		echo -e "See skimmingLoci --help \n"
		exit 1
	fi

	LABS=$(sed -e "s/.bam//g" bam.files)
	TOTAL_PP=$(ls $BAM_DIR/*.bam | wc -l)
	COUNTER=1

	### SNP call
	echo -e "\n SNP calling...\n"

	mkdir -p $VCF_DIR
	progress_bar 0

	for i in ${LABS}
		do
		### bcftools call 
		samtools mpileup -A -O -s -aa -E -t DP -aa -q30 -Q30 -uf $REFERENCE $BAM_DIR/$i.bam 2>/dev/null > $VCF_DIR/$i.temp.vcf
		bcftools call -c -M -Ov -o $VCF_DIR/$i.vcf $VCF_DIR/$i.temp.vcf 2>/dev/null
		rm $VCF_DIR/$i.temp.vcf
		vcftools --geno-depth --vcf $VCF_DIR/$i.vcf --out $VCF_DIR/$i 2>/dev/null
		
		### depth stats
		if [ "$STATS" = 1 ]; then
			echo "setwd('$VCF_DIR'); 
			read.table('$i.gdepth', header = T) -> dat; d <- dat[which(dat[,3] > 0),3]; median(d) -> m; sd(d) -> s.d.; length(d) -> bp; d[which(d < $DEPTH_MAX)] -> d; 
			jpeg('$i.depth.jpg', width = 400, height = 400); layout(matrix(c(1,1,1,2))); par(mar=c(3,3,3,1)); 
			hist(d, xlab='Depth', main='$i'); legend('topright', legend=paste('bp =', bp, '\ndepth median =', round(m,1), '\ndepth sd =',round(s.d.,1)), box.col = NA); 
			par(mar=c(3,3,1,1)); boxplot(d, horizontal = T, axes=F); axis(1); dev.off(); 
			c(bp,round(m,2),round(s.d.,2)) -> s; names(s) <- c('bp', 'depth median', 'depth s.d.'); t(as.data.frame(s)) -> s; rownames(s) <- '$i'; 
			write.table(s, '$i.idepth', quote=F)" > $WD/$i.R
			R CMD BATCH --no-save $WD/$i.R
			rm $WD/$i.R
			#rm $WD/$i.Rout
		fi
		
		(( PP = 100* $COUNTER/$TOTAL_PP ))
		(( COUNTER = $COUNTER+1 ))
		progress_bar "$PP"
		done
	rm bam.files
	if [ "$STATS" = 1 ]; then
		cd $VCF_DIR
		echo -e "Terminal bp_total depth_median depth_sd" > $WD/stats.samples.depth.txt
		tail -q -n 1 *.idepth >> $WD/stats.samples.depth.txt
		cd $WD
		echo "setwd('$VCF_DIR');
		list.files(pattern='.gdepth') -> files; unlist(lapply(sapply(files, read.table, header=T, simplify = F), '[', 3)) -> dat; dat[which(dat < $DEPTH_MAX)] -> dat;
		setwd('$WD');
		jpeg('stats.samples.depth.jpg'); layout(matrix(c(1,1,1,2))); par(mar=c(3,3,3,1));
		hist(dat, xlab='Depth', main='All samples'); legend('topright', legend=paste('depth median =', round(median(dat,2)), '\ndepth sd =', round(sd(dat), 2)), box.col = NA);
		par(mar=c(3,3,1,1)); boxplot(dat, horizontal = T, axes=F); axis(1); dev.off();" > $WD/all_dep.R
		R CMD BATCH --no-save $WD/all_dep.R
		rm $WD/all_dep.R
		rm $WD/all_dep.Rout
	fi
	echo -e "\n"
	
fi


# ------------------------------------------------------------------ #
# Consensus
# ------------------------------------------------------------------ #

if [ "$CONSENSUS" = 1 ]; then
	
	### Listing files

	if [ "$(ls $VCF_DIR/*.vcf 2>/dev/null | wc -l)" != "0" ]; then
		ls $VCF_DIR/*.vcf | awk -F/ '{print $NF}' > vcf.files
	else 
		echo "Could not find the *.vcf files."
		echo "See skimmingLoci --help"
		exit 1
	fi

	LABS=$(sed -e "s/.vcf//g" vcf.files)
	TOTAL_PP=$(ls $VCF_DIR/*.vcf | wc -l)
	COUNTER=1

	### Consensus	
	echo -e "\n Generating consensus...\n"

	mkdir -p $CONSENSUS_DIR
	progress_bar 0

	cat $REFERENCE | sed -e '/^>/! s/n//g' | sed -e '/^>/! s/N//g' | sed -e '/^>/! s/-//g' | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | tail -n +2 > $CONSENSUS_DIR/reference.lengths
			
	for i in ${LABS}
		do
		### vcfutils
		vcfutils.pl vcf2fq -D $DEPTH_MAX -d $DEPTH_MIN -Q $Q_MIN $VCF_DIR/$i.vcf > $CONSENSUS_DIR/$i.fastq 
		seqtk seq -a $CONSENSUS_DIR/$i.fastq > $CONSENSUS_DIR/$i.temp.fas
		sed -e '/^>/! s/[[:lower:]]/N/g' $CONSENSUS_DIR/$i.temp.fas > $CONSENSUS_DIR/$i.temp2.fas
		sed "s/>.*/&:@${i}/" $CONSENSUS_DIR/$i.temp2.fas > $CONSENSUS_DIR/$i.fas
		cat $CONSENSUS_DIR/$i.fas | sed -e '/^>/! s/n//g' | sed -e '/^>/! s/N//g' | sed -e '/^>/! s/-//g' | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | tail -n +2 > $CONSENSUS_DIR/$i.lengths
		rm $CONSENSUS_DIR/$i.temp.fas
		rm $CONSENSUS_DIR/$i.temp2.fas
		(( PP = 100* $COUNTER/$TOTAL_PP ))
		(( COUNTER = $COUNTER+1 ))
		progress_bar "$PP"
		done
	rm vcf.files
	echo -e "\n"
	
	### coverage stats
	if [ "$STATS" = 1 ]; then
		echo "setwd('$CONSENSUS_DIR'); list.files(pattern='.lengths') -> files; read.table(files[grep('reference.lengths', files)]) -> ref; files[-grep('reference.lengths', files)] -> files; 
		do.call(rbind, sapply(files, read.table, stringsAsFactors=F, simplify = F)) -> dat; 
		cbind(dat, unlist(lapply(strsplit(dat[,1], '@'), '[', 2))) -> dat; 
		cbind(dat, unlist(lapply(strsplit(dat[,1], ':'), '[', 1))) -> dat; 
		dat[which(dat[,2] > 0),-1] -> dat; unique(dat[,2]) -> spp; stats.all <- vector(); all <- vector();
		for (i in 1:length(spp)) {;
			dat[which(dat[,2] == spp[i]),] -> d0; 
			cbind(d0, ref[match(d0[,3], ref[,1]),]) -> d0;
			data.frame(locus=d0[,3], bp=d0[,1], coverage=round(d0[,1]/d0[,5],2)) -> d0;
			median(d0[,2]) -> bp.m; round(median(d0[,3]),2) -> cov.m; round(sd(d0[,3]),2) -> cov.sd;
			jpeg(paste(spp[i], '.coverage.jpg', sep=''));
			layout(matrix(c(1,1,1,2)));par(mar=c(3,3,3,1));
			hist(d0[,3], xlab='Coverage', main=spp[i]);
			legend('topright', legend=paste('bp median =', bp.m, '\ncoverage median =', cov.m, '\ncoverage sd =',cov.sd), box.col = NA);
			par(mar=c(3,3,1,1));boxplot(d0[,3], horizontal = T, axes=F); axis(1); dev.off(); 
			rbind(stats.all, c(spp[i], bp.m,cov.m,cov.sd)) -> stats.all; rbind(all, d0) -> all};
		setwd('$WD');
		jpeg('stats.samples.coverage.jpg');
		layout(matrix(c(1,1,1,2))); par(mar=c(3,3,3,1));
		hist(all[,3], xlab='Coverage', main='All samples');
		legend('topright', legend=paste('bp median =', round(median(all[,2]),2), '\ncoverage median =', round(median(all[,3]),2), '\ncoverage sd =', round(sd(all[,3]), 2)), box.col = NA);
		par(mar=c(3,3,1,1)); boxplot(all[,3], horizontal = T, axes=F); axis(1); dev.off();
		colnames(stats.all) <- c('Terminal', 'bp_locus_median', 'coverage_median', 'coverage_sd'); write.table(stats.all, file='stats.samples.coverage.txt', sep='\t', quote = F, row.names = F)" > cov.R
		R CMD BATCH --no-save cov.R
		rm $CONSENSUS_DIR/*.lengths
		rm cov.R
		rm cov.Rout
	fi
fi


# ------------------------------------------------------------------ #
# Aligning
# ------------------------------------------------------------------ #

if [ "$ALIGNING" = 1 ]; then 
	
	### Listing files
	if [ "$(ls $CONSENSUS_DIR/*.fas 2>/dev/null | wc -l)" != "0" ]; then
		ls $CONSENSUS_DIR/*.fas | awk -F/ '{print $NF}' > fasta.files
	else 
		echo "Could not find the *.fas files."
		echo "First run the mapping and snp calling steps."
		echo -e "See skimmingLoci --help \n"
		exit 1
	fi
	LABS=$(sed -e "s/.fas//g" fasta.files)
	TOTAL_PP=$(cat $REFERENCE.fai | wc -l)
	COUNTER=1

	### Listing loci
	LOCI=$(awk '{print $1}' $REFERENCE.fai)
	
	### Alignments 
	echo -e "\n Generating alignments...\n"	

	### Remove sequences with only missing data
	cat $CONSENSUS_DIR/*.fas > $CONSENSUS_DIR/cat.fasta
	sed -e '/^>/! s/N//g' $CONSENSUS_DIR/cat.fasta | awk 'NF' | awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' | grep "^>" | cut -c 1 --complement > $CONSENSUS_DIR/keep.list
	#seqtk subseq $CONSENSUS_DIR/cat.fasta $CONSENSUS_DIR/keep.list | sed -e '/^>/! s/N/-/g' > $CONSENSUS_DIR/cat.kept.fasta
	seqtk subseq $CONSENSUS_DIR/cat.fasta $CONSENSUS_DIR/keep.list > $CONSENSUS_DIR/cat.kept.fasta

	### Include reference
	sed "s/>.*/&:@reference/" $REFERENCE > $CONSENSUS_DIR/ref_temp.fas
	cat $CONSENSUS_DIR/ref_temp.fas $CONSENSUS_DIR/cat.kept.fasta > $CONSENSUS_DIR/cat.fasta
	rm $CONSENSUS_DIR/cat.kept.fasta
	rm $CONSENSUS_DIR/ref_temp.fas	

	### Split and align loci
	mkdir -p $ALIGNS_DIR
	mkdir -p $ALIGNS_F_DIR
	echo > nnn
	progress_bar 0
	for i in ${LOCI}
		do
		### Split
		grep ">$i:" $CONSENSUS_DIR/cat.fasta | cut -c 1 --complement > $ALIGNS_DIR/keep
		seqtk subseq $CONSENSUS_DIR/cat.fasta $ALIGNS_DIR/keep | cut -f2 -d ":" | sed "s/@/>/g" > $ALIGNS_DIR/$i.t0.fas
		### Filter
		cat $ALIGNS_DIR/$i.t0.fas | sed -e '/^>/! s/n//g' | sed -e '/^>/! s/N//g' | sed -e '/^>/! s/-//g' | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | tail -n +2 > xxx
		awk 'FNR==NR{max=($2+0>max)?$2:max;next} {print $1,$2, $2/max}' xxx xxx > $ALIGNS_DIR/$i.stats
		if  (( $(echo "$MIN_COV > 1" | bc -l) )); then
			awk -v n1=$MIN_COV '$2 > n1 {print $1}' $ALIGNS_DIR/$i.stats > keep.s
		else 
			awk -v n1=$MIN_COV '$3 > n1 {print $1}' $ALIGNS_DIR/$i.stats > keep.s
		fi
		seqtk subseq $ALIGNS_DIR/$i.t0.fas keep.s > $ALIGNS_DIR/$i.t1.fas
		### Align
		NSEQ=$(grep -c "^>" $ALIGNS_DIR/$i.t1.fas)
		if [ "$NSEQ" -gt "$MIN_SEQ_NUM" ]; then
			mafft --localpair --maxiterate 1000 --quiet $ALIGNS_DIR/$i.t1.fas > $ALIGNS_DIR/$i.fas
			grep ">" $ALIGNS_DIR/$i.fas | grep -v ">reference" | sed 's/>//g' > keep.f
			seqtk subseq $ALIGNS_DIR/$i.fas keep.f > $ALIGNS_F_DIR/$i.fas
			rm keep.f
			### trim and fill edges
			if [ "$TRIM" = 1 ]; then
				echo "library(skimmingLociR); library(ape); setwd('$WD'); read.dna('$ALIGNS_F_DIR/$i.fas', 'fasta') -> a0; trimAligns(a0, min.missing = 1-$MIN_SEQ_NUM/nrow(a0), quiet=T, edges.only = F) -> a0;  fillAligns(a0) -> a0; write.dna(a0, '$ALIGNS_F_DIR/$i.fas', 'fasta')" > fill.R
			else 
				echo "library(skimmingLociR); library(ape); setwd('$WD'); read.dna('$ALIGNS_F_DIR/$i.fas', 'fasta') -> a0; trimAligns(a0, min.missing = 1-$MIN_SEQ_NUM/nrow(a0), quiet=T, edges.only = T) -> a0;  fillAligns(a0) -> a0; write.dna(a0, '$ALIGNS_F_DIR/$i.fas', 'fasta')" > fill.R
			fi
			R CMD BATCH --no-save fill.R
			rm fill.R
			rm fill.Rout
		fi
		NSEQ=$((NSEQ - 1))
		echo -e "$i \t $NSEQ" >> nnn
		rm keep.s
		rm xxx
		rm $ALIGNS_DIR/keep
		rm $ALIGNS_DIR/$i.t0.fas
		rm $ALIGNS_DIR/$i.t1.fas
		(( PP = 100* $COUNTER/$TOTAL_PP ))
		(( COUNTER = $COUNTER+1 ))
		progress_bar "$PP"
		done
	tail -n +2 nnn > stats.loci.nseq.txt
	rm nnn
	rm fasta.files
	rm $CONSENSUS_DIR/cat.fasta
	rm $CONSENSUS_DIR/keep.list
	echo -e "\n"

	
fi


# ------------------------------------------------------------------ #
# Concatenate
# ------------------------------------------------------------------ #

### falta fazer um teste para ver se os alignment finals

if [ "$CONCATENATE" = 1 ]; then 
	
	mkdir -p $CONC_DIR
	### Concatenate
	echo -e "\n Concatenating alignments...\n"	
	echo "library(skimmingLociR); setwd('$WD'); fas.dir = '$ALIGNS_F_DIR'; conc.dir = '$CONC_DIR'; skimmingLociR:::skimming.concatenate(fas.dir, conc.dir, file = 'Concatenated.phy', trim.edges = T, trim.min = $MIN_SEQ_NUM)" > c.R
	R CMD BATCH --no-save c.R
	rm c.R
	rm c.Rout

	echo -e "\n"	

fi

# ------------------------------------------------------------------ #
# Stats
# ------------------------------------------------------------------ #

### falta fazer um teste para ver se tem o concatenado

if [ "$STATS" = 1 ]; then 
	
	### Stats
	echo -e "\n Generating basic statistics...\n"
	echo "library(skimmingLociR); setwd('$WD'); conc.dir = '$CONC_DIR'; reference = '$REFERENCE'; ref.dir = '$WD_IN'; skimmingLociR:::skimming.stats(conc.dir = conc.dir, reference = reference, ref.dir = ref.dir) -> dstats; 
	dstats\$loci.stats -> loci.stats; dstats\$spp.stats -> spp.stats; write.table(loci.stats, file = 'stats.loci.final.txt', sep='\t', quote=F); write.table(spp.stats, file = 'stats.samples.loci.txt', sep='\t', quote=F); 
	spp.stats[order(spp.stats\$coverage_median, decreasing = T),] -> mcov.spp; spp.stats[order(spp.stats\$loci_n, decreasing = T),] -> nloci.spp; pdf('stats.samples.loci.pdf'); par(mar=c(3,9,2,2)); 
	barplot(nloci.spp\$loci_n, horiz=T, names.arg = rownames(nloci.spp), las=2, space=2, cex.names=0.5, main='Number of loci across samples'); 
	barplot(mcov.spp\$coverage_median, horiz=T, names.arg = rownames(mcov.spp), las=2, space=2, cex.names=0.5, main='Loci median coverage across samples'); dev.off()" > stats.R
	R CMD BATCH --no-save stats.R
	rm stats.R
	rm stats.Rout

fi

echo -e "Analysis finished at: $(date) \n"
