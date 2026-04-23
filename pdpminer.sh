#!/bin/bash

####################################################
#--------- P D P - M I N E R --- v0.1.1 -----------#
#--------- (first GitHub implementation) ----------#
# Phage tail depolymerase mining tool in 3 steps:
# - Phage gene annotation with PHAROKKA
# - Depolymerase prediction with DEPP_CLI
# - Protein domain annotation with PFAMSCAN
#
# By Jeff Gauthier - Levesque Lab, Laval U.
####################################################

#---------------------------------------------------
# STEP 0a - LOADING ENVIRONMENT AND DEPENDENCIES
#---------------------------------------------------

# load conda; crash if uninstalled
if [ -z "$CONDA_PREFIX" ]; then
        echo "--------------------------------------------------"
        echo "ERROR: Miniconda3 environment not found!"
        echo "--------------------------------------------------"
	    echo "Download and install Miniconda3 for Linux;"
        echo "Then, run install.sh prior to running pdpminer.sh"
        echo "--------------------------------------------------"
        exit 1
else
	# load default conda base environment where dependencies
	# should be installed (see snippet below).
	source $CONDA_PREFIX/etc/profile.d/conda.sh
fi

# check that pdpminer conda env is activated
if [[ $CONDA_PREFIX == "*pdpminer_env*" ]; then
	echo "--------------------------------------------------"
	echo "ERROR: Environment not activated"
	echo "--------------------------------------------------"
	echo "Please follow install.sh instructions before" 
	echo "running pdpminer.sh"
	echo "--------------------------------------------------"
	exit 1
fi


# logo
echo "

░▒▓███████▓▒░░▒▓███████▓▒░░▒▓███████▓▒░░▒▓██████████████▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓████████▓▒░▒▓███████▓▒░  
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓██████▓▒░ ░▒▓███████▓▒░  
░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓█▓▒░      ░▒▓███████▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░ 
                                                                                                          
                                                                                                          
"


##############################################################################################


#---------------------------------------------------
# STEP 0b - ARGUMENT PARSER AND GREETER
# Here I use GETOPTS to make a nice and clean greeter
# and argument parser. I also made a heredoc containing
# usage info to be shown if help is called or if
# an error occurs.
#---------------------------------------------------

# usage information
function usage() {
echo -e "
####################################################
#--------- P D P - M I N E R --- v0.1.1 -----------#
  Phage tail depolymerase mining tool in 3 steps:
  - Phage gene annotation with PHAROKKA
  - Depolymerase prediction with DEPP_CLI
  - Protein domain annotation with PFAMSCAN

  By Jeff Gauthier - Levesque Lab, Laval U.
####################################################

USAGE
./pdpminer.sh -i PATH/TO/GENOME/FASTA
              -c Number of threads per Pharokka job
              -p minimum probability threshold (default 0.5)
              -m max number of concurrent jobs per step (default 5)
              -s Use SLURM grid
"
}

# SLURM GRID COMMANDS (default unset unless -s is provided)
unset SCMD
unset SCMD2

# option string for getopts
OPTSTRING="si:c:p:m:"

# the argument parser is here
while getopts ${OPTSTRING} opt; do

	case ${opt} in

		i)
			FASTA_DIR=$(realpath -e ${OPTARG})
			echo "Input genome FASTA folder:"
			echo $FASTA_DIR
			echo ""
			;;

		p)
			export MINPROB=${OPTARG}
			echo "Minimum probability score: $MINPROB"
			echo ""
			;;

		c)
			CPUS=${OPTARG}
			echo "Number of CPUs per Pharokka job: $CPUS"
			echo ""
			;;

		s)
			SCMD="srun -c $CPUS --mem 20G"
			SCMD2="srun -c 1 --mem 10G"
			echo "Using SLURM grid!"
			echo ""
			;;

		m)
			LIMIT=${OPTARG}
			echo "Max number of concurrent jobs per step: $LIMIT"
			echo ""
			;;

		:)
			echo "Option -${OPTARG} requires an argument."
			usage
			exit 1
			;;

		?)
			echo "Invalid option: -${OPTARG}."
			usage
			exit 1
			;;

	esac

done


# check for missing mandatory args
if [ -z $FASTA_DIR ] || [ -z $CPUS ]; then
	usage
	echo ""
        echo "--------------------------------------------------"
        echo "$(date) -- ERROR: MISSING ARGUMENTS -i or -c"
        echo "--------------------------------------------------"
	echo ""
	exit 1
fi

# set defaults for unset optional args

if [ -z $MINPROB ]; then
	export MINPROB=0.5
	echo "Minimum probability score set to default (0.5)"
	echo ""
fi

if [ -z $LIMIT ]; then
        export LIMIT=5
        echo "Max number of concurrent jobs set to default (5)"
        echo ""
fi


# #######################################################################################


# #---------------------------------------------------
# # STEP 1 - ANNOTATE PHAGE GENES
# # This step produces Prokka-like annotations for each
# # genome FASTA in the user input directory.
# # Annotations are returned in the pharokka_out folder,
# # with each individual genome annotations in separate subfolders.
# # Logfiles are kept for each genome FASTA annotation.
# #---------------------------------------------------

# output directory is same as input dir plus suffix
export OUTDIR="$(basename $FASTA_DIR)_out"

# remove previous output; make output directories
rm -rf $OUTDIR
mkdir $OUTDIR
mkdir $OUTDIR/temp
mkdir $OUTDIR/temp/pharokka_out
mkdir $OUTDIR/temp/pharokka_fix
mkdir $OUTDIR/temp/depp

# run pharokka separately on each genome FASTA with SLURM
COUNTER=1; for FASTA in $FASTA_DIR/*.fna; do

	# public announcement
	echo "$(date) -- Annotating genome $FASTA with Pharokka..."

	# the fasta file's name
	PREFIX=$(basename $FASTA .fna)

	# the actual pharokka command
	$SCMD pharokka.py \
		-i $FASTA \
		-o $OUTDIR/temp/pharokka_out/$PREFIX \
		-p $PREFIX \
		-l "CDS" \
		-d db \
		-t $CPUS -g prodigal --fast \
		> $OUTDIR/temp/pharokka_out/$PREFIX\_pharokka.log 2>&1 && echo "$(date) -- $FASTA annotation done!" &

	# limit nb simultaneous iterations
	#---------------------------------------------------------------------------------
	if (( "$COUNTER" >= "$LIMIT" )) ; then

		# wait for jobs to finish before continuing loop
		echo "$(date) -- Waiting for other jobs to finish..." && wait && COUNTER=1

	else
		# otherwise increment counter by one
        	COUNTER=$((COUNTER+1))

	fi
	#---------------------------------------------------------------------------------

done
wait && unset COUNTER

###################################################################################


#---------------------------------------------------
# STEP 2 - COMPUTE DEPP PROBABILITY SCORES
# This step uses Deoplymerase-Predictor to assign
# probability values to each protein FASTA
# generated with Pharokka.
# We used Pharokka indeed to avoid feeding non-phage
# proteins uselessly to DePP.
#---------------------------------------------------

# quick fix to accomodate DePP (it does not accept X's in amino acid sequences)
# CAUTION! If locus tags are fixed, they will mismatch those from DePP...
echo "$(date) -- Fixing missing amino acids (X -> M) to accomodate DePP..."
for FAA in $OUTDIR/temp/pharokka_out/*/prodigal.faa; do
	cp $FAA $OUTDIR/temp/pharokka_fix/$(basename $(dirname $FAA))_pharokka.faa
	sed -i 's/X/M/g' $OUTDIR/temp/pharokka_fix/$(basename $(dirname $FAA))_pharokka.faa
done
wait

# run DePP-CLI (now loop is throttled by the user)
COUNTER=1; for FAA in $OUTDIR/temp/pharokka_fix/*_pharokka.faa; do

	# public announcement
	echo "$(date) -- Finding candidate PDPs in $FAA ..."

	# the actual command
	DEPP_DIR="$CONDA_PREFIX/Depolymerase-Predictor/DePP_CLI"
	$SCMD2 $DEPP_DIR/depp_cli.py \
		-i $FAA \
		-t $DEPP_DIR/TrainingSet/TrainingSet.csv \
		-o $OUTDIR/temp/depp/$(basename $FAA _pharokka.faa)_depp.csv \
		> $OUTDIR/temp/depp/$(basename $FAA _pharokka.faa)_depp.log 2>&1 && echo "$(date) -- Done for $FAA" &

        # limit nb simultaneous iterations
        #---------------------------------------------------------------------------------
        if (( "$COUNTER" >= "$LIMIT" )) ; then

                # wait for jobs to finish before continuing loop
                echo "$(date) -- Waiting for other jobs to finish..." && wait && COUNTER=1

        else
                # otherwise increment counter by one
                COUNTER=$((COUNTER+1))

        fi
        #---------------------------------------------------------------------------------

done
wait && unset COUNTER


###############################################################################################


#---------------------------------------------------
# STEP 3 - MERGE PHAROKKA AND DEPP ANNOTATIONS
# DePP probability scores are merged to the Pharokka
# gene tables. These tables will subsequently be
# fitered out to select "tail proteins" sensu Pharokka.
#---------------------------------------------------

# merge probabilities with pharokka annotations
echo "$(date) -- Merging probabilities with Pharokka annotations..."

# function for merging
function mergeprobs() {

	# function argument
	DEPP=$1

	# capture genome name
	NAME=$(basename $DEPP _depp.csv)

	# swap columns and order by highest to lowest probability
	awk -F ',' '{ print $2, $1 }' $DEPP | sort -r > $OUTDIR/temp/depp/$NAME\_depp_swap.csv

	# add gene annotation to table
	echo -e "prob\tgene\tlocus\tstart\tstop\tframe\tcontig\tscore\tmmseqs_phrog\tmmseqs_alnScore\tmmseqs_seqIdentity\tmmseqs_eVal\tmmseqs_top_hit\tpyhmmer_phrog\tpyhmmer_bitscore\tpyhmmer_evalue\tcustom_hmm_id\tcustom_hmm_bitscore\tcustom_hmm_evalue\tphrog\tMethod\tRegion\tcolor\tannot\tcategory\tvfdb_hit\tvfdb_alnScore\tvfdb_seqIdentity\tvfdb_eVal\tvfdb_species\tvfdb_short_name\tvfdb_description\tCARD_hit\tCARD_alnScore\tCARD_seqIdentity\tCARD_eVal\tCARD_species\tARO_Accession\tCARD_short_name\tProtein_Accession\tDNA_Accession\tAMR_Gene_Family\tDrug_Class\tResistance_Mechanism\ttransl_table" > $OUTDIR/temp/depp/$NAME\_depp_annotated.csv
	while read LINE; do

		# add locus tag and DePP probability
		GENE=$(cut -d ' ' -f2 <<< $LINE)
		PROB=$(cut -d ' ' -f1 <<< $LINE)

		# this is the Pharokka annotation for the corresponding gene
		ANNOT=$(grep "$GENE" $OUTDIR/temp/pharokka_out/$NAME/$NAME\_cds_final_merged_output.tsv)

		# if the line is not empty + contains a TAIL PROTEIN annotation + probability above threshold;
		# then write merged line to final output file
		# (I had to use BC because bash has no builtin float arithmetic)
		if [[ "$ANNOT" =~ "tail protein" ]] && (( $(echo "$PROB > $MINPROB" | bc) )) ; then
			echo -e "$PROB\t$GENE\t$ANNOT" >> $OUTDIR/temp/depp/$NAME\_depp_annotated.csv && \
				echo "$(date) -- *** Found $GENE in $NAME with $PROB probability"
		fi

	done <<< $(tail -n +2 $OUTDIR/temp/depp/$NAME\_depp_swap.csv)

}

# run function for each DePP score table (now throttled by the user)
export -f mergeprobs
COUNTER=1; for CSV in $OUTDIR/temp/depp/*_depp.csv; do

	# public announcement
	echo "$(date) -- Producing final output for $(basename $CSV _depp.csv) ..."

	# launch merge function
	$SCMD2 bash -c "mergeprobs $CSV" && echo "$(date) -- Final output file done for $(basename $CSV _depp.csv)" &

	# limit nb simultaneous iterations
        #---------------------------------------------------------------------------------
        if (( "$COUNTER" >= "$LIMIT" )) ; then

                # wait for jobs to finish before continuing loop
                echo "$(date) -- Waiting for other jobs to finish..." && wait && COUNTER=1

        else
                # otherwise increment counter by one
                COUNTER=$((COUNTER+1))

        fi
        #---------------------------------------------------------------------------------

done
wait && unset COUNTER

#######################################################################################################

#---------------------------------------------------
# STEP 4 - RUN PFAMSCAN ON CANDIDATE PDPs
# pfam_scan.pl version 1.6 is run on a multiFASTA file
# containing all selected PDPs found in the pre-final results table.
# 1) for each pre-final file, list CDS found;
# 2) subsample them with seqtk;
# 3) run pfam_scan.pl on this particular file (not all the FAA, that's too long!)
#---------------------------------------------------

# public announcement
echo "$(date) -- Scanning protein domains within detected candidate PDPs..."

# make temp output dir for pfamscan and pre-requisites
mkdir $OUTDIR/temp/pfamscan

# main code block, throttled by counter
COUNTER=1; for DEPP in $OUTDIR/temp/depp/*_depp_annotated.csv; do

	# capture base name
	NAME=$(basename $DEPP _depp_annotated.csv)

	# list ORFs in DEPP results table:
	# load file without first line, keep the 2nd field (remove duplicates if any) and save to list
	tail -n +2 $DEPP | cut -d $'\t' -f2 | uniq > $OUTDIR/temp/pfamscan/$NAME\_proteins.txt

	# subset proteins of interest with seqtk, store to file
	echo "$(date) -- *** Subsetting candidate PDPs from all-protein Pharokka file for $NAME ..."
	seqtk subseq \
		$OUTDIR/temp/pharokka_fix/$NAME\_pharokka.faa \
		$OUTDIR/temp/pfamscan/$NAME\_proteins.txt \
		> $OUTDIR/temp/pfamscan/$NAME\_proteins.faa

	# run pfamscan on proteins of interest only ("-as" = predict active sites)
	echo "$(date) -- *** Running PFAM_SCAN.PL against candidate PDPs from $NAME ..."
	$SCMD pfam_scan.pl -as \
		-fasta $OUTDIR/temp/pfamscan/$NAME\_proteins.faa \
		-dir $CONDA_PREFIX/pfam_db/ \
		-cpu $CPUS \
		-outfile $OUTDIR/temp/pfamscan/$NAME\_proteins.pfamscan.txt && \
		echo "$(date) -- *** DONE PFAM_SCAN.PL against candidate PDPs from $NAME" &

        # limit nb simultaneous iterations
        #---------------------------------------------------------------------------------
        if (( "$COUNTER" >= "$LIMIT" )) ; then

                # wait for jobs to finish before continuing loop
                echo "$(date) -- Waiting for other jobs to finish..." && wait && COUNTER=1

        else
                # otherwise increment counter by one
                COUNTER=$((COUNTER+1))

        fi
        #---------------------------------------------------------------------------------

done
wait && unset COUNTER



#######################################################################################################

#---------------------------------------------------
# STEP 5 - MERGING PFAM ANNOTATIONS WITH FINAL OUTPUT
# pfam_scan.pl version 1.6 is run on a multiFASTA file
# containing all selected PDPs found in the pre-final results table.
#---------------------------------------------------

# main code block, throttled by counter
COUNTER=1; for DEPP in $OUTDIR/temp/depp/*_depp_annotated.csv; do

	# capture base name
	NAME=$(basename $DEPP _depp_annotated.csv)

	# final output header
	echo -e "domains\tprob\tgene\tlocus\tstart\tstop\tframe\tcontig\tscore\tmmseqs_phrog\tmmseqs_alnScore\tmmseqs_seqIdentity\tmmseqs_eVal\tmmseqs_top_hit\tpyhmmer_phrog\tpyhmmer_bitscore\tpyhmmer_evalue\tcustom_hmm_id\tcustom_hmm_bitscore\tcustom_hmm_evalue\tphrog\tMethod\tRegion\tcolor\tannot\tcategory\tvfdb_hit\tvfdb_alnScore\tvfdb_seqIdentity\tvfdb_eVal\tvfdb_species\tvfdb_short_name\tvfdb_description\tCARD_hit\tCARD_alnScore\tCARD_seqIdentity\tCARD_eVal\tCARD_species\tARO_Accession\tCARD_short_name\tProtein_Accession\tDNA_Accession\tAMR_Gene_Family\tDrug_Class\tResistance_Mechanism\ttransl_table" > $OUTDIR/$NAME\_final_output.csv

	# scan thru lines; match gene names to domains found by pfamscan
	while read LINE; do

		# capture gene ID
		GENE=$(cut -d $'\t' -f2 <<< $LINE)

		# capture all domain names found for GENE, store them in a separate file
		grep $GENE $OUTDIR/temp/pfamscan/$NAME\_proteins.pfamscan.txt >$OUTDIR/temp/pfamscan/$NAME\_proteins.pfamscan.$GENE.txt

		# store domains as (PFNUMBER__NAME), SEPARATE WITH SEMICOLON, then remove all whitespace including newline, THEN REMOVE TRAILING SEMICOLON
		DOMAINS=$(grep -oP "PF.*?Domain|PF.*?Family" $OUTDIR/temp/pfamscan/$NAME\_proteins.pfamscan.$GENE.txt | sed -E 's/ *Domain.*| *Family.*/;/g' | sed 's/  */:/g' | tr -d '[:space:]' | sed -E 's/\;$//')

		# if no domains then write none
		if [ -z "$DOMAINS" ]; then
			DOMAINS="none_found"
		else
			echo "$(date) -- FOUND domains $DOMAINS in gene $GENE from $NAME !!"
		fi

		# write lines to final output
		echo -e "$DOMAINS\t$LINE" >> $OUTDIR/$NAME\_final_output.csv

	done <<< $(tail -n +2 $DEPP)


done
wait && unset COUNTER

# Pipeline is over
echo "$(date) -- Back to shell"
