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
# STEP 1 - CREATE CONDA ENVIRONMENT
#---------------------------------------------------
# This installation requires Miniconda v3 or newer.
# The installer will crash if it can not find the
# default Conda environment path and suggest a link
# to help the user retrieve it first.
#---------------------------------------------------

# crash if a step fails
set -euo pipefail

# load conda base environment if available, otherwise exit
#if true; then
if [ -z "$CONDA_PREFIX" ]; then

	# Conda not installed, quit
	echo "--------------------------------------------------------------------------"
	echo "ERROR: Could not find a Conda base environment!"
	echo "Miniconda3 is required for installing PDP-Miner"
	echo "and its dependencies on Linux 64-bit systems."
	echo ""
	echo "To download Miniconda3 for Linux 64-bit, run this command:"
	echo "wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
	echo ""
	echo "After installing Miniconda3, please restart your shell."
	echo "You should see '(base)' to the left of your cursor".
	echo "--------------------------------------------------------------------------"
	exit 1

else

	# load Miniconda3 base environment
	source $CONDA_PREFIX/etc/profile.d/conda.sh

	# set flexible channel priority
	conda config --set channel_priority flexible

	# create conda env for all dependencies
	conda create -y -n pdpminer_env -c bioconda -c conda-forge \
		pharokka=1.7.1 \
		pfamscan>=1.6 \
		python>=3.9 \
		biopython>=1.77 \
		numpy>=1.22 \
		pandas>=1.4 \
		scikit-learn>=1.1 \
		seqtk
	conda activate pdpminer_env

	# add Pharokka db and test installation
	install_databases.py -o $CONDA_PREFIX/pharokka_db
	pharokka.py -h

	# test pfam_scan
	pfam_scan.pl -h

	# test seqtk
	seqtk

	# install pfam databases
	PFAM_URL='http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release'
	wget $PFAM_URL/Pfam-A.hmm.dat.gz
	wget $PFAM_URL/Pfam-A.hmm.gz
	wget $PFAM_URL/active_site.dat.gz
	mkdir $CONDA_PREFIX/pfam_db
	gunzip -cv Pfam-A.hmm.dat.gz > $CONDA_PREFIX/pfam_db/Pfam-A.hmm.dat
	gunzip -cv Pfam-A.hmm.gz > $CONDA_PREFIX/pfam_db/Pfam-A.hmm
	gunzip -cv active_site.dat.gz > $CONDA_PREFIX/pfam_db/active_site.dat
	rm -v Pfam-A.hmm.gz Pfam-A.hmm.dat.gz active_site.dat.gz
	hmmpress $CONDA_PREFIX/pfam_db/Pfam-A.hmm

	# install Depolymerase-Predictor in conda env
	cd $CONDA_PREFIX
	git clone https://github.com/DamianJM/Depolymerase-Predictor.git

	# test installation
	APP="$CONDA_PREFIX/Depolymerase-Predictor/DePP_CLI/depp_cli.py"
    $APP -h
	

	# go back to conda base env
	conda deactivate
	echo "INFO: Installation complete."

fi
