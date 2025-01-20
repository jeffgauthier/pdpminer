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

	# set strict channel priority
	conda config --set channel_priority flexible

	# remove previous installation
	conda env remove -y -n depp_cli; rm -rf Depolymerase-Predictor

	# download DePP (Depolymerase-Predictor)
	git clone https://github.com/DamianJM/Depolymerase-Predictor.git

	# go to directory
	cd Depolymerase-Predictor/DePP_CLI/

	# create environment for DePP (Depolymerase-Predictor)
	conda env create -y -f ./environment.yml
	conda activate depp_cli

	# go back
	cd ../../

	# test DePP installation
	APP="Depolymerase-Predictor/DePP_CLI/depp_cli.py"
	$APP -h

	# add phage annotator (pharokka)
	conda env remove -y -n pharokka_1.7.4
	conda create -y -n pharokka_1.7.4 -c bioconda -c conda-forge pharokka=1.7.4

	# test Pharokka installation
	conda activate pharokka_1.7.4
	install_databases.py -o db
	pharokka.py -h

	# add protein domain annotator (pfamscan) and seqtk to subset candidate proteins
	conda env remove -y -n pfamscan_1.6
	conda create -y -n pfamscan_1.6 -c bioconda -c conda-forge pfam_scan=1.6 seqtk

	# test pfam_scan
	conda activate pfamscan_1.6
	pfam_scan.pl -h
	seqtk

	# install pfam databases
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
	wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz
	mkdir pfamdb
	gunzip -c Pfam-A.hmm.dat.gz > pfamdb/Pfam-A.hmm.dat
	gunzip -c Pfam-A.hmm.gz > pfamdb/Pfam-A.hmm
	gunzip -c active_site.dat.gz > pfamdb/active_site.dat
	rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz active_site.dat.gz
	hmmpress pfamdb/Pfam-A.hmm

	# go back to conda base env
	conda deactivate
	echo "INFO: Installation complete."

fi
