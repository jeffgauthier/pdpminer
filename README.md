# PDP-Miner: 
A wrapper for DePolymerase-Predictor for searching phage tail depolymerases directly from genomic FASTA sequences. Annotates phage tail proteins first (with Pharokka), then runs DePP on this subset. Pfam domains are then annotated to corroborate DePP predictions for each gene candidate.

Implemented in Bash, supported on native Linux or WSL and supports submitting subtasks to a SLURM workload queue. This software is free and open source under the GNU General Public License v3.0.

---

# Installing PDP-Miner

## Requirements
Requires Miniconda3 to install main software dependencies (Pharokka, DePP and PfamScan). Miniconda3 can be found here: https://docs.anaconda.com/miniconda/install/. 

Tested on both Ubuntu Linux v24.04 (via WSL) and OpenSUSE Linux Leap v15.5.

## Installation
Clone this repository then enter directory. Otherwise download the ZIP archive on GitHub and extract its contents in a new folder. Then, run the script `install.sh` to install all dependencies and databases required. At the moment this install script creates separate conda environments for PDP-Miner's main dependencies (Pharokka, DePP and PfamScan) to avoid dependency conflicts which remain to be resolved. This will be fixed in the future.

---

# Running PDP-Miner

## Usage

```
./pdpminer.sh -i PATH/TO/GENOME/FASTA
              -c Number of threads per Pharokka job
              -p minimum probability threshold (default 0.5)
              -m max number of concurrent jobs per step (default 5)
              -s DO NOT use SLURM grid
```


