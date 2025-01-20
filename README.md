# PDP-Miner: an AI/ML tool to detect prophage tail proteins with depolymerase domains across thousands of bacterial genomes
A wrapper for DePolymerase-Predictor that allows searching phage or prophage tail depolymerases directly from whole-genome FASTA sequences. Annotates phage tail proteins first (with Pharokka), then runs DePP on this subset. Pfam domains are then annotated to corroborate DePP predictions for each gene candidate.

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

NOTE: When running PDP-Miner from a desktop computer, use option -s to avoid using the SLURM job scheduler when not available.

PDP-Miner automatically activates the installed environments as it needs them. Therefore, remain in `(base)` while running PDP-Miner.

```
./pdpminer.sh -i PATH/TO/GENOME/FASTA
              -c Number of threads per Pharokka job
              -p minimum probability threshold (default 0.5)
              -m max number of concurrent jobs per step (default 5)
              -s DO NOT use SLURM grid
```

Therefore, to analyse the test data provided in this repo on a desktop computer, the following need to be done to launch PDP-Miner with 4 threads per step, a minimum DePP score threshold of 0.8 and a max number of concurrent jobs of 2, without SLURM scheduling (which would be uncommon outside from a compute server):

```
git clone jeffgauthier/pdpminer.git
cd pdpminer
bash install.sh
bash pdpminer.sh -i testdata -c 4 -p 0.8 -m 2 -s
```

---

# Bugs/issues
I consider this software to still be in the Alpha stage. Though functional, it requires several improvements to ensure its reliability and usefulness. Please report any issues/bugs or quality-of-life requests in the Issues tab or by email: jeff.gauthier.1@ulaval.ca

---

# Acknowledgements
Thanks to Magill and Skvortsov for developing Depolymerase-Predictor (https://github.com/DamianJM/Depolymerase-Predictor), which served as the main engine of this genome mining workflow.

---

# How to cite
The manuscript supporting this work is currently under review; a preprint will be added shortly. Nevertheless this work can either be cited with this repo's URL, or via this temporary citation:

* Gauthier J, Kukavica-Ibrulj I, Levesque RC. (2025) PDP-Miner: an AI/ML tool to detect prophage tail proteins with depolymerase domains across thousands of bacterial genomes (in prep.)

