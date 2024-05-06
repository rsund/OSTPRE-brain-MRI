#!/bin/bash
# Jussi Tohka University of Eastern Finland, for runs in the UEF's "sampo" cluster
#SBATCH --ntasks 1              # Number of taskS
#SBATCH --time 72:00:00         # Runtime
#SBATCH --mem-per-cpu=10000      # Reserve 10 GB RAM for a task, this is the memory required for a single subject
#SBATCH --partition serial      # Partition to submit
#SBATCH --job-name synthseg_create
#SBATCH --output /research/users/justoh/slurm_logs/ostpre_synthseg_logs/ostpre_synth_create1998_%A_%a.err    # Standard output goes to here
#SBATCH --error  /research/users/justoh/slurm_logs/ostpre_synthseg_logs/ostpre_synth_create1998_%A_%a.out    # Standard error goes to here
#SBATCH --array=1-920%1     # Array range and number of simultanous jobs, note that these indexes refer to
                               # lines of the subjidses.txt file in the DATADIR, important to run only 1 job at a time as 
                               # this is for file generation 
                               # subjid.txt should be kept unchanged and contain all the subject and session indexes 
                               # note that now we need also sessions for better parallelization
                               # this allows for keeping track what subjects have been already processed based on log-files
                               # ,i.e., a log file always has an index corresponding to the line number in the subjid.txt file.
                               

FSDIR=/research/users/justoh/scripts/OSTPRE_sbatch
DATADIR=/research/groups/ostpre/ostpre_bids/bids
ROWINDEX=$((SLURM_ARRAY_TASK_ID+1998)) # to bypass MaxArraySize limit
SUBJSES=$(sed -n "$ROWINDEX"p ${DATADIR}/subjidses.txt)

# module load freesurfer/7.4.1 # load modules

bash synthseg_batch_bids_ostpre_createfiles.sh "${DATADIR}/${SUBJSES}"

