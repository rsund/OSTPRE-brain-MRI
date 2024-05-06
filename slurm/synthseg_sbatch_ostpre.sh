#!/bin/bash
# Jussi Tohka University of Eastern Finland, for runs in the UEF's "sampo" cluster
#SBATCH --ntasks 8              # Number of taskS
#SBATCH --nodes 1   #single cpu
#SBATCH --time 72:00:00         # Runtime
#SBATCH --mem-per-cpu=25000      # Reserve 10 GB RAM for a task, this is the memory required for a single subject
#SBATCH --partition serial      # Partition to submit
#SBATCH --job-name ss_ostpre
#SBATCH --output /research/users/justoh/slurm_logs/ostpre_synthseg_logs/ostpre_synth_k_%A.err    # Standard output goes to here
#SBATCH --error  /research/users/justoh/slurm_logs/ostpre_synthseg_logs/ostpre_synth_k_%A.out    # Standard error goes to here

FSDIR=/research/users/justoh/scripts/OSTPRE_sbatch
DATADIR=/research/groups/ostpre/ostpre_bids/bids

t1filename="${DATADIR}/synthseg_t1file.txt"
segfilename="${DATADIR}/synthseg_segfile.txt"
csvfilename="${DATADIR}/synthseg_csvfile.txt"
qcfilename="${DATADIR}/synthseg_qcfile.txt"
resamplefilename="${DATADIR}/synthseg_resamplefile.txt"

module load freesurfer/7.4.1 # load modules

mri_synthseg --i ${t1filename} --o ${segfilename} --vol ${csvfilename}  --qc ${qcfilename} --resample ${resamplefilename} --parc --threads 8 --cpu
