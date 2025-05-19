#!/bin/bash
#
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=48G
#SBATCH --job-name=MT_NFA2
#SBATCH --account=aces
#SBATCH --partition=aces
#
# End of embedded SBATCH options
#

# Define directories
BASE_DIR=~/mtme-gsi4wheat
SCRIPT_DIR=$BASE_DIR
JOBS_DIR=$BASE_DIR/Job
OUTPUT_DIR=$BASE_DIR/Out

# Job name
JOB_NAME="MTME_NFA2"

{
  echo "$JOB_NAME.sh started on $(hostname) at $(date)"
  echo ""
  echo "========== $JOB_NAME.R =========="
  cat $SCRIPT_DIR/$JOB_NAME.R
} | mail -s "$JOB_NAME Started" lucasb4@illinois.edu

# Run R script
module purge
module load R/4.4.2
R CMD BATCH $SCRIPT_DIR/$JOB_NAME.R $OUTPUT_DIR/$JOB_NAME.out

{
  echo "$JOB_NAME.sh finished on $(hostname) at $(date)"
  echo ""
  echo "========== $JOB_NAME.out =========="
  cat $OUTPUT_DIR/$JOB_NAME.out
} | mail -s "$JOB_NAME Finished" lucasb4@illinois.edu