#!/bin/bash
# #SBATCH --array=1-324
#SBATCH --array=1-1
#SBATCH -o logs/scan_%A_%a.out
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH --overcommit
#SBATCH --ntasks-per-core=2
#SBATCH -t 119:59:00
#SBATCH --mem=4000

OFFSET=-1
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)
rundir='v_a_vs_deltaL_and_c0'

outdir="../Data/Raw/${rundir}"
mkdir -p $outdir
outfile="$outdir/out_$LINE_NUM"
setup_file="$outdir/to_run.csv"
params_file="$outdir/params.json"

echo "Line $LINE_NUM ; infile $infile ; outfile $outfile"

python ode_solve.py -o $outfile -d $params_file -n $LINE_NUM -s $setup_file
