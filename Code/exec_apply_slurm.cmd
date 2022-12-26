#!/bin/bash
# #SBATCH --array=1-1
#SBATCH --array=1-200
#SBATCH -o logs/scan_%A_%a.out
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH --overcommit
#SBATCH --ntasks-per-core=2
#SBATCH -t 119:59:00
#SBATCH --mem=4000

OFFSET=-1
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)
echo "$LINE_NUM"
rundir='v_a_vs_deltaL_and_c0'
E_dir = 'E_1'

outdir="../Data/Raw/${rundir}/${E_dir}"
mkdir -p $outdir
outfile="$outdir/out_$LINE_NUM"
setup_file="$outdir/to_run.tsv"
default_params_file="../Data/Raw/${rundir}/default_params.json"

echo "Line $LINE_NUM ; params_file $params_file ; outfile $outfile"
python3 ./ode_solver.py -o $outfile -d $default_params_file -n $LINE_NUM -s $setup_file
