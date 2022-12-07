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

OFFSET=1
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)
rundir='ode_species_density'

outdir="../Data/Raw/${rundir}"
line=$(sed -n "$LINE_NUM"p $outdir/to_run.csv)

echo "Line $LINE_NUM ; infile $infile ; outfile $outfile"


matlab -r "load('$outdir/params.mat'); \
           params.log10c0 = $log10c0; \
           params.P = [$p1; 1-$p1]; \
           params.E = $E'; \
           params.log10delta = $log10delta'; \
           disp('An adaptors serial dilution simulation is starting...'); \
           output = serial_adaptors_odesolver(params,'$outfile'); \
           disp('Finished!'); \
           quit();"