import os

NUM_LINES = 2
ARRAY_START = 0
ARRAY_END = 999

for offset in range(0, NUM_LINES, ARRAY_END + 1):
    os.system(f'sbatch --array={ARRAY_START}-{ARRAY_END} exec_apply_slurm.cmd {offset}')
