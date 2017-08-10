#!/bin/bash
#SBATCH --job-name=NRKp41
#SBATCH -o NRKp41.out
#SBATCH -p normal 
#SBATCH -t 24:00:00

matlab -r get_base
