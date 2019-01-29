#!/bin/bash
#SBATCH -c 16
#SBATCH -N 1
#SBATCH -p general
#SBATCH -t 7-00:00
#SBATCH --mail-user=maximiliankasy@fas.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -J CalibratedSimulations
#SBATCH -o CalibratedSimulations.out
#SBATCH -e CalibratedSimulations.err
#SBATCH --mem=80000


echo "Print the current working directory..."
pwd

echo "Load R.."
module load R #Load R module

echo "Run R file..."

R CMD BATCH --quiet --no-restore --no-save CalibratedSimulationsServer.R CalibratedSimulationsServer.out


