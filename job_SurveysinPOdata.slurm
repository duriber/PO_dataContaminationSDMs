#!/bin/bash
#
#
#SBATCH --account=OD-230073
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80gb
#SBATCH --mail-user=David.Uriberivera@csiro.au
#SBATCH --mail-type=TIME_LIMIT_90

i=$1

module purge

module load gdal/3.7.2
module load proj/9.3.0
module load geos/3.12.0
module load udunits/2.2.28
module load R/4.3.1

# job submission directory;
cd /home/uri003/ReRun_SystematicRecordsinPOdata

Rscript --vanilla analysis_HPCBatch.R $i
