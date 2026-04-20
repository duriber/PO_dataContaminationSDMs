Code and data needed to replicate the workflow described in the manuscript "Mixing survey-derived and opportunistic records in presence-only 
data undermines model fit and predictive performance of species distribution models".

It is designed to be partially implemented on a unix cluster machine.
Repo contains both R and slurm code to do so, and the required structure to be run.
Results are stored in the '/out' folder.

After running the R script 'analysis_HPCBatch.R' as a batch job using 'batch_surveysinPOdata.slurm', the workflow continues by merging the 
output files into compilatory csvs with either predictions and performance metrics, or into coefficient estimates. 
To merge individual outputs, use 'merging_csvs.R' script in your local machine.

Finally, to recreate the figures, the R script 'Figures.R' has been provided.

Brisbane, Australia,
20 April 2026.
