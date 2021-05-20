# cell-line-signatures
Code base used to generate copy number signatures from pan-cancer cell line SNP6.0 array data
## Reproducing cell line data
### Setup
Clone this git repository and `cd` into the directory
```
git clone https://github.com/Phil9S/cell-line-signatures.git
cd cell-line-signatures
```
Install the conda environment
```
conda create --name cell-line-signatures --file conda_env.txt
conda activate cell-line-signatures
```
Run the setup bash script to install remote resources and reference files
```
./scripts/setup_dir.sh
```
### Download external data
Edit the line specifying the access credentials for the restricted access dataset [EGAD00010000644](https://ega-archive.org/datasets/EGAD00010000644). Remove the line `EGA_ACCESS="/home/smith10/ega_access"` and replace the filepath with a path for the file containing your own ega access credentials. See the details on setting an credentials file in the ega-download-client github repo [documentation](https://github.com/EGA-archive/ega-download-client#defining-credentials).
```
./scripts/download_cel.sh
```
### Running Affytools
_slurm-specific_
```
sbatch sbatch_affy_run
```
or
```
./script/affy_SNP6_pipeline.sh
```
### Running modified ASCAT pipeline
_slurm-specific (implemented as a SLURM array)_
```
sbatch sbatch_ascat_array_fixed_purity
```
This utilises the row number (after header) of the `cellLine_CEL_file_mapping.tsv` file to submit a number of array jobs. This submission works by sumbitting the array job index as a commandline variable to the ASCAT.sc script, where `Rscript scripts/ASCAT.on.celllines_fixed_purity.R ${SLURM_ARRAY_TASK_ID}` is run for each sample. The slurm job submission will need to be edited to match the number of rows (excluding header) in this file by editing the line `#SBATCH --array=1-2230%200`. 

To run this section without slurm, please see the documentation on array jobs for your job management software to provide the correct index value to the `Rscript scripts/ASCAT.on.celllines_fixed_purity.R ${SLURM_ARRAY_TASK_ID}`  or implement a parallelised looping script. 
#### PBS
PBS-torque can be used by setting up a PBS job submission script matching the job requirements specified in the slurm job submission script and setting the array configuration line to `#PBS -t 1-2230%200` and replacing the `${SLURM_ARRAY_TASK_ID}` variable with `${PBS_ARRAYID}`.
#### LSF
LSF can be used by setting up a LSF job submission script matching the job requirements specified in the slurm job submission script and setting the array configuration line to `#BSUB -J arrayJob[1-2230:200]` and replacing the `${SLURM_ARRAY_TASK_ID}` variable with `${LSB_JOBINDEX}`.
#### Local
The ASCAT.sc scatter-gather array can be run without a job management software. The slow solution is to implement the ASCAT.sc array job as a for loop;
```sh
for i in `seq 1 2230`; do
  Rscript scripts/ASCAT.on.celllines_fixed_purity.R $i;
done
```
A faster implementation, if xargs is available, is to run jobs in parallel though you should monitor resource usage carefully as this will start multiple R instances;
```sh
seq 1 2230 | xargs -n 1 -P 100 -I {} bash -c 'scripts/ASCAT.on.celllines_fixed_purity.R {}'
```
Where `-P N` is the number of concurrent jobs.
### Quality control of cell line fits
#### Generate QC data
```
Rscript scripts/get_ploidy_purity.R
Rscript scripts/get_sc_ploidy_purity.R
```
or 
_slurm-specific_
```
sbatch sbatch_get_ploidy_purity
sbatch sbatch_get_sc_ploidy_purity
```
#### Generate unrounded copy number segments
```
Rscript scripts/generate_unrounded_sc_fixedPurity_tCN_file.R
```
#### Plot QC data
```
Rscript scripts/plot_summarised_segData.R
```
_slurm-specific_
```
sbatch sbatch_plot_summarised_segData
```
### Perform fit selection and filtering

Inspect the copy number profile plots and statistics generated in the previous steps and update the `cell_fit_qc_table.tsv` file to include or exclude samples using a `TRUE` or `FALSE` boolean in the `use` column.

### Re-generate unrounded copy number segments

Regenerate the unrounded segment file which now uses the `cell_fit_qc_table.tsv` to remove unwanted samples.
```
Rscript scripts/generate_unrounded_sc_fixedPurity_tCN_file.R
```
### Generate copy number signatures
_slurm-specific_
```
sbatch sbatch_genPanCanSigs
```
or
```
Rscript scripts/0_Signature_activities_from_CN_profiles.R
```
### Finalise copy number signatures
```
Rscript scripts/finalise_signatures.R
```

The final output is a `RDS object` called `cellLine_signature_data.rds` which contains the copy nubmer segment data(`copy_number`), threshold corrected signatures (`signatures`), raw signature exposures (`signatures.raw`), and z-score normalised and threshold corrected signatures (`signatures.zscore`). 
