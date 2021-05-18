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
conda env create -f conda.yaml
conda activate cell-line-signatures
```
Run the setup bash script to install remote resources and reference files
```
./scripts/setup_dir.sh
```
### Download external data
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
This utilises the row number (after header) of the `cellLine_CEL_file_mapping.tsv` file to submit a number of array jobs.
The slurm job submission will need to be edited to match the number of rows (excluding header) in this file.
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
### Perform fit selection and filtering

Inspect the copy number profile plots and statistics generated in the previous steps and update the `cell_fit_qc_table.tsv` file to include or exclude samples using a `TRUE` or `FALSE` boolean in the use column.

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
