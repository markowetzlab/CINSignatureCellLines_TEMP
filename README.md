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
_implemented as a SLURM array_
```
sbatch sbatch_ascat_array_fixed_purity
```
### Generate unrounded copy number segments
```
Rscript scripts/generate_unrounded_sc_fixedPurity_tCN_file.R
```
### Quality control of cell line fits
```
./scripts/get_ploidy_purity.R
./scripts/get_sc_ploidy_purity.R
```
then
```
Rscript scripts/plot_summarised_segData.R
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
