# MasterThesis_Butterflies
Modelling landscape effects on gene flow and genomic erosion in grassland butterfly *Cyaniris semiargus*

## Structure
- `00-scripts/` analysis and preprocessing scripts
- `01-data/` input data (monitoring scheme, reference genome, etc.)
- `02-landscapes/` landscape analyses and processed tables
- `03-model_input/` all processed files need for the Simulation
- `04-model_development/` SLiM model develoment steps and script versions
- `05-burnin/` final burnin script and simulation output
- `06-simulation/` simulation scripts and output (historical and future simulation)
- `07-analysis/` comparison of model outputs to empirical data and downstream analysis

## Conda environment
```sh
conda activate butterflies
```
- samtools 1.13
- SLiM and SLiMgui 5.1


