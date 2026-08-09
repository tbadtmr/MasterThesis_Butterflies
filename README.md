# MasterThesis_Butterflies
Modelling landscape effects on gene flow and genomic erosion in grassland butterfly *Cyaniris semiargus*.

## Structure
- `00-scripts/` analysis and preprocessing scripts
- `01-data/` input data (monitoring scheme, reference genome, etc.)
- `02-landscapes/` landscape analyses and processed tables
- `03-model_input/` all processed files need for the Simulation
- `04-model_development/` SLiM model develoment steps and script versions
- `05-burnin/` final burnin script and simulation output
- `06-simulation/` simulation scripts and output (historical and future simulation)
- `07-analysis/` comparison of model outputs to empirical data and downstream analysis
- `logs`
- `plots`
- `summaries`

## Conda environment

(*add yml file)
```sh
conda activate butterflies
```
- samtools 1.13
- seqkit v2.13.0
- dustmasker 1.0.0
- SLiM and SLiMgui 5.1

## 1. Reference genome and chromosome selection
he reference genome of *Cyaniris semiargus* was obtained from NCBI:

- Assembly: **GCA_905187585.1 (ilCyaSemi1.1)**
- Source: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_905187585.1/
The genome FASTA file is not included in this repository due to file size limitations.

### Download
Download assembly GCA_905187585.1 from NCBI and place the archive in `01-data/`.

Optional (using the NCBI Datasets CLI):
```sh
cd 01-data
datasets download genome accession GCA_905187585.1
```
Unzip the dataset:
```sh
unzip ncbi_dataset.zip
```
To index the genome:
```sh
samtools faidx GCA_905187585.1_ilCyaSemi1.1_genomic.fna
```

### Chromosome selection
To select a suitable chromosome for the SLiM simulations, we computed per-chromosome quality metrics from the reference genome, including chromosome length, GC content, N content, overall low-complexity sequence content and the maximum proportion of low-complexity sequence within any 100 kb window.

```sh
./00-scripts/01_qc_summary_table.sh
```
Output: `summaries/per_chr_qc.tsv`

Most chromosomes showed highly consistent assembly characteristics, with low N content and similar GC content (~36–37%). Chromosome 9 (LR994555.1) was selected for the simulations because it has an intermediate chromosome size (19,361,588 bp; ~19.4 Mb), contains no N bases, and does not exhibit unusually high levels of low-complexity sequence. As a chromosome of representative size and sequence composition, it provides a realistic genomic scale for modelling while remaining computationally tractable.

### Generating the genomic coordination file

The genomic coordination file used by the SLiM model was generated from the chromosome 9 annotation. The script extracts gene and exon features from `genes.gff3` and converts them into a continuous segmentation of chromosome 9 consisting of exon, intron, and non_coding regions. Coordinates in the final file are stored as 1-based inclusive intervals. Consecutive regions are non-overlapping and together span the entire chromosome.

```sh
./00-scripts/02_make_coord_file.sh
```
Output: `03-model_input/coord_chr9.txt`
Additional chromosome-specific annotation output: `03-model_input/coord_chr9.gff3`
Two additional windows from the middle of the chromosome were generated for testing and debugging the model before running simulations on the full chromsome in `03-model_input/testing/` (200 kB, 2 Mb window).

## 2. Landscape analysis

### Landscape data

#### HILDA+ Land Use Data
This project uses the HILDA+ Global Land Use Change dataset (extended version; 1899–2019) as the spatial basis for both the simulation and landscape analyses (Winkler et al., 2021; Winkler et al., 2020).

Download `hildap_vGLOB-1.0_geotiff_eckert4.zip` from:
https://doi.pangaea.de/10.1594/PANGAEA.921846?format=html#download
Extract all `.tif` files to:
`01-data/HILDA/`

#### CORINE Land Cover 2018
We also use the CORINE Land Cover 2018 raster dataset (100 m resolution) to identify habitat classes and calibrate habitat suitability maps used in the landscape analyses.

Download the CORINE Land Cover 2018 raster dataset (100 m) from:
https://land.copernicus.eu/en/products/corine-land-cover/clc2018
Dataset DOI:
https://doi.org/10.2909/960998c1-1870-4e82-8051-6485205ebbac
Extract the downloaded raster files to:
`01-data/CORINE/`

### Analysis of monitoring and occurrences sites

MONITORING SCHEME
presence/absence site analysis
density analysis
forest / buffer / edge analysis

GBIF


### Generating landscape suitability maps

HILDA+ land-cover classes are reclassified into habitat suitability values for _C. semiargus_ and written as annual maps for 1899–2019. Suitability ranges from 0 (unsuitable) to 1 (primary breeding habitat):

| HILDA+ class | Code | Suitability |
|---|---|---|
| Pasture / rangeland | 33 | 1.0 |
| Unmanaged grass / shrubland | 55 | 1.0 |
| Cropland | 22 | 0.5 |
| Forest | 44 | 0.3 |
| Urban | 11 | 0.1 |
| Ocean / no data, sparse / no vegetation, water | 0, 66, 77 | 0.0 |

Since there is no semi-natural grassland class available, classes 33 and 55 are both treated as primary breeding habitat. Rasters are reprojected onto hilda_final_clip_20km_south.tif (340 × 277 cells, ~1110 × 1105 m, cell area 1.23 km²) and masked to Sweden; cells outside the mask become 0.

```sh
Rscript 00-scripts/R-make_suitability_maps.R
````
Outputs, per year, to `03-model_input/suit_maps/`:
- `hilda_suitability_<year>.tif` — suitability raster
- `hilda_suitability_<year>.txt` — flat row-major values, read directly by the SLiM model
- `hilda_suitability_<year>_dims.txt` — grid geometry

and to summaries/:
- `class_areas_by_year.tsv` — cells and km² per land-cover class per year
- `suitability_value_counts.tsv` — cells and km² per suitability value per year

**Habitat change through time**
Summarises the maps into the habitat-loss figure and the numbers quoted in the results.

```sh
Rscript 00-scripts/R-habitat_change_summary.R
```
Reads `summaries/class_areas_by_year.tsv`. Outputs to `plots/` as .png, .pdf and .svg: (A) weighted suitable habitat area and (B) grassland area split into pasture / rangeland and unmanaged grass / shrubland. Also prints the 1899 and 2019 values, the percentage changes, the year of minimum grassland area, and forest's share of the weighted total.

## References
Lohse, K. et al. (2023). Genome assembly of *Cyaniris semiargus* (Mazarine Blue). Wellcome Open Research / Darwin Tree of Life Project.

Winkler, K., Fuchs, R., Rounsevell, M.D.A., & Herold, M. (2021). Global land use changes are four times greater than previously estimated. Nature Communications, 12, 2501. https://doi.org/10.1038/s41467-021-22702-2

Winkler, K., Fuchs, R., Rounsevell, M.D.A., & Herold, M. (2020). HILDA+ Global Land Use Change between 1960 and 2019 (vGLOB-1.0). PANGAEA. https://doi.org/10.1594/PANGAEA.921846


