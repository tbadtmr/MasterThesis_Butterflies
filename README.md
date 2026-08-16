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
conda env create -f environment.yml
conda activate butterflies
```
- SLiM and SLiMgui 5.1
- samtools 1.13
- seqkit v2.13.0
- dustmasker 1.0.0
- BCFtools / RoH 1.20
- PLINK 1.9
- pixy
- R

## 1. Reference genome and chromosome selection
The reference genome of *Cyaniris semiargus* was obtained from NCBI:

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
The genome can then be indexed with:
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

The genomic coordination file used by the SLiM model was generated from the chromosome 9 annotation. The script extracts gene and exon features from `genes.gff3` and converts them into a continuous segmentation of chromosome 9 consisting of exon, intron and non_coding regions. Coordinates in the final file are stored as 1-based inclusive intervals. Consecutive regions are non-overlapping and together span the entire chromosome.

```sh
./00-scripts/02_make_coord_file.sh
```
Output: `03-model_input/coord_chr9.txt`
Additional chromosome-specific annotation output: `03-model_input/coord_chr9.gff3`
Two additional windows from the middle of the chromosome were generated for testing and debugging the model before running simulations on the full chromosome in `03-model_input/testing/` (200 kB, 2 Mb window).

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

### Analysis of monitoring and occurrence sites

Habitat associations of *Cyaniris semiargus* were evaluated using data from the Swedish Butterfly Monitoring Scheme (SeBMS), CORINE Land Cover 2018 and GBIF occurrence records. These analyses were used to inform the relative habitat-suitability values later assigned to the HILDA+ land-cover classes.

The monitoring input files are located in:

`01-data/monitoring/`

The analyses are run sequentially using scripts `03a`–`03e`.

#### Monitoring-site habitat extraction

CORINE 2018 and HILDA+ 2019 habitat classes were extracted at SeBMS monitoring sites. CORINE landscape composition was additionally calculated within 0.5, 1, 2, 5 and 10 km buffers around each site.

```sh
Rscript 00-scripts/03a_extract_monitoring_landscape.R
```

Inputs:
- `01-data/monitoring/presence_absence_master.tsv`
- `01-data/monitoring/corine_index_to_clc.tsv`
- CORINE Land Cover 2018 raster
- HILDA+ 2019 state raster

Outputs:
- `02-landscapes/derived/sites_with_classes.tsv`
- `02-landscapes/derived/site_buffer_composition_long.tsv`

The analysis includes 1,429 monitoring sites, of which 646 fall within the simulation study region. The CORINE classes extracted directly from the raster showed 100% agreement with the existing site classifications.

#### Habitat associations at monitoring sites

Recorded presence was compared among broad CORINE habitat classes. Pasture/rangeland and unmanaged grass/shrubland were pooled into a single Grassland category for statistical analyses, while Cropland, Forest, Urban and Other were retained as separate classes. Logistic regression models additionally included northing to account for the geographic distribution of the species.

```sh
Rscript 00-scripts/03b_monitoring_analysis.R
```

Outputs:
- `02-landscapes/tables/TABLE_S1_habitat_occurrence.tsv`
- `02-landscapes/tables/TABLE_S2_logistic_regression.tsv`
- `02-landscapes/figures/FIG_S1_habitat_association.*`
- `02-landscapes/figures/FIG_S2_CORINE_HILDA_correspondence.*`

This analysis provides the monitoring-site habitat associations and the comparison between the fine-scale CORINE and coarser HILDA+ classifications used in the thesis.

#### Surrounding landscape composition

To test whether recorded presence was related to the landscape surrounding monitoring sites, CORINE habitat proportions were calculated relative to terrestrial area within each buffer. Buffers containing less than 50% terrestrial area were excluded.

The primary analysis used a 2 km radius, with additional analyses at 0.5, 1, 5 and 10 km to assess sensitivity to spatial scale.

```sh
Rscript 00-scripts/03c_monitoring_buffer_analysis.R
```

Outputs:
- `02-landscapes/derived/site_buffer_composition_wide.tsv`
- `02-landscapes/tables/TABLE_S3_landscape_composition_2km.tsv`
- `02-landscapes/tables/grassland_association_by_radius.tsv`
- `02-landscapes/figures/grassland_association_by_radius.*`

The final 2 km analysis contained 628 monitoring sites within the study region. The positive association with grassland cover was retained across all tested buffer radii from 0.5 to 10 km.

#### Monitoring-based abundance estimates

SeBMS transect, transect-segment and point-count data were summarized to evaluate the plausible order of magnitude of local adult abundance used during demographic model calibration.

Transect counts were standardized using the 5 m wide Pollard-walk survey belt. Point counts conducted within a 25 m radius were converted to an equivalent survey area, allowing all three survey types to be expressed as area-standardized observed adults per km².

```sh
Rscript 00-scripts/03d_density_summary.R
```

Outputs:
- `02-landscapes/derived/density_records_standardized.tsv`
- `02-landscapes/tables/density_summary_by_method.tsv`
- `02-landscapes/tables/density_summary_by_habitat.tsv`
- `02-landscapes/figures/density_by_habitat.*`

These values are based only on records in which *C. semiargus* was observed and are not corrected for imperfect detection. They are therefore used to assess the plausible order of magnitude of local abundance rather than as estimates of absolute population density.

### GBIF occurrence data

Occurrence records for *Cyaniris semiargus* were downloaded from GBIF and used as an independent assessment of habitat associations.

**Download:** GBIF.org (17 March 2026), GBIF Occurrence Download  
**DOI:** https://doi.org/10.15468/dl.9kxv2n  
**Download key:** 0043200-260226173443078  
**Records:** 11,962

The original GBIF download contained human observations from Sweden with coordinates, present occurrence status, no reported geospatial issues, coordinate uncertainty up to 1,060 m and records from 2000–2026.

Place the downloaded Darwin Core Archive contents in:

`01-data/gbif/`

For the analysis, records were further restricted to:
- coordinate uncertainty ≤1,000 m
- years 2000–2025
- the May–August flight period
- unique combinations of coordinates and observation date

To reduce repeated sampling of frequently visited locations, occurrence records were spatially thinned to a maximum of one record per 1 km grid cell.

Habitat use was compared with:
1. exact habitat availability across all terrestrial CORINE cells within the simulation study region; and
2. the distribution of habitat classes among SeBMS monitoring sites, providing a second comparison that partially accounts for observer accessibility.

```sh
Rscript 00-scripts/03e_gbif_habitat_analysis.R
```

Outputs:
- `02-landscapes/tables/TABLE_S4_gbif_habitat_selection.tsv`
- `02-landscapes/tables/TABLE_S4_gbif_habitat_selection_numeric.tsv`
- `02-landscapes/tables/gbif_record_filtering.tsv`

After filtering and spatial thinning, 339 terrestrial GBIF occurrence records remained within the simulation study region. Grassland was the only habitat class consistently overrepresented relative to both landscape availability and the distribution of SeBMS monitoring sites. Together with the monitoring analyses and published ecological information, these results were used to inform the habitat-suitability values assigned in the simulation model.


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


