# MasterThesis_Butterflies
Modelling landscape effects on gene flow and genomic erosion in grassland butterfly *Cyaniris semiargus*.

## Structure
- `00-scripts/` analysis and preprocessing scripts
- `01-data/` input data (monitoring scheme, reference genome, etc.)
- `02-landscapes/` landscape analyses and processed tables
- `03-model_input/` all processed files need for the Simulation
- `04-model_development/` SLiM model develoment steps and script versions
- `05-burnin/` final burn-in candidates and selected ancestral burn-in model
- `06-simulation/` final historical and future simulation workflows
- `07-model_comparison/` comparison of model outputs to empirical data
- `08-analysis/` downstream analysis
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


## 3. Generating landscape suitability maps

HILDA+ land-cover classes are reclassified into relative habitat-suitability values for *C. semiargus* and converted into annual simulation landscapes for 1899–2019.

| HILDA+ class | Code | Suitability |
|---|---:|---:|
| Pasture / rangeland | 33 | 1.0 |
| Unmanaged grass / shrubland | 55 | 1.0 |
| Cropland | 22 | 0.5 |
| Forest | 44 | 0.3 |
| Urban | 11 | 0.1 |
| Ocean / no data, sparse / no vegetation, water | 0, 66, 77 | 0.0 |

Because HILDA+ does not contain a separate semi-natural grassland class, pasture/rangeland and unmanaged grass/shrubland are both treated as the highest-suitability grassland habitat.

The annual HILDA+ state rasters were first reprojected to WGS84 using nearest-neighbour resampling. These preprocessed rasters are stored locally in:

`01-data/HILDA/preprocessed_wgs84/`

They are not included in the repository because of file size. The final suitability maps are projected onto the simulation template `hilda_final_clip_20km_south.tif`, which defines the 340 × 277 cell landscape (~1110 × 1105 m per cell; ~1.23 km²), and are masked to Sweden. Cells outside Sweden or without suitable land cover are assigned suitability 0.

```sh
Rscript 00-scripts/03f_make_suitability_maps.R
```

Outputs for each year are written to:

`03-model_input/suit_maps/`

- `hilda_suitability_<year>.tif` — suitability raster
- `hilda_suitability_<year>.txt` — flat row-major suitability values read directly by the SLiM model
- `hilda_suitability_<year>_dims.txt` — grid dimensions and spatial extent

Landscape summaries are written to:

`02-landscapes/derived/`

- `class_areas_by_year.tsv` — number of cells and area (km²) per HILDA+ class and year
- `suitability_value_counts.tsv` — number of cells and area (km²) per suitability value and year

The regenerated suitability maps were checked against the maps used in the original simulations and matched cell-for-cell for tested historical and contemporary years.

#### Habitat change through time

Historical changes in habitat availability are summarized from `class_areas_by_year.tsv`.

```sh
Rscript 00-scripts/03g_habitat_change_summary.R
```

The script calculates weighted suitable habitat area, total grassland area and the individual trajectories of pasture/rangeland and unmanaged grass/shrubland.

Outputs:

`02-landscapes/figures/`
- `FIG_5_habitat_change_1899_2019.png`
- `FIG_5_habitat_change_1899_2019.pdf`
- `FIG_5_habitat_change_1899_2019.svg`

`02-landscapes/tables/`
- `TABLE_S5_land_cover_change.tsv`
- `TABLE_S5_land_cover_change_numeric.tsv`
- `habitat_change_summary.tsv`

Across the study period, weighted suitable habitat declined from approximately 20,830 to 18,710 weighted km², while total grassland declined from approximately 2,950 to 785 km².

## 4. Model development

This directory documents the main stages in the development of the spatially
explicit non-Wright-Fisher model for *Cyaniris semiargus*.

The scripts are retained as model-development provenance rather than as the
production workflow used to generate the final thesis results. Intermediate
debugging and parameter-testing versions that did not introduce a substantial
new model component have been omitted.

| Stage | Script | Main development |
|------:|--------|------------------|
| 1 | `01_basic_ecology_WF.slim` | Initial annual demographic proof-of-concept using a Wright-Fisher model |
| 2 | `02_basic_nonWF_spatial.slim` | Individual-based nonWF life cycle with explicit reproduction, continuous space and dispersal |
| 3 | `03_habitat_map.slim` | Habitat-suitability raster, habitat-restricted founder placement and offspring settlement |
| 4 | `04_suitability_relK_proxy.slim` | Monitoring-derived relative-capacity layer combined with habitat suitability |
| 5 | `05_local_density_dependence.slim` | Local competition, density-dependent recruit survival, spatial mate choice and rare long-distance dispersal |
| 6 | `06_local_carrying_capacity.slim` | Explicit local carrying capacity derived from the relative-capacity landscape |
| 7 | `07_suit_based_density_dependence.slim` | Simplified suitability-only carrying capacity and suitability-based density dependence |
| 8 | `08_genetics_added_2Mb.slim` | Addition of chromosome-based genetic architecture using a reduced 2 Mb chromosome-9 test region |

### Transition to the production model

The first seven stages focus on development of the ecological and spatial
framework. Stage 8 introduces the genomic component using a reduced 2 Mb
window of chromosome 9 to test the integration of genetic and ecological
processes while keeping computational requirements manageable.

The genomic model includes neutral and deleterious mutations, annotated
genomic element types, mutation and sex-specific recombination, variable
dominance effects and genetic summary statistics.

After these development stages, the model was expanded to the full chromosome
9 and several production-scale burn-in parameterisations were run. These are
documented in [`05-burnin/`](05-burnin/).

The burn-in parameterisation that reached the predefined nucleotide-diversity
convergence criterion was subsequently used to initialise the historical and
future simulations in [`06-simulation/`](06-simulation/).

## 5. Burn-in simulations

The full chromosome 9 model was run under the constant 1899 landscape to generate
an ancestral population for the historical simulations. Four long burn-in
parameterisations were explored during final model calibration.

| Version | Main differences | Use |
|---|---|---|
| `01_burnin2_selected` | Intermediate dispersal and density scaling (`DISP_SD = 2`, `DENS_MAX = 120`) | **Selected final burn-in** |
| `02_burnin1_narrow_dispersal` | Narrower local dispersal and higher density scaling (`DISP_SD = 1`, `DENS_MAX = 150`) | Alternative |
| `03_burnin4_broad_dispersal` | Broader dispersal and lower density scaling (`DISP_SD = 4`, `DENS_MAX = 100`) | Alternative |
| `04_burnin_dens_low_capacity` | Broad dispersal with the lowest density scaling (`DISP_SD = 4`, `DENS_MAX = 80`) | Alternative |

The parameterisations were explored iteratively rather than as a full factorial
parameter search. The first candidate to reach the predefined nucleotide-diversity
convergence criterion while maintaining stable, biologically plausible demographic
behaviour was `01_burnin2_selected`. Neutral diversity converged at generation
115,000, and the final burn-in state was written at generation 117,000. This state
was used to initialise all downstream historical simulations.

The alternative scripts represent additional long burn-in parameterisations that were initiated during model calibration but could not be completed within the timeframe of the thesis. They are retained as starting points for future simulation runs.

### Running a burn-in

Scripts use repository-relative paths and should be run from the repository root.
For example, the selected burn-in can be started in the background with:

```sh
nohup slim 05-burnin/01_burnin2_selected/burnin2_selected.slim \
  > 05-burnin/01_burnin2_selected/output/burnin2_selected.log 2>&1 &

# Progress can be followed with:
tail -f 05-burnin/01_burnin2_selected/output/burnin2_selected.log
```

## 6. Historical and future simulations

The selected generation-117000 burn-in state was used to initialize the final
historical and future simulations. All production runs use the same ecological
and genetic model and the final parameterisation (`DENS_MAX = 120`,
`K_EXPONENT = 1.3`, `K_GLOBAL_CAP = 120000`).

Each simulation starts from the selected burn-in state and undergoes a
100-generation stabilisation period on the constant 1899 landscape before the
historical landscape sequence begins. Annual HILDA+ suitability maps are then
used from 1899 to 2019. The 2019 landscape is retained as the common 2020
baseline, after which the future scenarios diverge from 2021 onward.

Four future scenarios were simulated:

- `status_quo` — the 2019 landscape is retained into the future
- `restore_2km` — restoration of historically suitable grassland within 2 km of current high-suitability habitat
- `restore_4km` — restoration within 4 km
- `restore_6km` — restoration within 6 km

Replicate IDs use deterministic random seeds (`100000 + replicate number`).
The same replicate number therefore uses the same seed across scenarios,
providing paired trajectories that are identical before the future scenarios
diverge in 2021.

### Simulation workflows

The final simulations are organised into three complementary workflows:

| Directory | Spatial extent | Replicates per scenario | Purpose |
|---|---|---:|---|
| `01-full-area-snapshots/` | Full study area | 50 | Main production simulations with complete population snapshots at selected historical, contemporary and future time points |
| `02-skane-summaries/` | Skåne | 50 | Computationally reduced simulations used to increase replication for demographic and genetic summary trajectories |
| `03-skane-structure/` | Skåne | 30 | Targeted full-population snapshots for downstream spatial sampling and population-structure analyses |

The full-area simulations represent the primary spatial model. Because complete
full-area simulations and population snapshots are computationally expensive,
additional simulations restricted to Skåne were subsequently used to increase
replication for analyses focused on the region containing the historical and
contemporary genomic samples. Restricting the spatial extent does not alter the
underlying demographic or genetic parameterisation.

The Skåne summary runs omit large full-population SLiM snapshots while retaining
the demographic and genetic summary calculations used to compare trajectories
among future scenarios.

The Skåne structure runs were further optimised for analyses requiring individual
genotypes. Repeated nucleotide-diversity, FROH and genetic-load calculations are
disabled in these runs. For the status-quo scenario, full population states are
written for 1900, 2020 and 2140. Because the historical trajectory is shared
among paired scenarios, restoration runs only require their scenario-specific
2140 state. These simulations terminate after the 2140 snapshot has been written.

### Full-area population snapshots

The main full-area simulations write complete SLiM population states at:

- 1900, 1951, 1956, 2020, 2140, 2260, 2380, 2500

The historical 1951 and 1956 states correspond to the time periods used for
comparison with historical genomic samples.

Some full-area simulations exceeded the HPC wall-time. The companion
`simulation_resume_from_snapshot.slim` and `run_resume_dardel.sh` scripts allow
an interrupted trajectory to continue from its most recent full-population
snapshot without duplicating output from the loaded generation.

### Running the simulations

The complete burn-in population state is not included in the repository because
of its size. Before running the simulations, paths to the selected generation-
117000 SLiM state and corresponding deleterious-mutation (`m2`) file must be
provided.

For example, a full-area status-quo simulation can be submitted on Dardel with:

```sh
export BURNIN_SNAPSHOT=/path/to/burnin_chr9_fulloutput_END_gen117000.slim
export BURNIN_M2_FILE=/path/to/burnin_chr9_m2_gen117000.txt

sbatch \
  --export=ALL,SCENARIO=status_quo \
  06-simulation/01-full-area-snapshots/run_dardel.sh
```
The other future scenarios are submitted by replacing `status quo` with `restore_2km`, `restore_4km` or `restore_6km`.
The Skåne summary simulations can similarly be submitted on LUNARC with:
```sh
sbatch \
  --export=ALL,SCENARIO=status_quo,BURNIN_SNAPSHOT=/path/to/burnin_chr9_fulloutput_END_gen117000.slim,BURNIN_M2_FILE=/path/to/burnin_chr9_m2_gen117000.txt \
  06-simulation/02-skane-summaries/run_lunarc.sh
```
and the targeted structure runs with:
```sh
sbatch \
  --export=ALL,SCENARIO=status_quo,BURNIN_SNAPSHOT=/path/to/burnin_chr9_fulloutput_END_gen117000.slim,BURNIN_M2_FILE=/path/to/burnin_chr9_m2_gen117000.txt \
  06-simulation/03-skane-structure/run_lunarc.sh
```
The provided SLURM launch scripts retain the cluster resources used during the
thesis. Account, partition and SLiM executable settings may need to be adjusted
when running on another HPC system.

## 7. Comparison with empirical genomic data

Simulation output was compared with historical and contemporary whole-genome
data from Skåne reported by Nolen et al. (2024). The aim was to evaluate whether
the model reproduced observed temporal and spatial genomic change under a
sampling design comparable to the empirical dataset.

The comparison used five independent full-area status-quo simulation replicates.
Simulated individuals were sampled around the corresponding empirical
localities using the same sample sizes as the empirical data:

- Western Skåne: 1951, n = 5; 2020, n = 8
- Eastern Skåne: 1956, n = 5; 2020, n = 8
- Southeastern Skåne: 2020, n = 7

Individuals within approximately 4 km of each empirical locality were eligible
for sampling. For each population and simulation replicate, 100 random samples
were drawn using the corresponding empirical sample size. Sampling was without
replacement within each draw, although individuals could occur in multiple
draws. Southeastern Skåne had no historical empirical sample and was therefore
compared with historical eastern Skåne.

The comparison included nucleotide diversity (π), individual heterozygosity
(H), runs of homozygosity (FROH) and pairwise genetic differentiation (FST).
Statistics from the 100 random samples were averaged within each simulation
replicate and replicate means were then summarized across the five independent
replicates.

### Analysis steps

The analysis scripts are located in `00-scripts/`, with matching SLURM launchers
for computationally intensive steps in `00-scripts/slurm/`. Script numbers
indicate the order of the analysis.

`04_export_positions.slim`  
Exports individual IDs and coordinates from the selected SLiM population
snapshots.

`05_check_sampling_pools.R`  
Identifies individuals within the sampling radius around each empirical
locality and checks whether the required sample size is available.

`06_make_validation_sampling_manifest.R`  
Generates 100 random samples for each locality and simulation replicate using
the corresponding empirical sample size.

`07_make_genotype_requests.R`  
Determines which simulated individuals are required for genotype export.

`08_vcf_from_validation_request.slim`  
Reloads the corresponding SLiM population snapshot and exports genotypes of the
requested individuals as VCF files.

`09_calculate_validation_pi.R` and `10_summarize_validation_pi.R`  
Calculate and summarize nucleotide diversity.

`11_calculate_validation_heterozygosity.R` and
`12_summarize_validation_heterozygosity.R`  
Calculate and summarize individual heterozygosity.

`13_parse_validation_froh.R` and `14_summarize_validation_froh.R`  
Parse BCFtools/RoH output and summarize FROH using runs longer than 100 kb.

`15_compare_froh_rates.R`  
Compares the tested recombination-rate settings for the ROH analysis.

`16_calculate_validation_fst.R` and `17_summarize_validation_fst.R`  
Calculate and summarize pairwise Hudson's FST.

### Running the analysis

The scripts are intended to be run sequentially because later steps use files
generated by earlier steps. SLURM launchers use the same number as their
corresponding analysis step.

For example, the initial steps can be submitted from the repository root with:

```sh
sbatch 00-scripts/slurm/04_export_positions.sh
sbatch 00-scripts/slurm/05_check_sampling_pools.sh
sbatch 00-scripts/slurm/06_make_validation_sampling_manifest.sh
sbatch 00-scripts/slurm/07_make_genotype_requests.sh
sbatch 00-scripts/slurm/08_export_validation_vcfs.sh
```
The subsequent metric-specific launchers follow the same numbering convention.
Large SLiM population states are not included in the repository, so
`SIM_ROOT` can be set to the location containing the required full-area
simulation snapshots when these steps are rerun on another system.

#### Outputs
Final validation summaries are stored in:
- `results/pi/`
- `results/heterozygosity/`
- `results/froh/`
- `results/fst/`
Small sampling, task-control and quality-control tables are stored directly in
`07-model_comparison/results/`. Large intermediate files, including exported
positions, complete sampling manifests, genotype requests, VCF files and HPC
logs, are reproducible from the supplied scripts and are therefore not tracked
by Git.

## 8. Downstream population structure analysis

Population structure was analysed using fixed regional sampling designs applied to the simulated population snapshots.

Two spatial designs were used:
- **Full landscape:** 7 regions × 3 fixed sites = 21 sampling sites (`08-population-analysis/full/design/final_21_regional_core_sites.tsv`)
- **Skåne-only:** 5 regions × 3 fixed sites = 15 sampling sites (`08-population-analysis/skane_only/design/final_15_regional_core_sites.tsv`)

The full design contains the five Skåne regions (NW, NE, W, E and SE), western Småland and Öland. The latter two are retained as `NEW02` and `NEW01`, respectively, in the internal design file.

At each site, up to 20 individuals within a radius of 4 model units were sampled. This corresponds to up to 60 individuals per region. Where fewer individuals were available, the smaller sample was retained while preserving the fixed-site design.

### Analysis steps

`18_export_population_positions.slim`  
Extracts individual simulation indices and spatial coordinates from the selected SLiM population snapshots.

`19_make_population_sampling_manifest.py`  
Assigns individuals to the fixed sampling sites and selects up to 20 individuals per site using deterministic SHA256 ranking. The resulting sampling manifest records the exact individuals selected for each site and snapshot.

`20_export_population_manifest_vcf.slim`  
Reloads each SLiM population snapshot and exports the individuals listed in the sampling manifest as VCF files.

`21_build_joint_population_vcf.py`  
Constructs a joint neutral VCF for each replicate using neutral variants shared across the compared time points and future scenarios.

`21_run_joint_population_pca_dardel.sh`  
Converts the joint VCF to PLINK format, performs LD pruning using

`--indep-pairwise 50 5 0.2`

and calculates 10 principal components.

For the future population-structure comparison, the joint PCA contains six states within each replicate:

- 1900
- 2020
- status quo 2140
- 2-km restoration 2140
- 4-km restoration 2140
- 6-km restoration 2140

Independent simulation replicates are analysed separately rather than pooled into a single PCA.

The position and genotype-export steps have matching SLURM launchers for Dardel and LUNARC in `00-scripts/slurm/`. The same scripts are used for the full and Skåne-only analyses by supplying the corresponding working directory, sampling design and task inventory.

### Pop structure




### Regional diversity and differentiation

`23_calculate_regional_pi_fst.py`  
Calculates regional nucleotide diversity (π) from neutral variants and pairwise Hudson's FST separately for each simulation replicate.

`24_summarize_regional_fst.R`  
Summarizes regional FST into the three comparisons reported in the thesis: within Skåne, Skåne–W Småland and Skåne–Öland. Pairwise values are first averaged within each simulation replicate and then summarized across replicates.

The regional π and FST analysis can be submitted on Dardel using:

```sh
sbatch 00-scripts/slurm/23_calculate_regional_pi_fst_dardel.sh
```
The grouped FST summary is then generated with:
```sh
Rscript 00-scripts/24_summarize_regional_fst.R
```
Final full-landscape results are stored in: `08-population-analysis/full/results/regional_genetics/`.


## References
Lohse, K. et al. (2023). Genome assembly of *Cyaniris semiargus* (Mazarine Blue). Wellcome Open Research / Darwin Tree of Life Project.

Winkler, K., Fuchs, R., Rounsevell, M.D.A., & Herold, M. (2021). Global land use changes are four times greater than previously estimated. Nature Communications, 12, 2501. https://doi.org/10.1038/s41467-021-22702-2

Winkler, K., Fuchs, R., Rounsevell, M.D.A., & Herold, M. (2020). HILDA+ Global Land Use Change between 1960 and 2019 (vGLOB-1.0). PANGAEA. https://doi.org/10.1594/PANGAEA.921846


