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

## Conda environment
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

#### Download
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

#### Chromosome selection
To select a suitable chromosome for the SLiM simulations, we computed per-chromosome quality metrics from the reference genome, including chromosome length, GC content, N content, overall low-complexity sequence content and the maximum proportion of low-complexity sequence within any 100 kb window.

```sh
./00-scripts/qc_summary_table.sh
```
Output: `summaries/per_chr_qc.tsv`
Most chromosomes showed highly consistent assembly characteristics, with low N content and similar GC content (~36–37%). Chromosome 9 was chosen as the initial test chromosome due to its medium chromosome size (~19.3 Mb) and no other abnormalities (no unusual N content or low-complexity regions) to provide a realistic genomic scale for testing the SLiM model while keeping simulations computationally managable.


## References
Lohse, K. et al. (2023). Genome assembly of *Cyaniris semiargus* (Mazarine Blue). Wellcome Open Research / Darwin Tree of Life Project.
