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
./00-scripts/01_qc_summary_table.sh
```
Output: `summaries/per_chr_qc.tsv`

Most chromosomes showed highly consistent assembly characteristics, with low N content and similar GC content (~36–37%). Chromosome 9 (LR994555.1) was selected for the simulations because it has an intermediate chromosome size (19,361,588 bp; ~19.4 Mb), contains no N bases, and does not exhibit unusually high levels of low-complexity sequence. As a chromosome of representative size and sequence composition, it provides a realistic genomic scale for modelling while remaining computationally tractable.

#### Generating the genomic coordination file

The genomic coordination file used by the SLiM model was generated from the chromosome 9 annotation. The script extracts gene and exon features from `genes.gff3` for chromosome 9 and converts them into a continuous segmentation of the chromosome into exon, intron and non_coding regions. Coordinates in the final file are 1-based inclusive.

```sh
./00-scripts/02_make_coord_file.sh
```
Output: `03-model_input/coord_chr9.txt`
Additional chromosome-specific annotation output: `03-model_input/coord_chr9.gff3`
Two additional windows from the middle of the chromosome were generated for testing and debugging the model before running simulations on the full chromsome in `03-model_input/testing/` (200 kB, 2 Mb window).




## References
Lohse, K. et al. (2023). Genome assembly of *Cyaniris semiargus* (Mazarine Blue). Wellcome Open Research / Darwin Tree of Life Project.
