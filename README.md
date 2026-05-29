# AMR Detection Pipeline in Elizabethkingia anophelis

Bioinformatics pipeline developed in Python for the identification and analysis of antimicrobial resistance (AMR) genes in public genomes of Elizabethkingia anophelis.

---

# Project Overview

This project aims to identify antimicrobial resistance genes using genome analysis, gene prediction, and AMR screening tools.

The pipeline performs:

- Genome parsing from FASTA files
- Basic genome statistics
- GC content analysis
- Gene prediction using Prodigal
- AMR gene detection using AMRFinderPlus
- Generation of summary outputs and visualizations

---

# Tools and Libraries

## Python Libraries

- pandas
- Biopython
- matplotlib

## Bioinformatics Tools

- Prodigal
- AMRFinderPlus

---

# Project Structure

```bash
project/
│
├── data/
│   └── Elizabethkingia_Anophelis_Ag1.fna
│
├── results/
│   ├── proteins.faa
│   ├── genes.fna
│   ├── genes.gff
│   └── amrfinder_results.tsv
│
├── proyecto_bioinformatico_mario_florez.py
│
└── README.md
```

---

# Pipeline Workflow

1. Load genome FASTA file
2. Compute genome statistics
3. Analyze nucleotide composition
4. Predict coding sequences with Prodigal
5. Generate protein and gene files
6. Detect AMR genes using AMRFinderPlus
7. Export results in TSV format

---

# Installation

## Create Conda Environment

```bash
conda create -n amr_env python=3.10
conda activate amr_env
```

## Install Python Dependencies

```bash
pip install biopython pandas matplotlib
```

## Install Prodigal

```bash
conda install -c bioconda prodigal
```

## Install AMRFinderPlus

```bash
conda install -c bioconda ncbi-amrfinderplus
```

---

# Running the Pipeline

```bash
python proyecto_bioinformatico_mario_florez.py
```

---

# Output Files

The pipeline generates:

- proteins.faa
- genes.fna
- genes.gff
- amrfinder_results.tsv

It also generates plots for:

- Genome nucleotide composition
- Candidate AMR gene length distribution

---

# Example AMR Genes Detected

- blaB
- blaGOB
- catB
- aadS

---

# Organism Analyzed

Elizabethkingia anophelis Ag1

---

# Author

Mario Florez

---

# License

MIT License
