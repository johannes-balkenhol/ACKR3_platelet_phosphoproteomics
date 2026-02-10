# ACKR3 Platelet Phosphoproteomics

Analysis scripts for the study of ACKR3 (CXCR7)-mediated signaling in human platelets using quantitative phosphoproteomics.

## Overview

This repository contains the computational pipeline for analyzing phosphoproteomic and proteomic data from ACKR3-stimulated human platelets. The study investigates time-resolved phosphorylation dynamics upon CXCR7/ACKR3 activation across two independent experimental cohorts (initial and validation), covering multiple time points (10s, 600s, 1800s) and conditions (ACKR3 agonist, DMSO vehicle, unstimulated control).

## Repository Structure

```
├── SubProjects/
│   ├── CXCR7_initial/                # Initial (discovery) dataset
│   │   ├── phosphoproteomics/
│   │   │   ├── scripts/R/            # Current phosphoproteomics pipeline
│   │   │   ├── phospho.a08.v2/       # Earlier pipeline version
│   │   │   └── scripts/perl/         # MaxQuant output parsers
│   │   └── proteomics/
│   │       └── scripts/R/            # Global proteomics analysis
│   │
│   └── CXCR7_validation/             # Validation dataset
│       ├── phosphoproteomics/
│       │   ├── scripts/R/            # Current phosphoproteomics pipeline
│       │   ├── scripts/Perl/         # MaxQuant output parsers
│       │   ├── scripts/SQL/          # Database queries for data integration
│       │   └── analysis/             # Pathway curation & network tools
│       └── proteomics/
│           └── scripts/R/            # Global proteomics analysis
│
├── scripts/                          # Cross-dataset integrated analyses
│   ├── R/                            # Comparative & reproducibility scripts
│   ├── integrated_network_analysis/  # Network diffusion (Jupyter/Python)
│   └── python/notebooks/             # Expression analysis & GSEA (Python)
│
├── share/scripts/                    # Shared/reference script copies
└── documents/supplemental_material/  # Publication-ready analysis code
```

## Analysis Pipeline

The phosphoproteomics pipeline follows a numbered step convention (`pl.phospho.N.*`):

| Step | Description |
|------|-------------|
| 0.1 | Remove modifications and miscleaved peptides |
| 1.0–1.1 | Preprocessing, QC, data integration |
| 2.0–2.1 | Empirical control site identification, filtering, normalization |
| 3.0–3.1 | Differential phosphorylation analysis (limma) |
| 4 | Volcano plots |
| 5 | Heatmaps |
| 6 | GSEA pathway enrichment |
| 7 | Reactome pathway enrichment (timeline, selected pathways, comparative) |
| 7.1 | Gene Ontology ORA |
| 8 | Signalosome & site-centric analysis |
| 9 | Network reconstruction (STRING/IntAct) |
| 10 | Network analysis (OmniPath, PSICQUIC, CausalPath, enrichment) |
| 11 | Kinase-substrate enrichment analysis (KSEA) |
| 12 | Kinase target heatmaps |
| 13 | KEGG pathway mapping |

The global proteomics pipeline covers preprocessing, normalization, differential expression, volcano plots, heatmaps, GSEA, and Reactome enrichment.

## Integrated Network Analysis

Python/Jupyter notebooks in `scripts/integrated_network_analysis/` implement network diffusion analysis using OmniPath protein interaction data to identify signaling modules downstream of ACKR3 activation.

## Languages & Key Dependencies

**R**: limma, PhosR, clusterProfiler, ReactomePA, KSEtool, OmnipathR, igraph, ComplexHeatmap, ggplot2

**Python**: pandas, networkx, omnipath, scipy

**Perl**: Custom parsers for MaxQuant PhosR/normalized output

**SQL**: Data integration queries for phosphosite tables

## Data

Raw mass spectrometry data and processed datasets are not included in this repository. Data files (`.RData`, `.csv`, `.xlsx`) are excluded via `.gitignore`.

## Authors

Johannes Balkenhol

University of Würzburg, Biozentrum

## License

This repository contains analysis code for an unpublished manuscript. Please contact the authors before reuse.
