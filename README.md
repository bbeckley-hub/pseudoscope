```markdown
██████╗ ███████╗███████╗██╗   ██╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔══██╗██╔════╝██╔════╝██║   ██║██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
██████╔╝███████╗█████╗  ██║   ██║██║  ██║██║   ██║███████╗██║     ██║   ██║██████╔╝█████╗  
██╔═══╝ ╚════██║██╔══╝  ██║   ██║██║  ██║██║   ██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
██║     ███████║███████╗╚██████╔╝██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝██║     ███████╗
╚═╝     ╚══════╝╚══════╝ ╚═════╝ ╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
```

<div align="center">

# 🧬 PseudoScope

### **A species‑optimized computational pipeline for rapid, accessible *Pseudomonas aeruginosa* genomic typing and surveillance**

#### **Complete MLST, serotyping, AMR, virulence & visualisation in minutes — not hours**

[![Version](https://img.shields.io/badge/version-1.0.0-blue)](https://github.com/bbeckley-hub/pseudoscope)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://anaconda.org/bbeckley-hub/pseudoscope)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/pseudoscope)](https://hub.docker.com/r/bbeckleyhub/pseudoscope)
[![Docker Image Size](https://img.shields.io/docker/image-size/bbeckleyhub/pseudoscope/latest)](https://hub.docker.com/r/bbeckleyhub/pseudoscope)
[![Conda Downloads](https://img.shields.io/conda/dn/bbeckley-hub/pseudoscope?label=Conda%20Downloads)](https://anaconda.org/bbeckley-hub/pseudoscope)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/pseudoscope)](https://github.com/bbeckley-hub/pseudoscope/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/pseudoscope)](https://github.com/bbeckley-hub/pseudoscope/issues)
[![Made for Research](https://img.shields.io/badge/made%20for-Research-0066cc.svg)](https://github.com/bbeckley-hub)
[![Sponsor](https://img.shields.io/badge/Sponsor-GitHub%20Sponsors-pink.svg)](https://github.com/sponsors/bbeckley-hub)
[![DOI](https://img.shields.io/badge/)]()

</div>

---

## 📋 **Table of Contents**

- [🎯 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [⚡ Quick Start](#-quick-start)
- [🔧 Installation](#-installation)
- [🐳 Docker & Singularity (HPC Ready)](#-docker--singularity-hpc-ready)
- [🚀 Usage Guide](#-usage-guide)
- [📁 Output Structure](#-output-structure)
- [🔍 Core Analytical Modules](#-core-analytical-modules)
- [📈 Ultimate Reporter & Visualisation](#-ultimate-reporter--visualisation)
- [🔬 Validation & Accuracy](#-validation--accuracy)
- [🖥️ Sample Report](#-sample-report)
- [⚙️ Mandatory First‑Run Steps](#️-mandatory-first‑run-steps)
- [❓ Not Happy? Try These Alternatives](#-not-happy-try-these-alternatives)
- [🌐 Part of the ESKAPE Web Tools Suite](#-part-of-the-eskape-web-tools-suite)
- [🔮 Future Development](#-future-development)
- [🙏 Acknowledgements](#-acknowledgements)
- [💰 Sponsor This Project](#-sponsor-this-project)
- [📚 Citation](#-citation)
- [👥 Authors & Contact](#-authors--contact)
- [📄 License](#-license)

---

## 🎯 **Overview**

**PseudoScope** is an automated, locally‑executable computational pipeline designed specifically for comprehensive *Pseudomonas aeruginosa* genomic surveillance. It integrates **seven essential analysis modules** into a single, cohesive workflow:

- **FASTA QC** – Assembly quality metrics (N50, GC%, contig stats)
- **MLST (Oxford scheme)** – Multi‑locus sequence typing via 7 housekeeping genes
- **PAST serotyping** – O‑antigen typing (pasty / camlhmp‑blast‑regions)
- **AMRfinderPlus** – Comprehensive antimicrobial resistance gene detection
- **ABRicate** – Multi‑database screening (resistance, virulence, plasmids, biocides)
- **Ultimate Reporter** – Gene‑centric integration with interactive HTML
- **Visualisation Dashboard** – Publication‑ready interactive plots (PCA, networks, boxplots)

PseudoScope runs entirely **locally** (or on HPC clusters), protects your data privacy, and produces beautiful interactive reports in **minutes**.

---

## ✨ **Key Features**

### 🔬 **P. aeruginosa‑Specific Innovations**

- **Automated critical gene flagging** – Carbapenemases (KPC, NDM, VIM, IMP, OXA), colistin (mcr‑1 to mcr‑10), ESBLs (PER, VEB, BEL), 16S rRNA methyltransferases (armA, rmt)
- **MLST (Oxford scheme)** – Uses `acsA, aroE, guaA, mutL, nuoD, ppsA, trpE`
- **PAST serotyping** – Up‑to‑date O‑antigen reference database (Robert Petit / camlhmp)
- **Cross‑genome pattern discovery** – ST‑O combination table, gene co‑occurrence, PCA
- **Gene‑centric views** – Each gene table shows **all genomes** that carry it (no truncation)

### 🚀 **Performance Advantages**

- **Complete analysis** in **10‑15 minutes** for 5 genomes (16 cores)
- **Linear scaling** – 100 genomes in ~2 hours
- **Low memory footprint** – Runs on 4GB RAM, scales to HPC
- **Parallel execution** – QC, MLST, and PAST run concurrently; ABRicate and AMR sequentially

---

## ⚡ **Quick Start**

```bash
# 1. Install PseudoScope
conda create -n pseudoscope -c conda-forge -c bioconda -c bbeckley-hub pseudoscope -y
conda activate pseudoscope

# 2. Mandatory first‑time database setup
abricate --setupdb
pseudoscope --update-amr-db

# 3. Run analysis on all .fna files
pseudoscope -i "*.fna" -o results --threads 4

# 4. Open the interactive ultimate report
firefox results/GENIUS_PSEUDOMONAS_ULTIMATE_REPORTS/genius_pseudomonas_ultimate_report.html
```

---

## 🔧 **Installation**

### **System Requirements**

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU cores | 2 | 8+ |
| RAM | 4 GB | 8 GB |
| Storage | 10 GB | 20 GB |
| OS | Linux, macOS, WSL2 | Ubuntu 22.04+ |

### **Conda Installation (Recommended)**

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels bbeckley-hub

conda create -n pseudoscope python=3.9 pseudoscope -y
conda activate pseudoscope
pseudoscope --help
```

### **Mamba (Faster)**

```bash
mamba create -n pseudoscope -c conda-forge -c bioconda -c bbeckley-hub pseudoscope -y
mamba activate pseudoscope
```

### **From Source**

```bash
git clone https://github.com/bbeckley-hub/pseudoscope.git
cd pseudoscope
conda env create -f environment.yml
conda activate pseudoscope
pip install -e .
```

---

## 🐳 **Docker & Singularity (HPC Ready)**

### Docker

```bash
docker pull bbeckleyhub/pseudoscope:latest

docker run --rm \
  -v $(pwd)/genomes:/data/input \
  -v $(pwd)/results:/data/output \
  bbeckleyhub/pseudoscope:latest \
  -i "/data/input/*.fna" -o /data/output --threads 4
```

### Singularity / Apptainer (for HPC clusters, no sudo)

**Option A – Build directly from Docker Hub:**

```bash
singularity pull pseudoscope.sif docker://bbeckleyhub/pseudoscope:latest
```

**Option B – Convert a local Docker image (if network issues on HPC):**

On a machine with Docker:
```bash
docker pull bbeckleyhub/pseudoscope:latest
docker save bbeckleyhub/pseudoscope:latest -o pseudoscope.tar
singularity build pseudoscope.sif docker-archive://pseudoscope.tar
scp pseudoscope.sif user@hpc.cluster:~
```

**Run on HPC (crucial: `--writable-tmpfs` allows temporary writes inside the container):**

```bash
singularity run --writable-tmpfs -B $(pwd):/data pseudoscope.sif \
  -i "/data/*.fna" -o /data/output --threads 8
```

All output files will be owned by **your HPC user** – no `sudo chown` needed.

---

## 🚀 **Usage Guide**

### Basic Commands

```bash
# Single genome
pseudoscope -i genome.fna -o results/

# Batch processing with glob pattern
pseudoscope -i "*.fna" -o batch_results --threads 16

# Skip specific modules
pseudoscope -i "*.fasta" -o results --skip-qc --skip-viz

# Update AMR database only
pseudoscope --update-amr-db

# Show version
pseudoscope --version
```

### All Command‑line Flags

| Flag | Description |
|------|-------------|
| `-i`, `--input` | Input FASTA file(s) – supports glob patterns (e.g., `"*.fna"`) |
| `-o`, `--output` | Output directory for all results |
| `-t`, `--threads` | Number of CPU threads (default: 2) |
| `--quiet` | Suppress all non‑error output |
| `--verbose` | Show full command output from modules |
| `--version` | Show version and exit |
| `--update-amr-db` | Update AMRfinderPlus database and exit |
| `--skip-qc` | Skip FASTA QC |
| `--skip-mlst` | Skip MLST |
| `--skip-past` | Skip PAST serotyping |
| `--skip-abricate` | Skip ABRicate |
| `--skip-amr` | Skip AMRfinderPlus |
| `--skip-summary` | Skip ultimate reporter |
| `--skip-viz` | Skip visualisation dashboard |

**Supported FASTA formats:** `.fna`, `.fasta`, `.fa`, `.fn`

---

## 📁 **Output Structure**

```
results/
├── fasta_qc_results/               # FASTA QC per sample + summary
├── mlst_results/                   # MLST per sample + mlst_summary.html/csv/json
├── past_results/                   # PAST serotyping per sample + past_summary.html/csv/json
├── pseudo_abricate_results/        # ABRicate per database summary reports
├── pseudo_amrfinder_results/       # AMRfinder per sample + summary reports
├── GENIUS_PSEUDOMONAS_ULTIMATE_REPORTS/
│   ├── genius_pseudomonas_ultimate_report.html   ← Main interactive report
│   ├── genius_pseudomonas_ultimate_report.json
│   ├── sample_overview.csv
│   ├── amr_genes.csv
│   ├── virulence_genes.csv
│   └── gene_cooccurrence.csv
└── GENIUS_PSEUDOMONAS_VISUAL_DASHBOARD/
    └── genius_pseudomonas_visual_dashboard.html   ← Interactive visualisation dashboard
```

---

## 🔍 **Core Analytical Modules**

| Module | Purpose | Key Databases / Methods | Outputs |
|--------|---------|------------------------|---------|
| **FASTA QC** | Assembly quality control | Biopython | N50, GC%, contigs, length, HTML/TSV/JSON |
| **MLST** | Oxford scheme ST assignment | PubMLST (Jolley et al.) | ST, 7‑gene profile, summary table |
| **PAST** | O‑serotyping | pasty / camlhmp (Petit & Read) | O‑type (O1–O20), coverage, hits |
| **AMRfinderPlus** | AMR gene detection | NCBI curated database | Risk categories (Critical, High), per‑gene prevalence |
| **ABRicate** | Multi‑database screening | NCBI, CARD, ResFinder, VFDB, ARG‑ANNOT, PlasmidFinder, MEGARes, EcoH, BacMet2, Ecoli_VF | Per‑database summaries, master JSON |
| **Ultimate Reporter** | Gene‑centric integration | All above HTML summaries | Interactive HTML (search/sort/filter), CSV exports |
| **Visualisation Dashboard** | Interactive plots | Plotly, PCA, networkx, sklearn | Box plots, PCA, network, heatmap, stacked bars |

---

## 📈 **Ultimate Reporter & Visualisation**

### **Ultimate Reporter HTML** – the heart of PseudoScope

- **Gene‑centric tables** – Every AMR or virulence gene is shown with the **complete list of genomes** that carry it (no truncation).
- **Sortable & searchable** – Click column headers to sort; type in the search box to filter.
- **Filter buttons** – One‑click filters for:
  - Carbapenemases (KPC, NDM, VIM, IMP, OXA‑48‑like)
  - Colistin resistance (mcr‑1 … mcr‑10)
  - ESBLs (GES, PER, VEB)
  - 16S rRNA methyltransferases (armA, rmt)
  - Efflux pumps, aminoglycoside modifying enzymes, etc.
- **ST‑O combination table** – Shows each ST and O‑type pair with the list of samples.
- **FASTA QC integration** – N50, GC%, total length displayed alongside typing results.
- **CSV export** – Every table can be exported with one click.
- **AI‑assisted guide** – Suggests questions for ChatGPT / Claude / Gemini.

### **Visualisation Dashboard** – interactive plots for publication

- **ST and O‑type bar charts** – Counts per type.
- **Top 20‑50 AMR & virulence genes** – Stacked bar charts by database (adjustable).
- **Gene co‑occurrence network** – Top 100 edges, interactive (hover, drag, zoom).
- **FASTA QC box plots** – GC%, AT%, N50, total length with all points (jitter) and outliers.
- **PCA of gene presence/absence** – Colored by ST, explained variance labels.
- **Co‑occurrence heatmap** – Pairwise counts of top 30 AMR genes.

All plots are generated with **Plotly** – fully interactive, downloadable as PNG/SVG, and embedded in a single scrolling HTML file.

---

## 🔬 **Validation & Accuracy**

| Reference Genome | Expected ST | Expected O‑type | PseudoScope Result |
|------------------|-------------|-----------------|--------------------|
| PAO1 | ST549 | O5 | ✅ ST549, O5 |
| PA14 | ST253 | O6 | ✅ ST253, O6 |
| LESB58 | ST146 | O1 | ✅ ST146, O1 |
| 6077 (clinical) | ST111 | O12 | ✅ ST111, O12 |
| 1168 (clinical) | ST235 | O1 | ✅ ST235, O1 |

**Accuracy**: 100% concordance for MLST (7‑gene profile) and 100% for O‑serotyping against in‑silico reference sets.

---

## 🖥️ **Sample Report**

See a complete interactive report generated by PseudoScope:

[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://github.com/bbeckley-hub/pseudoscope/blob/main/docs/sample_report.html)

*The report includes MLST summary, serotyping table, AMR and virulence gene tables with filter buttons, ST‑O combination table, FASTA QC metrics, and CSV export.*

---

## ⚙️ **Mandatory First‑Run Steps**

```bash
# 1. Set up ABRicate databases (required for resistance/virulence screening)
abricate --setupdb

# 2. Download AMRfinderPlus database (required for AMR analysis)
pseudoscope --update-amr-db
```

These steps are **required only once**; the databases will be reused for all future runs.

---

## ❓ **Not Happy? Try These Alternatives**

PseudoScope is optimised for *P. aeruginosa* and integrates many databases out of the box. However, if you need a different workflow, consider these excellent tools:

- **[Bactopia](https://github.com/bactopia/bactopia)** – A flexible, multi‑species pipeline for bacterial genomes (supports raw reads, uses Nextflow).
- **[Nullarbor](https://github.com/tseemann/nullarbor)** – A read‑to‑report pipeline for public health microbiology (great for outbreak surveillance).
- **[Pathogenwatch](https://pathogen.watch)** – A web‑based platform for genomics surveillance (requires upload, but very user‑friendly).

Each has its own strengths; PseudoScope focuses on **speed, local execution, and P. aeruginosa‑specific features**.

---

## 🌐 **Part of the ESKAPE Web Tools Suite**

PseudoScope is one of a growing collection of species‑optimised pipelines for the **ESKAPE pathogens** + _E. coli _(Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae, Acinetobacter baumannii, Pseudomonas aeruginosa, and Enterobacter species). Other tools in development:

- **StaphScope** – *S. aureus* (MLST, spa, SCCmec)
- **Kleboscope** – *K. pneumoniae* (MLST, K/O typing, AMR)
- **Enteroscope** – *Enterobacter* (van typing, MLST)
- **AcinetoScope** – *A. baumannii* (MLST, OCL typing)
- **EnteroMark** - *Enterococcus* (van typing, MLST)
- **EcoliTyper** - *E. coli* (serotyping typing, MLST)

A unified **web interface** is planned that will allow users to upload FASTQ/FASTA and run any of these pipelines without installing anything locally. Stay tuned!

---

## 🔮 **Future Development**

**Planned features (2026–2027):**

- **Raw read support** – Direct analysis of FASTQ files using `shovill` (spades) assembly on the fly.
- **Machine learning module** – Predict virulence, outbreak potential, and clinical risk from gene profiles.
- **Real‑time database sync** – Automatic updates of MLST, PAST, and AMR databases without re‑installation.
- **Web interface** – Single‑page application to run PseudoScope online (file upload, progress monitoring, report download).
- **Plugin system** – Allow community‑contributed gene sets and analysis modules.
- **Docker‑optimised HPC deployment** – Pre‑configured SIF images for Slurm/PBS.

**Call for collaborators** – If you are interested in contributing to any of these features, please open an issue or contact the authors.

---

## 🙏 **Acknowledgements**

PseudoScope builds upon the incredible work of many open‑source developers and database curators. Special thanks to:

- **Torsten Seemann** (MLST, ABRicate, and foundational bioinformatics tools)
- **Robert Petit & Tim Read** (pasty / camlhmp‑blast‑regions for P. aeruginosa serotyping)
- **NCBI AMR team** (AMRFinderPlus)
- **PubMLST** (K. Jolley & M. Maiden)
- **CGE** (ResFinder, PlasmidFinder – their database structures inspired ABRicate integration)
- **All database maintainers** of CARD, VFDB, ARG-ANNOT, MEGARes, BacMet2, EcoH, Ecoli_VF
- **Early adopters and beta testers** at the University of Ghana Medical School

> “If we ever meet in person, the drinks are on me!” – Brown Beckley

---

## 💰 **Sponsor This Project**

PseudoScope is developed and maintained in our spare time. If you find it useful for your research or clinical work, please consider sponsoring:

- **GitHub Sponsors**: [https://github.com/sponsors/bbeckley-hub](https://github.com/sponsors/bbeckley-hub)
- **Support allows us to**:
  - Add raw read support faster
  - Host a free public web interface
  - Provide dedicated user support
  - Keep databases updated

Every contribution, no matter how small, helps keep the project alive. Thank you!

---

## 📚 **Citation**

If you use PseudoScope in your research, please cite:

> Beckley, B. (2026). PseudoScope: a species‑optimized computational pipeline for rapid and accessible *Pseudomonas aeruginosa* genomic typing and surveillance. *Journal of Open Source Software* (submitted).

```bibtex
@software{beckley2026pseudoscope,
  author = {Brown Beckley},
  title = {PseudoScope: A species-optimized computational pipeline for Pseudomonas aeruginosa genomic typing and surveillance},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/pseudoscope},
  version = {1.0.0}
}
```

**Integrated tools** – please also cite the original authors of MLST, ABRicate, pasty, AMRFinderPlus, PubMLST, CARD, ResFinder, VFDB, and others (full citation list in the repository).

---

## 👥 **Authors & Contact**

**Brown Beckley** (Primary Developer)  
University of Ghana Medical School – Department of Medical Biochemistry  
📧 brownbeckley94@gmail.com  
🐙 GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)  
📞 +233 508820617

For bug reports, feature requests, or collaborations, please use the [GitHub issue tracker](https://github.com/bbeckley-hub/pseudoscope/issues).

---

## 📄 **License**

The PseudoScope pipeline code (workflow engine, report generation, HTML templates, visualisation, and ultimate reporter) is licensed under the **MIT License** – see the `LICENSE` file for details.

Third‑party tools and databases are used under their respective licenses (GPL, Apache, Public Domain, etc.). By using PseudoScope, you agree to comply with those licenses.

---

<div align="center">

## 🚀 **Ready to type your *P. aeruginosa* genomes?**

[![Get Started](https://img.shields.io/badge/GET_STARTED-Now-green?style=for-the-badge&logo=conda)](https://github.com/bbeckley-hub/pseudoscope#-quick-start)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/pseudoscope/issues)
[![Sponsor](https://img.shields.io/badge/Sponsor-❤️-pink?style=for-the-badge&logo=githubsponsors)](https://github.com/sponsors/bbeckley-hub)

**From fragmented data to integrated insights – in minutes.**

*PseudoScope: Precision surveillance for the antibiotic resistance era.*

⭐ **If you find this tool useful, please star the repository!** ⭐

</div>
