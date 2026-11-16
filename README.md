
<p align="center">
  <img src="https://raw.githubusercontent.com/bbeckley-hub/pseudoscope/main/pseudoscope.png" alt="PseudoScope Banner" width="100%">
</p>

<div align="center">

# 🧬 PseudoScope

### **A species-optimized computational pipeline for rapid, accessible *Pseudomonas aeruginosa* genomic surveillance**

#### **Complete MLST, serotyping, AMR, virulence, MGE, and visualisation in minutes, not hours**

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
![Version](https://anaconda.org/bioconda/pseudoscope/badges/version.svg)
![Last Updated](https://anaconda.org/bioconda/pseudoscope/badges/latest_release_relative_date.svg)
![Platforms](https://anaconda.org/bbeckley-hub/pseudoscope/badges/platforms.svg)
![License](https://anaconda.org/bbeckley-hub/pseudoscope/badges/license.svg)

[![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/pseudoscope)](https://hub.docker.com/r/bbeckleyhub/pseudoscope)
[![Docker Image Size](https://img.shields.io/docker/image-size/bbeckleyhub/pseudoscope/latest)](https://hub.docker.com/r/bbeckleyhub/pseudoscope)
[![Docker Version](https://img.shields.io/docker/v/bbeckleyhub/pseudoscope?sort=semver)](https://hub.docker.com/r/bbeckleyhub/pseudoscope)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](#)
[![Stage](https://img.shields.io/badge/status-active-brightgreen)](#)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Profile-0A66C2?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/brown-beckley-190315319)
![Conda Downloads](https://img.shields.io/conda/dn/bioconda/pseudoscope?label=Conda%20Downloads)

[![Documentation](https://img.shields.io/badge/docs-mkdocs-526CFE?logo=materialformkdocs)](https://bbeckley-hub.github.io/pseudoscope)
[![Last Commit](https://img.shields.io/github/last-commit/bbeckley-hub/pseudoscope)](https://github.com/bbeckley-hub/pseudoscope/commits)
[![Contributors](https://img.shields.io/github/contributors/bbeckley-hub/pseudoscope)](https://github.com/bbeckley-hub/pseudoscope/graphs/contributors)
[![Security: bandit](https://img.shields.io/badge/security-bandit-yellow.svg)](https://github.com/PyCQA/bandit)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![CI](https://img.shields.io/github/actions/workflow/status/bbeckley-hub/pseudoscope/ci.yml?branch=main&label=CI)](https://github.com/bbeckley-hub/pseudoscope/actions)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)](https://github.com/bbeckley-hub/pseudoscope/tests)

[![Try Web Version](https://img.shields.io/badge/🌐-Try%20PseudoScope%20Web-4285F4?logo=googlechrome&logoColor=white)](https://eskape.bio)
[![Carbapenemase](https://img.shields.io/badge/Detects-OXA%20%7C%20IMP%20%7C%20VIM%20%7C%20NDM%20%7C%20KPC-00BCD4.svg)](#)
[![Serotyping](https://img.shields.io/badge/O--Antigen-PAST%20serotyping-673AB7.svg)](#)
[![AI Ready](https://img.shields.io/badge/AI-Ready%20Reports-00C853.svg)](#-ai-disclosure)

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/en/latest/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub Issues](https://img.shields.io/badge/issues-Open-orange)](https://github.com/bbeckley-hub/pseudoscope/issues)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/pseudoscope)](https://github.com/bbeckley-hub/pseudoscope/stargazers)
[![Google Scholar](https://img.shields.io/badge/Google%20Scholar-Profile-4285F4?style=flat&logo=googlescholar&logoColor=white)](https://scholar.google.com/citations?user=CYNOsqIAAAAJ&hl=en)

**Two ways to use PseudoScope:**
🖥️ **Command-line tool** for high-throughput, local analysis
🌐 **PseudoScope Web** for non-bioinformaticians, at [https://eskape.bio](https://eskape.bio)

</div>

---

## 📋 **Table of Contents**

- [🎉 What's New in v1.3.0](#-whats-new-in-v130)
- [🎯 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [📊 Sample Reports](#-sample-reports)
- [🔥 Feature Deep-Dives](#-feature-deep-dives)
- [🌐 PseudoScope Web Platform](#-pseudoscope-web-platform)
- [⚡ Quick Start](#-quick-start)
- [🔧 Installation](#-installation)
- [🐳 Docker & Singularity](#-docker--singularity)
- [🚀 Usage Guide](#-usage-guide)
- [🔍 Core Analytical Modules](#-core-analytical-modules)
- [⚙️ Database Management](#️-database-management)
- [❓ Not Happy? Try These Alternatives](#-not-happy-try-these-alternatives)
- [🌐 Part of the ESKAPE Web Tools Suite](#-part-of-the-eskape-web-tools-suite)
- [🔮 Future Development](#-future-development)
- [❓ FAQ](#-faq)
- [🐛 Troubleshooting](#-troubleshooting)
- [🙏 Acknowledgements](#-acknowledgements)
- [💰 Sponsor This Project](#-sponsor-this-project)
- [📚 Citation](#-citation)
- [🤖 AI Disclosure](#-ai-disclosure)
- [📚 Third-Party Tool Citations](#-third-party-tool-citations)
- [👥 Authors & Contact](#-authors--contact)
- [📄 License](#-license)

---

## 🎉 **What's New in v1.3.0**

### PSEUDA is in development

**PSEUDA (Pseudomonas Unified Epidemiological Database)** is our answer to the "how do I find the genomes I actually need?" problem in *P. aeruginosa* surveillance. It follows the same model as [SAURUS](https://bbeckley-hub.github.io/saurus/) for *S. aureus*.

The common query: all ST235 bloodstream isolates from Asia collected after 2020, carrying *bla*KPC and exoU. Today that means downloading multiple compressed tables, writing join scripts, and cross-referencing typing against metadata. PSEUDA will turn it into a filter form. Click, filter, export an accession list ready for `datasets download genome accession --inputfile list.txt`.

The interface will run entirely over HTTPS on GitHub Pages. No installation, no account, no API key.

**Status:** under active development. Progress at [bbeckley-hub.github.io/pseuda](https://bbeckley-hub.github.io/pseuda/). Meanwhile, the CLI already runs species confirmation, typing, resistance, virulence, and MGE profiling end to end on any *P. aeruginosa* assembly you provide.

### v1.3.0 improvements

- **Report-level control.** New `--report-level` flag with `final`, `summary`, and `full` modes. Prevents disk exhaustion on 300+ sample cohorts by pruning intermediate files. `summary` is the default.
- **MLST database management.** New `--pull-mlst-db` (download pre-built scheme from GitHub, no API key) and `--update-mlst-db` (fetch the current Pasteur scheme from PubMLST via mlstdb, API key required).
- **MGE profiling with mobileOG-db.** Full module using Prodigal and DIAMOND against the 350 MB mobileOG-db catalogue. Baked into Docker and Singularity images. Conda users run `--pull-mge-db` once.
- **Environment-friendly orchestrator.** Every module runs in an isolated temporary directory, honouring `--tmpdir`, `PSEUDOSCOPE_TMPDIR`, and `TMPDIR`. `--keep-temp` preserves work dirs for debugging. Report-level pruning runs at the end. Tested on HPC login nodes, Docker, Singularity, and shared filesystems.
- **Virulence expander.** Click the count in the Sample Overview table to expand the full VFDB gene list for that isolate, without leaving the report.
- **Compare tab.** Pick any two isolates to get a similarity verdict, a typing-by-typing diff, and gene content split three ways (shared, only in A, only in B).
- **Cluster tab.** Automatic pairwise similarity matrix, colour-coded heatmap, union-find cluster detection, and per-cluster summaries.
- **Consistent report styling.** Unified footers across all modules, fixed ASCII banner rendering, and consistent colour palette throughout every HTML output.

---

## 🎯 **Overview**

**PseudoScope** is an automated pipeline for comprehensive *Pseudomonas aeruginosa* genomic surveillance. It integrates **eight essential analysis modules** into a single cohesive workflow, runs entirely on your own machine, and produces interactive HTML reports in minutes.

### The problem

- **Fragmented bioinformatics.** Typing, resistance detection, virulence profiling, and mobile genetic element analysis each need their own tool, their own dependencies, and their own quirks.
- **P. aeruginosa is not a generic Gram-negative.** Carbapenemase profiling matters more than KPC alone. Efflux systems matter more.
- **Web services raise data privacy questions** and usually have upload caps.
- **Generalist platforms take hours.** Outbreak investigations do not.

### Our solution

- Single-command installation via Conda.
- Complete analysis in tens of minutes for a typical batch.
- 100% local execution. Genomes stay on your machine.
- Interactive HTML reports with typing, resistance, virulence, and MGE context.
- Automatic carbapenemase and ESBL classification with dual-caller consensus.
- Report-level control to keep disk usage flat on large batches.
- A web interface for users who do not live in a terminal.

**Built for:** clinical labs, reference centres, outbreak investigations, genomic surveillance networks, and any group doing routine *P. aeruginosa* typing.

---

## ✨ **Key Features**

### 🔬 P. aeruginosa-specific innovations

- **Carbapenemase detection across all major families.** OXA-2, OXA-10, OXA-14, OXA-17, OXA-19, OXA-28, OXA-35, OXA-45, OXA-50, OXA-198, plus IMP, VIM, NDM, KPC, GES, SPM, AIM, DIM. Cross-checked between AMRFinderPlus and ABRicate so discordances are visible.
- **ESBL coverage tuned for P. aeruginosa.** *bla*PER, *bla*VEB, *bla*BEL, *bla*GES, plus the standard *bla*TEM, *bla*SHV, *bla*CTX-M.
- **Chromosomal colistin resistance.** *pmrA/B*, *phoP/Q*, *mgrB*, *lpxA/C/D*, *arnT*, *eptA/B*, *basS/R*, alongside the *mcr* family.
- **Efflux pump profiling.** *mexA/B/C/D/E/F/X/Y*, *oprM/N/J*.
- **Aminoglycoside methyltransferases.** *armA*, *rmtA* through *rmtH*, *npmA*, plus the common *aac/ant/aph* families.
- **T3SS effector typing.** exoU, exoS, exoT, exoY, cross-tabulated in the dashboard.
- **O-antigen serotyping** via PAST (pasty). Not a re-skinned spa typer.
- **MLST (Pasteur scheme)** for standard 7-gene sequence typing.
- **MGE profiling** via mobileOG-db, with 10 functional categories and per-contig composition.
- **Dual-caller consensus** for carbapenemases and other markers where AMRFinderPlus and ABRicate disagree.
- **Report-level control** so large cohorts never fill the disk.
- **Interactive compare and cluster tabs** for outbreak pair and cluster analysis.

### 🚀 Performance

- Complete analysis in 2 to 5 minutes for a single genome on a modern workstation.
- Roughly 15 to 30 minutes for 24 genomes on 16 cores.
- MGE dominates runtime; `--skip-mge` halves it when MGE is not the question.
- Runs on 4 GB RAM and scales to HPC.

---

## 📊 **Sample Reports**

Every report is a self-contained HTML file. No server, no dependencies, works offline. **[Open the live reports →](https://bbeckley-hub.github.io/pseudoscope/)**

| Report | What it shows | Direct link |
|---|---|---|
| 🧬 **Gene-centric report** | Every AMR, virulence, BacMet, and MGE gene, with the genomes that carry it | [View sample](https://bbeckley-hub.github.io/pseudoscope/genius_pseudomonas_ultimate_gene_centric_report.html) |
| 📋 **Sample-centric report** | Per-isolate cards with virulence expander, Compare, and Cluster tabs | [View sample](https://bbeckley-hub.github.io/pseudoscope/genius_pseudomonas_sample_centric_report.html) |
| 📈 **Interactive dashboard** | Plotly dashboard: Overview, Typing, QC, AMR, Virulence, MGE, Alerts, Story, Compare | [View sample](https://bbeckley-hub.github.io/pseudoscope/genius_pseudomonas_visual_dashboard.html) |
| 🔬 **MLST summary** | Sequence types and allele profiles across the cohort | [View sample](https://bbeckley-hub.github.io/pseudoscope/mlst_summary.html) |
| 🧬 **PAST summary** | O-antigen serotypes across the cohort | [View sample](https://bbeckley-hub.github.io/pseudoscope/past_summary.html) |
| 📏 **FASTA QC** | N50, GC%, contig stats plus fastANI species confirmation | [View sample](https://bbeckley-hub.github.io/pseudoscope/PseudoScope_FASTA_QC_summary.html) |
| 💊 **AMR summary** | Resistance genes and risk classification per isolate | [View sample](https://bbeckley-hub.github.io/pseudoscope/pseudo_amrfinder_summary_report.html) |
| 📱 **MGE profile** | mobileOG-db signatures across 10 functional categories | [View sample](https://bbeckley-hub.github.io/pseudoscope/pseudoscope_mge_summary.html) |

> Reports may take a few seconds to load from GitHub Pages. For the best experience, download the HTML files and open them locally. Every report is fully offline-capable.

---

## 🔥 **Feature Deep-Dives**

<details>
<summary><b>⚖️ The Compare tab. Because "I think they're the same bug" isn't good enough (click to expand)</b></summary>
<br>

Every outbreak investigation eventually reaches the same moment. Someone has two genomes. They open two reports. They scroll up and down, mentally diffing MLST, serotype, carbapenemase, AMR, virulence, MGE. Then they say "yeah, I think they match."

That's an eyeball, not a methodology. The Compare tab replaces it.

Pick any two isolates and you get:

1. **A verdict banner.** Red for near-identical (≥95% similarity), orange for high (≥85%), yellow for moderate (≥60%), blue for distinct.
2. **A similarity gauge.** Animated coloured ring showing overall match percentage.
3. **A typing table.** Every field side by side with green or red pills for match or diff.
4. **Gene content split three ways.** Shared, only in A, only in B. Per category: AMR, virulence, BacMet, MGE.
5. **Full profile expanders.** Click any category to see both isolates' complete gene lists.

**Why it matters in practice:**

- **Outbreak confirmation.** Two isolates from the same ward with identical MLST, O-type, carbapenemase profile, and T3SS effectors are a plausible transmission pair. Infection control will act on that.
- **Discordance analysis.** An ancestor and a resistant descendant look identical in every field except one. That one field tells you which element moved.
- **QC checks.** Two genomes from the same sequencing run coming back identical. Better to catch a duplicate now than publish it.
- **Teaching.** Explaining "why isn't this susceptible?" is faster with two profiles side by side than with two open tabs.

**A caveat we state openly.** "Identical typing = transmission pair" is a hypothesis generator, not a confirmation. Two isolates can share MLST, O-type, and resistance profile and still be unrelated. Those markers evolve too slowly for hospital-level resolution. The gold standard remains core-genome SNP distance, which PseudoScope does not compute and should not at this scope. We use the phrase *"possible transmission pair"* deliberately. If we ever lose that hedge, shout at us.

</details>

<details>
<summary><b>🧩 Cluster mode. Because 40 isolates is not one answer (click to expand)</b></summary>
<br>

When you have more than two samples, the Compare tab's pairwise view becomes a chore. Switch to Cluster mode.

- **Automatic pairwise similarity matrix** across every isolate pair.
- **Colour-coded heatmap** so the eye catches clusters without prompting.
- **Union-find cluster detection** on your own typing data. No external reference, no database lookup. It groups isolates whose typing overlap is high enough that clustering them is defensible.
- **Per-cluster summary cards** showing mean intra-cluster similarity, representative isolate, and lineage markers.

The output reads like a briefing: "we have 40 isolates in 6 clusters, 3 of which warrant follow-up." That is what infection control can actually use.

</details>

<details>
<summary><b>🦠 The virulence expander. Because 400 gene names don't fit in a cell (click to expand)</b></summary>
<br>

Before v1.3.0, if you wanted to know which virulence genes a specific isolate carried, you opened the Virulence tab, found the row for each gene, and mentally ticked off which genomes were listed. For 40 samples, you did that once and never wanted to do it again.

Now, the Sample Overview table has a Virulence column with a click-to-expand count. Click the number and a panel slides open showing every VFDB gene that isolate carries, colour-coded to match the rest of the report.

**Why it matters:**

- **Case-by-case drill-down.** When a clinician asks what an isolate carries, you have the answer in one click.
- **Teaching.** Explaining PVL, exotoxin A, or T3SS effector combinations is easier when the gene list is right there.
- **Screening.** Sort by the count column and the highest-virulence isolates float to the top.
- **No page reload.** This is a `details/summary` element. Still one file. Still offline. No JavaScript framework.

</details>

<details>
<summary><b>📱 The MGE tab. Where resistance gets a postal address (click to expand)</b></summary>
<br>

You have AMR genes. You have virulence genes. You know what a genome carries. But how did it get there?

Mobile genetic elements are the answer. Plasmids, prophages, transposons, integrative conjugative elements. These are the vehicles that move resistance and virulence genes between unrelated bacteria. Without understanding the MGE landscape you're looking at the mail without seeing the postal system.

**What the MGE tab gives you:**

- **Per-sample counts** across ten functional categories: Integrase, Transfer, Stability, Phage, Replication, IS-associated, ICE-associated, Plasmid-associated, Phage-associated, Key MGE signatures.
- **Aggregate profiles by typing.** Mean MGE load per ST, per serotype, per carbapenemase status.
- **Top-N views.** The five isolates with the most IS-associated proteins, or the most plasmid proteins, or the highest key-signature count.
- **Per-contig composition.** Which contigs carry which types of MGE proteins.

**What it helps with:**

- **Resistance dissemination.** A plasmid-associated gene cluster in two unrelated STs is a signal for horizontal transfer.
- **Virulence acquisition.** Exotoxin and T3SS genes can be prophage- or PAPI-borne.
- **Persistent infections.** Strong stability systems let plasmids survive without selection.
- **Outlier detection.** An isolate with 2x the cohort's MGE load is either something new or something broken. Either way, look.

**A word of honesty.** mobileOG-db counts protein signatures, not reconstructed elements. Multiple proteins from the same element produce multiple hits, and some categories also match chromosomal genes. That's why we say "signatures" everywhere. Element-level claims need co-localization analysis, which is a bioinformatics project, not a column in a dashboard.

</details>

<details>
<summary><b>🛡️ The Carbapenemase tab. Because one caller is never enough (click to expand)</b></summary>
<br>

Carbapenem resistance in *P. aeruginosa* is not one gene. It is OXA variants, metallo-beta-lactamases, and chromosomal AmpC upregulation, all coexisting. Treating it as "did we find blaKPC?" misses most of the picture.

**The tab shows:**

- **OXA-type distribution** across the cohort. OXA-2, OXA-10, OXA-14, OXA-17, OXA-19, OXA-28, OXA-35, OXA-45, OXA-50, OXA-198. Cross-tabbed with ST, serotype, and MGE load.
- **MBL distribution.** IMP, VIM, NDM, SPM, AIM, DIM.
- **Serine carbapenemase distribution.** KPC, GES.
- **Dual-caller agreement.** AMRFinderPlus and ABRicate each report carbapenemases. Where both agree, high confidence. Where they disagree, a signal to look closer.
- **Source class breakdown.** Chromosomal, plasmid-borne, integron-associated.

**Why two callers matter:**

- AMRFinderPlus is curated and conservative.
- ABRicate sweeps ten databases and can flag homologues AMRFinderPlus misses.
- When both call the same variant, you have a defensible answer for a clinical report.
- When only one calls it, you have a specific thing to verify.

</details>

<details>
<summary><b>💊 O-antigen serotyping. Because Gram-negatives have coats too (click to expand)</b></summary>
<br>

*P. aeruginosa* wears an O-antigen coat on its lipopolysaccharide. The serotype matters for vaccine coverage, for outbreak tracing, and for understanding persistence in cystic fibrosis.

**Why serotype matters:**

- **Vaccine relevance.** O-antigen is a target for investigational *P. aeruginosa* vaccines. If a candidate targets O5 and your cohort is O11-dominant, the coverage story changes.
- **Lineage correlation.** Specific O-types cluster with specific clonal complexes. PAO1 is ST549/O5, PA14 is ST253/O10, LESB58 is ST146/O6. Divergence is informative.
- **Cross-tabulation.** The tab cross-tabs serotype against MLST, carbapenemase profile, and T3SS effector type.
- **Confidence score.** PAST reports coverage and hits per candidate O-type, so you can see whether the call is clean or contested.

**How it's called.** Serotyping uses pasty against the PAST O-antigen locus catalogue. Two thresholds are configurable (`--past-min-pident`, `--past-min-coverage`, default 95% for both) so you can tighten or loosen stringency depending on assembly quality.

</details>

---

## 🌐 **PseudoScope Web Platform**

PseudoScope is part of the **ESKAPE AMR Platform**, a unified web suite for species-optimized genomic surveillance.

- Drag-and-drop upload (single FASTA, multiple files, or a ZIP).
- Module selection. Pick which analyses to run.
- Live progress with real-time logs.
- Download all results as a single ZIP.
- No installation required. Any modern browser works.

> 🌐 **Try PseudoScope Web:** [https://eskape.bio](https://eskape.bio)
> 📦 **Web repository:** [https://github.com/bbeckley-hub/eskape-web-platform](https://github.com/bbeckley-hub/eskape-web-platform)

*Note: the web version caps uploads at 10 files per job. Larger datasets need the CLI.*

---

## ⚡ **Quick Start**

```bash
# 1. Install PseudoScope
conda create -n pseudoscope -c conda-forge -c bioconda pseudoscope -y
conda activate pseudoscope

# 2. Run your first analysis
pseudoscope -i "*.fna" -o results --threads 4

# 3. Open the sample-centric report
xdg-open results/PSEUDOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS/*.html
```

### Database setup (conda users only, once)

```bash
abricate --setupdb           # ABRicate databases
pseudoscope --update-amr-db  # AMRFinderPlus database
pseudoscope --pull-mlst-db   # Pasteur MLST scheme from GitHub (no API key)
pseudoscope --pull-mge-db    # mobileOG-db payload (~300 MB)
```

Docker and Singularity users skip this block. Every database is baked into the image.

---

## 🔧 **Installation**

### System requirements

| Resource | Minimum | Recommended |
|---|---|---|
| CPU cores | 2 | 8+ |
| RAM | 4 GB | 16 GB |
| Storage | 10 GB | 20 GB |
| OS | Linux, macOS, WSL2 | Ubuntu 22.04+ |

<details>
<summary><b>Installation methods (click to expand)</b></summary>

**Conda (recommended)**

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n pseudoscope python=3.9 pseudoscope -y
conda activate pseudoscope
```

**Mamba (faster)**

```bash
mamba create -n pseudoscope -c conda-forge -c bioconda pseudoscope -y
mamba activate pseudoscope
```

**From source**

```bash
git clone https://github.com/bbeckley-hub/pseudoscope.git
cd pseudoscope
conda env create -f environment.yml
conda activate pseudoscope
pip install -e .
```

</details>

---

## 🐳 **Docker & Singularity**

Docker and Singularity users get everything pre-installed. No first-run setup, no `--pull-mlst-db`, no `--update-amr-db`, no `abricate --setupdb`. Pull and run.

> **Conda users:** run the four database commands in Quick Start once after install. Docker and Singularity users skip that step entirely.

---

### 🐳 **Docker (laptops, workstations, cloud VMs)**

Docker runs containers as root by default. Files written to a bind mount come out owned by `root:root`, so you end up running `sudo chown` after every analysis. Two fixes, choose one.

<details>
<summary><b>Option A: pre-create the output directory yourself (click to expand)</b></summary>

```bash
mkdir -p "$PWD/results"
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$PWD/genomes:/data:ro" \
    -v "$PWD/results:/results" \
    bbeckleyhub/pseudoscope:1.3.0 \
    -i "/data/*.fna" -o /results -t 4 --report-level summary
```

</details>

<details>
<summary><b>Option B: wrapper script (click to expand)</b></summary>

```bash
cat > ~/bin/pseudoscope-run << 'EOF'
#!/bin/bash
set -e
mkdir -p "$PWD/results"
exec docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$PWD:/work" \
    -w /work \
    bbeckleyhub/pseudoscope:latest \
    "$@"
EOF
chmod +x ~/bin/pseudoscope-run
```

Then:

```bash
cd ~/PA/test_run
pseudoscope-run -i "/work/genomes/*.fna" -o /work/results -t 4 --report-level summary
```

Files land owned by you, on any host UID. No `sudo chown` required.

</details>

<details>
<summary><b>Verify the image and open a debugging shell (click to expand)</b></summary>

```bash
docker run --rm bbeckleyhub/pseudoscope:1.3.0 --help
docker run --rm bbeckleyhub/pseudoscope:1.3.0 abricate --list
docker run --rm -it bbeckleyhub/pseudoscope:1.3.0 bash
```

Image size: about 1.2 GB compressed, about 5.2 GB on disk.

</details>

---

### 🔬 **Singularity / Apptainer (HPC clusters)**

HPC login and compute nodes do not run the Docker daemon. Singularity (or Apptainer) is the standard. Same SIF format, rootless by default, no UID juggling.

<details>
<summary><b>1. Load the runtime (click to expand)</b></summary>

```bash
module load apptainer         # newer clusters
# or
module load containers/singularity/3.6.4

module load squashfs-tools/4.4  # only needed for the pull
singularity --version
```
</details>

<details>
<summary><b>2. Pull the image (click to expand)</b></summary>

```bash
singularity pull pseudoscope_1.3.0.sif docker://bbeckleyhub/pseudoscope:1.3.0
```

Creates a ~1.2 GB SIF in the current directory. `squashfs-tools` is only needed for the pull, not for running.
</details>

<details>
<summary><b>3. Run the pipeline (click to expand)</b></summary>

```bash
singularity run --bind "$PWD:/data" \
    pseudoscope_1.3.0.sif \
    -i "/data/*.fna" -o /data/output -t 16
```

- No `-u` needed. Singularity runs as your real UID.
- No `-e HOME=` needed. MLST and AMR caches persist under your real `$HOME`.

> ⚠️ Use `singularity run`, not `singularity exec`. The image entrypoint sets `PATH` and `LD_LIBRARY_PATH` for the conda env. `exec` skips it and breaks AMRFinderPlus.
</details>

<details>
<summary><b>4. SLURM job script (click to expand)</b></summary>

```bash
#!/bin/bash
#SBATCH --job-name=pseudoscope
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=pseudoscope_%j.log

module load containers/singularity/3.6.4
cd $SLURM_SUBMIT_DIR

singularity run --bind "$PWD:/data" \
    /path/to/pseudoscope_1.3.0.sif \
    -i "/data/*.fna" \
    -o /data/output \
    -t ${SLURM_CPUS_PER_TASK}
```
</details>

<details>
<summary><b>5. HPC notes (click to expand)</b></summary>

- Pull once, run everywhere. The SIF is a single self-contained file. Share it on scratch and every compute node can use it.
- Small home quota? Redirect the build cache to scratch:
  ```bash
  export SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
  ```
- MLST and AMR caches live under `$HOME/.local/share/pseudoscope/`. On tight quotas:
  ```bash
  export PSEUDOSCOPE_MLST_DB=/scratch/$USER/pseudoscope/mlst_db
  export PSEUDOSCOPE_MGE_DATA_DIR=/scratch/$USER/pseudoscope/mobileOG-db
  ```
</details>

<details>
<summary><b>6. Troubleshooting (click to expand)</b></summary>

| Error | Fix |
|---|---|
| `mksquashfs: executable file not found` | `module load squashfs-tools/4.4` before `singularity pull` |
| `manifest unknown` | Tag not pushed or name wrong. Check https://hub.docker.com/r/bbeckleyhub/pseudoscope/tags |
| `singularity: command not found` | `module avail 2>&1 \| grep -i singularity` |
| `Bundled AMRfinderPlus not properly installed` | You used `singularity exec`. Use `run`. |
| `Permission denied: /results/...` in Docker | `mkdir -p results` as yourself before `docker run` |
</details>

---

### 🤔 **Which should you use?**

| Platform | Recommendation | Output ownership |
|---|---|---|
| Laptop, workstation | Docker | You (with `--user` + `mkdir`) |
| HPC login and compute nodes | **Singularity / Apptainer** | You (automatic) |
| Cloud VM | Either. Docker is simpler | You |
| Shared cluster without Docker | Singularity | You (automatic) |

The rule: Docker where you control the daemon, Singularity everywhere else.

---

## 🚀 **Usage Guide**

### Basic commands

```bash
# Single genome
pseudoscope -i genome.fna -o results

# Batch with glob pattern
pseudoscope -i "*.fna" -o batch_results --threads 16

# Skip specific modules
pseudoscope -i "*.fasta" -o results --skip-mge --skip-viz

# Report-level control for large cohorts
pseudoscope -i "*.fna" -o results --report-level final
```

<details>
<summary><b>All command-line flags (click to expand)</b></summary>

| Flag | Description |
|---|---|
| `-i`, `--input` | Input FASTA file(s). Supports glob patterns. |
| `-o`, `--output` | Output directory for all results. |
| `-t`, `--threads` | Number of CPU threads (default: 2). |
| `--quiet` | Suppress all non-error output. |
| `--keep-temp` | Preserve temporary directories (debugging). |
| `--tmpdir PATH` | Override temporary directory root (HPC, Docker). |
| `--clean-output` | Delete output directory before analysis. |
| `--report-level {final,summary,full}` | Controls which files are kept after the run. |
| `--version` | Show version and exit. |
| `--update-amr-db` | Update AMRFinderPlus database (incremental). |
| `--force-update-amr-db` | Force complete AMR database overwrite. |
| `--pull-mlst-db` | Download pre-built Pasteur MLST scheme from GitHub. |
| `--update-mlst-db` | Update MLST scheme from PubMLST via mlstdb. |
| `--pull-mge-db` | Download the mobileOG-db payload (~300 MB). |
| `--skip-qc` | Skip FASTA QC. |
| `--skip-mlst` | Skip MLST. |
| `--skip-past` | Skip PAST serotyping. |
| `--skip-abricate` | Skip ABRicate. |
| `--skip-amr` | Skip AMRFinderPlus. |
| `--skip-mge` | Skip MGE profiling. |
| `--skip-summary` | Skip gene-centric ultimate summary. |
| `--skip-sample-centric` | Skip sample-centric interactive report. |
| `--skip-visualization` | Skip visualisation dashboard. |
| `--abricate-minid` | ABRicate minimum % identity (default: 80). |
| `--abricate-mincov` | ABRicate minimum % coverage (default: 80). |
| `--past-min-pident` | PAST minimum % identity (default: 95). |
| `--past-min-coverage` | PAST minimum % coverage (default: 95). |
| `--amr-min-identity` | AMR minimum identity, 0 to 1. |
| `--amr-min-coverage` | AMR minimum coverage, 0 to 1. |
| `--amr-skip-mutations` | Disable point mutation reporting. |
| `--amr-force-update` | Force AMR database update before analysis. |

**Supported FASTA formats:** `.fna`, `.fasta`, `.fa`, `.fsa`

</details>

<details>
<summary><b>Report levels (click to expand)</b></summary>

`--report-level` controls what survives a run. Choose based on batch size and how much you want to archive.

| Level | What is kept | Best for |
|---|---|---|
| **final** | Only `PseudoScope_final_report/`. All module directories are deleted. | Archival, very large batches. |
| **summary** | `PseudoScope_final_report/` plus batch TSVs, HTMLs, and JSONs in each module directory. Per-sample folders are removed. | Default. Most workflows. |
| **full** | Everything, including all per-sample HTML reports. | Small batches, debugging. |

</details>

<details>
<summary><b>Output structure (click to expand)</b></summary>

```
results/
├── PseudoScope_final_report/
│   ├── PSEUDOSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS/
│   ├── PSEUDOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS/
│   ├── PSEUDOSCOPE_VISUAL_DASHBOARD/
│   ├── pseudoscope_mge_summary.html
│   ├── pseudoscope_mge_hits.tsv
│   └── pseudoscope_mge_summary.json
├── fasta_qc_results/
├── mlst_results/
├── past_results/
├── abricate_results/
├── amrfinder_results/
├── mge_results/
└── pseudoscope_run.log
```

Per-sample folders inside each module directory contain the raw and HTML reports for that isolate. They survive only when `--report-level full` is used.

</details>

---

## 🔍 **Core Analytical Modules**

| Module | Purpose | Key databases or methods | Outputs |
|---|---|---|---|
| **FASTA QC** | Assembly quality and species confirmation | Biopython, fastANI | N50, GC%, contig stats, ANI %, species pass/fail |
| **MLST** | Pasteur scheme typing | mlst, PubMLST | ST, 7-gene allele profile |
| **PAST** | O-antigen serotyping | pasty, camlhmp | O-type, coverage, hits per candidate |
| **ABRicate** | Multi-database screening | NCBI, CARD, ResFinder, VFDB, ARG-ANNOT, PlasmidFinder, MEGARes, EcoH, BacMet2, Ecoli_VF | Per-database summaries, master JSON |
| **AMRFinderPlus** | AMR genes plus point mutations | NCBI AMRFinderPlus | Resistance genes, point mutations, mechanism classification |
| **MGE** | Mobile genetic element profiling | mobileOG-db, Prodigal, DIAMOND | 10 functional categories, per-contig composition, key MGE signatures |
| **Gene-centric reporter** | Cross-genome integration | All modules' summaries | Interactive HTML, CSV exports |
| **Sample-centric reporter** | Isolate-level view with Compare and Cluster | All modules' summaries | Per-isolate cards, similarity matrix, cluster detection |
| **Visualisation** | Interactive plots | Plotly, PCA, networkx, scikit-learn | Box plots, PCA, co-occurrence network, heatmap, stacked bars |

---

## ⚙️ **Database Management**

PseudoScope ships with the databases it needs but they can be refreshed.

<details>
<summary><b>ABRicate</b></summary>

```bash
abricate --setupdb
```

Ten databases: NCBI, CARD, ResFinder, VFDB, ARG-ANNOT, PlasmidFinder, MEGARes, EcoH, BacMet2, Ecoli_VF.

</details>

<details>
<summary><b>AMRFinderPlus</b></summary>

```bash
pseudoscope --update-amr-db          # incremental
pseudoscope --force-update-amr-db    # full overwrite
```

Requires outbound HTTPS to `ftp.ncbi.nlm.nih.gov`.

</details>

<details>
<summary><b>MLST</b></summary>

Two commands. Different sources, different keys.

```bash
# Download pre-built Pasteur scheme from GitHub (no API key)
pseudoscope --pull-mlst-db

# Fetch the live Pasteur scheme from PubMLST (requires mlstdb API key)
pseudoscope --update-mlst-db
```

If you don't have `mlstdb` configured, one-time setup:

```bash
mlstdb connect --db pubmlst --api-key
# Paste your API key when prompted
# Register at https://pubmlst.org/site-accounts and generate a key in your profile
```

For most users `--pull-mlst-db` is sufficient. Use `--update-mlst-db` if you're typing an isolate with a novel allele that the GitHub snapshot doesn't know about yet.

</details>

<details>
<summary><b>MGE (mobileOG-db)</b></summary>

```bash
pseudoscope --pull-mge-db
```

Downloads about 300 MB from the GitHub release. One-time setup for conda users. Docker and Singularity images already have it baked in.

</details>

---

## ❓ **Not Happy? Try These Alternatives**

PseudoScope is optimized for *P. aeruginosa* and integrates many databases out of the box. If you need a different workflow:

- **[Bactopia](https://github.com/bactopia/bactopia)** is a flexible multi-species pipeline for bacterial genomes. Supports raw reads. Uses Nextflow.
- **[Nullarbor](https://github.com/tseemann/nullarbor)** is a read-to-report pipeline for public health microbiology. Excellent for outbreak surveillance.
- **[Pathogenwatch](https://pathogen.watch)** is a web-based platform for genomic surveillance. Upload-based, very user-friendly.

Each has its own strengths. PseudoScope focuses on speed, local execution, and *P. aeruginosa*-specific features.

---

## 🌐 **Part of the ESKAPE Web Tools Suite**

PseudoScope is one of a growing collection of species-optimized pipelines for the ESKAPE pathogens plus *E. coli*: *Enterococcus faecium*, *Staphylococcus aureus*, *Klebsiella pneumoniae*, *Acinetobacter baumannii*, *Pseudomonas aeruginosa*, and *Enterobacter* species. Other tools in development:

- **StaphScope** for *S. aureus* (MLST, spa, SCCmec, agr, capsule)
- **Kleboscope** for *K. pneumoniae* (MLST, K/O typing, AMR)
- **AcinetoScope** for *A. baumannii* (MLST, OCL typing)
- **EnteroMark** for *Enterococcus* (van typing, MLST)
- **EcoliTyper** for *E. coli* (serotyping, MLST)

A unified web interface is live at [https://eskape.bio](https://eskape.bio). Upload a FASTA, pick a species, run the analysis without installing anything locally.

---

## 🔮 **Future Development**

Planned features for 2026 and 2027:

- **PSEUDA** (browser-based genome retrieval for *P. aeruginosa*). Filter by MLST, O-type, carbapenemase profile, AMR genes, geography, isolation source. Click, filter, export. [bbeckley-hub.github.io/pseuda](https://bbeckley-hub.github.io/pseuda/)
- **Raw read support.** Direct FASTQ analysis with integrated assembly.
- **Lineage database.** Curated reference panel with HA-MRSA/CA-MRSA equivalents for *P. aeruginosa* epidemic lineages.
- **Machine learning module.** Outbreak prediction and phenotype inference.
- **Plugin system.** Community-contributed analysis modules.
- **Expanded ESKAPE coverage.** PseudoScope's architecture ported to other members of the ESKAPE group.

**Call for collaborators.** If you're interested in contributing to any of these, open an issue or contact the authors.

---

## ❓ **FAQ**

**Is PseudoScope free?**
Yes. MIT licensed. Free for academic, clinical, and commercial use.

**What makes PseudoScope different from generalist tools?**
Species-specific gene sets. Carbapenemase detection tuned for *P. aeruginosa*. O-antigen serotyping via PAST. Efflux pump and aminoglycoside methyltransferase detection built in. Report-level control for large batches.

**Can I use PseudoScope for clinical diagnosis?**
No. It is a research tool. Results should be confirmed with orthogonal methods before any clinical decision.

**Which version should I use, CLI or Web?**
Web for convenience and small batches. CLI for large batches, pipeline integration, or when genomes must stay strictly local.

**Does PseudoScope work with short-read and long-read assemblies?**
Yes. Any FASTA works. Assembly quality is assessed by the QC module and shown in every report.

**Why is MGE the slowest module?**
mobileOG-db is 350 MB of protein families. DIAMOND is fast, but the database is large. Use `--skip-mge` if MGE is not the question.

**What is the difference between `--pull-mlst-db` and `--update-mlst-db`?**
`--pull-mlst-db` downloads a pre-built snapshot from GitHub. No API key. Fast.
`--update-mlst-db` fetches the current Pasteur scheme from PubMLST via mlstdb. Requires an API key. Slower, but guarantees you have the newest allele variants.

---

## 🐛 **Troubleshooting**

```bash
# AMR database missing or outdated
pseudoscope --update-amr-db
pseudoscope --force-update-amr-db   # full overwrite

# ABRicate databases not found
abricate --setupdb

# MLST scheme missing
pseudoscope --pull-mlst-db

# Refresh MLST from PubMLST (API key required)
pseudoscope --update-mlst-db

# mobileOG-db missing (conda users)
pseudoscope --pull-mge-db

# Docker permission errors
mkdir -p "$PWD/results"   # create as your user first
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" ...

# "master TSV not found" in visualisation
# The pipeline must run through the gene-centric reporter. Do not use --skip-summary.
```

### Getting help

1. Check existing issues: [GitHub Issues](https://github.com/bbeckley-hub/pseudoscope/issues)
2. Open a new issue with: the full error, `conda list` output, the failing command, and `pseudoscope_run.log`
3. Email: brownbeckley94@gmail.com

---

## 🙏 **Acknowledgements**

PseudoScope is built on tools that other people spent years on.

- **Torsten Seemann** for MLST, ABRicate, and the wider ecosystem.
- **Robert Petit and Tim Read** for pasty and camlhmp.
- **The NCBI AMR team** for AMRFinderPlus.
- **The PubMLST team** (Keith Jolley and Martin Maiden) for the Pasteur scheme.
- **The CGE team** for ResFinder and PlasmidFinder database structures.
- **All database maintainers** of CARD, VFDB, ARG-ANNOT, MEGARes, BacMet2, EcoH, Ecoli_VF.
- **ParBLiSS** for fastANI.
- **Connor L. Brown and James I. Mullet** for mobileOG-db and the parsing logic.
- **Doug Hyatt** for Prodigal.
- **Benjamin Buchfink** for DIAMOND.
- **The Chart.js and Plotly teams** for interactive visualisation support.
- **Early adopters and beta testers** at the University of Ghana Medical School.

> "If we ever meet in person, the drinks are on me." — Brown Beckley

---

## 💰 **Sponsor This Project**

PseudoScope is developed and maintained in our spare time. If you find it useful for your research or clinical work, please consider sponsoring.

- **GitHub Sponsors:** [https://github.com/sponsors/bbeckley-hub](https://github.com/sponsors/bbeckley-hub)

Sponsorships go towards:

- Adding raw read support faster
- Hosting a free public web interface
- Providing dedicated user support
- Keeping databases updated

Every contribution, no matter how small, helps keep the project alive. Thank you.

---

## 📚 **Citation**

A peer-reviewed publication for PseudoScope is in preparation. Until then, please cite the repository:

```bibtex
@software{beckley2026pseudoscope,
  author = {Brown Beckley and Vincent Amarh},
  title = {PseudoScope: a species-specific bioinformatics suite for rapid and accessible Pseudomonas aeruginosa genomic analysis},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/pseudoscope},
  version = {1.3.0}
}
```

**Integrated tools.** Please also cite the original authors of MLST, ABRicate, pasty, AMRFinderPlus, PubMLST, CARD, ResFinder, VFDB, and others. Full list below.

---

## 🤖 **AI Disclosure**

PseudoScope is developed with the assistance of AI tools, and we are transparent about it.

**How we use AI:**

- **Code drafting and refactoring.** AI assistants help generate initial implementations, suggest refactors, and catch patterns we might miss.
- **Documentation.** Some prose in this README, module docstrings, and HTML report text has been drafted with AI assistance and edited by human authors.
- **Testing scaffolds.** AI helps write boilerplate tests and mock fixtures.

**What we do not do:**

- **We do not ship code we don't understand.** Every line that lands in the repository is reviewed by at least one human author.
- **We do not ship code we can't debug.** If a bug appears in code you can't reason about, that's a bug in the process, not a feature.
- **We do not use AI for scientific claims.** Biological interpretations, thresholds, and validation calls are made by human authors based on peer-reviewed literature and reference strains.
- **We do not train models on your data.** PseudoScope is a local tool. Nothing runs through a third-party service. Your genomes stay on your machine.

**Responsibility:** The authors of PseudoScope take full responsibility for the correctness, behaviour, and safety of the released code. AI is a productivity tool in our hands, not an author on our papers. If something in PseudoScope is wrong, the fault is ours.

---

## 📚 **Third-Party Tool Citations**

PseudoScope integrates several third-party tools and databases. If you use PseudoScope in your research, please also cite the originals.

<details>
<summary><b>Full citation list (click to expand)</b></summary>

#### MLST
```bibtex
@software{seemann_mlst_2018,
  author = {Seemann, T.},
  title = {MLST: Scan contig files against traditional PubMLST typing schemes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/mlst}
}
```

#### PubMLST
```bibtex
@article{jolley_pubmlst_2018,
  author = {Jolley, K. A. and Bray, J. E. and Maiden, M. C. J.},
  title = {Open-access bacterial population genomics: {BIGSdb} software, the {PubMLST.org} website and their applications},
  journal = {Wellcome Open Research},
  volume = {3},
  pages = {124},
  year = {2018},
  doi = {10.12688/wellcomeopenres.14826.1}
}
```

#### ABRicate
```bibtex
@software{seemann_abricate_2018,
  author = {Seemann, T.},
  title = {ABRicate: Mass screening of contigs for antimicrobial resistance and virulence genes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/abricate}
}
```

#### AMRFinderPlus
```bibtex
@article{feldgarden_amrfinderplus_2021,
  author = {Feldgarden, M. et al.},
  title = {AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence},
  journal = {Scientific Reports},
  volume = {11},
  pages = {12728},
  year = {2021},
  doi = {10.1038/s41598-021-91456-0}
}
```

#### pasty
```bibtex
@software{petit_pasty_2024,
  author = {Petit, R. A.},
  title = {PASTY: In silico serogrouping of Pseudomonas aeruginosa},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/rpetit3/pasty}
}
```

#### camlhmp
```bibtex
@software{petit_camlhmp_2024,
  author = {Petit, R. A.},
  title = {camlhmp: Classification through yAML Heuristic Mapping Protocol},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/rpetit3/camlhmp}
}
```

#### fastANI
```bibtex
@article{jain_fastani_2018,
  author = {Jain, C. and Rodriguez-R, L. M. and Phillippy, A. M. and Konstantinidis, K. T. and Aluru, S.},
  title = {High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries},
  journal = {Nature Communications},
  volume = {9},
  number = {1},
  pages = {5114},
  year = {2018},
  doi = {10.1038/s41467-018-07641-9}
}
```

#### mobileOG-db
```bibtex
@article{brown_mobileogdb_2022,
  author = {Brown, C. L. and Mullet, J. and Hindi, F. and Stoll, J. E. and Gupta, S. and Choi, M. and Keenum, I. and Vikesland, P. and Pruden, A. and Zhang, L.},
  title = {mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements},
  journal = {Applied and Environmental Microbiology},
  volume = {88},
  number = {18},
  pages = {e00991-22},
  year = {2022},
  doi = {10.1128/aem.00991-22}
}
```

#### Prodigal
```bibtex
@article{hyatt_prodigal_2010,
  author = {Hyatt, D. and Chen, G.-L. and LoCascio, P. F. and Land, M. L. and Larimer, F. W. and Hauser, L. J.},
  title = {Prodigal: prokaryotic gene recognition and translation initiation site identification},
  journal = {BMC Bioinformatics},
  volume = {11},
  pages = {119},
  year = {2010},
  doi = {10.1186/1471-2105-11-119}
}
```

#### DIAMOND
```bibtex
@article{buchfink_diamond_2021,
  author = {Buchfink, B. and Reuter, K. and Drost, H.-G.},
  title = {Sensitive protein alignments at tree-of-life scale using DIAMOND},
  journal = {Nature Methods},
  volume = {18},
  pages = {366-368},
  year = {2021},
  doi = {10.1038/s41592-021-01101-x}
}
```

#### CARD
```bibtex
@article{alcock_card_2023,
  author = {Alcock, B. P. et al.},
  title = {CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database},
  journal = {Nucleic Acids Research},
  volume = {51},
  number = {D1},
  pages = {D690-D699},
  year = {2023},
  doi = {10.1093/nar/gkac920}
}
```

#### ResFinder
```bibtex
@article{bortolaia_resfinder_2020,
  author = {Bortolaia, V. et al.},
  title = {ResFinder 4.0 for predictions of phenotypes from genotypes},
  journal = {Journal of Antimicrobial Chemotherapy},
  volume = {75},
  number = {12},
  pages = {3491-3500},
  year = {2020},
  doi = {10.1093/jac/dkaa345}
}
```

#### VFDB
```bibtex
@article{chen_vfdb_2016,
  author = {Chen, L. et al.},
  title = {VFDB 2016: hierarchical and refined dataset for big data analysis - 10 years on},
  journal = {Nucleic Acids Research},
  volume = {44},
  number = {D1},
  pages = {D694-D697},
  year = {2016},
  doi = {10.1093/nar/gkv1239}
}
```

#### PlasmidFinder
```bibtex
@article{carattoli_plasmidfinder_2014,
  author = {Carattoli, A. et al.},
  title = {In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing},
  journal = {Antimicrobial Agents and Chemotherapy},
  volume = {58},
  number = {7},
  pages = {3895-3903},
  year = {2014},
  doi = {10.1128/AAC.02412-14}
}
```

#### BacMet2
```bibtex
@article{pal_bacmet_2014,
  author = {Pal, C. et al.},
  title = {BacMet: antibacterial biocide and metal resistance genes database},
  journal = {Nucleic Acids Research},
  volume = {42},
  pages = {D737-D743},
  year = {2014},
  doi = {10.1093/nar/gkt1252}
}
```

#### MEGARes
```bibtex
@article{doster_megares_2020,
  author = {Doster, E. et al.},
  title = {MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data},
  journal = {Nucleic Acids Research},
  volume = {48},
  pages = {D561-D569},
  year = {2020},
  doi = {10.1093/nar/gkz1010}
}
```

#### Biopython
```bibtex
@article{cock_biopython_2009,
  author = {Cock, P. J. A. et al.},
  title = {Biopython: freely available Python tools for computational molecular biology and bioinformatics},
  journal = {Bioinformatics},
  volume = {25},
  number = {11},
  pages = {1422-1423},
  year = {2009},
  doi = {10.1093/bioinformatics/btp163}
}
```

</details>

---

## 👥 **Authors & Contact**

**Brown Beckley** (Primary Developer)
University of Ghana Medical School, Department of Medical Biochemistry
📧 brownbeckley94@gmail.com
🐙 GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)
LinkedIn: [@brownbeckley](https://www.linkedin.com/in/brown-beckley-190315319/)

**Vincent Amarh** (Co-author, Advisor)
University of Ghana Medical School, Department of Medical Biochemistry

For bug reports, feature requests, or collaborations, please use the [GitHub issue tracker](https://github.com/bbeckley-hub/pseudoscope/issues).

---

## 📄 **License**

### Core PseudoScope code

The PseudoScope pipeline code (orchestrator, module scripts, HTML templates, and Python modules) is MIT licensed. See [LICENSE](LICENSE).

### Third-party tools

PseudoScope invokes external tools installed as Conda dependencies. Each tool is the property of its respective developers and is used under its own licence. By using PseudoScope you agree to comply with those licences.

---

<div align="center">

## 🚀 **Ready to type your P. aeruginosa genomes?**

| **Choose your platform** | |
|---|---|
| 🖥️ **Command line** | For high-throughput, local analysis |
| 🌐 **PseudoScope Web** | For non-bioinformaticians at [https://eskape.bio](https://eskape.bio) |
| 🧬 **PSEUDA** (in development) | Genome retrieval for *P. aeruginosa* at [bbeckley-hub.github.io/pseuda](https://bbeckley-hub.github.io/pseuda/) |

[![Get Started](https://img.shields.io/badge/GET_STARTED-Now-green?style=for-the-badge&logo=conda)](https://github.com/bbeckley-hub/pseudoscope#-quick-start)
[![Try Web Version](https://img.shields.io/badge/TRY_WEB_VERSION-Here-blue?style=for-the-badge&logo=html5)](https://eskape.bio)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/pseudoscope/issues)
[![Sponsor](https://img.shields.io/badge/Sponsor-❤️-pink?style=for-the-badge&logo=githubsponsors)](https://github.com/sponsors/bbeckley-hub)

**From days to minutes. From fragmented to integrated. From raw assemblies to actionable surveillance.**

*PseudoScope: precision surveillance for the carbapenemase era.*

⭐ **If you find this tool useful, please star the repository.** ⭐

*P. aeruginosa is one of the toughest pathogens in the ESKAPE group. Carbapenem resistance in this species is not one gene or one mechanism. It is a set of overlapping adaptations that differ by lineage and by geography. We invite researchers, clinical labs, and public health teams to collaborate on validating the pipeline, sharing regional data, and expanding the PSEUDA dataset as it comes online.*

**Together, we can make P. aeruginosa surveillance routine.**

</div>
