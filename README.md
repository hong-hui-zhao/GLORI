

---

````markdown
# GLORI: Reproducing the Analysis Pipeline

This repository aims to **reproduce and understand the analysis pipeline of GLORI** (Glyoxal and Nitrite-based RNA sequencing), a transcriptome-wide technology for **absolute quantification of m⁶A at single-base resolution**.  
The reproduction is performed by re-implementing the computational steps based on the official [`GLORI-tools`](https://github.com/liucongcas/GLORI-tools)and the descriptions in the original paper.

---

## 📘 Project Objective

The goal of this project is to:
1. Reproduce the data processing and analysis steps described in the GLORI paper.
2. Write and document custom code to perform each step manually, instead of directly running the provided pipeline.
3. Understand the computational logic behind **m⁶A site calling**, **A-rate calculation**, and **filtering of false positives**.

This serves both as a **learning project** and a **validation of reproducibility** for the GLORI method.

---

## 🧬 Background

**GLORI (Glyoxal and Nitrite Reaction-based Sequencing)** is a chemical RNA-seq approach that converts unmethylated adenosines to inosines (read as G), while methylated adenosines (m⁶A) remain as A.  
By comparing the **A/G ratio** at each position, one can infer the **m⁶A methylation level** at single-base resolution.

**Reference paper:**  
> Wang et al., *Nature Biotechnology*, 2023.  
> “Absolute quantification of m6A at single-base resolution using GLORI-seq.”

---

## ⚙️ Pipeline Overview

Below is the step-by-step workflow that will be reproduced in this repository:

1. **Preprocessing raw reads**
   - Adapter trimming and quality filtering using **Trim Galore** (v0.6.6).
   - Removal of PCR duplicates based on UMI sequences.
   - Stripping UMIs from the 5′ end to obtain clean reads.

2. **Read alignment**
   - Align reads to the reference genome (e.g., `hg38`) using **STAR**.
   - Generate BAM files for downstream analysis.

3. **Site calling and A-rate calculation**
   - For each A site in the transcriptome, calculate:
     \[
     A\text{-rate} = \frac{\text{Number of A reads}}{\text{Total reads covering the site}}
     \]
   - Compare A-rates between treated and control samples.

4. **Filtering false positives**
   - Remove sites affected by sequencing or reverse transcription bias.
   - Apply thresholds for coverage depth and background noise.

5. **m⁶A level estimation and visualization**
   - Derive the methylation level from the A-rate.
   - Visualize results (metagene plots, distribution of A-rates, site overlap with known m⁶A peaks).

---

## 🧰 Requirements

| Tool | Version | Description |
|------|----------|-------------|
| Trim Galore | ≥0.6.6 | Adapter and quality trimming |
| STAR | ≥2.7 | Read alignment |
| Samtools | ≥1.9 | BAM file processing |
| Bedtools | ≥2.29 | Genomic coordinate operations |
| Python | ≥3.8 | Data analysis scripts |
| R | ≥4.1 | Statistical analysis and visualization |

Dependencies will be managed through `conda`:
```bash
conda env create -f environment.yml
conda activate glori
````

---

## 🧪 Data

The datasets are obtained from the **GLORI-tools** repository and corresponding GEO entries.
For example:

```bash
wget https://ftp.ncbi.nlm.nih.gov/.../GCF_000001405.39_GRCh38.p13_assembly_report.txt
```

> Note: Due to data size and copyright constraints, raw FASTQ files are not included here.

---

## 🧠 Reproduction Strategy

The analysis is implemented **step-by-step**, each with detailed comments and explanations:

| Step | Script                | Description                                     |
| ---- | --------------------- | ----------------------------------------------- |
| 01   | `preprocess_reads.sh` | Trimming, deduplication, and UMI removal        |
| 02   | `align_reads.sh`      | Read mapping with STAR                          |
| 03   | `calculate_A_rate.py` | Parse BAM to compute A-rate                     |
| 04   | `filter_sites.R`      | Apply filters and generate candidate m⁶A sites  |
| 05   | `visualize_results.R` | Plot site distribution and methylation profiles |

Each script is self-contained, allowing readers to trace the logic behind the GLORI analysis rather than simply rerunning existing tools.

---

## 📊 Expected Outputs

* Cleaned reads (FASTQ)
* Aligned reads (BAM)
* Site-level A-rate tables (CSV/TSV)
* Filtered m⁶A candidate list
* Visualization plots:

  * A-rate distribution
  * m⁶A site density across transcripts
  * Example locus-level coverage plots

---

## 🧩 Repository Structure

```
GLORI-reproduce/
│
├── data/                # input data, genome files, annotation
├── scripts/             # shell and python scripts for each step
├── results/             # output tables and figures
├── environment.yml      # conda environment file
└── README.md            # this document
```

---

## 📚 References

* Wang et al., *Nature Biotechnology*, 2023 — GLORI original paper
* [GLORI-tools (official repo)](https://github.com/sxjscience/GLORI-tools)
* Trim Galore documentation: [https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* STAR manual: [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)

---

## ✍️ Author Notes

This repository is maintained by **[Your Name]**, as part of a project to deeply understand and reproduce the computational methodology of GLORI.
All scripts are independently written while cross-referencing the official implementation to ensure transparency and reproducibility.

---

## 📄 License

This project is open-sourced under the **MIT License**.

---

## 🚧 Status

**In progress.**
Currently reproducing the preprocessing and alignment steps (Trim Galore → STAR).
Subsequent scripts for A-rate calculation and filtering will be added progressively.

---

> 🧠 *Goal:* “Not just to reproduce results — but to understand every line that makes them possible.”

```

---

Would you like me to make this **fit your actual repo structure** (e.g., you already have `/scripts`, `/data`, `/results` folders and some `.py` or `.sh` files)? I can tailor the `README` to your current setup and fill in command examples automatically.
```
