# SingleCell_PRISM_Daratumumab

# SingleCell-PRISM: Primary AL Amyloidosis and Daratumumab Response

## ⚠️ Disclaimer

The accompanying manuscript is a preprint and has not yet undergone peer review. The repository is provided for research and reproducibility purposes only.

The code and results should not be used independently to guide clinical diagnosis, prognosis, or treatment decisions.

> Single-cell profiling of plasma-cell evolution and the bone marrow immune microenvironment in primary light-chain amyloidosis during daratumumab-based treatment

This repository contains the analysis code accompanying:

**Single-Cell Analysis Reveals Inflammatory–Immunosuppressive Niches in Daratumumab-Resistant Primary AL Amyloidosis**

The study integrates single-cell transcriptomics, paired B-cell receptor sequencing, plasma-cell clonal analysis, gene-expression program inference, immune-cell state characterization, and cell–cell communication analysis to investigate the biological mechanisms associated with suboptimal response to daratumumab-based treatment in primary light-chain amyloidosis.

## ✨ Overview

Primary light-chain amyloidosis, also known as primary AL amyloidosis or pAL, is driven by clonal plasma cells that produce misfolded immunoglobulin light chains. Although daratumumab-based therapy has improved treatment outcomes, a subset of patients does not achieve an optimal hematologic response.

This study constructs a single-cell bone marrow atlas from patients with primary AL amyloidosis treated with daratumumab–bortezomib–dexamethasone. The analytical workflow examines:

1. **Single-cell data preprocessing and quality control**
2. **Amyloidogenic plasma-cell identification**
3. **Plasma-cell gene-expression programs**
4. **B-cell receptor clonality and clonal evolution**
5. **Copy-number alterations in plasma cells**
6. **Bone marrow immune-cell composition and functional states**
7. **Myeloid, T-cell, NK-cell, and B-cell subpopulations**
8. **Cell–cell communication within the bone marrow microenvironment**

The analyses identify plasma-cell-intrinsic programs and inflammatory–immunosuppressive interactions associated with suboptimal response to daratumumab-based therapy.

## 🔬 Analysis Modules

The notebooks are organized into five principal analytical modules.

### 1. Data preprocessing and quality control

The initial workflow performs sample-level quality control, filtering, normalization, dimensionality reduction, integration, clustering, and broad cell-type annotation across the single-cell dataset.

### 2. Plasma-cell programs and functional states

Plasma-cell analyses include:

* Identification of amyloidogenic clonal plasma cells using paired B-cell receptor information
* Consensus non-negative matrix factorization
* Projection and quantification of plasma-cell gene-expression programs
* Differential-expression and functional-enrichment analyses
* Evaluation of transcriptional signatures in external bulk datasets

### 3. Plasma-cell clonal evolution

Longitudinal plasma-cell analyses include:

* Reconstruction of clonal relationships with CoSpar
* B-cell receptor repertoire analysis with Dandelion
* Validation of clonotypes using light-chain mass-spectrometry information
* Inference of plasma-cell copy-number alterations

### 4. Bone marrow immune microenvironment

The immune-microenvironment workflow characterizes:

* Broad bone marrow cell populations
* Myeloid-cell states
* T-cell states and functional programs
* NK-cell states
* B-cell states
* Treatment- and response-associated differences in immune composition and activity

### 5. Cell–cell communication

Customized CellChat analyses are used to infer ligand–receptor interactions among plasma cells, myeloid cells, T cells, NK cells, and other bone marrow populations.

These analyses focus particularly on inflammatory and immunosuppressive signaling associated with treatment response.

## 🗂 Repository Structure

```text
SingleCell_PRISM_pALDaratumumab/
├── README.md
│
├── 1.1_QC_standardpipeline_78samples.ipynb
│
├── 2.1_plasmacell_78samples-cloneBCR.ipynb
├── 2.2_cNMF_plasmacell_k=13.ipynb
├── 2.5_Functional-string-DEG.ipynb
├── 2.6_Functional-cNMF.ipynb
├── 2.7_bulk_signature_scores.ipynb
│
├── 3.1_CoSpar-ALI9.ipynb
├── 3.2_CoSpar-ALk3.ipynb
├── 3.3_CoSpar-ALl5.ipynb
├── 3.4_dandelion-allclones.ipynb
├── 3.5_dandelion-LCMSvalidation.ipynb
├── 3.6_infercnvpy_plasma.ipynb
│
├── 4.1_BM-TME-singleR.ipynb
├── 4.2_BM-TME-singleR-results.ipynb
├── 4.3_BM-TME-azimuth-myeloid.ipynb
├── 4.4_BM-TME-T cells.ipynb
├── 4.5_starCAT_TCAT_Tcells.ipynb
├── 4.6_BM-TME-NK cells.ipynb
├── 4.7_BM-TME-B cells.ipynb
│
└── 5.1_cellchat-custom_78samples.ipynb
```

## 📁 Notebook Description

| Notebook                                  | Description                                                                                                |
| ----------------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `1.1_QC_standardpipeline_78samples.ipynb` | Quality control, preprocessing, integration, clustering, and initial annotation of the single-cell dataset |
| `2.1_plasmacell_78samples-cloneBCR.ipynb` | Identification and characterization of plasma cells using transcriptomic and B-cell receptor information   |
| `2.2_cNMF_plasmacell_k=13.ipynb`          | Inference of plasma-cell gene-expression programs using consensus non-negative matrix factorization        |
| `2.5_Functional-string-DEG.ipynb`         | Differential-expression analysis and functional interpretation using STRING-associated resources           |
| `2.6_Functional-cNMF.ipynb`               | Functional characterization of plasma-cell cNMF programs                                                   |
| `2.7_bulk_signature_scores.ipynb`         | Evaluation of single-cell-derived transcriptional signatures in bulk-expression datasets                   |
| `3.1_CoSpar-ALI9.ipynb`                   | CoSpar-based plasma-cell state-transition analysis for patient ALI9                                        |
| `3.2_CoSpar-ALk3.ipynb`                   | CoSpar-based plasma-cell state-transition analysis for patient ALk3                                        |
| `3.3_CoSpar-ALl5.ipynb`                   | CoSpar-based plasma-cell state-transition analysis for patient ALl5                                        |
| `3.4_dandelion-allclones.ipynb`           | B-cell receptor clonotype and repertoire analysis using Dandelion                                          |
| `3.5_dandelion-LCMSvalidation.ipynb`      | Validation of plasma-cell clonotypes using light-chain mass-spectrometry information                       |
| `3.6_infercnvpy_plasma.ipynb`             | Inference and visualization of copy-number alterations in plasma cells                                     |
| `4.1_BM-TME-singleR.ipynb`                | Reference-based annotation of the bone marrow tumor microenvironment                                       |
| `4.2_BM-TME-singleR-results.ipynb`        | Downstream analysis and visualization of SingleR annotation results                                        |
| `4.3_BM-TME-azimuth-myeloid.ipynb`        | Detailed annotation and characterization of myeloid-cell populations                                       |
| `4.4_BM-TME-T cells.ipynb`                | Identification and characterization of T-cell subpopulations                                               |
| `4.5_starCAT_TCAT_Tcells.ipynb`           | T-cell functional-state analysis using starCAT and TCAT-associated programs                                |
| `4.6_BM-TME-NK cells.ipynb`               | Identification and characterization of NK-cell subpopulations                                              |
| `4.7_BM-TME-B cells.ipynb`                | Identification and characterization of non-plasma B-cell populations                                       |
| `5.1_cellchat-custom_78samples.ipynb`     | Customized cell–cell communication analysis using CellChat                                                 |

## 🚀 Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/Selecton98/SingleCell_PRISM_pALDaratumumab.git
cd SingleCell_PRISM_pALDaratumumab
```

### 2. Create a computational environment

The notebooks use both Python- and R-based single-cell analysis tools. Because different analytical modules may require different environments, users should install the packages required by the specific notebook they intend to run.

Frequently used frameworks include:

* Scanpy
* AnnData
* Seurat
* SingleR
* Azimuth
* cNMF
* CoSpar
* Dandelion
* infercnvpy
* CellChat
* starCAT
* TCAT
* pandas
* NumPy
* SciPy
* Matplotlib
* seaborn

A Python environment can be initialized with:

```bash
conda create -n singlecell-prism python=3.10
conda activate singlecell-prism

pip install jupyterlab scanpy anndata pandas numpy scipy matplotlib seaborn
```

Additional packages should be installed according to the corresponding notebook documentation and software requirements.

### 3. Start JupyterLab

```bash
jupyter lab
```

The notebooks are numbered according to their approximate position in the analysis workflow. They should generally be reviewed in numerical order, although some patient-level and cell-type-specific analyses can be run independently after the required processed objects have been generated.

## 💾 Data Availability

The repository contains analysis code but does not include patient-level raw sequencing data.

Access to controlled human sequencing data, processed objects, and relevant metadata should follow the data-availability statement and institutional requirements described in the associated manuscript.

Users reproducing the analyses should update the input and output paths in each notebook to reflect their local directory structure.

## 📊 Reproducibility

The notebooks document the principal computational analyses used in the accompanying study.

Exact reproduction may require:

* The corresponding raw or processed single-cell datasets
* Sample-level clinical and treatment-response metadata
* Paired single-cell B-cell receptor sequencing data
* External reference datasets
* Gene-signature and pathway-resource files
* Software versions consistent with the original computational environment

Because some notebooks were developed within the original institutional computing environment, users may need to modify file paths, environment-specific settings, and package versions before execution.

## 📚 Related Publications

This repository accompanies the following publications and presentations:

1. **Wang X, Xiong X, Han H, et al.** Single-Cell Analysis Reveals Inflammatory–Immunosuppressive Niches in Daratumumab-Resistant Primary AL Amyloidosis. *medRxiv*. 2026.
   DOI: [10.64898/2026.03.28.26349317](https://doi.org/10.64898/2026.03.28.26349317)

2. **Single-cell analysis of plasma cell–centered inflammatory–immunosuppressive niches and response to daratumumab in primary light-chain amyloidosis.** *Journal of Clinical Oncology*. 2026;44(suppl 16):7553.

3. **Single-cell analysis of BCR and transcriptome uncovers plasma cell clonal dynamics and immune exhaustion linked to daratumumab response in primary light-chain amyloidosis.** *Blood*. 2025;146:2158–2159.

4. **MM-913: Defining Gene Expression Programs of Amyloidogenic Clonal Plasma Cells Characterized by Single-Cell B-Cell Receptor and Transcriptomic Profiling.** *Clinical Lymphoma, Myeloma & Leukemia*. 2025;25(suppl 1):S951.

## 📖 Citation

When using this repository or its analytical workflow, please cite the associated preprint:

```bibtex
@article{wang2026singlecell,
  title   = {Single-Cell Analysis Reveals Inflammatory--Immunosuppressive Niches in Daratumumab-Resistant Primary AL Amyloidosis},
  author  = {Wang, Xuezhu and Xiong, Xinyi and Han, Hongxiao and Guan, Ai and Gao, Yajuan and Yan, Qi and Shen, Kaini and Li, Jian},
  journal = {medRxiv},
  year    = {2026},
  doi     = {10.64898/2026.03.28.26349317}
}
```

Because the manuscript is currently a preprint, users should check for an updated peer-reviewed citation before publication.

## 📄 License

This project is released under the **MIT License**. See the [`LICENSE`](LICENSE) file for details.

