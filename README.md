# SingleCell-PRISM: AL Amyloidosis and Daratumumab Response

## ✨ Overview

Systemic light-chain amyloidosis (AL) is driven by clonal plasma cells that produce misfolded immunoglobulin light chains. Although daratumumab-based therapy has improved treatment outcomes, a subset of patients does not achieve an optimal hematologic response.

The analyses identify plasma-cell-intrinsic programs and inflammatory–immunosuppressive interactions associated with suboptimal response to daratumumab-based therapy. The analytical workflow examines:

1. **Single-cell cohort assembly, preprocessing, and quality control**
2. **Amyloidogenic plasma-cell identification**
3. **Plasma-cell gene-expression programs and functional states**
4. **B-cell receptor clonality and longitudinal clonal evolution**
5. **Plasma-cell developmental potential and copy-number alterations**
6. **Bone-marrow immune-cell composition and functional states**
7. **Cell–cell communication and niche-associated signaling**
8. **Validation in independent single-cell and bulk-expression cohorts**

<img width="191" height="191" alt="image" src="https://github.com/user-attachments/assets/8fef46f0-c1e6-433a-bf15-4e7a50ae077f" />

## 📁 Notebook Descriptions

Assembly of sample-level 10x matrices into cohort objects is omitted because the source data include protected health information.

| Notebook | Description |
|---|---|
| `1.1_doubletsambientRNA.ipynb` | Quality control and preprocessing by doublets and ambient RNA removal. |
| `1.2_standardpipelines.ipynb` | Adjusted embedding, clustering, and initial annotation with marker genes. |
| `2.1_plasmacellcloneVDJ.ipynb` | Identification of amyloidogenic and polyclonal plasma-cell clones using transcriptomic and V(D)J information. |
| `2.2_cNMF_plasmacellsGEP.ipynb` | Inference of plasma-cell gene-expression programs using cNMF. |
| `2.3_starCAT_plasmacellsGEP.ipynb` | Projection of discrete and continuous plasma-cell programs with starCAT. |
| `2.4_plasmacellclone-starCAT-part1.ipynb` | Sample-level analysis of starCAT-projected plasma-cell programs: post/pre-treatment distributions and pretreatment Good/Poor response comparisons across all patients and t(11;14) strata. |
| `2.4_plasmacellclone-starCAT-part2.ipynb` | Clinical annotation of sample-level starCAT usages, cytogenetic-subgroup response comparisons, and export of subgroup figures. |
| `2.5_Functional-string.ipynb` | Recurrent differential-expression analysis and functional interpretation of amyloidogenic plasma-cell genes. |
| `2.6_Functional-cNMF.ipynb` | Functional enrichment and interpretation of cNMF programs. |
| `2.7_bulk_signature_GSVA_updated.ipynb` | Evaluation of single-cell-derived signatures in an external bulk RNA-seq cohort using GSVA. |
| `3.1_cytotrace_cloneGEP5.ipynb` | CytoTRACE analysis of plasma-cell developmental potential and gene-expression programs. |
| `3.2_ALl9_ALk3_ALl5-clonetrack.ipynb` | Clone2vec-based analysis of pre/post-treatment plasma-cell clonal states. |
| `3.3_dandelion-LCMS.ipynb` | Dandelion repertoire analysis and light-chain mass-spectrometry validation of clonotypes. |
| `3.4_dandelion-allclones.ipynb` | Cohort-wide B-cell receptor clonotype and repertoire analysis. |
| `3.5_infercnvpy_plasma.ipynb` | Inference and visualization of plasma-cell copy-number alterations. |
| `4.0_BM-TME-singleR.ipynb` | Selection and reference-based annotation of the non-plasma-cell bone-marrow microenvironment. |
| `4.1_BM-TME-singleR-percelltype.ipynb` | Per-lineage refinement and visualization of SingleR annotation results. |
| `4.2_BM-TME-azimuth-myeloid.ipynb` | Detailed annotation and characterization of myeloid-cell populations. |
| `4.3_BM-TME-T cells.ipynb` | Identification and characterization of T-cell subpopulations. |
| `4.4_starCAT_Tcells.ipynb` | Projection of continuous T-cell functional programs with starCAT. |
| `4.5_BM-TME-NK cells.ipynb` | Identification and characterization of NK-cell subpopulations. |
| `4.6_BM-TME-B cells.ipynb` | Identification and characterization of non-plasma B-cell populations. |
| `4.7_upsetplot_interferon.ipynb` | Comparison of interferon-response gene sets across immune-cell populations. |
| `5.1_cellchat-78samples-trial3_annotations.ipynb` | Customized, sample-aware cell–cell communication analysis using CellChat. |
| `5.2_cellchat-ALNicheGraph.ipynb` | Construction and querying of a signaling-to-phenotype AL niche graph from saved CellChat and enrichment results. |
| `5.3_CD38_publication_figures.ipynb` | Generation of CD38-centered publication figures and cell-type comparisons. |
| `6.0_GSE175385_Alameda_preprocess.ipynb` | Preparation of the GSE175385 Alameda validation dataset. |
| `6.1_plasmacell-Alameda.ipynb` | Validation of disease-associated plasma-cell programs in the Alameda cohort. |
| `6.2_starCAT_plasmacells-Alameda.ipynb` | Projection of reference plasma-cell states in the Alameda cohort. |
| `7.0_GSE292189_Gort_preprocessing.ipynb` | Assembly and harmonization of the GSE292189 validation dataset. |
| `7.1_plasmacell_Gort.ipynb` | Validation of amyloidogenic and polyclonal plasma-cell states in the GSE292189 cohort. |
| `7.2_BM-TME-azimuth-myleoid_Gort.ipynb` | Validation of myeloid-cell annotations in the GSE292189 cohort. |
| `7.3_starCAT_plasmacells-Gort.ipynb` | Projection of reference plasma-cell states in the GSE292189 cohort. |
| `8.1_BM-TME-azimuth-myeloid_Safina.ipynb` | Validation of myeloid-cell states in an additional external cohort. |

### Plasma-cell program workflow (2.2–2.4)

1. Run `2.2_cNMF_plasmacellsGEP.ipynb` to infer the reference programs.
2. Run `2.3_starCAT_plasmacellsGEP.ipynb` to project those reference programs onto query cells with starCAT and export `plasma_cnmf_3_usage_norm_withoutIG_k=13_allplasmacell.csv`.
3. Run either or both `2.4_plasmacellclone-starCAT-part1.ipynb` and `2.4_plasmacellclone-starCAT-part2.ipynb` for downstream R analyses. Each notebook prepares its own sample summaries and clinical response table; neither requires objects left in memory by the other notebook.

The `cnmf` substring in the usage CSV is the historical reference-program filename. The usages consumed by both Part 1 and Part 2 come from starCAT in notebook 2.3.

Both Part 1 and Part 2 require the cohort Seurat object `250817-AL-all78umap-singletonly-umap.rds`, the additional plasma-cell barcode list `sce1_ALl19clustermeta.data.csv`, the starCAT usage CSV above, and `sampleannotation-all.csv` with at least `samplename`, `bestorlast.response`, and `mut`. Supply these files in the working directory and configure the reticulate Python/Conda paths for the local environment. These cohort and clinical input files are not included in this repository.

- **Part 1 outputs:** statistical CSVs and summary PDFs in `GEP_prepost_simplified/` and `GEP_pretreatment_GoodPoor/`. Tests compare sample-level distributions with unpaired Wilcoxon tests and Benjamini–Hochberg FDR within each analysis stratum.
- **Part 2 outputs:** annotated heatmaps and nine-program subgroup comparison plots inline, plus `PlotB_*`, `PlotC_*`, and `PlotD_*` PDFs for the seven-program subset (excluding GEP3 and GEP7) in `usage_plots_by_mut/`. The exploratory subgroup figures display nominal Wilcoxon p-values with at least three samples per comparison arm and Cliff's delta for pre/post comparisons when estimable.

## 🧰 Major Software Dependencies

The repository contains both R- and Python-based analyses. Package requirements vary by notebook, but the major dependencies are summarized below.

### R workflows

- **Single-cell analysis:** Seurat, SingleCellExperiment, DropletUtils, SoupX, DoubletFinder/scDblFinder, celda, harmony, and scCustomize
- **R/Python interoperability:** reticulate
- **Reference-based cell annotation:** SingleR, celldex, Azimuth, and sceasy
- **Cell–cell communication and networks:** CellChat and igraph
- **Gene-set and differential-expression analysis:** GSVA, clusterProfiler, msigdbr, limma, edgeR, speckle, org.Hs.eg.db, and AnnotationDbi
- **Data manipulation and visualization:** dplyr, tidyr, readr, data.table, stringr, purrr, ggplot2, patchwork, pheatmap, ComplexHeatmap, ggpubr, scales, broom, effsize, and ComplexUpset

### Python workflows

- **Single-cell data structures and analysis:** Scanpy and AnnData
- **Plasma-cell programs and state projection:** cNMF and starCAT
- **Clonotype and repertoire analysis:** Dandelion and clone2vec
- **Developmental-potential and trajectory analysis:** CytoTRACE/CellRank and scVelo
- **Copy-number analysis:** infercnvpy
- **Scientific computing and visualization:** pandas, NumPy, SciPy, statsmodels, Matplotlib, and seaborn

## 📊 Reproducibility

Exact reproduction may require:

- Raw or processed single-cell datasets
- Sample-level clinical and treatment-response metadata
- Paired single-cell VDJ sequencing data
- External validation datasets
- Gene-signature and pathway-resource files
- Software and package versions consistent with the original computational environment

## 📚 Related Publications

When using this repository or its analytical workflow, please cite the associated preprint:

1. **Wang X, Xiong X, Xu C, Han H, Guan A, Gao Y, Shen K, Li J.** Single-Cell Analysis Reveals Inflammatory–Immunosuppressive Niches Associated with Suboptimal Response to Daratumumab-Based Therapy in AL Amyloidosis. *Nature Communications* (Accepted, In press) *medRxiv*. 2026.(https://www.medrxiv.org/content/10.64898/2026.03.28.26349317v2))
2. **Single-cell analysis of plasma cell–centered inflammatory–immunosuppressive niches and response to daratumumab in primary light-chain amyloidosis.** *Journal of Clinical Oncology*. 2026;44(suppl 16):7553. (Poster, American Society of Clinical Oncology Annual Meeting 2026)
3. **Single-cell analysis of BCR and transcriptome uncovers plasma cell clonal dynamics and immune exhaustion linked to daratumumab response in primary light-chain amyloidosis.** *Blood*. 2025;146:2158–2159. (Poster, American Society of Hematology Annual Meeting 2025; Young Investigator Award by CAHON)
4. **MM-913: Defining Gene Expression Programs of Amyloidogenic Clonal Plasma Cells Characterized by Single-Cell B-Cell Receptor and Transcriptomic Profiling.** *Clinical Lymphoma, Myeloma & Leukemia*. 2025;25(suppl 1):S951. (Poster, Society of Hematologic Oncology 2025)


## 📄 License

This project is licensed under the [MIT License](LICENSE).
