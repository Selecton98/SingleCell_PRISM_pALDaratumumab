# Dependencies and recorded software versions

**Audit date:** 2026-11-01. **Source snapshot:** [2fbafac](https://github.com/Selecton98/SingleCell_PRISM_pALDaratumumab/tree/2fbafacab68426f0d52f457f6161cf45904acf7c).

This document inventories all **42 files** in that snapshot: **34 Jupyter notebooks, 4 R helper scripts, and 4 supporting files**. It collects versions from notebook kernel metadata, retained output, startup banners, source comments, and script headers. It also lists package imports and local helper dependencies for every analysis file.

The analysis notebooks and their outputs were read without re-executing them. Versions from the machine used to prepare this document were not substituted for repository evidence.

## How to interpret the records

- **Recorded runtime/package version:** printed by `sessionInfo()`, `reticulate::py_config()`, a package startup banner, or an explicit runtime diagnostic.
- **Kernel metadata:** the saved `metadata.language_info.version`; this is a notebook record, not proof that every retained cell output came from that same execution environment.
- **Source-comment version:** documented in code, but not independently confirmed by retained runtime output.
- **NR:** no exact version is recorded for that file/package. A version recorded in another notebook is not automatically assigned to it.
- Build warnings and API-change notices are separated below because their numbers are not necessarily installed package versions. Kernel display names such as “Python (cellpose)” and Conda environment names are not package version records.
- Cell references such as **c3** use **one-based notebook cell order, including Markdown cells**, rather than the mutable `In[...]` execution counter. **F01–F42** link to the per-file audit.

The repository records multiple environments. This inventory is **not a single tested environment lockfile**.

## R and Python runtime records

| Runtime | Version | Evidence | Files |
|---|---|---|---|
| Python | `3.10.0` | Notebook kernel metadata | [F14](#file-f14) |
| Python | `3.10.18` | Notebook kernel metadata | [F12](#file-f12), [F13](#file-f13) |
| Python | `3.8.20` | Notebook kernel metadata | [F04](#file-f04), [F05](#file-f05), [F11](#file-f11), [F15](#file-f15), [F20](#file-f20), [F29](#file-f29), [F33](#file-f33) |
| R | `4.2.1` | Notebook kernel metadata | [F10](#file-f10), [F24](#file-f24) |
| R | `4.3.3` | Notebook kernel metadata | [F01](#file-f01), [F02](#file-f02), [F03](#file-f03), [F06](#file-f06), [F07](#file-f07), [F08](#file-f08), [F09](#file-f09), [F16](#file-f16), [F17](#file-f17), [F18](#file-f18), [F19](#file-f19), [F21](#file-f21), [F22](#file-f22), [F23](#file-f23), [F25](#file-f25), [F26](#file-f26), [F27](#file-f27), [F28](#file-f28), [F30](#file-f30), [F31](#file-f31), [F32](#file-f32), [F34](#file-f34) |
| Python called from R | `3.8.19` | Retained `py_config()` output | [F03](#file-f03) c4, [F06](#file-f06) c4, [F07](#file-f07) c4, [F19](#file-f19) c4, [F28](#file-f28) c4, [F31](#file-f31) c4 |

R `4.3.3` is additionally confirmed by the full `sessionInfo()` output in [7.0_GSE292189_Gort_preprocessing.ipynb](7.0_GSE292189_Gort_preprocessing.ipynb) c3. The R `4.2.1` notebooks are the bulk GSVA workflow and the main CellChat workflow; their versions should not be silently replaced by the more common R `4.3.3` record.

## Additional package versions recorded outside the full session snapshot

| Software | Version | Evidence type | Sources |
|---|---|---|---|
| `circlize` | `0.4.16` | startup banner | [F03](#file-f03) c597 |
| `clusterProfiler` | `4.10.1` | startup banner | [F09](#file-f09) c5 |
| `ComplexHeatmap` | `2.18.0` | startup banner | [F03](#file-f03) c597 |
| `dendextend` | `1.18.1` | startup banner | [F17](#file-f17) c19 |
| `dplyr` | `1.1.4` | tidyverse startup banner | [F24](#file-f24) c6, [F26](#file-f26) c6 |
| `forcats` | `1.0.0` | tidyverse startup banner | [F03](#file-f03) c575, [F24](#file-f24) c6, [F26](#file-f26) c6 |
| `lubridate` | `1.9.4` | tidyverse startup banner | [F03](#file-f03) c575, [F24](#file-f24) c6, [F26](#file-f26) c6 |
| `matrixStats` | `1.3.0` | source comment (not independently verified) | [F10](#file-f10) c3 |
| `Modularity Optimizer` | `1.3.0` | startup banner | [F01](#file-f01) c15, [F01](#file-f01) c32, [F02](#file-f02) c4, [F03](#file-f03) c15, [F06](#file-f06) c13, [F07](#file-f07) c13, [F16](#file-f16) c15, [F17](#file-f17) c31, [F17](#file-f17) c37, [F17](#file-f17) c41, [F17](#file-f17) c47, [F17](#file-f17) c57, [F17](#file-f17) c61, [F17](#file-f17) c65, [F17](#file-f17) c75, [F17](#file-f17) c78, [F17](#file-f17) c89, [F17](#file-f17) c93, [F18](#file-f18) c5, [F28](#file-f28) c15, [F31](#file-f31) c20 |
| `NumPy` | `1.24.4` | reticulate runtime | [F03](#file-f03) c4, [F06](#file-f06) c4, [F07](#file-f07) c4, [F19](#file-f19) c4, [F28](#file-f28) c4, [F31](#file-f31) c4 |
| `purrr` | `1.1.0` | tidyverse startup banner | [F03](#file-f03) c575, [F26](#file-f26) c6 |
| `purrr` | `1.2.0` | tidyverse startup banner | [F24](#file-f24) c6 |
| `readr` | `2.1.5` | tidyverse startup banner | [F03](#file-f03) c575, [F24](#file-f24) c6, [F26](#file-f26) c6 |
| `scCustomize` | `2.1.2` | startup banner | [F01](#file-f01) c17, [F02](#file-f02) c6, [F03](#file-f03) c10, [F06](#file-f06) c7, [F07](#file-f07) c7, [F16](#file-f16) c3, [F17](#file-f17) c3, [F18](#file-f18) c3, [F19](#file-f19) c3, [F21](#file-f21) c3, [F22](#file-f22) c3, [F24](#file-f24) c8, [F26](#file-f26) c8, [F28](#file-f28) c6, [F31](#file-f31) c11, [F34](#file-f34) c3 |
| `scikit-learn` | `1.0.2` | serialized-model version | [F20](#file-f20) c12 |
| `scikit-learn` | `1.3.2` | runtime warning reports current version | [F20](#file-f20) c12 |
| `stringr` | `1.6.0` | tidyverse startup banner | [F03](#file-f03) c575, [F24](#file-f24) c6, [F26](#file-f26) c6 |
| `tibble` | `3.3.0` | tidyverse startup banner | [F24](#file-f24) c6 |
| `tibble` | `3.3.1` | tidyverse startup banner | [F03](#file-f03) c575, [F26](#file-f26) c6 |
| `tidyr` | `1.3.1` | tidyverse startup banner | [F24](#file-f24) c6 |
| `tidyr` | `1.3.2` | tidyverse startup banner | [F26](#file-f26) c6 |
| `tidyverse` | `2.0.0` | tidyverse startup banner | [F03](#file-f03) c575, [F24](#file-f24) c6, [F26](#file-f26) c6 |

`matrixStats 1.3.0` is stated only in a source comment in the bulk GSVA notebook. `matrixStats 1.5.0` is separately printed in the Gort preprocessing session below. Both records are retained with their own provenance.

The `scikit-learn` warning in the T-cell starCAT notebook explicitly distinguishes the saved estimator’s originating version (`1.0.2`) from the loading runtime (`1.3.2`). The saved-estimator version is not treated as the current runtime version.

## Full R package snapshot

**Source:** [7.0_GSE292189_Gort_preprocessing.ipynb](7.0_GSE292189_Gort_preprocessing.ipynb) c3, retained `sessionInfo()` output. These **118 package/version records apply to that saved session**. They include dependencies loaded indirectly through namespaces; they are not claimed for every R notebook.

| Package | Recorded version | Session status |
|---|---|---|
| `abind` | `1.4-8` | Namespace |
| `base64enc` | `0.1-6` | Namespace |
| `cli` | `3.6.3` | Namespace |
| `cluster` | `2.1.7` | Namespace |
| `codetools` | `0.2-19` | Namespace |
| `colorspace` | `2.1-1` | Namespace |
| `compiler` | `4.3.3` | Namespace; bundled with R |
| `cowplot` | `1.2.0` | Namespace |
| `crayon` | `1.5.3` | Namespace |
| `data.table` | `1.18.2.1` | Attached |
| `deldir` | `2.0-4` | Namespace |
| `dichromat` | `2.0-0.1` | Namespace |
| `digest` | `0.6.39` | Namespace |
| `dotCall64` | `1.2` | Namespace |
| `dplyr` | `1.1.4` | Attached |
| `evaluate` | `1.0.5` | Namespace |
| `farver` | `2.1.2` | Namespace |
| `fastDummies` | `1.7.5` | Namespace |
| `fastmap` | `1.2.0` | Namespace |
| `fitdistrplus` | `1.2-6` | Namespace |
| `future` | `1.69.0` | Namespace |
| `future.apply` | `1.20.1` | Namespace |
| `generics` | `0.1.4` | Namespace |
| `ggplot2` | `3.5.2` | Attached |
| `ggrepel` | `0.9.6` | Namespace |
| `ggridges` | `0.5.7` | Namespace |
| `globals` | `0.19.0` | Namespace |
| `glue` | `1.8.0` | Namespace |
| `goftest` | `1.2-3` | Namespace |
| `grid` | `4.3.3` | Namespace; bundled with R |
| `gridExtra` | `2.3` | Namespace |
| `gtable` | `0.3.6` | Namespace |
| `htmltools` | `0.5.9` | Namespace |
| `htmlwidgets` | `1.6.4` | Namespace |
| `httpuv` | `1.6.16` | Namespace |
| `httr` | `1.4.7` | Namespace |
| `ica` | `1.0-3` | Namespace |
| `igraph` | `2.1.2` | Namespace |
| `IRdisplay` | `1.1` | Namespace |
| `IRkernel` | `1.3.2` | Namespace |
| `irlba` | `2.3.7` | Namespace |
| `jsonlite` | `2.0.0` | Namespace |
| `KernSmooth` | `2.23-22` | Namespace |
| `later` | `1.4.5` | Namespace |
| `lattice` | `0.22-6` | Namespace |
| `lazyeval` | `0.2.2` | Namespace |
| `leiden` | `0.4.3.1` | Namespace |
| `lifecycle` | `1.0.5` | Namespace |
| `listenv` | `0.10.0` | Namespace |
| `lmtest` | `0.9-40` | Namespace |
| `magrittr` | `2.0.4` | Namespace |
| `MASS` | `7.3-60.0.1` | Namespace |
| `Matrix` | `1.6-5` | Attached |
| `matrixStats` | `1.5.0` | Namespace |
| `mime` | `0.13` | Namespace |
| `miniUI` | `0.1.2` | Namespace |
| `nlme` | `3.1-164` | Namespace |
| `parallel` | `4.3.3` | Namespace; bundled with R |
| `parallelly` | `1.46.1` | Namespace |
| `patchwork` | `1.3.2` | Attached |
| `pbapply` | `1.7-4` | Namespace |
| `pbdZMQ` | `0.3-11` | Namespace |
| `pillar` | `1.11.1` | Namespace |
| `pkgconfig` | `2.0.3` | Namespace |
| `plotly` | `4.11.0` | Namespace |
| `plyr` | `1.8.9` | Namespace |
| `png` | `0.1-8` | Namespace |
| `polyclip` | `1.10-7` | Namespace |
| `progressr` | `0.18.0` | Namespace |
| `promises` | `1.3.3` | Namespace |
| `purrr` | `1.1.0` | Attached |
| `R6` | `2.6.1` | Namespace |
| `ragg` | `1.3.3` | Namespace |
| `RANN` | `2.6.2` | Namespace |
| `RColorBrewer` | `1.1-3` | Namespace |
| `Rcpp` | `1.1.1` | Namespace |
| `RcppAnnoy` | `0.0.23` | Namespace |
| `RcppHNSW` | `0.6.0` | Namespace |
| `repr` | `1.1.7` | Namespace |
| `reshape2` | `1.4.5` | Namespace |
| `reticulate` | `1.37.0` | Namespace |
| `rlang` | `1.1.4` | Namespace |
| `ROCR` | `1.0-12` | Namespace |
| `RSpectra` | `0.16-2` | Namespace |
| `Rtsne` | `0.17` | Namespace |
| `scales` | `1.4.0` | Namespace |
| `scattermore` | `1.2` | Namespace |
| `sctransform` | `0.4.1` | Namespace |
| `Seurat` | `5.1.0` | Attached |
| `SeuratObject` | `5.0.2` | Attached |
| `shiny` | `1.11.1` | Namespace |
| `sp` | `2.1-4` | Attached |
| `spam` | `2.11-3` | Namespace |
| `spatstat.data` | `3.1-9` | Namespace |
| `spatstat.explore` | `3.3-3` | Namespace |
| `spatstat.geom` | `3.3-4` | Namespace |
| `spatstat.random` | `3.3-2` | Namespace |
| `spatstat.sparse` | `3.1-0` | Namespace |
| `spatstat.univar` | `3.1-1` | Namespace |
| `spatstat.utils` | `3.1-5` | Namespace |
| `splines` | `4.3.3` | Namespace; bundled with R |
| `stringi` | `1.8.4` | Namespace |
| `stringr` | `1.6.0` | Attached |
| `survival` | `3.5-8` | Namespace |
| `systemfonts` | `1.2.3` | Namespace |
| `tensor` | `1.5.1` | Namespace |
| `textshaping` | `0.4.1` | Namespace |
| `tibble` | `3.3.1` | Attached |
| `tidyr` | `1.3.2` | Attached |
| `tidyselect` | `1.2.1` | Namespace |
| `tools` | `4.3.3` | Namespace; bundled with R |
| `uuid` | `1.2-1` | Namespace |
| `uwot` | `0.2.4` | Namespace |
| `vctrs` | `0.6.5` | Namespace |
| `viridisLite` | `0.4.3` | Namespace |
| `withr` | `3.0.2` | Namespace |
| `xtable` | `1.8-4` | Namespace |
| `zoo` | `1.8-15` | Namespace |

The unversioned “attached base packages” line in the same session lists `stats`, `graphics`, `grDevices`, `utils`, `datasets`, `methods`, and `base`; these belong to the recorded R `4.3.3` runtime.

## Platform and compiler records

| Component | Version | Meaning | Sources |
|---|---|---|---|
| GCC | `12.3.0` | Python build compiler | [F03](#file-f03) c4, [F06](#file-f06) c4, [F07](#file-f07) c4, [F19](#file-f19) c4, [F28](#file-f28) c4, [F31](#file-f31) c4 |
| LAPACK | `3.12.0` | sessionInfo platform | [F30](#file-f30) c3 |
| OpenBLAS | `0.3.26` | sessionInfo platform | [F30](#file-f30) c3 |
| Rocky Linux | `8.10` | sessionInfo platform | [F30](#file-f30) c3 |

The full R session reports `x86_64-conda-linux-gnu` (64-bit). The GCC number is part of the Python build string, not a newly verified system compiler requirement.

## Version numbers that must not be mistaken for installed package versions

| Software/context | Recorded number(s) | Interpretation | Sources |
|---|---|---|---|
| R package build warnings | R `4.2.2`: AnnotationDbi, S4Vectors, edgeR, limma; R `4.2.3`: data.table, dplyr, stringr, tidyr, ggplot2 | R used to build those packages; the package versions themselves are NR in this notebook. | [F10](#file-f10) c3/c6 |
| SeuratObject API notices | `5.0.0` | Deprecation introduction, not a reported installed version. | [F01](#file-f01) c5/c12; [F18](#file-f18) c61 |
| dplyr API notice | `1.1.0` | Deprecation introduction for `across(...)`, not a reported installed version. | F03 c597 |
| ggplot2 API notices | `3.0.0`, `3.4.0`; source comment `3.4+` | API-change thresholds, not runtime pins. | F06 c7; F23 c4; F24 c31 |
| matrixStats diagnostic | `>= 1.2.0` | Diagnostic/deprecation threshold; not an exact package-version printout. | F10 c4 |
| nxviz API notice | `0.7.4`, future `1.0` | Warning recommends an old-API compatibility pin and mentions a future release; installed version is NR. | F13/F14 c3 |
| seaborn API notice | future `0.14.0` | Future removal target; installed version is NR. | F29 c28 |

## Dependencies without an exact recorded version

The core analysis tools **cNMF (`cnmf`), starCAT (`starcat`), CellChat, GSVA, Scanpy (`scanpy`), Dandelion (`dandelion`), clone2vec, CellRank (`cellrank`), scVelo (`scvelo`), and infercnvpy** are referenced/imported, but their exact package versions are not printed anywhere in the audited snapshot.

The following lists mean **no exact installed package version anywhere in this snapshot**. Per-file imports below remain unpinned even when another file records a version for the same package.

- **R imports:** `AnnotationDbi`, `Azimuth`, `BiocParallel`, `broom`, `Cairo`, `celda`, `CellChat`, `celldex`, `ComplexUpset`, `DoubletFinder`, `DropletUtils`, `edgeR`, `effsize`, `enrichplot`, `ggdendro`, `ggpubr`, `GSVA`, `harmony`, `limma`, `msigdbr`, `org.Hs.eg.db`, `pheatmap`, `Polychrome`, `progeny`, `rstatix`, `scDblFinder`, `sceasy`, `SingleCellExperiment`, `SingleR`, `SoupX`, `speckle`, `VennDiagram`, `viridis`.
- **Python import names:** `cellrank`, `clone2vec`, `cnmf`, `dandelion`, `infercnvpy`, `IPython`, `matplotlib`, `pandas`, `partipy`, `scanpy`, `scipy`, `scvelo`, `seaborn`, `starcat`, `statsmodels`.
- **NumPy:** `1.24.4` is recorded only in the six R-to-Python bridge configurations; the ten Python-kernel notebooks do not print their NumPy versions.
- **README declarations:** AnnData and CytoTRACE are named in the repository overview without exact versions. The development-potential notebook directly imports CellRank and scVelo; a README label or analysis name does not establish a separately installed CytoTRACE package version.

Python standard-library modules are supplied by Python and are omitted from third-party package lists. R base/recommended runtime components are identified separately where relevant. Import names are kept as written in the code rather than guessing PyPI/Conda distribution names.

## Per-file audit

Each original repository file is listed below. Declared dependencies are extracted from executable R `library`/`require`/`requireNamespace` and namespace calls, Python imports, and literal R `source()` calls. Conditional dependencies are included. Comment-only imports are excluded. These are dependency declarations, not proof that every branch was executed.

<a id="file-f01"></a>

### F01 · [1.1_doubletsambientRNA.ipynb](1.1_doubletsambientRNA.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `BiocParallel`, `DoubletFinder`, `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SoupX`, `celda`, `ggplot2`, `limma`, `patchwork`, `scCustomize`, `scDblFinder`, `sceasy`, `speckle`.
- **Version evidence, cell 15 output 2 text:** Modularity Optimizer `1.3.0`.
- **Version evidence, cell 17 output 1 text:** scCustomize `2.1.2`.

<a id="file-f02"></a>

### F02 · [1.2_standardpipelines.ipynb](1.2_standardpipelines.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DoubletFinder`, `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SoupX`, `celda`, `ggplot2`, `limma`, `patchwork`, `scCustomize`, `sceasy`, `speckle`.
- **Version evidence, cell 4 output 2 text:** Modularity Optimizer `1.3.0`.
- **Version evidence, cell 6 output 1 text:** scCustomize `2.1.2`.

<a id="file-f03"></a>

### F03 · [2.1_plasmacellcloneVDJ.ipynb](2.1_plasmacellcloneVDJ.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `ComplexHeatmap`, `DropletUtils`, `Matrix`, `Polychrome`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `VennDiagram`, `celda`, `celldex`, `circlize`, `dplyr`, `ggplot2`, `ggpubr`, `limma`, `patchwork`, `pheatmap`, `readr`, `reticulate`, `rstatix`, `scCustomize`, `scales`, `sceasy`, `speckle`, `stringr`, `tibble`, `tidyr`, `tidyverse`.
- **Version evidence, cell 4 output 1 text/plain:** Python `3.8.19`; NumPy `1.24.4`.
- **Version evidence, cell 10 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 15 output 2 text:** Modularity Optimizer `1.3.0`.
- **Version evidence, cell 575 output 1 text:** tidyverse `2.0.0`; forcats `1.0.0`; readr `2.1.5`; lubridate `1.9.4`; stringr `1.6.0`; purrr `1.1.0`; tibble `3.3.1`.
- **Version evidence, cell 597 output 1 text:** ComplexHeatmap `2.18.0`; circlize `0.4.16`.

<a id="file-f04"></a>

### F04 · [2.2_cNMF_plasmacellsGEP.ipynb](2.2_cNMF_plasmacellsGEP.ipynb)

- **Kernel runtime:** Python `3.8.20` (saved metadata).
- **Third-party Python imports:** `IPython`, `cnmf`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `scipy`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f05"></a>

### F05 · [2.3_starCAT_plasmacellsGEP.ipynb](2.3_starCAT_plasmacellsGEP.ipynb)

- **Kernel runtime:** Python `3.8.20` (saved metadata).
- **Third-party Python imports:** `cnmf`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `starcat`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f06"></a>

### F06 · [2.4_plasmacellclone-starCAT-part1.ipynb](2.4_plasmacellclone-starCAT-part1.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `celda`, `celldex`, `dplyr`, `ggplot2`, `limma`, `patchwork`, `pheatmap`, `reticulate`, `scCustomize`, `scales`, `speckle`, `stringr`, `tibble`, `tidyr`.
- **Version evidence, cell 4 output 1 text/plain:** Python `3.8.19`; NumPy `1.24.4`.
- **Version evidence, cell 7 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 13 output 2 text:** Modularity Optimizer `1.3.0`.

<a id="file-f07"></a>

### F07 · [2.4_plasmacellclone-starCAT-part2.ipynb](2.4_plasmacellclone-starCAT-part2.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `broom`, `celda`, `celldex`, `dplyr`, `effsize`, `ggplot2`, `ggpubr`, `limma`, `patchwork`, `pheatmap`, `purrr`, `reticulate`, `scCustomize`, `scales`, `speckle`, `tibble`, `tidyr`.
- **Version evidence, cell 4 output 1 text/plain:** Python `3.8.19`; NumPy `1.24.4`.
- **Version evidence, cell 7 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 13 output 2 text:** Modularity Optimizer `1.3.0`.

<a id="file-f08"></a>

### F08 · [2.5_Functional-string.ipynb](2.5_Functional-string.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `dplyr`, `ggplot2`, `readr`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f09"></a>

### F09 · [2.6_Functional-cNMF.ipynb](2.6_Functional-cNMF.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `clusterProfiler`, `dplyr`, `enrichplot`, `forcats`, `ggplot2`, `org.Hs.eg.db`, `readr`, `stringr`, `tidyr`, `viridis`.
- **Version evidence, cell 5 output 1 text:** clusterProfiler `4.10.1`.

<a id="file-f10"></a>

### F10 · [2.7_bulk_signature_GSVA_updated.ipynb](2.7_bulk_signature_GSVA_updated.ipynb)

- **Kernel runtime:** R `4.2.1` (saved metadata).
- **R package declarations:** `AnnotationDbi`, `GSVA`, `data.table`, `dplyr`, `edgeR`, `ggplot2`, `ggpubr`, `org.Hs.eg.db`, `pheatmap`, `rstatix`, `stringr`, `tibble`, `tidyr`.
- **Version evidence, cell 3 source line 1:** matrixStats `1.3.0` (source comment only).
- **Additional records:** package build warnings mention R `4.2.2`/`4.2.3`; the matrixStats warning mentions the threshold `>= 1.2.0`. c19 writes `audit/sessionInfo.txt` and `audit/org.Hs.eg.db_metadata.txt`, but neither output file is committed in this snapshot.

<a id="file-f11"></a>

### F11 · [3.1_cytotrace_cloneGEP5.ipynb](3.1_cytotrace_cloneGEP5.ipynb)

- **Kernel runtime:** Python `3.8.20` (saved metadata).
- **Third-party Python imports:** `cellrank`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `scipy`, `scvelo`, `statsmodels`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f12"></a>

### F12 · [3.2_ALl9_ALk3_ALl5-clonetrack.ipynb](3.2_ALl9_ALk3_ALl5-clonetrack.ipynb)

- **Kernel runtime:** Python `3.10.18` (saved metadata).
- **Third-party Python imports:** `clone2vec`, `matplotlib`, `numpy`, `pandas`, `partipy`, `scanpy`, `scipy`, `seaborn`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f13"></a>

### F13 · [3.3_dandelion-LCMS.ipynb](3.3_dandelion-LCMS.ipynb)

- **Kernel runtime:** Python `3.10.18` (saved metadata).
- **Third-party Python imports:** `IPython`, `dandelion`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `scipy`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.
- **Indirect dependency notice:** c3 mentions nxviz API releases `0.7.4` and `1.0`; the installed nxviz version is NR.

<a id="file-f14"></a>

### F14 · [3.4_dandelion-allclones.ipynb](3.4_dandelion-allclones.ipynb)

- **Kernel runtime:** Python `3.10.0` (saved metadata).
- **Third-party Python imports:** `IPython`, `dandelion`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `scipy`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.
- **Indirect dependency notice:** c3 mentions nxviz API releases `0.7.4` and `1.0`; the installed nxviz version is NR.

<a id="file-f15"></a>

### F15 · [3.5_infercnvpy_plasma.ipynb](3.5_infercnvpy_plasma.ipynb)

- **Kernel runtime:** Python `3.8.20` (saved metadata).
- **Third-party Python imports:** `infercnvpy`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `scipy`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f16"></a>

### F16 · [4.0_BM-TME-singleR.ipynb](4.0_BM-TME-singleR.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `celda`, `celldex`, `ggplot2`, `limma`, `patchwork`, `scCustomize`, `sceasy`, `speckle`.
- **Version evidence, cell 3 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 15 output 2 text:** Modularity Optimizer `1.3.0`.

<a id="file-f17"></a>

### F17 · [4.1_BM-TME-singleR-percelltype.ipynb](4.1_BM-TME-singleR-percelltype.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `Seurat`, `SingleCellExperiment`, `SingleR`, `celda`, `celldex`, `dendextend`, `dplyr`, `ggdendro`, `ggplot2`, `harmony`, `limma`, `patchwork`, `reshape2`, `scCustomize`, `scales`, `speckle`, `tidyr`.
- **Version evidence, cell 3 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 19 output 1 text:** dendextend `1.18.1`.
- **Version evidence, cell 31 output 2 text:** Modularity Optimizer `1.3.0`.

<a id="file-f18"></a>

### F18 · [4.2_BM-TME-azimuth-myeloid.ipynb](4.2_BM-TME-azimuth-myeloid.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `Azimuth`, `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `VennDiagram`, `celda`, `celldex`, `dplyr`, `forcats`, `ggplot2`, `ggrepel`, `limma`, `msigdbr`, `patchwork`, `pheatmap`, `purrr`, `readr`, `reshape2`, `scCustomize`, `scales`, `speckle`, `tibble`, `tidyr`, `viridis`.
- **Local helper dependencies:** [enrichment_utils_C5.R](enrichment_utils_C5.R), [enrichment_utils_C5_custom_IFN.R](enrichment_utils_C5_custom_IFN.R), [enrichment_utils_C5_custom_IFNGR.R](enrichment_utils_C5_custom_IFNGR.R), [enrichment_utils_C5_custom_PGE.R](enrichment_utils_C5_custom_PGE.R). Their package dependencies and header revisions are listed under F39–F42.
- **Version evidence, cell 3 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 5 output 1 text:** Modularity Optimizer `1.3.0`.

<a id="file-f19"></a>

### F19 · [4.3_BM-TME-T cells.ipynb](4.3_BM-TME-T%20cells.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `broom`, `celda`, `celldex`, `dplyr`, `effsize`, `forcats`, `ggplot2`, `ggpubr`, `limma`, `patchwork`, `pheatmap`, `purrr`, `readr`, `reshape2`, `reticulate`, `scCustomize`, `scales`, `sceasy`, `speckle`, `stringr`, `tibble`, `tidyr`, `viridis`, `viridisLite`.
- **Local helper dependencies:** [enrichment_utils_C5.R](enrichment_utils_C5.R), [enrichment_utils_C5_custom_IFN.R](enrichment_utils_C5_custom_IFN.R), [enrichment_utils_C5_custom_IFNGR.R](enrichment_utils_C5_custom_IFNGR.R), [enrichment_utils_C5_custom_PGE.R](enrichment_utils_C5_custom_PGE.R). Their package dependencies and header revisions are listed under F39–F42.
- **Version evidence, cell 3 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 4 output 1 text/plain:** Python `3.8.19`; NumPy `1.24.4`.

<a id="file-f20"></a>

### F20 · [4.4_starCAT_Tcells.ipynb](4.4_starCAT_Tcells.ipynb)

- **Kernel runtime:** Python `3.8.20` (saved metadata).
- **Third-party Python imports:** `matplotlib`, `numpy`, `pandas`, `scanpy`, `seaborn`, `starcat`.
- **Version evidence, cell 12 output 3 text:** scikit-learn `1.0.2` (saved estimator); scikit-learn `1.3.2` (loading runtime).

<a id="file-f21"></a>

### F21 · [4.5_BM-TME-NK cells.ipynb](4.5_BM-TME-NK%20cells.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `broom`, `celda`, `celldex`, `dplyr`, `effsize`, `forcats`, `ggplot2`, `ggpubr`, `limma`, `patchwork`, `pheatmap`, `purrr`, `reshape2`, `scCustomize`, `scales`, `speckle`, `stringr`, `tibble`, `tidyr`, `viridis`, `viridisLite`.
- **Local helper dependencies:** [enrichment_utils_C5.R](enrichment_utils_C5.R), [enrichment_utils_C5_custom_IFN.R](enrichment_utils_C5_custom_IFN.R), [enrichment_utils_C5_custom_IFNGR.R](enrichment_utils_C5_custom_IFNGR.R), [enrichment_utils_C5_custom_PGE.R](enrichment_utils_C5_custom_PGE.R). Their package dependencies and header revisions are listed under F39–F42.
- **Version evidence, cell 3 output 1 text:** scCustomize `2.1.2`.

<a id="file-f22"></a>

### F22 · [4.6_BM-TME-B cells.ipynb](4.6_BM-TME-B%20cells.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `celda`, `celldex`, `dplyr`, `ggplot2`, `limma`, `patchwork`, `pheatmap`, `purrr`, `reshape2`, `scCustomize`, `speckle`, `tidyr`.
- **Version evidence, cell 3 output 1 text:** scCustomize `2.1.2`.

<a id="file-f23"></a>

### F23 · [4.7_upsetplot_interferon.ipynb](4.7_upsetplot_interferon.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `ComplexUpset`, `dplyr`, `ggplot2`, `msigdbr`, `readr`, `tidyr`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f24"></a>

### F24 · [5.1_cellchat-78samples-trial3_annotations.ipynb](5.1_cellchat-78samples-trial3_annotations.ipynb)

- **Kernel runtime:** R `4.2.1` (saved metadata).
- **R package declarations:** `Cairo`, `CellChat`, `DoubletFinder`, `DropletUtils`, `Matrix`, `Seurat`, `SingleCellExperiment`, `SoupX`, `celda`, `dplyr`, `forcats`, `ggplot2`, `limma`, `patchwork`, `purrr`, `readr`, `scCustomize`, `speckle`, `stringr`, `tibble`, `tidyr`, `tidyverse`.
- **Version evidence, cell 6 output 1 text:** tidyverse `2.0.0`; dplyr `1.1.4`; readr `2.1.5`; forcats `1.0.0`; stringr `1.6.0`; lubridate `1.9.4`; tibble `3.3.0`; purrr `1.2.0`; tidyr `1.3.1`.
- **Version evidence, cell 8 output 1 text:** scCustomize `2.1.2`.

<a id="file-f25"></a>

### F25 · [5.2_cellchat-ALNicheGraph.ipynb](5.2_cellchat-ALNicheGraph.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `CellChat`, `RColorBrewer`, `dplyr`, `forcats`, `ggplot2`, `igraph`, `purrr`, `readr`, `scales`, `stringr`, `tibble`, `tidyr`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f26"></a>

### F26 · [5.3_CD38_publication_figures.ipynb](5.3_CD38_publication_figures.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `ComplexHeatmap`, `DoubletFinder`, `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SoupX`, `celda`, `circlize`, `dplyr`, `ggplot2`, `ggpubr`, `limma`, `patchwork`, `readr`, `scCustomize`, `scales`, `speckle`, `stringr`, `tibble`, `tidyr`, `tidyverse`.
- **Version evidence, cell 6 output 1 text:** tidyverse `2.0.0`; dplyr `1.1.4`; readr `2.1.5`; forcats `1.0.0`; stringr `1.6.0`; lubridate `1.9.4`; tibble `3.3.1`; purrr `1.1.0`; tidyr `1.3.2`.
- **Version evidence, cell 8 output 1 text:** scCustomize `2.1.2`.

<a id="file-f27"></a>

### F27 · [6.0_GSE175385_Alameda_preprocess.ipynb](6.0_GSE175385_Alameda_preprocess.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `Seurat`, `dplyr`, `ggplot2`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f28"></a>

### F28 · [6.1_plasmacell-Alameda.ipynb](6.1_plasmacell-Alameda.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `Seurat`, `SingleCellExperiment`, `SingleR`, `celda`, `celldex`, `dplyr`, `ggplot2`, `limma`, `patchwork`, `purrr`, `reticulate`, `scCustomize`, `sceasy`, `speckle`, `stringr`, `tibble`.
- **Version evidence, cell 4 output 1 text/plain:** Python `3.8.19`; NumPy `1.24.4`.
- **Version evidence, cell 6 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 15 output 2 text:** Modularity Optimizer `1.3.0`.

<a id="file-f29"></a>

### F29 · [6.2_starCAT_plasmacells-Alameda.ipynb](6.2_starCAT_plasmacells-Alameda.ipynb)

- **Kernel runtime:** Python `3.8.20` (saved metadata).
- **Third-party Python imports:** `cnmf`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `scipy`, `seaborn`, `starcat`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.
- **Diagnostic-only reference:** c28 mentions a future seaborn `0.14.0` removal; this does not identify the installed seaborn version.

<a id="file-f30"></a>

### F30 · [7.0_GSE292189_Gort_preprocessing.ipynb](7.0_GSE292189_Gort_preprocessing.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `Matrix`, `Seurat`, `data.table`, `dplyr`, `ggplot2`, `patchwork`, `purrr`, `stringr`, `tibble`, `tidyr`.
- **Version evidence:** c3 prints R `4.3.3`, the complete 118-package session snapshot above, Rocky Linux `8.10`, OpenBLAS `0.3.26`, and LAPACK `3.12.0`.

<a id="file-f31"></a>

### F31 · [7.1_plasmacell_Gort.ipynb](7.1_plasmacell_Gort.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `DropletUtils`, `Matrix`, `Polychrome`, `Seurat`, `SingleCellExperiment`, `SingleR`, `celda`, `celldex`, `dplyr`, `forcats`, `ggplot2`, `igraph`, `limma`, `patchwork`, `purrr`, `reticulate`, `scCustomize`, `scales`, `sceasy`, `speckle`, `stringr`, `tibble`, `tidyr`.
- **Version evidence, cell 4 output 1 text/plain:** Python `3.8.19`; NumPy `1.24.4`.
- **Version evidence, cell 11 output 1 text:** scCustomize `2.1.2`.
- **Version evidence, cell 20 output 2 text:** Modularity Optimizer `1.3.0`.

<a id="file-f32"></a>

### F32 · [7.2_BM-TME-azimuth-myleoid_Gort.ipynb](7.2_BM-TME-azimuth-myleoid_Gort.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `Azimuth`, `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `VennDiagram`, `celda`, `celldex`, `dplyr`, `ggplot2`, `limma`, `msigdbr`, `patchwork`, `pheatmap`, `purrr`, `readr`, `reshape2`, `scCustomize`, `speckle`, `tibble`, `tidyr`, `viridis`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f33"></a>

### F33 · [7.3_starCAT_plasmacells-Gort.ipynb](7.3_starCAT_plasmacells-Gort.ipynb)

- **Kernel runtime:** Python `3.8.20` (saved metadata).
- **Third-party Python imports:** `cnmf`, `matplotlib`, `numpy`, `pandas`, `scanpy`, `scipy`, `seaborn`, `starcat`.
- **Exact package versions:** NR in this file; imports or diagnostic notices alone do not pin them.

<a id="file-f34"></a>

### F34 · [8.1_BM-TME-azimuth-myeloid_Safina.ipynb](8.1_BM-TME-azimuth-myeloid_Safina.ipynb)

- **Kernel runtime:** R `4.3.3` (saved metadata).
- **R package declarations:** `Azimuth`, `DropletUtils`, `Matrix`, `RColorBrewer`, `Seurat`, `SingleCellExperiment`, `SingleR`, `VennDiagram`, `celda`, `celldex`, `dplyr`, `ggplot2`, `limma`, `msigdbr`, `patchwork`, `pheatmap`, `purrr`, `readr`, `reshape2`, `scCustomize`, `speckle`, `tibble`, `tidyr`, `viridis`.
- **Version evidence, cell 3 output 1 text:** scCustomize `2.1.2`.

<a id="file-f35"></a>

### F35 · [LICENSE](LICENSE)

MIT license text; no software version or executable dependency declarations.

<a id="file-f36"></a>

### F36 · [README.md](README.md)

Overview and major software dependency names; no exact software versions. See the per-file imports for actual code references.

<a id="file-f37"></a>

### F37 · [Reproducibility](Reproducibility)

Program, sample, monocyte, and NK-cell naming crosswalks; no software versions.

<a id="file-f38"></a>

### F38 · [chr_pq.txt](chr_pq.txt)

Chromosome p/q-arm annotation data; no software version or executable dependency declarations.

<a id="file-f39"></a>

### F39 · [enrichment_utils_C5.R](enrichment_utils_C5.R)

- **Runtime:** R; exact R version NR in this script.
- **R package declarations:** `GSVA`, `Matrix`, `Seurat`, `dplyr`, `forcats`, `ggplot2`, `msigdbr`, `progeny`, `readr`, `rlang`, `scales`, `stringr`, `tibble`, `tidyr`.
- **Script header revision:** `2026-06-05-C5-support`; package versions NR.

<a id="file-f40"></a>

### F40 · [enrichment_utils_C5_custom_IFN.R](enrichment_utils_C5_custom_IFN.R)

- **Runtime:** R; exact R version NR in this script.
- **R package declarations:** `GSVA`, `Matrix`, `Seurat`, `dplyr`, `forcats`, `ggplot2`, `msigdbr`, `progeny`, `readr`, `rlang`, `scales`, `stringr`, `tibble`, `tidyr`.
- **Script header revision:** `2026-06-05-C5-custom-IFN-support`; package versions NR.

<a id="file-f41"></a>

### F41 · [enrichment_utils_C5_custom_IFNGR.R](enrichment_utils_C5_custom_IFNGR.R)

- **Runtime:** R; exact R version NR in this script.
- **R package declarations:** `GSVA`, `Matrix`, `Seurat`, `dplyr`, `forcats`, `ggplot2`, `msigdbr`, `progeny`, `readr`, `rlang`, `scales`, `stringr`, `tibble`, `tidyr`.
- **Script header revision:** `2026-07-09-C5-IFNGR-no-custom-panels-no-HPO-WP`; package versions NR.

<a id="file-f42"></a>

### F42 · [enrichment_utils_C5_custom_PGE.R](enrichment_utils_C5_custom_PGE.R)

- **Runtime:** R; exact R version NR in this script.
- **R package declarations:** `GSVA`, `Matrix`, `Seurat`, `dplyr`, `forcats`, `ggplot2`, `msigdbr`, `progeny`, `readr`, `rlang`, `scales`, `stringr`, `tibble`, `tidyr`.
- **Script header revision:** `2026-06-05-C5-custom-IFN-support`; package versions NR.

## Reproducibility coverage

- Only one complete R `sessionInfo()` snapshot is retained: [F30](#file-f30) c3. No complete `pip freeze`, `conda list`, Conda environment YAML, `renv.lock`, or equivalent environment lockfile is committed in this source snapshot.
- Source helper revision strings are copied literally. In particular, `enrichment_utils_C5_custom_PGE.R` retains the header revision `2026-06-05-C5-custom-IFN-support`; the filename does not supply a different version.
- The tables report historical repository evidence. They do not claim that the packages were freshly installed, that external reference models were regenerated, or that all notebooks were rerun together successfully.
