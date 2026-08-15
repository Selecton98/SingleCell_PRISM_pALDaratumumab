# ============================================================
# enrichment_utils.R
# Version: 2026-07-09-C5-IFNGR-no-custom-panels-no-HPO-WP
# Reusable utilities for Seurat-based PROGENy + ssGSEA analyses
#
# Main functions:
#
#   1) run_receptor_stratified_ssgsea()
#      - receptor-expression grouping, e.g. PTGER2/PTGER4
#      - ssGSEA phenotype analysis across receptor groups
#      - receptor 2D plot, proportion plot, heatmap, dotplot
#
# Usage in notebook:
#   source("enrichment_utils.R")
#   res <- run_combined_progeny_ssgsea(obj = sce1_tcell, focus_regex = "INTERFERON")
# ============================================================

# -----------------------------
# Package loading
# -----------------------------
.load_enrichment_packages <- function(include_progeny = TRUE) {
  pkgs <- c(
    "Seurat", "GSVA", "msigdbr", "dplyr", "tidyr", "tibble",
    "readr", "stringr", "forcats", "scales", "ggplot2", "Matrix"
  )
  if (include_progeny) pkgs <- c("progeny", pkgs)

  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      "\nInstall them before running this utility script.",
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(Seurat)
    if (include_progeny) library(progeny)
    library(GSVA)
    library(msigdbr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(readr)
    library(stringr)
    library(forcats)
    library(scales)
    library(ggplot2)
    library(Matrix)
  })
}

# -----------------------------
# General helpers
# -----------------------------
map_species_for_progeny <- function(x) {
  if (grepl("human|sapiens", x, ignore.case = TRUE)) "Human" else "Mouse"
}

get_assay_data_compat <- function(object, assay = NULL, layer_or_slot = "data") {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  if ("layer" %in% names(formals(Seurat::GetAssayData))) {
    Seurat::GetAssayData(object, assay = assay, layer = layer_or_slot)
  } else {
    Seurat::GetAssayData(object, assay = assay, slot = layer_or_slot)
  }
}

assay_layer_exists <- function(object, assay = NULL, layer_or_slot = "data") {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  x <- tryCatch(
    get_assay_data_compat(object, assay = assay, layer_or_slot = layer_or_slot),
    error = function(e) NULL
  )
  !is.null(x) && nrow(x) > 0 && ncol(x) > 0
}

zscore_safe <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

star_from_p <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  )
}

make_resp_group <- function(x, good_set, poor_set) {
  dplyr::case_when(
    x %in% good_set ~ "Good",
    x %in% poor_set ~ "Poor",
    TRUE ~ NA_character_
  )
}

fetch_data_compat <- function(object, vars, assay = NULL, layer_or_slot = "data") {
  args <- list(object = object, vars = vars)
  if (!is.null(assay)) args$assay <- assay
  if ("layer" %in% names(formals(Seurat::FetchData))) {
    args$layer <- layer_or_slot
  } else if ("slot" %in% names(formals(Seurat::FetchData))) {
    args$slot <- layer_or_slot
  }
  do.call(Seurat::FetchData, args)
}

run_ssgsea_compat <- function(expr_mat, gene_sets, kcdf = "Gaussian", abs_ranking = TRUE, seed = 1) {
  set.seed(seed)
  out <- tryCatch({
    GSVA::gsva(
      as.matrix(expr_mat),
      gene_sets,
      method = "ssgsea",
      kcdf = kcdf,
      abs.ranking = abs_ranking,
      parallel.sz = 1
    )
  }, error = function(e) {
    message("Old gsva() interface failed. Trying GSVA >= 1.50 ssgseaParam() interface.")
    param <- GSVA::ssgseaParam(
      exprData = as.matrix(expr_mat),
      geneSets = gene_sets,
      assay = NA,
      normalize = TRUE
    )
    GSVA::gsva(param)
  })
  out
}

# -----------------------------
# MSigDB gene set helper
# -----------------------------
get_msigdb_gene_sets <- function(
    species = "Homo sapiens",
    focus_only = TRUE,
    focus_regex = "INTERFERON",
    include_C2_C7 = FALSE,
    include_C5 = FALSE,
    C2_subcats = c("CP:REACTOME", "CP:KEGG"),
    C2_regex = "",
    C5_subcats = c("GO:BP", "GO:CC", "GO:MF"),
    C5_regex = "",
    C7_regex = "",
    min_genes_per_set = 5
) {
  msig_H <- msigdbr::msigdbr(species = species, category = "H") %>%
    dplyr::select(gs_name, gene_symbol, gs_cat)

  if (isTRUE(focus_only) && !is.null(focus_regex) && nzchar(focus_regex)) {
    msig_H <- msig_H %>%
      dplyr::filter(grepl(focus_regex, gs_name, ignore.case = TRUE))
  }

  msig_C2 <- tibble::tibble(gs_name = character(), gene_symbol = character(), gs_cat = character())
  msig_C5 <- tibble::tibble(gs_name = character(), gene_symbol = character(), gs_cat = character())
  msig_C7 <- tibble::tibble(gs_name = character(), gene_symbol = character(), gs_cat = character())

  if (isTRUE(include_C2_C7)) {
    if (!is.null(C2_regex) && nzchar(C2_regex)) {
      msig_C2 <- msigdbr::msigdbr(species = species, category = "C2") %>%
        dplyr::filter(gs_subcat %in% C2_subcats) %>%
        dplyr::filter(grepl(C2_regex, gs_name, ignore.case = TRUE)) %>%
        dplyr::select(gs_name, gene_symbol, gs_cat)
    }

    if (!is.null(C7_regex) && nzchar(C7_regex)) {
      msig_C7 <- msigdbr::msigdbr(species = species, category = "C7") %>%
        dplyr::filter(grepl(C7_regex, gs_name, ignore.case = TRUE)) %>%
        dplyr::select(gs_name, gene_symbol, gs_cat)
    }
  }

  if (isTRUE(include_C5) && !is.null(C5_regex) && nzchar(C5_regex)) {
    msig_C5 <- msigdbr::msigdbr(species = species, category = "C5") %>%
      dplyr::filter(gs_subcat %in% C5_subcats) %>%
      dplyr::filter(grepl(C5_regex, gs_name, ignore.case = TRUE)) %>%
      dplyr::select(gs_name, gene_symbol, gs_cat)
  }

  msig_all <- dplyr::bind_rows(msig_H, msig_C2, msig_C5, msig_C7) %>%
    dplyr::distinct(gs_name, gene_symbol)

  pathways_list <- split(msig_all$gene_symbol, msig_all$gs_name)
  pathways_list <- lapply(pathways_list, unique)
  pathways_list <- pathways_list[lengths(pathways_list) >= min_genes_per_set]

  if (length(pathways_list) == 0) {
    stop(
      "No MSigDB gene sets selected. Check focus_regex/C2_regex/C5_regex/C7_regex, include_C5/include_C2_C7, or set focus_only = FALSE.",
      call. = FALSE
    )
  }

  list(
    gene_sets = pathways_list,
    msig_H = msig_H,
    msig_C2 = msig_C2,
    msig_C5 = msig_C5,
    msig_C7 = msig_C7,
    msig_all = msig_all
  )
}

write_gene_sets_txt <- function(pathways_list, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    unlist(lapply(names(pathways_list), function(nm) {
      c(
        paste0("### ", nm),
        paste(sort(unique(pathways_list[[nm]])), collapse = ", "),
        ""
      )
    })),
    con = file
  )
}

# -----------------------------
# Custom gene-set helpers
# -----------------------------
get_recommended_ifn_gene_sets <- function() {
list(

  # Pan-cell type-I IFNα/β response
  IFNAB_CORE = c(
    "MX1", "ISG15", "IFIT3", "OAS1", "OASL",
    "IRF7", "IFI44", "IFI44L", "DDX60", "HERC6",
    "USP18", "LY6E", "BST2", "IFITM1", "PLSCR1", "IFI35"
  ),

  # Pan-cell type-II IFNγ response
  IFNG_CORE = c(
    "IRF1", "STAT1", "CXCL9", "CXCL10", "CXCL11",
    "IDO1", "GBP4", "PSMB8", "PSMB9", "PSMB10",
    "TAP1", "NLRC5", "LAP3", "WARS", "SOCS1"
  ),

  # More IFNγ-discriminating genes
  IFNG_DISCRIMINATING = c(
    "CXCL9", "CXCL11", "IDO1", "GBP4",
    "PSMB8", "PSMB9", "PSMB10",
    "TAP1", "NLRC5", "SOCS1", "IRF1"
  ),

  # Myeloid/TAM-specific type-I IFNα/β response
  MYELOID_IFNAB_SPECIFIC = c(
    "IFI30", "LPAR6", "VSIG4", "IL10",
    "IL1RN", "C1QB", "SIGLEC1"
  ),

  # Myeloid/TAM-specific type-II IFNγ response
  MYELOID_IFNG_SPECIFIC = c(
    "CXCL9", "FPR1", "JAK2", "TNFAIP2", "PTGS2",
    "IL15", "LAMP3", "NFKB1", "SLAMF8", "PSTPIP2", "FCGR1B"
  ),

  # Myeloid shared IFN activation / inhibitory-state genes
  MYELOID_IFN_SHARED = c(
    "CD86", "FCGR1A", "MYD88", "SIGLEC10", "LILRB1",
    "SECTM1", "CXCL10", "GCH1", "MX2", "TNFSF13B",
    "FPR2", "ANKRD22"
  ),

  # T/NK IFN-associated genes
  TNK_IFN_SHARED = c(
    "CCL5", "GZMA", "IL2RB", "NLRC5", "OASL", "STAT4"
  )

)
}

add_custom_gene_sets_to_list <- function(
    pathways_list,
    custom_gene_sets = NULL,
    custom_gene_set_prefix = "CUSTOM",
    genes_present = NULL,
    min_genes_per_custom_set = 1
) {
  if (is.null(custom_gene_sets)) {
    return(pathways_list)
  }

  if (!is.list(custom_gene_sets)) {
    stop("custom_gene_sets must be a named list, e.g. list(MY_SET = c('GENE1', 'GENE2')).", call. = FALSE)
  }

  if (is.null(names(custom_gene_sets)) || any(names(custom_gene_sets) == "")) {
    stop("custom_gene_sets must be a named list with non-empty names.", call. = FALSE)
  }

  custom_gene_sets <- lapply(custom_gene_sets, function(x) {
    unique(as.character(x[!is.na(x) & x != ""]))
  })

  if (!is.null(genes_present)) {
    custom_gene_sets <- lapply(custom_gene_sets, function(x) intersect(x, genes_present))
  }

  custom_gene_sets <- custom_gene_sets[lengths(custom_gene_sets) >= min_genes_per_custom_set]

  if (length(custom_gene_sets) == 0) {
    warning("No custom gene sets retained after intersecting with genes present in the object.", call. = FALSE)
    return(pathways_list)
  }

  if (!is.null(custom_gene_set_prefix) && nzchar(custom_gene_set_prefix)) {
    names(custom_gene_sets) <- paste0(custom_gene_set_prefix, "_", names(custom_gene_sets))
  }

  duplicated_names <- intersect(names(pathways_list), names(custom_gene_sets))
  if (length(duplicated_names) > 0) {
    stop("Custom gene-set names duplicate existing gene-set names: ", paste(duplicated_names, collapse = ", "), call. = FALSE)
  }

  message("Added custom gene sets: ", paste(names(custom_gene_sets), collapse = ", "))
  c(pathways_list, custom_gene_sets)
}

.filter_pathways_after_intersection <- function(
    pathways_list,
    custom_gene_set_prefix = "CUSTOM",
    min_genes_per_set = 5,
    min_genes_per_custom_set = 1
) {
  if (length(pathways_list) == 0) return(pathways_list)

  if (!is.null(custom_gene_set_prefix) && nzchar(custom_gene_set_prefix)) {
    is_custom <- grepl(paste0("^", custom_gene_set_prefix, "_"), names(pathways_list))
  } else {
    is_custom <- rep(FALSE, length(pathways_list))
  }

  keep <- ifelse(is_custom, lengths(pathways_list) >= min_genes_per_custom_set, lengths(pathways_list) >= min_genes_per_set)
  pathways_list[keep]
}

# -----------------------------
# Delta/stat helpers
# -----------------------------
make_post_vs_pre_table <- function(long_df, method_name, min_cells_per_group = 3) {
  long_df %>%
    dplyr::group_by(Annotation, resp_group, pathway) %>%
    dplyr::group_modify(~{
      dat <- .x
      n_pre  <- sum(dat$treatment == "Pre", na.rm = TRUE)
      n_post <- sum(dat$treatment == "Post", na.rm = TRUE)

      mean_pre <- ifelse(n_pre > 0, mean(dat$score_z[dat$treatment == "Pre"], na.rm = TRUE), NA_real_)
      mean_post <- ifelse(n_post > 0, mean(dat$score_z[dat$treatment == "Post"], na.rm = TRUE), NA_real_)

      p_val <- NA_real_
      if (n_pre >= min_cells_per_group && n_post >= min_cells_per_group) {
        p_val <- suppressWarnings(stats::wilcox.test(score_z ~ treatment, data = dat, exact = FALSE)$p.value)
      }

      tibble::tibble(
        contrast = "Post_vs_Pre_within_response",
        n_pre = n_pre,
        n_post = n_post,
        mean_pre = mean_pre,
        mean_post = mean_post,
        deltaNES = ifelse(n_pre > 0 && n_post > 0, mean_post - mean_pre, NA_real_),
        p_value = p_val
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pathway, Annotation) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(star = star_from_p(padj), Method = method_name) %>%
    dplyr::filter(!is.na(deltaNES)) %>%
    dplyr::rename(response_group = resp_group) %>%
    dplyr::select(
      Method, contrast, Annotation, response_group, pathway,
      n_pre, n_post, mean_pre, mean_post, deltaNES, p_value, padj, star
    )
}

make_poor_vs_good_table <- function(long_df, method_name, min_cells_per_group = 3) {
  long_df %>%
    dplyr::group_by(Annotation, treatment, pathway) %>%
    dplyr::group_modify(~{
      dat <- .x
      n_good <- sum(dat$resp_group == "Good", na.rm = TRUE)
      n_poor <- sum(dat$resp_group == "Poor", na.rm = TRUE)

      mean_good <- ifelse(n_good > 0, mean(dat$score_z[dat$resp_group == "Good"], na.rm = TRUE), NA_real_)
      mean_poor <- ifelse(n_poor > 0, mean(dat$score_z[dat$resp_group == "Poor"], na.rm = TRUE), NA_real_)

      p_val <- NA_real_
      if (n_good >= min_cells_per_group && n_poor >= min_cells_per_group) {
        p_val <- suppressWarnings(stats::wilcox.test(score_z ~ resp_group, data = dat, exact = FALSE)$p.value)
      }

      tibble::tibble(
        contrast = "Poor_vs_Good_within_treatment",
        n_good = n_good,
        n_poor = n_poor,
        mean_good = mean_good,
        mean_poor = mean_poor,
        deltaNES = ifelse(n_good > 0 && n_poor > 0, mean_poor - mean_good, NA_real_),
        p_value = p_val
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pathway, Annotation) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(star = star_from_p(padj), Method = method_name) %>%
    dplyr::filter(!is.na(deltaNES)) %>%
    dplyr::rename(treatment_group = treatment) %>%
    dplyr::select(
      Method, contrast, Annotation, treatment_group, pathway,
      n_good, n_poor, mean_good, mean_poor, deltaNES, p_value, padj, star
    )
}

# -----------------------------
# Heatmap helper
# -----------------------------
plot_combined_delta_heatmap <- function(
    delta_df,
    facet_col,
    fill_label,
    out_pdf = NULL,
    out_png = NULL,
    celltype_order = NULL,
    width = 4.8,
    height = 2.8,
    base_size = 5.5,
    x_text_size = 4.8,
    y_text_size = 4.2,
    star_size = 1.4,
    low_color = "#2166AC",
    mid_color = "white",
    high_color = "#B2182B"
) {
  facet_col <- rlang::ensym(facet_col)

  df_all <- delta_df %>%
    dplyr::mutate(
      Method = factor(Method, levels = c("PROGENy", "ssGSEA")),
      pathway_method = paste(Method, pathway, sep = " • "),
      facet_value = as.character(!!facet_col)
    )

  if (all(c("Good", "Poor") %in% unique(df_all$facet_value))) {
    df_all <- df_all %>% dplyr::mutate(facet_value = factor(facet_value, levels = c("Good", "Poor")))
  } else if (all(c("Pre", "Post") %in% unique(df_all$facet_value))) {
    df_all <- df_all %>% dplyr::mutate(facet_value = factor(facet_value, levels = c("Pre", "Post")))
  } else {
    df_all <- df_all %>% dplyr::mutate(facet_value = factor(facet_value))
  }

  obs_ann <- df_all %>% dplyr::distinct(Annotation) %>% dplyr::pull(Annotation)
  if (!is.null(celltype_order) && length(celltype_order) > 0) {
    ann_levels <- c(intersect(celltype_order, obs_ann), setdiff(sort(obs_ann), celltype_order))
  } else {
    ann_levels <- df_all %>%
      dplyr::group_by(Annotation) %>%
      dplyr::summarise(m = mean(abs(deltaNES), na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(m)) %>%
      dplyr::pull(Annotation)
  }

  row_order <- df_all %>%
    dplyr::group_by(Method, pathway) %>%
    dplyr::summarise(m = mean(deltaNES, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(Method) %>%
    dplyr::arrange(dplyr::desc(m), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pathway_method = paste(Method, pathway, sep = " • ")) %>%
    dplyr::arrange(factor(Method, levels = c("PROGENy", "ssGSEA")), dplyr::desc(m)) %>%
    dplyr::pull(pathway_method)

  df_all <- df_all %>%
    dplyr::mutate(
      Annotation = factor(Annotation, levels = ann_levels),
      pathway_method = factor(pathway_method, levels = row_order)
    )

  lim <- max(abs(df_all$deltaNES), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1e-6

  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = Annotation, y = pathway_method, fill = deltaNES)) +
    ggplot2::geom_tile(linewidth = 0.05) +
    ggplot2::geom_text(ggplot2::aes(label = star), size = star_size) +
    ggplot2::scale_fill_gradient2(
      name = fill_label,
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish,
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(0.22, "cm"),
        barheight = grid::unit(1.4, "cm")
      )
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::facet_grid(. ~ facet_value, scales = "free_y", space = "free_y") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.5, "mm"),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_text_size, margin = ggplot2::margin(t = 0)),
      axis.text.y = ggplot2::element_text(size = y_text_size, margin = ggplot2::margin(r = 0)),
      axis.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(size = base_size, face = "bold", margin = ggplot2::margin(1, 1, 1, 1)),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = base_size - 0.5),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(1, 1, 1, 1)
    )

  if (!is.null(out_pdf)) ggplot2::ggsave(out_pdf, plot = p, width = width, height = height, limitsize = FALSE, useDingbats = FALSE)
  if (!is.null(out_png)) ggplot2::ggsave(out_png, plot = p, width = width, height = height, dpi = 300, limitsize = FALSE)

  p
}


# ============================================================
# Simplified receptor-stratified ssGSEA + pseudobulk DEG
# Patched:
#   1) Default comparison = High receptor vs Low receptor
#   2) ssGSEA gene sets restricted by EP2/EP4-focused regex:
#      cAMP / PI3K / ERK / IL2-STAT5 / MYC / NF-kB / TCR /
#      inflammatory / AP-1 / OXPHOS / RP / interferon response
# ============================================================

run_receptor_stratified_ssgsea <- function(
    obj,
    assay_use = NULL,
    layer_use = "data",
    deg_layer_use = layer_use,

    celltype_col = "Annotation",
    sample_col = NULL,          # strongly recommended for real pseudobulk
    treat_col = "treatment",
    resp_col = "response",

    receptor_genes = c("PTGER2", "PTGER4"),
    receptor_group_name = NULL,
    celltypes_keep = NULL,

    species = "Homo sapiens",
    good_set = c("AL_CR", "AL_VGPR"),
    poor_set = c("AL_PR", "AL_NR/SD"),

    focus_only = FALSE,
    focus_regex = "",
    include_C2_C7 = TRUE,
    include_C5 = TRUE,
    C2_subcats = c("CP:REACTOME", "CP:KEGG", "CP:BIOCARTA"),
    C2_regex = "",
    C5_subcats = c("GO:BP", "GO:CC", "GO:MF"),
    C5_regex = "",
    C7_regex = "",
    custom_gene_sets = NULL,
    custom_gene_set_prefix = "CUSTOM",
    exclude_ssgsea_regex = "GSE|GSM",

    # New: focused enrichment regex.
    # If NULL, a built-in EP2/EP4 downstream + immune-program regex is used.
    enrichment_focus_regex = NULL,

    # New default:
    # contrast_group1 becomes Low PTGER2/PTGER4
    # contrast_group2 becomes High PTGER2/PTGER4
    contrast_group1 = NULL,
    contrast_group2 = NULL,

    min_cells_per_group = 10,
    min_cells_for_deg = 20,
    min_cells_per_pseudobulk = 5,
    min_genes_per_set = 5,
    min_genes_per_custom_set = 1,

    top_n_pathways = 20,
    top_n_deg = 10,
    exclude_receptor_genes_from_deg = TRUE,

    out_dir = "./results/receptor_simple_ssGSEA_DEG",
    save_outputs = TRUE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(forcats)
    library(stringr)
  })

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------
  # Local helpers
  # -----------------------------
  .get_mat <- function(obj, assay, layer_or_slot) {
    if (exists("get_assay_data_compat", mode = "function")) {
      return(get_assay_data_compat(obj, assay = assay, layer_or_slot = layer_or_slot))
    }

    out <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = layer_or_slot),
      error = function(e) {
        Seurat::GetAssayData(obj, assay = assay, slot = layer_or_slot)
      }
    )
    out
  }

  .row_means_safe <- function(x) {
    if (inherits(x, "Matrix")) {
      Matrix::rowMeans(x)
    } else {
      rowMeans(x)
    }
  }

  .zscore_safe <- function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    m <- mean(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - m) / s
  }

  .star_from_p <- function(p) {
    dplyr::case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  }

  .clean_pathway <- function(x) {
    x %>%
      stringr::str_replace("^HALLMARK_", "") %>%
      stringr::str_replace("^REACTOME_", "REACTOME ") %>%
      stringr::str_replace("^KEGG_", "KEGG ") %>%
      stringr::str_replace_all("_", " ")
  }

  .write_gene_sets_txt_simple <- function(pathways_list, file) {
    con <- file(file, open = "wt")
    on.exit(close(con), add = TRUE)

    for (nm in names(pathways_list)) {
      writeLines(paste0(">", nm), con)
      writeLines(paste(pathways_list[[nm]], collapse = "\t"), con)
    }
  }

  make_ep2_ep4_focus_regex <- function() {
    paste(
      c(
        # ----------------------------------------------------
        # EP2 / EP4 downstream signaling
        # ----------------------------------------------------
        "CAMP",
        "CYCLIC_AMP",
        "ADENYLATE_CYCLASE",
        "ADENYLYL_CYCLASE",
        "ADCY",
        "PKA",
        "PRKA",
        "CREB",
        "CREM",

        "PI3K",
        "PHOSPHATIDYLINOSITOL",
        "PIP3",
        "AKT",
        "MTOR",

        "ERK",
        "ERKS",
        "MAPK",
        "MAP_KINASE",
        "MEK",
        "RAF",

        # ----------------------------------------------------
        # HALLMARK / T cell activation / IL2 / STAT5
        # ----------------------------------------------------
        "HALLMARK_IL2_STAT5_SIGNALING",
        "IL2_STAT5",
        "IL_2",
        "INTERLEUKIN_2",
        "IL2",
        "STAT5",
        "TCR",
        "T_CELL_RECEPTOR",
        "T_CELL_ACTIVATION",
        "LYMPHOCYTE_ACTIVATION",
        "LEUKOCYTE_ACTIVATION",
        "IMMUNE_ACTIVATION",
        "ACTIVATION",

        # ----------------------------------------------------
        # NF-kB / TNF / inflammatory / AP-1
        # ----------------------------------------------------
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "TNFA_SIGNALING_VIA_NFKB",
        "TNF",
        "TNFA",
        "NFKB",
        "NF_KB",
        "NF-?KB",
        "NKFB",       # typo-tolerant
        "INFLAMMATORY_RESPONSE",
        "INFLAMMATION",
        "AP_1",
        "AP1",
        "JUN",
        "FOS",

        # ----------------------------------------------------
        # MYC
        # ----------------------------------------------------
        "HALLMARK_MYC_TARGETS",
        "MYC",

        # ----------------------------------------------------
        # OXPHOS / oxidative phosphorylation / mitochondria
        # ----------------------------------------------------
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "OXIDATIVE_PHOSPHORYLATION",
        "OXITATIVE_PHOSPHORYLATION",   # typo-tolerant
        "OXIDATIVE_PHOSPORYLATION",    # typo-tolerant
        "OXPHOS",
        "RESPIRATORY_CHAIN",
        "ELECTRON_TRANSPORT",
        "MITOCHONDRIAL_RESPIRATION",

        # ----------------------------------------------------
        # Ribosome / RP / translation
        # ----------------------------------------------------
        "RIBOSOME",
        "RIBOSOMAL",
        "TRANSLATION",
        "PROTEIN_TRANSLATION",
        "TRANSLATIONAL",
        "(^|_)RP($|_)",
        "(^|_)RPL",
        "(^|_)RPS",

        # ----------------------------------------------------
        # Interferon response
        # ----------------------------------------------------
        "INTERFERON",
        "IFN",
        "IFNA",
        "IFNB",
        "IFNG",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE"
      ),
      collapse = "|"
    )
  }

  # -----------------------------
  # Labels and colors
  # -----------------------------
  receptor_label <- paste(receptor_genes, collapse = "/")
  low_label  <- paste0("Low ", receptor_label)
  mid_label  <- paste0("Mid ", receptor_label)
  high_label <- paste0("High ", receptor_label)

  if (is.null(contrast_group1)) contrast_group1 <- low_label
  if (is.null(contrast_group2)) contrast_group2 <- high_label

  if (is.null(receptor_group_name)) {
    receptor_group_name <- paste0(paste(receptor_genes, collapse = "_"), "_group")
  }

  receptor_group_levels <- c("Undetectable", low_label, mid_label, high_label)

  # Python matplotlib tab10-like colors
  receptor_group_colors <- c(
    "Undetectable" = "#1f77b4",
    setNames("#ff7f0e", low_label),
    setNames("#2ca02c", mid_label),
    setNames("#d62728", high_label)
  )

  direction_colors <- c(
    setNames("#1f77b4", paste0("Up in ", contrast_group1)),
    setNames("#d62728", paste0("Up in ", contrast_group2))
  )

  # -----------------------------
  # Assay and metadata
  # -----------------------------
  if (is.null(assay_use)) assay_use <- Seurat::DefaultAssay(obj)
  Seurat::DefaultAssay(obj) <- assay_use

  E <- .get_mat(obj, assay = assay_use, layer_or_slot = layer_use)

  missing_genes <- setdiff(receptor_genes, rownames(E))
  if (length(missing_genes) > 0) {
    stop(
      "Missing receptor genes in assay/layer: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  md <- obj@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id")

  if (!(celltype_col %in% colnames(md))) {
    stop("Missing metadata column: ", celltype_col, call. = FALSE)
  }

  meta <- tibble::tibble(
    cell_id = md$cell_id,
    Annotation = as.character(md[[celltype_col]]),
    sample = if (!is.null(sample_col) && sample_col %in% colnames(md)) {
      as.character(md[[sample_col]])
    } else {
      "all_cells"
    },
    treatment = if (treat_col %in% colnames(md)) as.character(md[[treat_col]]) else NA_character_,
    response = if (resp_col %in% colnames(md)) as.character(md[[resp_col]]) else NA_character_
  ) %>%
    dplyr::mutate(
      resp_group = dplyr::case_when(
        response %in% good_set ~ "Good",
        response %in% poor_set ~ "Poor",
        TRUE ~ NA_character_
      )
    )

  if (!is.null(celltypes_keep)) {
    meta <- meta %>% dplyr::filter(Annotation %in% celltypes_keep)
    obj <- subset(obj, cells = meta$cell_id)
    E <- E[, colnames(obj), drop = FALSE]
  }

  # -----------------------------
  # 1. Strict receptor grouping
  # Undetectable = all receptor genes <= 0
  # Detectable cells only are split into Low/Mid/High tertiles within cell type
  # -----------------------------
  receptor_expr <- as.data.frame(
    t(as.matrix(E[receptor_genes, meta$cell_id, drop = FALSE])),
    check.names = FALSE
  ) %>%
    tibble::rownames_to_column("cell_id")

  names(receptor_expr)[names(receptor_expr) %in% receptor_genes] <-
    paste0(receptor_genes, "_expr")

  expr_cols <- paste0(receptor_genes, "_expr")
  detect_cols <- paste0(receptor_genes, "_detected")

  meta_rec <- meta %>%
    dplyr::left_join(receptor_expr, by = "cell_id")

  for (i in seq_along(receptor_genes)) {
    meta_rec[[detect_cols[i]]] <- meta_rec[[expr_cols[i]]] > 0
  }

  meta_rec <- meta_rec %>%
    dplyr::mutate(
      any_receptor_detected =
        rowSums(as.data.frame(dplyr::across(dplyr::all_of(detect_cols))), na.rm = TRUE) > 0,
      receptor_score =
        rowMeans(as.data.frame(dplyr::across(dplyr::all_of(expr_cols))), na.rm = TRUE)
    ) %>%
    dplyr::group_by(Annotation) %>%
    dplyr::mutate(
      receptor_tertile_num = {
        out <- rep(NA_integer_, dplyr::n())
        idx <- which(any_receptor_detected & !is.na(receptor_score))

        if (length(idx) > 0) {
          out[idx] <- dplyr::ntile(receptor_score[idx], 3)
        }

        out
      },
      receptor_group = dplyr::case_when(
        !any_receptor_detected ~ "Undetectable",
        receptor_tertile_num == 1 ~ low_label,
        receptor_tertile_num == 2 ~ mid_label,
        receptor_tertile_num == 3 ~ high_label,
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      receptor_group = factor(receptor_group, levels = receptor_group_levels)
    )

  # Individual receptor bins for co-expression dot plot
  for (i in seq_along(receptor_genes)) {
    g <- receptor_genes[i]
    expr_col <- paste0(g, "_expr")
    bin_col <- paste0(g, "_bin")

    meta_rec <- meta_rec %>%
      dplyr::group_by(Annotation) %>%
      dplyr::mutate(
        !!bin_col := {
          out <- rep("Undetectable", dplyr::n())
          idx <- which(.data[[expr_col]] > 0 & !is.na(.data[[expr_col]]))

          if (length(idx) > 0) {
            tile <- dplyr::ntile(.data[[expr_col]][idx], 3)
            out[idx] <- c("Low", "Mid", "High")[tile]
          }

          factor(out, levels = c("Undetectable", "Low", "Mid", "High"))
        }
      ) %>%
      dplyr::ungroup()
  }

  obj[[receptor_group_name]] <-
    meta_rec$receptor_group[match(colnames(obj), meta_rec$cell_id)]

  for (g in receptor_genes) {
    obj[[paste0(g, "_expr")]] <-
      meta_rec[[paste0(g, "_expr")]][match(colnames(obj), meta_rec$cell_id)]
  }

  # -----------------------------
  # Plot 1: stacked receptor-group proportion
  # -----------------------------
  prop_df <- meta_rec %>%
    dplyr::filter(!is.na(Annotation), !is.na(receptor_group)) %>%
    dplyr::count(Annotation, receptor_group, name = "n") %>%
    dplyr::group_by(Annotation) %>%
    dplyr::mutate(freq = n / sum(n)) %>%
    dplyr::ungroup()

  p_stack <- prop_df %>%
    ggplot2::ggplot(ggplot2::aes(x = Annotation, y = freq, fill = receptor_group)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = receptor_group_colors, drop = FALSE) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      x = NULL,
      y = "Fraction of cells",
      fill = "Receptor group",
      title = paste0(receptor_label, " receptor-state composition")
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  # -----------------------------
  # Plot 2: receptor co-expression dot plot
  # x = receptor 1 tertile, y = receptor 2 tertile
  # size = cell fraction, color = mean combined receptor expression
  # -----------------------------
  p_coexp <- NULL
  coexp_df <- NULL

  if (length(receptor_genes) >= 2) {
    x_bin <- paste0(receptor_genes[1], "_bin")
    y_bin <- paste0(receptor_genes[2], "_bin")

    coexp_df <- meta_rec %>%
      dplyr::filter(!is.na(Annotation), !is.na(.data[[x_bin]]), !is.na(.data[[y_bin]])) %>%
      dplyr::mutate(
        receptor1_bin = .data[[x_bin]],
        receptor2_bin = .data[[y_bin]]
      ) %>%
      dplyr::group_by(Annotation, receptor1_bin, receptor2_bin) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean_receptor_score = mean(receptor_score, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(Annotation) %>%
      dplyr::mutate(freq = n / sum(n)) %>%
      dplyr::ungroup()

    p_coexp <- coexp_df %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = receptor1_bin,
          y = receptor2_bin,
          size = freq,
          color = mean_receptor_score
        )
      ) +
      ggplot2::geom_point(alpha = 0.9) +
      ggplot2::facet_wrap(~ Annotation, ncol = 4) +
      ggplot2::scale_color_gradient(low = "#1f77b4", high = "#d62728") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(
        x = paste0(receptor_genes[1], " bin"),
        y = paste0(receptor_genes[2], " bin"),
        size = "Fraction",
        color = "Mean receptor score",
        title = paste0(receptor_genes[1], " / ", receptor_genes[2], " co-expression")
      ) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(linewidth = 0.2),
        legend.position = "bottom",
        strip.text = ggplot2::element_text(size = 8)
      )
  }

  # -----------------------------
  # 3. ssGSEA differential score plot
  # -----------------------------
  if (exists(".load_enrichment_packages", mode = "function")) {
    .load_enrichment_packages(include_progeny = FALSE)
  }

  if (!exists("get_msigdb_gene_sets", mode = "function")) {
    stop("Missing helper function: get_msigdb_gene_sets()", call. = FALSE)
  }

  if (!exists("run_ssgsea_compat", mode = "function")) {
    stop("Missing helper function: run_ssgsea_compat()", call. = FALSE)
  }

  gs <- get_msigdb_gene_sets(
    species = species,
    focus_only = focus_only,
    focus_regex = focus_regex,
    include_C2_C7 = include_C2_C7,
    include_C5 = include_C5,
    C2_subcats = C2_subcats,
    C2_regex = C2_regex,
    C5_subcats = C5_subcats,
    C5_regex = C5_regex,
    C7_regex = C7_regex,
    min_genes_per_set = min_genes_per_set
  )

  pathways_list <- gs$gene_sets

  if (!is.null(custom_gene_sets)) {
    custom_gene_sets <- custom_gene_sets[!is.na(names(custom_gene_sets))]

    for (nm in names(custom_gene_sets)) {
      pathways_list[[paste0(custom_gene_set_prefix, "_", nm)]] <-
        intersect(custom_gene_sets[[nm]], rownames(E))
    }
  }

  # Strict EP2/EP4-focused pathway filtering
  if (is.null(enrichment_focus_regex) || !nzchar(enrichment_focus_regex)) {
    enrichment_focus_regex <- make_ep2_ep4_focus_regex()
  }

  pathways_list <- pathways_list[
    grepl(enrichment_focus_regex, names(pathways_list), ignore.case = TRUE)
  ]

  if (length(pathways_list) == 0) {
    stop(
      "No gene sets matched enrichment_focus_regex. ",
      "Relax the regex or inspect names(gs$gene_sets).",
      call. = FALSE
    )
  }

  # Exclude unwanted signatures before ssGSEA
  if (!is.null(exclude_ssgsea_regex) && nzchar(exclude_ssgsea_regex)) {
    pathways_list <- pathways_list[
      !grepl(exclude_ssgsea_regex, names(pathways_list), ignore.case = TRUE)
    ]
  }

  if (length(pathways_list) == 0) {
    stop(
      "No gene sets remained after exclude_ssgsea_regex filtering.",
      call. = FALSE
    )
  }

  genes_keep <- intersect(rownames(E), unique(unlist(pathways_list)))
  E_ssgsea <- E[genes_keep, , drop = FALSE]

  pathways_list <- lapply(pathways_list, function(v) intersect(v, rownames(E_ssgsea)))

  is_custom <- startsWith(names(pathways_list), custom_gene_set_prefix)

  keep_sets <- ifelse(
    is_custom,
    lengths(pathways_list) >= min_genes_per_custom_set,
    lengths(pathways_list) >= min_genes_per_set
  )

  pathways_list <- pathways_list[keep_sets]

  if (length(pathways_list) == 0) {
    stop(
      "No focused gene set has enough genes after intersection with this object.",
      call. = FALSE
    )
  }

  message("Focused gene sets retained: ", length(pathways_list))
  message("Examples: ", paste(head(names(pathways_list), 20), collapse = " | "))

  ssm <- run_ssgsea_compat(E_ssgsea, pathways_list)

  ssgsea_long <- as.data.frame(t(ssm), check.names = FALSE) %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::left_join(meta_rec, by = "cell_id") %>%
    dplyr::filter(!is.na(Annotation), !is.na(receptor_group)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(rownames(ssm)),
      names_to = "pathway",
      values_to = "score_raw"
    ) %>%
    dplyr::group_by(pathway) %>%
    dplyr::mutate(score_z = .zscore_safe(score_raw)) %>%
    dplyr::ungroup()

  ssgsea_stats <- ssgsea_long %>%
    dplyr::filter(receptor_group %in% c(contrast_group1, contrast_group2)) %>%
    dplyr::group_by(Annotation, pathway) %>%
    dplyr::group_modify(~{
      dat <- .x

      n_g1 <- sum(dat$receptor_group == contrast_group1, na.rm = TRUE)
      n_g2 <- sum(dat$receptor_group == contrast_group2, na.rm = TRUE)

      mean_g1 <- mean(dat$score_z[dat$receptor_group == contrast_group1], na.rm = TRUE)
      mean_g2 <- mean(dat$score_z[dat$receptor_group == contrast_group2], na.rm = TRUE)

      p_value <- NA_real_

      if (n_g1 >= min_cells_per_group && n_g2 >= min_cells_per_group) {
        p_value <- suppressWarnings(
          stats::wilcox.test(score_z ~ receptor_group, data = dat, exact = FALSE)$p.value
        )
      }

      tibble::tibble(
        group1 = contrast_group1,
        group2 = contrast_group2,
        n_group1 = n_g1,
        n_group2 = n_g2,
        mean_group1 = mean_g1,
        mean_group2 = mean_g2,
        delta_score = mean_g2 - mean_g1,
        p_value = p_value
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Annotation) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      star = .star_from_p(padj),
      pathway_clean = .clean_pathway(pathway),
      direction = ifelse(
        delta_score >= 0,
        paste0("Up in ", contrast_group2),
        paste0("Up in ", contrast_group1)
      )
    ) %>%
    dplyr::arrange(padj, dplyr::desc(abs(delta_score)))

  top_pathways <- ssgsea_stats %>%
    dplyr::filter(!is.na(delta_score), is.finite(delta_score)) %>%
    dplyr::group_by(pathway_clean) %>%
    dplyr::summarise(
      max_abs_delta = max(abs(delta_score), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(max_abs_delta)) %>%
    dplyr::slice_head(n = top_n_pathways) %>%
    dplyr::pull(pathway_clean)

  ssgsea_plot_df <- ssgsea_stats %>%
    dplyr::filter(pathway_clean %in% top_pathways) %>%
    dplyr::mutate(
      pathway_clean = factor(pathway_clean, levels = rev(top_pathways))
    )

  if (nrow(ssgsea_plot_df) > 0) {
    p_ssgsea_diff <- ssgsea_plot_df %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = delta_score,
          y = pathway_clean,
          fill = direction
        )
      ) +
      ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = star,
          hjust = ifelse(delta_score >= 0, -0.1, 1.1)
        ),
        size = 3
      ) +
      ggplot2::facet_wrap(~ Annotation, scales = "free_x") +
      ggplot2::scale_fill_manual(values = direction_colors) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(
        x = paste0(contrast_group2, " - ", contrast_group1, " ssGSEA z-score"),
        y = NULL,
        fill = NULL,
        title = "Focused ssGSEA pathway score differential"
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        strip.text = ggplot2::element_text(size = 8)
      )
  } else {
    p_ssgsea_diff <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No ssGSEA pathway difference available.\nCheck receptor groups, min_cells_per_group, and enrichment_focus_regex."
      )
  }

  # -----------------------------
  # 4. Pseudobulk DEG top up/down bar plot
  # -----------------------------
  E_deg <- .get_mat(obj, assay = assay_use, layer_or_slot = deg_layer_use)
  E_deg <- E_deg[, meta_rec$cell_id, drop = FALSE]

  run_simple_pseudobulk_deg <- function(E_deg, meta_rec) {
    out <- list()
    anns <- unique(meta_rec$Annotation)

    for (ct in anns) {
      dat_meta <- meta_rec %>%
        dplyr::filter(
          Annotation == ct,
          receptor_group %in% c(contrast_group1, contrast_group2),
          cell_id %in% colnames(E_deg)
        )

      if (nrow(dat_meta) == 0) next

      total_counts <- dat_meta %>%
        dplyr::count(receptor_group, name = "n_total")

      if (!all(c(contrast_group1, contrast_group2) %in% total_counts$receptor_group)) next
      if (any(total_counts$n_total < min_cells_for_deg)) next

      pb_meta <- dat_meta %>%
        dplyr::count(sample, receptor_group, name = "n_cells") %>%
        dplyr::filter(n_cells >= min_cells_per_pseudobulk) %>%
        dplyr::mutate(
          Annotation = ct,
          pb_id = paste(ct, receptor_group, sample, sep = "||")
        )

      if (!all(c(contrast_group1, contrast_group2) %in% pb_meta$receptor_group)) next

      pb_list <- vector("list", nrow(pb_meta))

      for (i in seq_len(nrow(pb_meta))) {
        cells_i <- dat_meta %>%
          dplyr::filter(
            sample == pb_meta$sample[i],
            receptor_group == pb_meta$receptor_group[i]
          ) %>%
          dplyr::pull(cell_id)

        pb_list[[i]] <- .row_means_safe(E_deg[, cells_i, drop = FALSE])
      }

      pb_mat <- do.call(cbind, pb_list)
      rownames(pb_mat) <- rownames(E_deg)
      colnames(pb_mat) <- pb_meta$pb_id

      g1_idx <- which(pb_meta$receptor_group == contrast_group1)
      g2_idx <- which(pb_meta$receptor_group == contrast_group2)

      avg_g1 <- if (length(g1_idx) == 1) {
        pb_mat[, g1_idx]
      } else {
        rowMeans(pb_mat[, g1_idx, drop = FALSE])
      }

      avg_g2 <- if (length(g2_idx) == 1) {
        pb_mat[, g2_idx]
      } else {
        rowMeans(pb_mat[, g2_idx, drop = FALSE])
      }

      p_value <- rep(NA_real_, nrow(pb_mat))

      if (length(g1_idx) >= 2 && length(g2_idx) >= 2) {
        p_value <- apply(pb_mat, 1, function(v) {
          suppressWarnings(
            stats::wilcox.test(v[g2_idx], v[g1_idx], exact = FALSE)$p.value
          )
        })
      }

      deg_ct <- tibble::tibble(
        Annotation = ct,
        gene = rownames(pb_mat),
        group1 = contrast_group1,
        group2 = contrast_group2,
        n_cells_group1 = total_counts$n_total[match(contrast_group1, total_counts$receptor_group)],
        n_cells_group2 = total_counts$n_total[match(contrast_group2, total_counts$receptor_group)],
        n_pseudobulk_group1 = length(g1_idx),
        n_pseudobulk_group2 = length(g2_idx),
        avg_group1 = as.numeric(avg_g1),
        avg_group2 = as.numeric(avg_g2),
        avg_logFC = as.numeric(avg_g2 - avg_g1),
        p_value = p_value,
        padj = stats::p.adjust(p_value, method = "BH")
      )

      out[[ct]] <- deg_ct
    }

    dplyr::bind_rows(out)
  }

  pseudobulk_deg <- run_simple_pseudobulk_deg(E_deg, meta_rec)

  if (exclude_receptor_genes_from_deg && nrow(pseudobulk_deg) > 0) {
    pseudobulk_deg <- pseudobulk_deg %>%
      dplyr::filter(!gene %in% receptor_genes)
  }

  deg_top <- pseudobulk_deg %>%
    dplyr::filter(!is.na(avg_logFC), is.finite(avg_logFC), avg_logFC != 0) %>%
    dplyr::group_by(Annotation) %>%
    dplyr::group_modify(~{
      up <- .x %>%
        dplyr::arrange(dplyr::desc(avg_logFC)) %>%
        dplyr::slice_head(n = top_n_deg)

      down <- .x %>%
        dplyr::arrange(avg_logFC) %>%
        dplyr::slice_head(n = top_n_deg)

      dplyr::bind_rows(up, down) %>%
        dplyr::distinct(gene, .keep_all = TRUE)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      direction = ifelse(
        avg_logFC >= 0,
        paste0("Up in ", contrast_group2),
        paste0("Up in ", contrast_group1)
      ),
      gene_plot = paste(gene, Annotation, sep = "___")
    )

  if (nrow(deg_top) > 0) {
    deg_top <- deg_top %>%
      dplyr::arrange(Annotation, avg_logFC) %>%
      dplyr::mutate(gene_plot = factor(gene_plot, levels = unique(gene_plot)))

    p_deg_bar <- deg_top %>%
      ggplot2::ggplot(ggplot2::aes(x = avg_logFC, y = gene_plot, fill = direction)) +
      ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::facet_wrap(~ Annotation, scales = "free_y") +
      ggplot2::scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
      ggplot2::scale_fill_manual(values = direction_colors) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(
        x = paste0(contrast_group2, " - ", contrast_group1, " pseudobulk average expression"),
        y = NULL,
        fill = NULL,
        title = "Top pseudobulk DEG: up/down genes"
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        strip.text = ggplot2::element_text(size = 8)
      )
  } else {
    p_deg_bar <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No pseudobulk DEG result passed filtering.\nCheck sample_col, min_cells_for_deg, and receptor groups."
      )
  }

  # -----------------------------
  # Output files
  # -----------------------------
  output_files <- list(
    meta_receptor_csv = file.path(out_dir, "cell_receptor_groups.csv.gz"),
    receptor_stack_pdf = file.path(out_dir, "receptor_group_stackplot_by_celltype.pdf"),
    receptor_stack_png = file.path(out_dir, "receptor_group_stackplot_by_celltype.png"),
    receptor_coexp_dot_pdf = file.path(out_dir, "receptor_coexpression_dotplot.pdf"),
    receptor_coexp_dot_png = file.path(out_dir, "receptor_coexpression_dotplot.png"),
    gene_sets_txt = file.path(out_dir, "focused_ssgsea_gene_sets_used.txt"),
    ssgsea_cell_scores_csv = file.path(out_dir, "single_cell_focused_ssGSEA_scores.csv.gz"),
    ssgsea_stats_csv = file.path(out_dir, "focused_ssGSEA_high_vs_low_stats.csv.gz"),
    ssgsea_diff_pdf = file.path(out_dir, "focused_ssGSEA_high_vs_low_barplot.pdf"),
    ssgsea_diff_png = file.path(out_dir, "focused_ssGSEA_high_vs_low_barplot.png"),
    pseudobulk_deg_csv = file.path(out_dir, "pseudobulk_DEG_high_vs_low_stats.csv.gz"),
    pseudobulk_deg_top_pdf = file.path(out_dir, "pseudobulk_DEG_high_vs_low_top_up_down_barplot.pdf"),
    pseudobulk_deg_top_png = file.path(out_dir, "pseudobulk_DEG_high_vs_low_top_up_down_barplot.png")
  )

  if (isTRUE(save_outputs)) {
    readr::write_csv(meta_rec, output_files$meta_receptor_csv)
    readr::write_csv(ssgsea_long, output_files$ssgsea_cell_scores_csv)
    readr::write_csv(ssgsea_stats, output_files$ssgsea_stats_csv)
    readr::write_csv(pseudobulk_deg, output_files$pseudobulk_deg_csv)

    .write_gene_sets_txt_simple(pathways_list, output_files$gene_sets_txt)

    ggplot2::ggsave(output_files$receptor_stack_pdf, p_stack, width = 7.5, height = 5.5)
    ggplot2::ggsave(output_files$receptor_stack_png, p_stack, width = 7.5, height = 5.5, dpi = 300)

    if (!is.null(p_coexp)) {
      ggplot2::ggsave(output_files$receptor_coexp_dot_pdf, p_coexp, width = 10, height = 7)
      ggplot2::ggsave(output_files$receptor_coexp_dot_png, p_coexp, width = 10, height = 7, dpi = 300)
    }

    ggplot2::ggsave(output_files$ssgsea_diff_pdf, p_ssgsea_diff, width = 10, height = 7)
    ggplot2::ggsave(output_files$ssgsea_diff_png, p_ssgsea_diff, width = 10, height = 7, dpi = 300)

    ggplot2::ggsave(output_files$pseudobulk_deg_top_pdf, p_deg_bar, width = 10, height = 7)
    ggplot2::ggsave(output_files$pseudobulk_deg_top_png, p_deg_bar, width = 10, height = 7, dpi = 300)
  }

  if (is.null(sample_col)) {
    message(
      "sample_col is NULL. Pseudobulk DEG is descriptive mean difference only; ",
      "p-values require biological/sample replicates."
    )
  }

  message("Comparison: ", contrast_group2, " - ", contrast_group1)
  message("Done. Output directory: ", out_dir)

  invisible(list(
    obj = obj,
    meta_rec = meta_rec,
    prop_df = prop_df,
    coexp_df = coexp_df,
    enrichment_focus_regex = enrichment_focus_regex,
    gene_sets = pathways_list,
    gene_set_sources = gs,
    ssgsea_long = ssgsea_long,
    ssgsea_stats = ssgsea_stats,
    pseudobulk_deg = pseudobulk_deg,
    deg_top = deg_top,
    plots = list(
      receptor_stack = p_stack,
      receptor_coexpression_dot = p_coexp,
      ssgsea_differential = p_ssgsea_diff,
      pseudobulk_deg_top = p_deg_bar
    ),
    output_files = output_files
  ))
}

# ============================================================
# IFNGR-focused adapter
# Added: IFNGR1/IFNGR2 receptor-stratified ssGSEA + pseudobulk DEG
# Focused pathways:
#   1) Interleukin responses: IL-6, IL-10, IL-17, IL-18,
#      and inflammatory cytokine-response programs
#   2) Interferon responses: IFN-gamma/STAT1 axis and broader IFN programs
#   3) Non-classical MHC class I / MHC-I-like antigen presentation
# ============================================================

make_ifngr_focus_regex <- function() {
  paste(
    c(
      # ----------------------------------------------------
      # Interleukin responses / inflammatory cytokines
      # ----------------------------------------------------
 #     "IL_6", "IL6", "INTERLEUKIN_6",
 #     "IL_10", "IL10", "INTERLEUKIN_10",
 #     "IL_17", "IL17", "INTERLEUKIN_17",
 #     "IL_18", "IL18", "INTERLEUKIN_18",

      # ----------------------------------------------------
      # IFN-gamma / STAT1 axis and broad interferon programs
      # ----------------------------------------------------
 #     "INTERFERON",
      "IFNG", "IFNA", "IFNB",
      "JAK_STAT",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",


      # ----------------------------------------------------
      # Non-classical MHC class I / MHC-I-like programs
      # ----------------------------------------------------
      "NON_CLASSICAL_MHC",
      "NONCLASSICAL_MHC",
      "HLA_E", "HLA_F", "HLA_G"
    ),
    collapse = "|"
  )
}



run_ifngr_stratified_ssgsea <- function(
    obj,
    assay_use = NULL,
    layer_use = "data",
    deg_layer_use = layer_use,

    celltype_col = "Annotation",
    sample_col = NULL,
    treat_col = "treatment",
    resp_col = "response",

    receptor_genes = c("IFNGR1", "IFNGR2"),
    receptor_group_name = "IFNGR1_IFNGR2_group",
    celltypes_keep = NULL,

    species = "Homo sapiens",
    good_set = c("AL_CR", "AL_VGPR"),
    poor_set = c("AL_PR", "AL_NR/SD"),

    enrichment_focus_regex = NULL,

    min_cells_per_group = 10,
    min_cells_for_deg = 20,
    min_cells_per_pseudobulk = 5,
    min_genes_per_set = 5,
    min_genes_per_custom_set = 1,

    top_n_pathways = 30,
    top_n_deg = 10,

    out_dir = "./results/IFNGR1_IFNGR2_High_vs_Low_focused_ssGSEA_DEG",
    save_outputs = TRUE
) {
  if (is.null(enrichment_focus_regex) || !nzchar(enrichment_focus_regex)) {
    enrichment_focus_regex <- make_ifngr_focus_regex()
  }

  # No custom gene-set panels are injected in the IFNGR wrapper.
  # Enrichment is restricted to MSigDB H/C2/C5/C7 candidates matching enrichment_focus_regex.

  run_receptor_stratified_ssgsea(
    obj = obj,
    assay_use = assay_use,
    layer_use = layer_use,
    deg_layer_use = deg_layer_use,

    celltype_col = celltype_col,
    sample_col = sample_col,
    treat_col = treat_col,
    resp_col = resp_col,

    receptor_genes = receptor_genes,
    receptor_group_name = receptor_group_name,
    celltypes_keep = celltypes_keep,

    species = species,
    good_set = good_set,
    poor_set = poor_set,

    # Pull broad HALLMARK + focused C2/C5/C7 candidates,
    # then the main function applies enrichment_focus_regex again strictly.
    focus_only = FALSE,
    focus_regex = "",
    include_C2_C7 = TRUE,
    include_C5 = TRUE,
    # Exclude HPO/HP and WikiPathways/WP from IFNGR enrichment.
    C2_subcats = c("CP:REACTOME", "CP:KEGG", "CP:BIOCARTA"),
    C5_subcats = c("GO:BP", "GO:CC", "GO:MF"),
    C2_regex = enrichment_focus_regex,
    C5_regex = enrichment_focus_regex,
    C7_regex = enrichment_focus_regex,

    custom_gene_sets = NULL,
    custom_gene_set_prefix = "CUSTOM_IFNGR",
    exclude_ssgsea_regex = "GSE|GSM",

    enrichment_focus_regex = enrichment_focus_regex,

    # Default inside main function becomes:
    # High IFNGR1/IFNGR2 - Low IFNGR1/IFNGR2
    contrast_group1 = NULL,
    contrast_group2 = NULL,

    min_cells_per_group = min_cells_per_group,
    min_cells_for_deg = min_cells_for_deg,
    min_cells_per_pseudobulk = min_cells_per_pseudobulk,
    min_genes_per_set = min_genes_per_set,
    min_genes_per_custom_set = min_genes_per_custom_set,

    top_n_pathways = top_n_pathways,
    top_n_deg = top_n_deg,
    exclude_receptor_genes_from_deg = TRUE,

    out_dir = out_dir,
    save_outputs = save_outputs
  )
}

# ------------------------------------------------------------
# Example call
# ------------------------------------------------------------
# res_ifngr <- run_ifngr_stratified_ssgsea(
#   obj = obj,
#   assay_use = assay_use,
#   layer_use = layer_use,
#   deg_layer_use = layer_use,
#   celltype_col = celltype_col,
#   sample_col = NULL,  # replace with patient/sample column if available
#   treat_col = treat_col,
#   resp_col = resp_col,
#   receptor_genes = c("IFNGR1", "IFNGR2"),
#   min_cells_per_group = 10,
#   min_cells_for_deg = 20,
#   top_n_pathways = 30,
#   top_n_deg = 10,
#   out_dir = "./results/IFNGR1_IFNGR2_High_vs_Low_focused_ssGSEA_DEG"
# )
#
# head(res_ifngr$ssgsea_stats)
# head(res_ifngr$pseudobulk_deg)
# res_ifngr$output_files
