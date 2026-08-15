# ============================================================
# enrichment_utils.R
# Version: 2026-06-05-C5-custom-IFN-support
# Reusable utilities for Seurat-based PROGENy + ssGSEA analyses
#
# Main functions:
#   1) run_combined_progeny_ssgsea()
#      - PROGENy + ssGSEA per-cell scores
#      - Post vs Pre within response group
#      - Poor vs Good within treatment group
#      - compact heatmaps

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
    C5_subcats = c("GO:BP", "GO:CC", "GO:MF", "HPO"),
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
# 1) Combined PROGENy + ssGSEA response-treatment analysis
# ============================================================
run_combined_progeny_ssgsea <- function(
    obj,
    assay_use = NULL,
    layer_use = "data",
    celltype_col = "Annotation",
    treat_col = "treatment",
    resp_col = "response",
    good_set = c("AL_CR", "AL_VGPR"),
    poor_set = c("AL_PR", "AL_NR/SD"),
    species = "Homo sapiens",
    focus_only = TRUE,
    focus_regex = "INTERFERON",
    include_C2_C7 = FALSE,
    include_C5 = FALSE,
    C2_subcats = c("CP:REACTOME", "CP:KEGG"),
    C2_regex = "",
    C5_subcats = c("GO:BP", "GO:CC", "GO:MF", "HPO"),
    C5_regex = "",
    C7_regex = "",
    custom_gene_sets = NULL,
    custom_gene_set_prefix = "CUSTOM",
    exclude_ssgsea_regex = "GSE|GSM",
    progeny_top = 500,
    exclude_progeny_pathways = c("NFkB","JAK-STAT","TNFa","Hypoxia","TGFb","Trail", "VEGF", "p53", "PI3K", "MAPK", "EGFR", "Estrogen", "Androgen", "WNT"),
    min_cells_per_group = 3,
    min_genes_per_set = 5,
    min_genes_per_custom_set = 1,
    celltype_order = NULL,
    out_dir = "./combined_PROGENy_ssGSEA_response_treatment",
    heatmap_width = 4.8,
    heatmap_height = 2.8,
    base_size = 5.5,
    x_text_size = 4.8,
    y_text_size = 4.2,
    star_size = 1.4,
    save_outputs = TRUE
) {
  .load_enrichment_packages(include_progeny = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(assay_use)) assay_use <- Seurat::DefaultAssay(obj)
  Seurat::DefaultAssay(obj) <- assay_use

  if (!assay_layer_exists(obj, assay_use, layer_use)) {
    message(sprintf("Requested layer/slot '%s' not found or empty; trying 'data'.", layer_use))
    layer_use <- "data"
  }

  md <- obj@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id")

  required_cols <- c(celltype_col, treat_col, resp_col)
  missing_cols <- setdiff(required_cols, colnames(md))
  if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

  meta <- md %>%
    dplyr::transmute(
      cell_id,
      Annotation = .data[[celltype_col]],
      treatment = .data[[treat_col]],
      response = .data[[resp_col]],
      resp_group = make_resp_group(.data[[resp_col]], good_set, poor_set)
    )

  if (is.null(celltype_order)) celltype_order <- unique(meta$Annotation)

  # -------- PROGENy --------
  organism_tag <- map_species_for_progeny(species)
  obj <- progeny::progeny(
    obj,
    scale = FALSE,
    organism = organism_tag,
    top = progeny_top,
    perm = 1,
    return_assay = TRUE,
    assay = assay_use,
    slot = layer_use
  )

  obj <- Seurat::ScaleData(obj, assay = "progeny", verbose = FALSE)
  m_progeny <- Seurat::GetAssayData(obj, assay = "progeny", slot = "scale.data")
  progeny_pathways <- rownames(m_progeny)

  progeny_long <- t(m_progeny) %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::left_join(meta, by = "cell_id") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(progeny_pathways), names_to = "pathway", values_to = "score_raw") %>%
    dplyr::filter(!is.na(resp_group), !is.na(treatment), !is.na(Annotation)) %>%
    dplyr::filter(if (length(exclude_progeny_pathways) > 0) !(pathway %in% exclude_progeny_pathways) else TRUE) %>%
    dplyr::group_by(pathway) %>%
    dplyr::mutate(score_z = zscore_safe(score_raw)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Method = "PROGENy")

  # -------- ssGSEA --------
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

  E <- get_assay_data_compat(obj, assay = assay_use, layer_or_slot = layer_use)

  pathways_list <- add_custom_gene_sets_to_list(
    pathways_list = pathways_list,
    custom_gene_sets = custom_gene_sets,
    custom_gene_set_prefix = custom_gene_set_prefix,
    genes_present = rownames(E),
    min_genes_per_custom_set = min_genes_per_custom_set
  )

  genes_keep <- intersect(rownames(E), unique(unlist(pathways_list)))
  E <- E[genes_keep, , drop = FALSE]

  pathways_list <- lapply(pathways_list, function(v) intersect(v, rownames(E)))
  pathways_list <- .filter_pathways_after_intersection(
    pathways_list = pathways_list,
    custom_gene_set_prefix = custom_gene_set_prefix,
    min_genes_per_set = min_genes_per_set,
    min_genes_per_custom_set = min_genes_per_custom_set
  )
  if (length(pathways_list) == 0) stop("No selected gene set has enough genes in this object.", call. = FALSE)

  ssm <- run_ssgsea_compat(E, pathways_list)

  ssgsea_long <- as.data.frame(t(ssm), check.names = FALSE) %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::left_join(meta, by = "cell_id") %>%
    dplyr::filter(!is.na(resp_group), !is.na(treatment), !is.na(Annotation)) %>%
    tidyr::pivot_longer(cols = -c(cell_id, Annotation, treatment, response, resp_group), names_to = "pathway", values_to = "score_raw") %>%
    dplyr::filter(if (!is.null(exclude_ssgsea_regex) && nzchar(exclude_ssgsea_regex)) !grepl(exclude_ssgsea_regex, pathway, ignore.case = TRUE) else TRUE) %>%
    dplyr::group_by(pathway) %>%
    dplyr::mutate(score_z = zscore_safe(score_raw)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Method = "ssGSEA")

  # -------- Tables --------
  all_cell_scores <- dplyr::bind_rows(
    progeny_long %>% dplyr::select(Method, cell_id, Annotation, treatment, response, resp_group, pathway, score_raw, score_z),
    ssgsea_long %>% dplyr::select(Method, cell_id, Annotation, treatment, response, resp_group, pathway, score_raw, score_z)
  )

  post_vs_pre_delta <- dplyr::bind_rows(
    make_post_vs_pre_table(progeny_long, "PROGENy", min_cells_per_group),
    make_post_vs_pre_table(ssgsea_long, "ssGSEA", min_cells_per_group)
  )

  poor_vs_good_delta <- dplyr::bind_rows(
    make_poor_vs_good_table(progeny_long, "PROGENy", min_cells_per_group),
    make_poor_vs_good_table(ssgsea_long, "ssGSEA", min_cells_per_group)
  )

  all_delta_scores <- dplyr::bind_rows(
    post_vs_pre_delta %>% dplyr::mutate(facet_group = response_group),
    poor_vs_good_delta %>% dplyr::mutate(facet_group = treatment_group)
  )

  # -------- Outputs --------
  output_files <- list(
    gene_sets_txt = file.path(out_dir, "ssgsea_gene_sets_used.txt"),
    all_cell_scores_csv = file.path(out_dir, "all_cell_scores_PROGENy_ssGSEA.csv.gz"),
    post_vs_pre_csv = file.path(out_dir, "deltaNES_Post_vs_Pre_within_GoodPoor.csv.gz"),
    poor_vs_good_csv = file.path(out_dir, "deltaNES_Poor_vs_Good_within_PrePost.csv.gz"),
    all_delta_scores_csv = file.path(out_dir, "all_delta_scores_two_contrasts_PROGENy_ssGSEA.csv.gz"),
    heatmap_post_vs_pre_pdf = file.path(out_dir, "heatmap1_Post_vs_Pre_grouped_by_GoodPoor_extra_compact.pdf"),
    heatmap_post_vs_pre_png = file.path(out_dir, "heatmap1_Post_vs_Pre_grouped_by_GoodPoor_extra_compact.png"),
    heatmap_poor_vs_good_pdf = file.path(out_dir, "heatmap2_Poor_vs_Good_grouped_by_PrePost_extra_compact.pdf"),
    heatmap_poor_vs_good_png = file.path(out_dir, "heatmap2_Poor_vs_Good_grouped_by_PrePost_extra_compact.png")
  )

  if (isTRUE(save_outputs)) {
    write_gene_sets_txt(pathways_list, output_files$gene_sets_txt)
    readr::write_csv(all_cell_scores, output_files$all_cell_scores_csv)
    readr::write_csv(post_vs_pre_delta, output_files$post_vs_pre_csv)
    readr::write_csv(poor_vs_good_delta, output_files$poor_vs_good_csv)
    readr::write_csv(all_delta_scores, output_files$all_delta_scores_csv)
  }

  p_post_vs_pre <- plot_combined_delta_heatmap(
    delta_df = post_vs_pre_delta,
    facet_col = response_group,
    fill_label = "ΔNES\nPost-Pre",
    out_pdf = if (save_outputs) output_files$heatmap_post_vs_pre_pdf else NULL,
    out_png = if (save_outputs) output_files$heatmap_post_vs_pre_png else NULL,
    celltype_order = celltype_order,
    width = heatmap_width,
    height = heatmap_height,
    base_size = base_size,
    x_text_size = x_text_size,
    y_text_size = y_text_size,
    star_size = star_size
  )

  p_poor_vs_good <- plot_combined_delta_heatmap(
    delta_df = poor_vs_good_delta,
    facet_col = treatment_group,
    fill_label = "ΔNES\nPoor-Good",
    out_pdf = if (save_outputs) output_files$heatmap_poor_vs_good_pdf else NULL,
    out_png = if (save_outputs) output_files$heatmap_poor_vs_good_png else NULL,
    celltype_order = celltype_order,
    width = heatmap_width,
    height = heatmap_height,
    base_size = base_size,
    x_text_size = x_text_size,
    y_text_size = y_text_size,
    star_size = star_size
  )

  message("Done. Output directory: ", out_dir)
  message("Gene set counts — H: ", length(unique(gs$msig_H$gs_name)),
          " | C2 selected: ", length(unique(gs$msig_C2$gs_name)),
          " | C5 selected: ", length(unique(gs$msig_C5$gs_name)),
          " | C7 selected: ", length(unique(gs$msig_C7$gs_name)))

  invisible(list(
    obj = obj,
    meta = meta,
    gene_sets = pathways_list,
    gene_set_sources = gs,
    progeny_long = progeny_long,
    ssgsea_long = ssgsea_long,
    all_cell_scores = all_cell_scores,
    post_vs_pre_delta = post_vs_pre_delta,
    poor_vs_good_delta = poor_vs_good_delta,
    all_delta_scores = all_delta_scores,
    plots = list(
      post_vs_pre_heatmap = p_post_vs_pre,
      poor_vs_good_heatmap = p_poor_vs_good
    ),
    output_files = output_files
  ))
}
