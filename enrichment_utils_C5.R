# ============================================================
# enrichment_utils.R
# Version: 2026-06-05-C5-support
# Reusable utilities for Seurat-based PROGENy + ssGSEA analyses
#
# Main functions:
#   1) run_combined_progeny_ssgsea()
#      - PROGENy + ssGSEA per-cell scores
#      - Post vs Pre within response group
#      - Poor vs Good within treatment group
#      - compact heatmaps
#
#   2) run_receptor_stratified_ssgsea()
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
    exclude_ssgsea_regex = "GSE|GSM",
    progeny_top = 500,
    exclude_progeny_pathways = c("Trail", "VEGF", "p53", "PI3K", "MAPK", "EGFR", "Estrogen", "Androgen", "WNT"),
    min_cells_per_group = 3,
    min_genes_per_set = 5,
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
  genes_keep <- intersect(rownames(E), unique(unlist(pathways_list)))
  E <- E[genes_keep, , drop = FALSE]

  pathways_list <- lapply(pathways_list, function(v) intersect(v, rownames(E)))
  pathways_list <- pathways_list[lengths(pathways_list) >= min_genes_per_set]
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

# ============================================================
# 2) Receptor-expression stratified ssGSEA analysis
# ============================================================
run_receptor_stratified_ssgsea <- function(
    obj,
    assay_use = NULL,
    layer_use = "data",
    celltype_col = "Annotation",
    treat_col = "treatment",
    resp_col = "response",
    receptor_genes = c("PTGER2", "PTGER4"),
    receptor_group_name = NULL,
    celltypes_keep = NULL,
    species = "Homo sapiens",
    good_set = c("AL_CR", "AL_VGPR"),
    poor_set = c("AL_PR", "AL_NR/SD"),
    focus_only = TRUE,
    focus_regex = "",
    include_C2_C7 = TRUE,
    include_C5 = FALSE,
    C2_subcats = c("CP:REACTOME", "CP:KEGG", "CP:BIOCARTA", "CP:WIKIPATHWAYS"),
    C2_regex = "",
    C5_subcats = c("GO:BP", "GO:CC", "GO:MF", "HPO"),
    C5_regex = "",
    C7_regex = "",
    exclude_ssgsea_regex = "GSE|GSM",
    contrast_group1 = "Undetectable",
    contrast_group2 = NULL,
    min_cells_per_group = 10,
    min_genes_per_set = 5,
    top_n_pathways = 30,
    out_dir = "./receptor_stratified_ssGSEA",
    save_outputs = TRUE
) {
  .load_enrichment_packages(include_progeny = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(assay_use)) assay_use <- Seurat::DefaultAssay(obj)
  Seurat::DefaultAssay(obj) <- assay_use

  if (is.null(receptor_group_name)) {
    receptor_group_name <- paste0(paste(receptor_genes, collapse = "_"), "_group")
  }
  if (is.null(contrast_group2)) {
    contrast_group2 <- paste0("High ", paste(receptor_genes, collapse = "/"))
  }

  if (!assay_layer_exists(obj, assay_use, layer_use)) {
    message(sprintf("Requested layer/slot '%s' not found or empty; trying 'data'.", layer_use))
    layer_use <- "data"
  }

  md <- obj@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell_id")
  if (!(celltype_col %in% colnames(md))) stop("Missing metadata column: ", celltype_col, call. = FALSE)

  missing_genes <- setdiff(receptor_genes, rownames(obj[[assay_use]]))
  if (length(missing_genes) > 0) stop("Missing receptor genes in assay '", assay_use, "': ", paste(missing_genes, collapse = ", "), call. = FALSE)

  meta <- md %>%
    dplyr::transmute(
      cell_id,
      Annotation = .data[[celltype_col]],
      treatment = if (treat_col %in% colnames(md)) .data[[treat_col]] else NA_character_,
      response = if (resp_col %in% colnames(md)) .data[[resp_col]] else NA_character_,
      resp_group = if (resp_col %in% colnames(md)) make_resp_group(.data[[resp_col]], good_set, poor_set) else NA_character_
    )

  if (!is.null(celltypes_keep)) {
    cells_keep <- meta %>% dplyr::filter(Annotation %in% celltypes_keep) %>% dplyr::pull(cell_id)
    obj <- subset(obj, cells = cells_keep)
    meta <- meta %>% dplyr::filter(cell_id %in% cells_keep)
  }

  receptor_expr <- fetch_data_compat(obj, vars = receptor_genes, assay = assay_use, layer_or_slot = layer_use) %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("cell_id")

  names(receptor_expr)[names(receptor_expr) %in% receptor_genes] <- paste0(receptor_genes, "_expr")
  expr_cols <- paste0(receptor_genes, "_expr")
  detect_cols <- paste0(receptor_genes, "_detected")
  z_cols <- paste0(receptor_genes, "_z")

  for (i in seq_along(expr_cols)) {
    receptor_expr[[detect_cols[i]]] <- receptor_expr[[expr_cols[i]]] > 0
  }

  meta_rec <- meta %>%
    dplyr::left_join(receptor_expr, by = "cell_id") %>%
    dplyr::mutate(any_receptor_detected = rowSums(as.data.frame(dplyr::across(dplyr::all_of(detect_cols)))) > 0) %>%
    dplyr::group_by(Annotation)

  for (i in seq_along(expr_cols)) {
    meta_rec <- meta_rec %>%
      dplyr::mutate(!!z_cols[i] := zscore_safe(.data[[expr_cols[i]]]))
  }

  meta_rec <- meta_rec %>%
    dplyr::mutate(receptor_score = rowMeans(as.data.frame(dplyr::across(dplyr::all_of(z_cols))), na.rm = TRUE)) %>%
    dplyr::mutate(
      receptor_tertile = dplyr::case_when(
        !any_receptor_detected ~ "Undetectable",
        TRUE ~ as.character(dplyr::ntile(receptor_score, 3))
      ),
      receptor_group = dplyr::case_when(
        receptor_tertile == "Undetectable" ~ "Undetectable",
        receptor_tertile == "1" ~ paste0("Low ", paste(receptor_genes, collapse = "/")),
        receptor_tertile == "2" ~ paste0("Mid ", paste(receptor_genes, collapse = "/")),
        receptor_tertile == "3" ~ paste0("High ", paste(receptor_genes, collapse = "/")),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      receptor_group = factor(
        receptor_group,
        levels = c(
          "Undetectable",
          paste0("Low ", paste(receptor_genes, collapse = "/")),
          paste0("Mid ", paste(receptor_genes, collapse = "/")),
          paste0("High ", paste(receptor_genes, collapse = "/"))
        )
      )
    )

  obj[[receptor_group_name]] <- meta_rec$receptor_group[match(colnames(obj), meta_rec$cell_id)]
  for (g in receptor_genes) {
    obj[[paste0(g, "_expr")]] <- meta_rec[[paste0(g, "_expr")]][match(colnames(obj), meta_rec$cell_id)]
  }

  # -------- Receptor plots --------
  output_files <- list(
    gene_sets_txt = file.path(out_dir, "ssgsea_gene_sets_used.txt"),
    cell_scores_csv = file.path(out_dir, "single_cell_receptor_ssGSEA_scores.csv.gz"),
    stats_csv = file.path(out_dir, "receptor_group_ssGSEA_stats.csv.gz"),
    stats_baseline_csv = file.path(out_dir, "receptor_group_ssGSEA_stats_baseline_Pre_only.csv.gz"),
    receptor_2d_pdf = file.path(out_dir, "receptor_2D_groups_by_celltype.pdf"),
    receptor_2d_png = file.path(out_dir, "receptor_2D_groups_by_celltype.png"),
    prop_pdf = file.path(out_dir, "receptor_group_proportion_by_celltype.pdf"),
    prop_png = file.path(out_dir, "receptor_group_proportion_by_celltype.png"),
    pathway_heat_pdf = file.path(out_dir, "receptor_group_ssGSEA_heatmap.pdf"),
    pathway_heat_png = file.path(out_dir, "receptor_group_ssGSEA_heatmap.png"),
    pathway_dot_pdf = file.path(out_dir, "receptor_group_ssGSEA_dotplot.pdf"),
    pathway_dot_png = file.path(out_dir, "receptor_group_ssGSEA_dotplot.png")
  )

  if (length(receptor_genes) >= 2) {
    xg <- paste0(receptor_genes[1], "_expr")
    yg <- paste0(receptor_genes[2], "_expr")
    p_receptor_2d <- meta_rec %>%
      dplyr::filter(!is.na(Annotation), !is.na(receptor_group)) %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[[xg]], y = .data[[yg]], color = receptor_group)) +
      ggplot2::geom_point(size = 0.35, alpha = 0.45) +
      ggplot2::facet_wrap(~ Annotation, scales = "free", ncol = 4) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(
        x = paste0(receptor_genes[1], " normalized expression"),
        y = paste0(receptor_genes[2], " normalized expression"),
        color = "Receptor group",
        title = paste0("Single cells grouped by ", paste(receptor_genes, collapse = "/"), " expression")
      ) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 8), panel.grid = ggplot2::element_blank(), legend.position = "bottom")
  } else {
    p_receptor_2d <- NULL
  }

  prop_df <- meta_rec %>%
    dplyr::filter(!is.na(Annotation), !is.na(receptor_group)) %>%
    dplyr::count(Annotation, receptor_group, name = "n") %>%
    dplyr::group_by(Annotation) %>%
    dplyr::mutate(freq = n / sum(n)) %>%
    dplyr::ungroup()

  p_prop <- prop_df %>%
    ggplot2::ggplot(ggplot2::aes(x = Annotation, y = freq, fill = receptor_group)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = NULL, y = "Fraction of cells", fill = "Receptor group", title = "Distribution of receptor-defined immune-cell states") +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "bottom")

  if (isTRUE(save_outputs)) {
    if (!is.null(p_receptor_2d)) {
      ggplot2::ggsave(output_files$receptor_2d_pdf, p_receptor_2d, width = 11, height = 8)
      ggplot2::ggsave(output_files$receptor_2d_png, p_receptor_2d, width = 11, height = 8, dpi = 300)
    }
    ggplot2::ggsave(output_files$prop_pdf, p_prop, width = 7.5, height = 5.5)
    ggplot2::ggsave(output_files$prop_png, p_prop, width = 7.5, height = 5.5, dpi = 300)
  }

  # -------- ssGSEA gene sets and scores --------
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
  genes_keep <- intersect(rownames(E), unique(unlist(pathways_list)))
  E <- E[genes_keep, , drop = FALSE]
  pathways_list <- lapply(pathways_list, function(v) intersect(v, rownames(E)))
  pathways_list <- pathways_list[lengths(pathways_list) >= min_genes_per_set]
  if (length(pathways_list) == 0) stop("No selected gene set has enough genes in this object.", call. = FALSE)

  ssm <- run_ssgsea_compat(E, pathways_list)

  ssgsea_long <- as.data.frame(t(ssm), check.names = FALSE) %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::left_join(meta_rec, by = "cell_id") %>%
    dplyr::filter(!is.na(Annotation), !is.na(receptor_group)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(rownames(ssm)), names_to = "pathway", values_to = "score_raw") %>%
    dplyr::filter(if (!is.null(exclude_ssgsea_regex) && nzchar(exclude_ssgsea_regex)) !grepl(exclude_ssgsea_regex, pathway, ignore.case = TRUE) else TRUE) %>%
    dplyr::group_by(pathway) %>%
    dplyr::mutate(score_z = zscore_safe(score_raw)) %>%
    dplyr::ungroup()

  make_receptor_stats <- function(df) {
    df %>%
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
          p_value <- suppressWarnings(stats::wilcox.test(score_z ~ receptor_group, data = dat, exact = FALSE)$p.value)
        }
        tibble::tibble(
          group1 = contrast_group1,
          group2 = contrast_group2,
          n_group1 = n_g1,
          n_group2 = n_g2,
          mean_group1 = mean_g1,
          mean_group2 = mean_g2,
          deltaNES = mean_g2 - mean_g1,
          p_value = p_value
        )
      }) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(pathway) %>%
      dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(star = star_from_p(padj)) %>%
      dplyr::arrange(padj, dplyr::desc(abs(deltaNES)))
  }

  ssgsea_stats <- make_receptor_stats(ssgsea_long)
  ssgsea_stats_baseline <- ssgsea_long %>% dplyr::filter(treatment == "Pre") %>% make_receptor_stats()

  if (isTRUE(save_outputs)) {
    write_gene_sets_txt(pathways_list, output_files$gene_sets_txt)
    readr::write_csv(ssgsea_long, output_files$cell_scores_csv)
    readr::write_csv(ssgsea_stats, output_files$stats_csv)
    readr::write_csv(ssgsea_stats_baseline, output_files$stats_baseline_csv)
  }

  # -------- Pathway plots --------
  plot_df <- ssgsea_stats %>%
    dplyr::filter(!is.na(deltaNES)) %>%
    dplyr::mutate(
      pathway_clean = pathway %>% stringr::str_replace("^HALLMARK_", "") %>% stringr::str_replace_all("_", " "),
      Annotation = factor(Annotation, levels = unique(meta_rec$Annotation))
    )

  top_pathways <- plot_df %>%
    dplyr::group_by(pathway_clean) %>%
    dplyr::summarise(max_abs_delta = max(abs(deltaNES), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(max_abs_delta)) %>%
    dplyr::slice_head(n = top_n_pathways) %>%
    dplyr::pull(pathway_clean)

  plot_df_top <- plot_df %>%
    dplyr::filter(pathway_clean %in% top_pathways) %>%
    dplyr::mutate(pathway_clean = forcats::fct_reorder(pathway_clean, abs(deltaNES), .fun = max))

  p_heat <- plot_df_top %>%
    ggplot2::ggplot(ggplot2::aes(x = Annotation, y = pathway_clean, fill = deltaNES)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = star), size = 3) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(x = NULL, y = NULL, fill = paste0(contrast_group2, " - ", contrast_group1, "\nssGSEA z-score"), title = "Pathway activity associated with high receptor expression") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1), panel.grid = ggplot2::element_blank())

  p_dot <- plot_df_top %>%
    ggplot2::ggplot(ggplot2::aes(x = Annotation, y = pathway_clean)) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(padj), color = deltaNES), alpha = 0.9) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(x = NULL, y = NULL, color = paste0(contrast_group2, " - ", contrast_group1), size = "-log10(BH-adjusted P)", title = "ssGSEA pathway differences by receptor state") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1), panel.grid.major = ggplot2::element_line(linewidth = 0.2))

  if (isTRUE(save_outputs)) {
    ggplot2::ggsave(output_files$pathway_heat_pdf, p_heat, width = 10, height = 7)
    ggplot2::ggsave(output_files$pathway_heat_png, p_heat, width = 10, height = 7, dpi = 300)
    ggplot2::ggsave(output_files$pathway_dot_pdf, p_dot, width = 10, height = 7)
    ggplot2::ggsave(output_files$pathway_dot_png, p_dot, width = 10, height = 7, dpi = 300)
  }

  message("Done. Output directory: ", out_dir)
  message("Gene set counts — H: ", length(unique(gs$msig_H$gs_name)),
          " | C2 selected: ", length(unique(gs$msig_C2$gs_name)),
          " | C5 selected: ", length(unique(gs$msig_C5$gs_name)),
          " | C7 selected: ", length(unique(gs$msig_C7$gs_name)))

  invisible(list(
    obj = obj,
    meta_rec = meta_rec,
    prop_df = prop_df,
    gene_sets = pathways_list,
    gene_set_sources = gs,
    ssgsea_long = ssgsea_long,
    ssgsea_stats = ssgsea_stats,
    ssgsea_stats_baseline = ssgsea_stats_baseline,
    plots = list(
      receptor_2d = p_receptor_2d,
      receptor_proportion = p_prop,
      pathway_heatmap = p_heat,
      pathway_dotplot = p_dot
    ),
    output_files = output_files
  ))
}
