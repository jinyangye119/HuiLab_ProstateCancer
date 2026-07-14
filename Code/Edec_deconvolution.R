suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(EDec)
})

project_dir <- "/Users/YUC199/Library/CloudStorage/OneDrive-JohnsHopkins/Personal_data/Post-doc/Qian Lab JHU/Zhang Hui Lab/Hongyi/"
ref_dir <- file.path(project_dir, "raw/raw/Hongyi/Methelaytion/Reference_list")
default_methylation_file <- file.path(project_dir, "raw/raw/Hongyi/Methelaytion/10_mvalues_annotated_fulldata-001.csv")
methylation_file_override <- Sys.getenv("EDEC_METH_FILE", default_methylation_file)
output_tag <- Sys.getenv("EDEC_OUTPUT_TAG", "default")
out_dir <- file.path("/Users/YUC199/Chatgpt", paste0("edec_parameter_sweep_output_", output_tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  methylation_file = methylation_file_override,
  metadata_file = file.path(project_dir, "Output/Meta_Mathylation.csv"),
  extended_metadata_file = file.path(project_dir, "raw/raw/Hongyi/extended meta data - 4_3_2026 - Updated by Hongyi.csv"),
  reference_sample_list = file.path(ref_dir, "reference_dataset_list_strict_luminal.csv"),
  reference_manifest = file.path(ref_dir, "reference_dataset_manifest_wide_strict_luminal.csv"),
  primary_tumor_matrix = file.path(ref_dir, "prostate_primary_tumor_highcellularity_reference_GSE112047_450k_beta_matrix.tsv"),
  adjacent_normal_matrix = file.path(ref_dir, "prostate_adjacent_normal_reference_GSE112047_450k_beta_matrix.tsv"),
  bulk_cache_file = file.path(out_dir, "bulk_beta_filtered_cache.rds"),
  prepared_cache_file = file.path(out_dir, "prepared_edec_inputs.rds"),
  tumor_labels = c("No.Tumor", "Low.Purity", "multifocal", "Sufficient.Purity"),
  normal_labels = c("Normal"),
  stage0_version = "one.vs.rest",
  seed = 1L
)

sample_like <- function(x) grepl("^[A-Z0-9]+\\.[A-Z0-9]+\\.[NT]$", x)
normalize_probe_ids <- function(x) sub("_.*$", "", as.character(x))

collapse_duplicate_rows <- function(mat) {
  ids <- rownames(mat)
  if (!anyDuplicated(ids)) return(mat)
  keep_first <- !duplicated(ids)
  dup_any <- !(keep_first & !duplicated(ids, fromLast = TRUE))
  nondup_mat <- mat[!dup_any, , drop = FALSE]
  dup_mat <- mat[dup_any, , drop = FALSE]
  dup_ids <- ids[dup_any]
  dup_zero <- dup_mat
  dup_zero[is.na(dup_zero)] <- 0
  dup_sums <- rowsum(dup_zero, group = dup_ids, reorder = FALSE)
  dup_counts <- rowsum(
    matrix(as.numeric(!is.na(dup_mat)), nrow = nrow(dup_mat), ncol = ncol(dup_mat)),
    group = dup_ids,
    reorder = FALSE
  )
  dup_collapsed <- dup_sums / dup_counts
  dup_collapsed[dup_counts == 0] <- NA_real_
  combined <- rbind(nondup_mat, dup_collapsed)
  combined <- combined[order(match(rownames(combined), ids)), , drop = FALSE]
  storage.mode(combined) <- "numeric"
  combined
}

read_reference_matrix <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- normalize_probe_ids(dt[[1]])
  collapse_duplicate_rows(mat)
}

read_bulk_methylation <- function(path, keep_probes = NULL) {
  dt <- fread(path, header = TRUE, data.table = FALSE)
  sample_cols <- names(dt)[vapply(names(dt), sample_like, logical(1))]
  if (length(sample_cols) == 0) stop("Could not detect methylation sample columns.")
  probe_col <- if ("Name" %in% names(dt)) {
    "Name"
  } else if ("CpG" %in% names(dt)) {
    "CpG"
  } else if ("ID_REF" %in% names(dt)) {
    "ID_REF"
  } else {
    names(dt)[1]
  }
  probe_ids <- normalize_probe_ids(dt[[probe_col]])
  if (!is.null(keep_probes)) {
    keep_idx <- probe_ids %in% keep_probes
    dt <- dt[keep_idx, , drop = FALSE]
    probe_ids <- probe_ids[keep_idx]
  }
  mat <- as.matrix(dt[, sample_cols, drop = FALSE])
  mat[mat == "None"] <- NA
  storage.mode(mat) <- "numeric"
  rownames(mat) <- probe_ids
  collapse_duplicate_rows(mat)
}

convert_to_beta_if_needed <- function(mat) {
  rng <- range(mat, na.rm = TRUE)
  if (rng[1] >= 0 && rng[2] <= 1) mat else 2^mat / (1 + 2^mat)
}

make_unique_cols <- function(mat, prefix = NULL) {
  if (!is.null(prefix)) colnames(mat) <- paste(prefix, colnames(mat), sep = "_")
  colnames(mat) <- make.unique(colnames(mat))
  mat
}

resolve_reference_path <- function(path) {
  if (file.exists(path)) return(path)
  candidate <- file.path(ref_dir, basename(path))
  if (file.exists(candidate)) return(candidate)
  stop("Could not resolve reference path: ", path)
}

edec_class_from_sample_list <- function(sample_list) {
  ifelse(
    sample_list$cell_group == "immune", "Immune",
    ifelse(
      sample_list$sub_cell_group == "fibroblast" | sample_list$cell_group %in% c("CAF", "stromal"),
      "Stroma",
      ifelse(
        sample_list$cell_group == "epithelial_cancer", "CancerEp",
        ifelse(sample_list$cell_group == "epithelial", "NormalEp", NA_character_)
      )
    )
  )
}

load_manifest_references <- function(manifest_path, sample_list_path, keep_classes) {
  manifest <- fread(manifest_path, data.table = FALSE)
  sample_list <- fread(sample_list_path, data.table = FALSE)
  sample_list$assigned_class <- edec_class_from_sample_list(sample_list)
  sample_list <- sample_list[sample_list$assigned_class %in% keep_classes, , drop = FALSE]
  out_mats <- list()
  out_samples <- list()
  out_idx <- 1L
  for (i in seq_len(nrow(manifest))) {
    matrix_file <- resolve_reference_path(manifest$matrix_file[i])
    mat <- read_reference_matrix(matrix_file)
    keep_cols <- intersect(colnames(mat), sample_list$accession)
    if (length(keep_cols) == 0) next
    mat <- mat[, keep_cols, drop = FALSE]
    mat <- make_unique_cols(mat, manifest$dataset_name[i])
    sample_block <- sample_list[match(keep_cols, sample_list$accession), , drop = FALSE]
    sample_block$sample <- colnames(mat)
    sample_block$source <- manifest$dataset_name[i]
    out_mats[[out_idx]] <- mat
    out_samples[[out_idx]] <- sample_block[, c("sample", "assigned_class", "source", "accession", "cell_type", "cell_group", "sub_cell_group", "tissue", "platform"), drop = FALSE]
    out_idx <- out_idx + 1L
  }
  if (length(out_mats) == 0) return(list(reference_meth = NULL, reference_samples = NULL))
  common_rows <- Reduce(intersect, lapply(out_mats, rownames))
  out_mats <- lapply(out_mats, function(x) x[common_rows, , drop = FALSE])
  list(reference_meth = do.call(cbind, out_mats), reference_samples = rbindlist(out_samples, fill = TRUE))
}

prepare_inputs <- function() {
  if (file.exists(cfg$prepared_cache_file)) return(readRDS(cfg$prepared_cache_file))

  meta <- fread(cfg$metadata_file, data.table = FALSE, header = TRUE)
  if (!("Meth_sample" %in% names(meta)) && all(names(meta)[seq_len(min(3, ncol(meta)))] %in% c("V1", "V2", "V3"))) {
    meta <- fread(cfg$metadata_file, data.table = FALSE, header = FALSE)
    names(meta)[seq_len(min(3, ncol(meta)))] <- c("row_id", "Meth_sample", "condition")
    if (nrow(meta) > 0 && identical(as.character(meta$Meth_sample[1]), "Meth_sample")) meta <- meta[-1, , drop = FALSE]
  }
  stopifnot(all(c("Meth_sample", "condition") %in% names(meta)))
  meta$sample_id <- meta$Meth_sample

  if (file.exists(cfg$extended_metadata_file)) {
    purity_meta <- fread(cfg$extended_metadata_file, data.table = FALSE)
    purity_meta <- purity_meta[, intersect(c("common_ID", "FirstCategory", "Purity"), names(purity_meta)), drop = FALSE]
    if (all(c("common_ID", "Purity") %in% names(purity_meta))) {
      purity_meta$Purity <- suppressWarnings(as.numeric(purity_meta$Purity))
      meta <- merge(meta, purity_meta, by.x = "sample_id", by.y = "common_ID", all.x = TRUE)
    }
  }

  primary_tumor <- make_unique_cols(read_reference_matrix(cfg$primary_tumor_matrix), "PrimaryTumor")
  adjacent_normal <- make_unique_cols(read_reference_matrix(cfg$adjacent_normal_matrix), "AdjacentNormal")
  extra_reference <- load_manifest_references(
    cfg$reference_manifest,
    cfg$reference_sample_list,
    keep_classes = c("Immune", "Stroma")
  )

  reference_blocks <- list(primary_tumor, adjacent_normal)
  sample_blocks <- list(
    data.frame(sample = colnames(primary_tumor), assigned_class = "CancerEp", source = "GSE112047_primary_high_cellularity_tumor", stringsAsFactors = FALSE),
    data.frame(sample = colnames(adjacent_normal), assigned_class = "NormalEp", source = "GSE112047_adjacent_normal_tissue", stringsAsFactors = FALSE)
  )
  if (!is.null(extra_reference$reference_meth)) {
    reference_blocks <- c(reference_blocks, list(extra_reference$reference_meth))
    sample_blocks <- c(sample_blocks, list(as.data.frame(extra_reference$reference_samples)))
  }

  common_ref <- Reduce(intersect, lapply(reference_blocks, rownames))
  reference_blocks <- lapply(reference_blocks, function(x) x[common_ref, , drop = FALSE])
  reference_meth <- do.call(cbind, reference_blocks)
  reference_samples <- rbindlist(sample_blocks, fill = TRUE)
  reference_samples <- reference_samples[match(colnames(reference_meth), reference_samples$sample), , drop = FALSE]
  reference_classes <- reference_samples$assigned_class

  if (file.exists(cfg$bulk_cache_file)) {
    bulk_beta <- readRDS(cfg$bulk_cache_file)
  } else {
    bulk_raw <- read_bulk_methylation(cfg$methylation_file, keep_probes = rownames(reference_meth))
    bulk_beta <- convert_to_beta_if_needed(bulk_raw)
    keep_rows <- rowMeans(is.na(bulk_beta)) <= 0.05
    bulk_beta <- bulk_beta[keep_rows, , drop = FALSE]
    row_med <- rowMedians(bulk_beta, na.rm = TRUE)
    na_idx <- which(is.na(bulk_beta), arr.ind = TRUE)
    if (nrow(na_idx) > 0) bulk_beta[na_idx] <- row_med[na_idx[, 1]]
    bulk_beta <- bulk_beta[rowVars(bulk_beta) > 0, , drop = FALSE]
    saveRDS(bulk_beta, cfg$bulk_cache_file)
  }

  matched_samples <- intersect(colnames(bulk_beta), meta$sample_id)
  bulk_beta <- bulk_beta[, matched_samples, drop = FALSE]
  meta <- meta[match(matched_samples, meta$sample_id), , drop = FALSE]
  common_cpgs <- Reduce(intersect, list(rownames(bulk_beta), rownames(reference_meth)))
  bulk_beta <- bulk_beta[common_cpgs, , drop = FALSE]
  reference_meth <- reference_meth[common_cpgs, , drop = FALSE]

  class_levels <- unique(reference_classes)
  keep_cpg_for_stage0 <- apply(reference_meth, 1, function(x) {
    all(vapply(class_levels, function(cl) sum(!is.na(x[reference_classes == cl])) >= 2, logical(1)))
  })
  reference_meth_stage0 <- reference_meth[keep_cpg_for_stage0, , drop = FALSE]

  primary_delta <- abs(
    rowMeans(reference_meth[, reference_classes == "CancerEp", drop = FALSE], na.rm = TRUE) -
      rowMeans(reference_meth[, reference_classes == "NormalEp", drop = FALSE], na.rm = TRUE)
  )
  normal_mask <- meta$condition %in% cfg$normal_labels
  tumor_mask <- meta$condition %in% cfg$tumor_labels
  bulk_delta <- abs(
    rowMeans(bulk_beta[, tumor_mask, drop = FALSE], na.rm = TRUE) -
      rowMeans(bulk_beta[, normal_mask, drop = FALSE], na.rm = TRUE)
  )

  prepared <- list(
    meta = meta,
    bulk_beta = bulk_beta,
    reference_meth = reference_meth,
    reference_meth_stage0 = reference_meth_stage0,
    reference_samples = as.data.frame(reference_samples),
    reference_classes = reference_classes,
    primary_delta = primary_delta,
    bulk_delta = bulk_delta
  )
  saveRDS(prepared, cfg$prepared_cache_file)
  prepared
}

rename_duplicate_components <- function(labels) {
  label_counts <- ave(seq_along(labels), labels, FUN = seq_along)
  label_totals <- ave(seq_along(labels), labels, FUN = length)
  ifelse(label_totals > 1, paste0(labels, "_", label_counts), labels)
}

score_one_run <- function(prepared, stage0_max_p_value, stage0_num_markers, primary_tumor_vs_normal_marker_n, bulk_contrast_marker_n, stage1_num_cell_types) {
  markers_stage0 <- EDec::run_edec_stage_0(
    reference_meth = prepared$reference_meth_stage0,
    reference_classes = prepared$reference_classes,
    max_p_value = stage0_max_p_value,
    num_markers = stage0_num_markers,
    version = cfg$stage0_version
  )
  primary_markers <- names(sort(prepared$primary_delta, decreasing = TRUE))[seq_len(min(primary_tumor_vs_normal_marker_n, length(prepared$primary_delta)))]
  bulk_markers <- names(sort(prepared$bulk_delta, decreasing = TRUE))[seq_len(min(bulk_contrast_marker_n, length(prepared$bulk_delta)))]
  markers <- unique(c(markers_stage0, primary_markers, bulk_markers))
  markers <- intersect(markers, rownames(prepared$bulk_beta))
  if (length(markers) < 100) stop("Too few markers selected.")

  set.seed(cfg$seed)
  stage1_result <- EDec::run_edec_stage_1(
    meth_bulk_samples = prepared$bulk_beta,
    informative_loci = markers,
    num_cell_types = stage1_num_cell_types
  )
  component_ids <- paste0("Component_", seq_len(ncol(stage1_result$methylation)))
  colnames(stage1_result$methylation) <- component_ids
  colnames(stage1_result$proportions) <- component_ids

  cors_deconv_refs <- cor(
    prepared$reference_meth[markers, , drop = FALSE],
    stage1_result$methylation[markers, , drop = FALSE],
    use = "pairwise.complete.obs"
  )
  colnames(cors_deconv_refs) <- component_ids
  best_ref_idx <- apply(cors_deconv_refs, 2, which.max)
  best_ref_cor <- apply(cors_deconv_refs, 2, max)
  best_ref_sample <- rownames(cors_deconv_refs)[best_ref_idx]
  best_ref_class <- prepared$reference_samples$assigned_class[match(best_ref_sample, prepared$reference_samples$sample)]

  component_annotation <- data.frame(
    component = component_ids,
    best_reference_sample = best_ref_sample,
    assigned_class = best_ref_class,
    best_reference_correlation = best_ref_cor,
    component_label = rename_duplicate_components(best_ref_class),
    stringsAsFactors = FALSE
  )

  cancer_cols <- component_annotation$component[component_annotation$assigned_class == "CancerEp"]
  tumor_fraction <- if (length(cancer_cols) > 0) {
    rowSums(stage1_result$proportions[, cancer_cols, drop = FALSE])
  } else {
    rep(0, nrow(stage1_result$proportions))
  }
  names(tumor_fraction) <- rownames(stage1_result$proportions)

  summary_df <- data.frame(
    sample_id = names(tumor_fraction),
    tumor_fraction = as.numeric(tumor_fraction),
    stringsAsFactors = FALSE
  )
  summary_df <- merge(summary_df, prepared$meta[, intersect(c("sample_id", "condition", "Purity"), names(prepared$meta)), drop = FALSE], by = "sample_id", all.x = TRUE)
  normal_vals <- summary_df$tumor_fraction[summary_df$condition == "Normal"]
  notumor_vals <- summary_df$tumor_fraction[summary_df$condition == "No.Tumor"]
  sufficient <- summary_df[summary_df$condition == "Sufficient.Purity" & !is.na(summary_df$Purity), , drop = FALSE]
  purity_cor <- if (nrow(sufficient) >= 4) suppressWarnings(cor(sufficient$tumor_fraction, as.numeric(sufficient$Purity), use = "complete.obs")) else NA_real_

  class_counts <- table(factor(component_annotation$assigned_class, levels = c("CancerEp", "NormalEp", "Immune", "Stroma")))
  result <- data.frame(
    stage0_max_p_value = stage0_max_p_value,
    stage0_num_markers = stage0_num_markers,
    primary_tumor_vs_normal_marker_n = primary_tumor_vs_normal_marker_n,
    bulk_contrast_marker_n = bulk_contrast_marker_n,
    stage1_num_cell_types = stage1_num_cell_types,
    selected_markers = length(markers),
    Normal_mean = mean(normal_vals, na.rm = TRUE),
    Normal_median = median(normal_vals, na.rm = TRUE),
    Normal_max = max(normal_vals, na.rm = TRUE),
    NoTumor_mean = mean(notumor_vals, na.rm = TRUE),
    NoTumor_median = median(notumor_vals, na.rm = TRUE),
    NoTumor_max = max(notumor_vals, na.rm = TRUE),
    NoTumor_minus_Normal_mean = mean(notumor_vals, na.rm = TRUE) - mean(normal_vals, na.rm = TRUE),
    SufficientPurity_cor = purity_cor,
    CancerEp_components = as.integer(class_counts[["CancerEp"]]),
    NormalEp_components = as.integer(class_counts[["NormalEp"]]),
    Immune_components = as.integer(class_counts[["Immune"]]),
    Stroma_components = as.integer(class_counts[["Stroma"]]),
    component_labels = paste(component_annotation$component_label, collapse = ";"),
    min_best_ref_cor = min(component_annotation$best_reference_correlation, na.rm = TRUE),
    mean_best_ref_cor = mean(component_annotation$best_reference_correlation, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  list(result = result, component_annotation = component_annotation, sample_summary = summary_df)
}

prepared <- prepare_inputs()
message("Prepared bulk samples: ", ncol(prepared$bulk_beta), "; CpGs: ", nrow(prepared$bulk_beta))
message("Reference classes: ", paste(names(table(prepared$reference_classes)), as.integer(table(prepared$reference_classes)), sep = "=", collapse = ", "))

sweep_mode <- Sys.getenv("EDEC_SWEEP_MODE", "focused")
if (identical(sweep_mode, "full")) {
  grid <- expand.grid(
    stage0_max_p_value = c(1e-3, 1e-4, 1e-5, 1e-6),
    stage0_num_markers = c(150L, 300L, 500L),
    primary_tumor_vs_normal_marker_n = c(100L, 300L, 500L),
    bulk_contrast_marker_n = c(100L, 300L, 500L),
    stage1_num_cell_types = 3:8,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  marker_configs <- data.frame(
    stage0_max_p_value = c(1e-3, 1e-5, 1e-5, 1e-4),
    stage0_num_markers = c(150L, 300L, 150L, 300L),
    primary_tumor_vs_normal_marker_n = c(100L, 300L, 500L, 500L),
    bulk_contrast_marker_n = c(100L, 300L, 100L, 500L)
  )
  grid <- merge(marker_configs, data.frame(stage1_num_cell_types = 3:8))
}

results_file <- file.path(out_dir, paste0("edec_parameter_sweep_results_", sweep_mode, ".tsv"))
done <- if (file.exists(results_file)) fread(results_file, data.table = FALSE) else data.frame()
if (nrow(done) > 0) {
  done$key <- with(done, paste(stage0_max_p_value, stage0_num_markers, primary_tumor_vs_normal_marker_n, bulk_contrast_marker_n, stage1_num_cell_types, sep = "|"))
} else {
  done$key <- character(0)
}
grid$key <- with(grid, paste(stage0_max_p_value, stage0_num_markers, primary_tumor_vs_normal_marker_n, bulk_contrast_marker_n, stage1_num_cell_types, sep = "|"))
grid <- grid[!grid$key %in% done$key, , drop = FALSE]

message("Runs remaining: ", nrow(grid))
for (i in seq_len(nrow(grid))) {
  row <- grid[i, , drop = FALSE]
  message("Run ", i, "/", nrow(grid), ": ", row$key)
  run_result <- tryCatch(
    score_one_run(
      prepared = prepared,
      stage0_max_p_value = row$stage0_max_p_value,
      stage0_num_markers = row$stage0_num_markers,
      primary_tumor_vs_normal_marker_n = row$primary_tumor_vs_normal_marker_n,
      bulk_contrast_marker_n = row$bulk_contrast_marker_n,
      stage1_num_cell_types = row$stage1_num_cell_types
    ),
    error = function(e) e
  )
  if (inherits(run_result, "error")) {
    result <- data.frame(
      stage0_max_p_value = row$stage0_max_p_value,
      stage0_num_markers = row$stage0_num_markers,
      primary_tumor_vs_normal_marker_n = row$primary_tumor_vs_normal_marker_n,
      bulk_contrast_marker_n = row$bulk_contrast_marker_n,
      stage1_num_cell_types = row$stage1_num_cell_types,
      error = conditionMessage(run_result),
      stringsAsFactors = FALSE
    )
  } else {
    result <- run_result$result
    detail_prefix <- paste("p", row$stage0_max_p_value, "s0", row$stage0_num_markers, "pt", row$primary_tumor_vs_normal_marker_n, "bulk", row$bulk_contrast_marker_n, "k", row$stage1_num_cell_types, sep = "_")
    saveRDS(run_result, file.path(out_dir, paste0("run_detail_", gsub("[^A-Za-z0-9_]", "", detail_prefix), ".rds")))
  }
  fwrite(result, results_file, sep = "\t", append = file.exists(results_file), col.names = !file.exists(results_file))
}

all_results <- fread(results_file, data.table = FALSE)
ok_results <- all_results[!is.na(all_results$Normal_mean), , drop = FALSE]
ok_results$rank_low_normal_high_delta <- rank(ok_results$Normal_mean, ties.method = "min") + rank(-ok_results$NoTumor_minus_Normal_mean, ties.method = "min")
ok_results <- ok_results[order(ok_results$rank_low_normal_high_delta, ok_results$Normal_mean, -ok_results$NoTumor_minus_Normal_mean), , drop = FALSE]
fwrite(ok_results, file.path(out_dir, "edec_parameter_sweep_results_ranked.tsv"), sep = "\t")
fwrite(head(ok_results, 50), file.path(out_dir, "edec_parameter_sweep_top50.tsv"), sep = "\t")
print(head(ok_results, 20))
