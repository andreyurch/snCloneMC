#!/usr/bin/env Rscript
# =============================================================================
# snCloneMC - single-nucleus Clones from Mutations and Copy number
# Step 3: clustering and clonal inference for targeted single-nucleus
# DNA-seq (snDNA-seq) data based on mutations and copy-number alterations
#
# Usage:   Rscript scripts/03_snclonemc.R config/my_config.R
# Author:  Andrey A. Yurchenko (Inserm U1361, Gustave Roussy)
# Licence: MIT (see LICENSE)
#
# Design notes (v2)
#   * All cell/feature objects are indexed BY NAME (no positional matching).
#   * Variant -> label mapping is built explicitly (no colnames(...) <- ids).
#   * Arm values = SUMMED raw counts per arm (stable), converted to fractions
#     AFTER bad arms are removed, divided by the spike-in median fraction per arm
#     and log2-transformed (median used for references).
#   * Data-driven amplicon QC (spike-in detection + spike-in noise) used only to
#     FILTER amplicons before summing.
#   * Cell-level QC (reads, MAD-based low/high read filter, amplicon detection,
#     VAF missingness) applied to the whole analysis.
#   * Robust normal-cluster call (never the noise cluster, size + altered-arm
#     criteria, falls back to spike-in-only reference when no normal found).
#   * NA-safe mutation filtering; minimum carrier cells for rescued variants.
#   * Block weighting of VAF vs CN features; NA VAFs imputed explicitly.
#   * Same UMAP/HDBSCAN settings for the full and tumour-only analysis.
#   * pheatmap(silent = TRUE) + safe device handling; ggsave() everywhere.
#   * Everything per-sample in run_sample(); shared files loaded once.
#   * Per-cell table + RDS + QC tables + run overview are saved.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(uwot)
  library(dbscan)
  library(pheatmap)
  library(rstatix)
})

`%nin%` <- Negate(`%in%`)

# =============================================================================
# 0. CONFIG
# Defaults for all parameters. Paths are NOT set here: they are read from a
# config file passed on the command line (see config/config_example.R), which
# defines a list `cfg_user`. Any entry of `cfg_user` overrides the default below.
# =============================================================================
cfg <- list(
  # ---- paths (set in the config file) ---------------------------------------
  samples_file   = NULL,   # text file, one sample name per line, no header
  data_dir       = NULL,   # directory with the per-sample input files
  spikein_file   = NULL,   # read counts (cells x amplicons) of trusted normal cells
  arm_annot_file = NULL,   # amplicon -> chromosome arm / focal gene annotation
  idt_file       = NULL,   # OncoKB-annotated list of known variants (rescue list)
  blacklist_file = NULL,   # recurrent artefactual variants (columns id, fdr, log)
  outdir         = NULL,   # output directory (created if needed)
  n_cores        = 1,
  seed           = 12345,

  # ---- spike-in normals -----------------------------------------------------
  n_spikein = 300,              # FIRST n rows of the spike-in file = trusted normal cells (no resampling)
  qc_filter_spikein = FALSE,    # TRUE = also apply cell QC to those spike-in cells (QC is always reported)

  # ---- cell-level QC (applied to tumour AND spike-in cells) -----------------
  qc_min_reads          = 2000, # absolute floor of total reads / cell (tune!)
  qc_mad_k              = 3.5,  # log10(total reads) outside median +/- k*MAD -> fail
  qc_detect_reads       = 1,    # amplicon "detected" if reads >= this
  qc_min_frac_detected  = 0.80, # min fraction of amplicons detected per cell
  max_na_vaf_frac       = 0.50, # drop cells ungenotyped for > this fraction of selected variants
  min_cells_sample      = 100,  # skip sample if fewer cells after QC

  # ---- amplicon QC (from spike-in cells) ------------------------------------
  amp_detect_reads    = 1,
  amp_min_detect_frac = 0.90,   # amplicon detected in >= 90% of spike-in cells
  amp_noise_mad_k     = 4,      # drop amplicons with spike-in noise > median + k*MAD
  amp_noise_filter    = TRUE,   # drop amplicons that are noisy in spike-ins (filter only, counts are then summed)
  pseudocount         = 0.5,    # only used for the amplicon noise QC
  min_ratio           = 0.05,   # floor of arm ratio before log2 (arm completely lost / dropout)
  min_amp_per_arm     = 3,
  min_amp_per_gene    = 1,      # features ending in "_CN" (focal amplifications)
  # Arms NOT listed here are excluded (this reproduces the old bad-arm list:
  # 18p, 4p, 20q, 9p, 16p, 3p, 2p are intentionally absent).
  arm_order = c("CCND1_CN","FGFR1_CN","MYC_CN","MDM2_CN","ERBB2_CN","1p","1q","2q","3q","4q",
                "5p","5q","6p","6q","7p","7q","8p","8q","9q","10p","10q","11p","11q","12p",
                "12q","13q","14q","15q","16q","17p","17q","18q","19p","19q","20p","21q","22q",
                "Xp","Xq"),
  lr_clip   = 2,                # clip log2 ratios to +/- this for clustering/heatmap
  heat_lim  = 1,                # heatmap colour range (log2): -1 ~ 1 copy loss, +1 ~ 4 copies

  # ---- mutation selection ---------------------------------------------------
  min_genotyped_frac = 0.70,
  min_AF             = 0.001,
  min_het_frac       = 0.20,
  germline_vaf_lo    = 40,      # drop variants with VAF_by_cell_count in [lo, hi]
  germline_vaf_hi    = 60,      # (set lo=hi=NA to disable)  <- check this is intended
  onc_classes        = c("Likely Oncogenic", "Oncogenic", "Resistance"),   # OncoKB classes treated as oncogenic
  min_rescue_cells   = 5,       # oncogenic/IDT-listed variants need >= this many carrier cells
  extra_blacklist_ids = c("ZNF292_S2049Vfs*11"),   # recurrent false positive
  vaf_weight         = 0.5,     # relative weight of the VAF block vs CN block in clustering

  # ---- clustering -----------------------------------------------------------
  n_neighbors = 50,
  min_pts     = 20,

  # ---- normal cluster call --------------------------------------------------
  # The normal cluster is the cluster with the FLATTEST median profile (lowest SD of the
  # per-arm median log2 ratios). This tolerates constant per-arm offsets between the
  # spike-in and the sample (batch / amplicon-efficiency effects), which are removed
  # afterwards by re-referencing to the called normal cluster.
  normal_min_cells         = 30,
  # The cluster with the LOWEST profile SD is ALWAYS called normal. The two cutoffs below only
  # label how clear the call is: 'lowest_sd_clear' if SD <= normal_max_sd OR SD <= normal_rel_sd *
  # median SD of the other clusters, otherwise 'lowest_sd_low_confidence' (still used for this
  # sample, but NOT used to learn the arm bias). Set normal_require_flat = TRUE to refuse low-confidence calls.
  normal_max_sd            = 0.30,
  normal_rel_sd            = 0.90,
  normal_require_flat      = FALSE,
  # Other clusters whose median profile differs from the called normal by <= this in EVERY arm
  # (log2) are also labelled normal (same CN as normal). 0 = disable.
  normal_merge_max_diff    = 0.35,
  normal_alt_threshold     = 0.30,  # reporting only: |median log2| above this = altered arm
  # Manual override for samples where the call is wrong/missing, e.g.
  # normal_override = c(SAMPLE_A = "5", SAMPLE_B = "2")   (value = cluster id from CLUSTER_STATS.tsv)
  normal_override          = character(0),
  exclude_unassigned_from_tumor = TRUE,
  min_cells_tumor_recluster = 100,

  # ---- systematic spike-in vs. sample arm bias --------------------------------
  # The spike-in is a different library, so many arms are shifted (e.g. CCND1_CN, 10q, 11p/q,
  # 19p, 20p) in ALL cells of a sample incl. normal ones. After a first run the pipeline learns
  # this bias from the called normal clusters of all samples and writes <outdir>/ARM_BIAS.tsv.
  # Put that path here and re-run: bias is subtracted before clustering / normal calling.
  arm_bias_file   = NULL,         # e.g. file.path(outdir, "ARM_BIAS.tsv")
  learn_arm_bias  = TRUE,
  min_bias_samples = 5,           # need normal calls from >= this many samples

  # ---- putative doublet clusters (normal + tumour cell) ----------------------
  detect_doublets          = TRUE,
  doublet_w_range          = c(0.2, 0.8), # fitted tumour fraction (CN and VAF) must lie in this range
  doublet_min_arm_delta    = 0.2,         # arm informative if |ratio_parent - 1| >= this
  doublet_min_arms         = 3,
  doublet_min_cn_cor       = 0.6,         # correlation of cluster vs parent tumour CN deviations
  doublet_min_vaf_delta    = 0.1,         # variant informative if |VAF_parent - VAF_normal| >= this (0-1 scale)
  # Reads: median reads of the cluster vs. median of (normal + parent) cells, one-sided Wilcoxon p below.
  # If BOTH CN and VAF fit a normal+tumour mixture, a small but significant increase is enough;
  # if only one of them can be evaluated, a clearer increase is required.
  doublet_min_read_ratio_both = 1.05,
  doublet_min_read_ratio   = 1.25,
  doublet_reads_p          = 0.01,        # one-sided Wilcoxon p for "more reads"

  # ---- significance of CN shifts vs. normal (written to tables only, used by nothing else) ----
  # shift = cluster median - normal median (log2). A gain/loss needs ALL of:
  #   (1) cell level: BH-adjusted Wilcoxon < sig_alpha and bootstrap CI excluding 0 (a weak gate:
  #       with thousands of cells everything passes);
  #   (2) amplicon level: the amplicons of the arm agree (one-sample t-test of their shifts,
  #       p < sig_amp_alpha, arm needs >= 3 amplicons; leave-one-amplicon-out mean keeps its sign);
  #   (3) |shift| >= max(sig_min_effect, pseudo-cluster null of that arm, sig_rel_effect * cluster scale).
  # Other labels: uncertain_few_amplicons / uncertain_amplicon_inconsistent /
  # uncertain_small_effect / uncertain_not_significant / neutral.
  sig_alpha              = 0.01,  # BH-adjusted cell-level Wilcoxon AND bootstrap CI excluding 0
  sig_amp_alpha          = 0.05,  # amplicon-level t-test (amplicons are the replicates)
  sig_min_effect         = 0.25,  # absolute minimum |shift| (log2) for a gain/loss call
  # Cluster-adaptive threshold: a call must also reach sig_rel_effect * (typical size of the cluster's
  # strongest events = median of its sig_scale_top_n largest |shifts|). The assay's dynamic range is
  # compressed and differs between clusters (clonal vs. diluted/doublet), and small consistent shifts
  # (~0.1-0.25) occur on many arms of heavily rearranged clusters (closure / tumour-vs-normal effects).
  sig_rel_effect         = 0.40,
  sig_scale_top_n        = 5,
  sig_pseudocluster_null = TRUE,  # per-arm threshold = max difference between pseudo-clusters
                                  # obtained by clustering the NORMAL cells alone (selection effect)
  sig_min_group          = 20,    # min cells of a pseudo-cluster used for that null
  sig_n_boot             = 200,
  sig_max_cells          = 1500   # per group (cluster / normal) subsample for bootstrap + Wilcoxon
)

# ---- read the user config ---------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 03_snclonemc.R <config.R>")
if (!file.exists(args[1])) stop("Config file not found: ", args[1])
cfg_env <- new.env()
sys.source(args[1], envir = cfg_env)
if (!exists("cfg_user", envir = cfg_env, inherits = FALSE))
  stop("The config file must define a list named `cfg_user`")
cfg <- modifyList(cfg, get("cfg_user", envir = cfg_env))

req_paths <- c("samples_file", "data_dir", "spikein_file", "arm_annot_file",
               "idt_file", "blacklist_file", "outdir")
miss <- req_paths[vapply(req_paths, function(k) is.null(cfg[[k]]), logical(1))]
if (length(miss)) stop("Missing in config: ", paste(miss, collapse = ", "))
for (k in setdiff(req_paths, "outdir"))
  if (!file.exists(cfg[[k]])) stop(k, " not found: ", cfg[[k]])

# =============================================================================
# 1. HELPERS
# =============================================================================

make_palette <- function(lv) {
  base <- c("#0072B5FF","#20854EFF","#631879FF","#6F99ADFF","#EFC000FF","#3B7A57FF",
            "#868686FF","#FDAF91FF","#8A9045FF","grey10","#D72631FF","pink","#9B59B6",
            "#F2AA4CFF","#E18727FF","#C16622FF")
  other <- setdiff(lv, "unassigned")
  c(unassigned = "grey80", setNames(rep_len(base, length(other)), other))
}

heat_palette <- colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(50)

# ---- cell-level QC -----------------------------------------------------------
qc_cells <- function(counts, cfg) {
  total <- rowSums(counts)
  frac  <- rowMeans(counts >= cfg$qc_detect_reads)
  lt    <- log10(pmax(total, 1))
  lo    <- max(log10(cfg$qc_min_reads), median(lt) - cfg$qc_mad_k * mad(lt))
  hi    <- median(lt) + cfg$qc_mad_k * mad(lt)
  reason <- rep("", length(total))
  reason[lt < lo] <- "low_reads"
  reason[reason == "" & lt > hi] <- "high_reads_possible_doublet"
  reason[reason == "" & frac < cfg$qc_min_frac_detected] <- "low_amplicon_detection"
  data.frame(cell = rownames(counts), total_reads = total, frac_detected = frac,
             qc_pass = reason == "",
             qc_fail_reason = ifelse(reason == "", "pass", reason),
             row.names = rownames(counts), stringsAsFactors = FALSE)
}

# ---- shared inputs (loaded once) ---------------------------------------------
load_shared <- function(cfg) {
  sp  <- as.matrix(read.delim(cfg$spikein_file, row.names = 1))
  sp  <- sp[seq_len(min(nrow(sp), cfg$n_spikein)), , drop = FALSE]    # first N cells are trusted normals
  rownames(sp) <- paste0("normal_", rownames(sp))
  qc_s  <- qc_cells(sp, cfg)
  sp_ok <- if (isTRUE(cfg$qc_filter_spikein)) sp[qc_s$qc_pass, , drop = FALSE] else sp
  message("Spike-in: first ", nrow(sp), " cells taken; ", sum(qc_s$qc_pass), " would pass QC; ", nrow(sp_ok), " used")

  ann <- read.delim(cfg$arm_annot_file)[, c(1, 8)]
  colnames(ann) <- c("amplicon", "arm")
  ann <- unique(ann[!is.na(ann$arm) & ann$arm != "", ])

  idt <- read.delim(cfg$idt_file)
  idt$maf_position <- paste(paste0("chr", idt$CHROMOSOME), idt$START_POSITION, idt$END_POSITION,
                            idt$REFERENCE_ALLELE, idt$TUMOR_SEQ_ALLELE2, sep = "_")
  idt <- idt[idt$ONCOGENIC %in% c("Inconclusive", cfg$onc_classes), ]
  idt$id <- paste(idt$ONCOKB_HUGO_SYMBOL, idt$ONCOKB_PROTEIN_CHANGE, sep = "_")

  black <- read.delim(cfg$blacklist_file)
  black <- black[which(black$fdr < 0.1 & black$log > 0), ]
  black <- black[black$id %nin% idt$id, ]

  bias <- NULL
  if (!is.null(cfg$arm_bias_file) && file.exists(cfg$arm_bias_file)) {
    b <- read.delim(cfg$arm_bias_file)
    bias <- setNames(b$bias, b$arm)
    message("Arm bias loaded from ", cfg$arm_bias_file, " (", length(bias), " arms)")
  }

  list(spike = sp_ok, spike_qc = qc_s, ann = ann, idt = idt, black = black, arm_bias = bias)
}

# ---- arm-level log2 ratios from SUMMED counts -------------------------------
compute_arm_lr <- function(tum, spk, ann, cfg) {
  amps0 <- Reduce(intersect, list(colnames(tum), colnames(spk), ann$amplicon))
  if (length(amps0) < 50) stop("Only ", length(amps0), " amplicons shared between tumour, spike-in and annotation")
  amp_qc <- data.frame(amplicon = amps0, stringsAsFactors = FALSE)

  # amplicon QC from spike-in cells (used ONLY to filter amplicons)
  det <- colMeans(spk[, amps0, drop = FALSE] >= cfg$amp_detect_reads)
  amp_qc$spike_detect_frac <- det
  amps <- amps0[det >= cfg$amp_min_detect_frac]
  amp_qc$spike_noise <- NA_real_
  if (isTRUE(cfg$amp_noise_filter) && length(amps) > 10) {
    sp  <- spk[, amps, drop = FALSE]
    spn <- sp / rowSums(sp) * median(rowSums(sp))
    ref <- apply(spn, 2, median)
    lr  <- log2(sweep(spn + cfg$pseudocount, 2, ref + cfg$pseudocount, "/"))
    noise <- apply(lr, 2, mad)
    amp_qc$spike_noise[match(names(noise), amp_qc$amplicon)] <- noise
    amps <- amps[noise <= median(noise) + cfg$amp_noise_mad_k * mad(noise)]
  }
  amp_qc$kept <- amp_qc$amplicon %in% amps

  # arms (only arms listed in arm_order; enough amplicons)
  a2   <- ann[ann$amplicon %in% amps, ]
  arms <- split(a2$amplicon, a2$arm)
  arms <- arms[intersect(cfg$arm_order, names(arms))]
  min_n <- ifelse(grepl("_CN$", names(arms)), cfg$min_amp_per_gene, cfg$min_amp_per_arm)
  arms <- arms[lengths(arms) >= min_n]
  if (length(arms) < 5) stop("Fewer than 5 arms/features survive amplicon QC")

  # 1. SUM raw counts per arm  2. fractions across the retained arms (after arm removal)
  all_counts <- rbind(spk[, amps, drop = FALSE], tum[, amps, drop = FALSE])
  S    <- do.call(cbind, lapply(arms, function(a) rowSums(all_counts[, a, drop = FALSE])))
  frac <- S / rowSums(S)
  # 3. divide by spike-in median fraction per arm, 4. log2
  is_sp <- rownames(frac) %in% rownames(spk)
  ref   <- apply(frac[is_sp, , drop = FALSE], 2, median)
  ok    <- ref > 0
  frac <- frac[, ok, drop = FALSE]; ref <- ref[ok]; arms <- arms[ok]
  ratio <- sweep(frac, 2, ref, "/")
  ratio[ratio < cfg$min_ratio] <- cfg$min_ratio
  lr <- log2(ratio)

  list(arm_lr = t(lr[!is_sp, , drop = FALSE]),      # arms x tumour cells
       n_amplicons = lengths(arms),
       arm_amplicons = arms,                      # list: arm -> retained amplicons
       amplicon_qc = amp_qc)
}

# ---- mutation selection (NA-safe) -------------------------------------------
select_variants <- function(mut, n_cells_file, sh, cfg) {
  num <- c("AF","AF_popmax","Genotyped_cells","HET_cells","WT_cells","HOM_cells","VAF_by_cell_count")
  for (cl in intersect(num, colnames(mut))) mut[[cl]] <- suppressWarnings(as.numeric(mut[[cl]]))
  mut$AF_popmax[is.na(mut$AF_popmax)] <- 0
  mut$AF[is.na(mut$AF)] <- 0

  is_onc <- mut$ONCOGENIC %in% cfg$onc_classes
  mut$id <- ifelse(is_onc, paste(mut$ONCOKB_HUGO_SYMBOL, mut$ONCOKB_PROTEIN_CHANGE, sep = "_"), mut$variant)
  mut$maf_position <- paste(mut$Chromosome, mut$Start_Position, mut$End_Position,
                            mut$Reference_Allele, mut$Tumor_Seq_Allele2, sep = "_")

  tot <- mut$WT_cells + mut$HET_cells + mut$HOM_cells
  bases <- c("A","C","G","T")
  germ <- if (is.na(cfg$germline_vaf_lo) || is.na(cfg$germline_vaf_hi)) rep(TRUE, nrow(mut)) else
    (mut$VAF_by_cell_count > cfg$germline_vaf_hi | mut$VAF_by_cell_count < cfg$germline_vaf_lo)

  pass_main <- mut$Genotyped_cells > n_cells_file * cfg$min_genotyped_frac &
    mut$AF > cfg$min_AF &
    (mut$HET_cells / tot) > cfg$min_het_frac & germ &
    mut$Reference_Allele %in% bases & mut$Tumor_Seq_Allele2 %in% bases
  pass_rescue <- (is_onc | mut$maf_position %in% sh$idt$maf_position) &
    (mut$HET_cells + mut$HOM_cells) >= cfg$min_rescue_cells
  pass_main[is.na(pass_main)] <- FALSE
  pass_rescue[is.na(pass_rescue)] <- FALSE

  mut <- mut[pass_main | pass_rescue, , drop = FALSE]
  mut <- mut[mut$id %nin% c(sh$black$id, cfg$extra_blacklist_ids), , drop = FALSE]
  mut <- mut[!duplicated(mut$variant), , drop = FALSE]
  data.frame(variant = mut$variant, id = mut$id,
             label = make.unique(as.character(mut$id)),
             oncogenic = mut$ONCOGENIC %in% cfg$onc_classes,
             stringsAsFactors = FALSE)
}

impute_median <- function(m) {
  for (j in seq_len(ncol(m))) {
    na <- is.na(m[, j])
    if (any(na)) {
      med <- if (all(na)) 0 else median(m[!na, j])
      m[na, j] <- med
    }
  }
  m
}

# ---- feature matrix: cells x (CN log2 + weighted VAF) ------------------------
build_features <- function(cn, vaf_imp, vaf_scale, cfg) {
  cn[cn >  cfg$lr_clip] <-  cfg$lr_clip
  cn[cn < -cfg$lr_clip] <- -cfg$lr_clip
  x <- t(cn)
  if (ncol(vaf_imp) > 0) {
    v <- vaf_imp[rownames(x), , drop = FALSE] / vaf_scale
    v <- v[, which(apply(v, 2, sd) > 0), drop = FALSE]
    if (ncol(v) > 0) x <- cbind(x, v * (cfg$vaf_weight * ncol(x) / ncol(v)))
  }
  x
}

# ---- UMAP (manhattan, min_dist = 0) + HDBSCAN --------------------------------
embed_cluster <- function(feat, cfg) {
  n  <- nrow(feat)
  nn <- max(5, min(cfg$n_neighbors, floor(n / 2)))
  mp <- max(5, min(cfg$min_pts, floor(n / 5)))
  set.seed(cfg$seed)
  emb <- uwot::umap(feat, n_neighbors = nn, metric = "manhattan", min_dist = 0, n_sgd_threads = 1)
  rownames(emb) <- rownames(feat); colnames(emb) <- c("UMAP1", "UMAP2")
  hd  <- dbscan::hdbscan(emb, minPts = mp)
  lab <- ifelse(hd$cluster == 0, "unassigned", as.character(hd$cluster))
  ids <- sort(as.numeric(setdiff(unique(lab), "unassigned")))
  lv  <- c(if ("unassigned" %in% lab) "unassigned", as.character(ids))
  cl  <- factor(lab, levels = lv); names(cl) <- rownames(feat)
  list(umap = as.data.frame(emb), cluster = cl)
}

# ---- normal-cluster call -----------------------------------------------------
call_normal <- function(cn, cl, cfg, override = NULL) {
  none <- function(st) list(cluster = NA_character_, members = character(0), stats = st,
                            ref = NULL, method = "none")
  lv <- setdiff(levels(droplevels(cl)), "unassigned")
  if (!length(lv)) return(none(NULL))
  profs <- setNames(lapply(lv, function(k) apply(cn[, names(cl)[cl == k], drop = FALSE], 1, median)), lv)
  st <- do.call(rbind, lapply(lv, function(k) {
    prof <- profs[[k]]
    data.frame(cluster = k, n_cells = sum(cl == k),
               profile_sd = sd(prof), mean_abs_log2 = mean(abs(prof)),
               n_altered_arms = sum(abs(prof) > cfg$normal_alt_threshold),
               stringsAsFactors = FALSE)
  }))
  st$candidate <- st$n_cells >= cfg$normal_min_cells

  # k = best normal cluster; other clusters with an (almost) identical CN profile are merged
  make <- function(k, how) {
    d   <- vapply(lv, function(j) max(abs(profs[[j]] - profs[[k]])), numeric(1))
    mem <- union(k, if (cfg$normal_merge_max_diff > 0) lv[d <= cfg$normal_merge_max_diff] else character(0))
    st2 <- st
    st2$max_diff_to_normal <- unname(d)
    st2$called_normal <- st2$cluster %in% mem
    list(cluster = k, members = mem, stats = st2, method = how,
         ref = apply(cn[, names(cl)[cl %in% mem], drop = FALSE], 1, median))
  }

  if (length(override) == 1 && !is.na(override) && nzchar(override)) {
    if (override %in% st$cluster) return(make(override, "manual_override"))
    warning("normal_override cluster '", override, "' not found; using automatic call")
  }
  cand <- st[st$candidate, , drop = FALSE]
  if (!nrow(cand)) return(none(st))
  best   <- cand[which.min(cand$profile_sd), ]
  others <- st$profile_sd[st$cluster != best$cluster]
  ratio  <- if (length(others)) best$profile_sd / median(others) else NA_real_
  clear  <- best$profile_sd <= cfg$normal_max_sd || (!is.na(ratio) && ratio <= cfg$normal_rel_sd)
  if (isTRUE(cfg$normal_require_flat) && !clear) return(none(st))
  res <- make(best$cluster, if (clear) "lowest_sd_clear" else "lowest_sd_low_confidence")
  res$stats$best_sd_ratio <- ratio          # SD of chosen cluster / median SD of the others
  res
}

# ---- putative doublet clusters ------------------------------------------------
# A doublet cluster (normal + tumour cell) should look like a MIXTURE of the normal
# cluster and a parent tumour cluster: intermediate CN deviations, intermediate VAF
# of tumour variants, and more reads per cell. Evaluated for every non-normal cluster
# against every larger-or-equal non-normal cluster; best parent is reported.
detect_doublet_clusters <- function(cn, cl, vaf_imp, vaf_scale, reads, normal, cfg) {
  if (!length(normal) || all(is.na(normal))) return(NULL)
  lv <- setdiff(levels(droplevels(cl)), c("unassigned", normal))
  if (length(lv) < 2) return(NULL)
  has_vaf <- ncol(vaf_imp) > 0
  idx <- setNames(lapply(lv, function(k) names(cl)[cl == k]), lv)
  idx[[".normal"]] <- names(cl)[cl %in% normal]      # all normal clusters pooled
  rr  <- lapply(idx, function(i) 2^apply(cn[, i, drop = FALSE], 1, median))        # linear ratio vs normal
  vv  <- if (has_vaf) lapply(idx, function(i) colMeans(vaf_imp[i, , drop = FALSE]) / vaf_scale) else NULL
  wr  <- cfg$doublet_w_range
  in_rng <- function(w) !is.na(w) && w >= wr[1] && w <= wr[2]

  out <- lapply(lv, function(k) {
    cand <- lapply(setdiff(lv, k), function(j) {
      if (length(idx[[j]]) < length(idx[[k]])) return(NULL)        # parent must be at least as large
      # CN
      dj <- rr[[j]] - 1; dk <- rr[[k]] - 1; inf <- abs(dj) >= cfg$doublet_min_arm_delta
      cn_eval <- sum(inf) >= cfg$doublet_min_arms
      cn_cor <- if (cn_eval) suppressWarnings(cor(dk[inf], dj[inf])) else NA_real_
      cn_w   <- if (cn_eval) sum(dk[inf] * dj[inf]) / sum(dj[inf]^2) else NA_real_
      cn_pass <- if (cn_eval) (!is.na(cn_cor) && cn_cor >= cfg$doublet_min_cn_cor && in_rng(cn_w)) else NA
      # VAF
      vaf_w <- NA_real_; n_inf_v <- 0; vaf_pass <- NA
      if (has_vaf) {
        dvj <- vv[[j]] - vv[[".normal"]]; dvk <- vv[[k]] - vv[[".normal"]]
        iv <- abs(dvj) >= cfg$doublet_min_vaf_delta; n_inf_v <- sum(iv)
        if (n_inf_v >= 1) { vaf_w <- sum(dvk[iv] * dvj[iv]) / sum(dvj[iv]^2); vaf_pass <- in_rng(vaf_w) }
      }
      # reads
      ref_cells <- c(idx[[".normal"]], idx[[j]])
      rk <- reads[idx[[k]]]; rf <- reads[ref_cells]
      read_ratio <- median(rk) / median(rf)
      read_p <- tryCatch(wilcox.test(rk, rf, alternative = "greater")$p.value, error = function(e) NA_real_)
      evid <- c(cn_pass, vaf_pass)
      need_ratio <- if (all(evid %in% TRUE)) cfg$doublet_min_read_ratio_both else cfg$doublet_min_read_ratio
      reads_pass <- !is.na(read_p) && read_ratio >= need_ratio && read_p < cfg$doublet_reads_p
      is_d <- reads_pass && any(!is.na(evid)) && !any(evid %in% FALSE)
      data.frame(cluster = k, n_cells = length(idx[[k]]), median_reads = median(reads[idx[[k]]]),
                 parent = j, cn_cor = cn_cor, cn_w = cn_w, vaf_w = vaf_w, n_informative_vaf = n_inf_v,
                 read_ratio = read_ratio, read_ratio_needed = need_ratio,
                 read_ratio_vs_normal = median(rk) / median(reads[idx[[".normal"]]]),
                 read_ratio_vs_parent = median(rk) / median(reads[idx[[j]]]),
                 read_p = read_p,
                 cn_pass = cn_pass, vaf_pass = vaf_pass, reads_pass = reads_pass,
                 is_doublet = is_d, stringsAsFactors = FALSE)
    })
    cand <- do.call(rbind, cand)
    if (is.null(cand)) return(NULL)
    cand[order(!cand$is_doublet, -ifelse(is.na(cand$cn_cor), 0, cand$cn_cor)), ][1, ]
  })
  do.call(rbind, out)
}

# ---- plotting ----------------------------------------------------------------
plot_umap <- function(df, pal, title = NULL) {
  ggscatter(df, x = "UMAP1", y = "UMAP2", color = "cluster", alpha = 0.7, size = 1.2, title = title) +
    scale_color_manual(name = "cluster", values = pal)
}

plot_heatmap <- function(cn, cl, vaf, pal, file, cfg, title = NULL, tags = NULL) {
  cells <- names(cl)[order(cl)]
  cl_o  <- cl[cells]
  mat <- cn[, cells, drop = FALSE]
  mat[mat >  cfg$heat_lim] <-  cfg$heat_lim
  mat[mat < -cfg$heat_lim] <- -cfg$heat_lim
  ann <- data.frame(cluster = cl_o, row.names = cells)
  if (ncol(vaf) > 0) ann <- data.frame(ann, vaf[cells, , drop = FALSE], check.names = FALSE)

  # cluster names printed under the middle of each cluster block (no gaps)
  runs   <- rle(as.character(cl_o))
  ends   <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  lab    <- rep("", length(cells))
  nm  <- runs$values
  if (length(tags)) nm <- ifelse(nm %in% names(tags), paste0(nm, " (", tags[nm], ")"), nm)
  big <- runs$lengths >= 0.01 * length(cells)        # tiny blocks (e.g. a few unassigned cells) would be clipped
  lab[floor((starts + ends) / 2)[big]] <- nm[big]

  ph <- pheatmap::pheatmap(mat, color = heat_palette,
                           breaks = seq(-cfg$heat_lim, cfg$heat_lim, length.out = length(heat_palette) + 1),
                           cluster_cols = FALSE, cluster_rows = FALSE,
                           show_colnames = TRUE, labels_col = lab, angle_col = 0, fontsize_col = 11,
                           annotation_col = ann, annotation_colors = list(cluster = pal),
                           border_color = NA, main = title, silent = TRUE)
  # page height grows with the number of annotation legends so none (incl. "cluster") is cut off
  h <- max(8, 2.5 + 1.1 * ncol(vaf))
  grDevices::pdf(file, width = 10, height = h)
  on.exit(grDevices::dev.off())
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  invisible(NULL)
}

plot_vaf_panels <- function(df_umap, vaf_raw, pal, vaf_scale, fp, ggs) {
  long <- data.frame(df_umap, vaf_raw[rownames(df_umap), , drop = FALSE], check.names = FALSE) %>%
    pivot_longer(all_of(colnames(vaf_raw)), names_to = "mutation", values_to = "vaf")
  nf <- ncol(vaf_raw); nc <- min(4, nf); nr <- ceiling(nf / nc)
  p1 <- ggscatter(long, x = "UMAP1", y = "UMAP2", color = "vaf", alpha = 0.8, size = 0.4,
                  facet.by = "mutation", ncol = nc) +
    scale_color_gradient(low = "grey90", high = "red", na.value = "white")
  ggs(p1, "VAF_ON_UMAP", 3 * nc + 1, 3 * nr + 1)

  med <- long %>% group_by(mutation, cluster) %>% summarise(m = median(vaf, na.rm = TRUE), .groups = "drop")
  med$lab <- format(round(med$m, ifelse(vaf_scale == 1, 2, 0)))
  p2 <- ggboxplot(long, x = "cluster", y = "vaf", fill = "cluster", outlier.shape = NA,
                  facet.by = "mutation", ncol = nc) +
    geom_jitter(position = position_jitter(0.1), alpha = 0.1, size = 0.4) +
    geom_text(data = med, aes(x = cluster, y = m, label = lab), vjust = -0.5, colour = "red",
              size = 3, inherit.aes = FALSE) +
    scale_fill_manual(name = "cluster", values = pal)
  ggs(p2, "VAF_BOXPLOT", 3 * nc + 1, 3 * nr + 1)
}

# =============================================================================
# 2. PER-SAMPLE ANALYSIS
# =============================================================================
run_sample <- function(name, cfg, sh) {
  t0  <- Sys.time()
  out <- file.path(cfg$outdir, name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  fp  <- function(suffix, ext = "pdf") file.path(out, paste0(name, "_", suffix, ".", ext))
  ggs <- function(p, suffix, w, h) ggsave(fp(suffix), p, width = w, height = h, limitsize = FALSE)
  message("[", name, "] start")

  # ---- 1. load -------------------------------------------------------------
  tum_raw <- as.matrix(read.delim(file.path(cfg$data_dir, paste0(name, ".dna.CNV.tsv")), row.names = 1))
  vaf_all <- read.delim(file.path(cfg$data_dir, paste0(name, ".dna.AF.tsv")), check.names = FALSE, row.names = 1)
  mut     <- read.delim(file.path(cfg$data_dir, paste0(name, ".dna.var_stat.tsv.vcf.hg19_multianno.txt.maf.oncoKB.txt")),
                        check.names = FALSE)

  # ---- 2. cell-level QC (reads, doublets, amplicon detection) -----------------
  qc_t <- qc_cells(tum_raw, cfg)
  write.table(rbind(cbind(group = "tumour", qc_t),
                    cbind(group = "spikein", sh$spike_qc)),
              fp("CELL_QC", "tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cells0 <- intersect(qc_t$cell[qc_t$qc_pass], rownames(vaf_all))
  message("[", name, "] cells: ", nrow(tum_raw), " in CNV file; ", sum(qc_t$qc_pass),
          " pass QC; ", length(cells0), " also present in AF file")
  if (length(cells0) < cfg$min_cells_sample) stop("Too few cells after QC: ", length(cells0))

  # ---- 3. amplicon-level log2 ratios -> arm level ----------------------------
  cnv <- compute_arm_lr(tum_raw[cells0, , drop = FALSE], sh$spike, sh$ann, cfg)
  write.table(cnv$amplicon_qc, fp("AMPLICON_QC", "tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cn_ss <- cnv$arm_lr                                         # arms x cells, spike-in referenced
  b_full <- setNames(rep(0, nrow(cn_ss)), rownames(cn_ss))    # learned spike-in/sample arm bias
  if (!is.null(sh$arm_bias)) {
    bb <- sh$arm_bias[rownames(cn_ss)]; bb[is.na(bb)] <- 0
    b_full[] <- bb
    cn_ss <- sweep(cn_ss, 1, b_full, "-")
  }

  # ---- 4. mutations ----------------------------------------------------------
  feats <- select_variants(mut, nrow(tum_raw), sh, cfg)
  feats <- feats[feats$variant %in% colnames(vaf_all), , drop = FALSE]
  vaf_raw <- as.matrix(vaf_all[cells0, feats$variant, drop = FALSE])
  colnames(vaf_raw) <- feats$label
  has_vaf <- ncol(vaf_raw) > 0
  if (has_vaf) {
    keep_cells <- rowMeans(is.na(vaf_raw)) <= cfg$max_na_vaf_frac
    message("[", name, "] ", sum(!keep_cells), " cells dropped (too many ungenotyped variants)")
    vaf_raw <- vaf_raw[keep_cells, , drop = FALSE]
    cn_ss   <- cn_ss[, rownames(vaf_raw), drop = FALSE]
  } else {
    message("[", name, "] no variants selected - clustering on CN only")
  }
  cells <- colnames(cn_ss)
  if (length(cells) < cfg$min_cells_sample) stop("Too few cells after VAF QC: ", length(cells))
  mx <- suppressWarnings(max(vaf_raw, na.rm = TRUE))
  vaf_scale <- if (is.finite(mx) && mx > 1) 100 else 1
  vaf_imp <- impute_median(vaf_raw)

  # ---- 5. UMAP + HDBSCAN (all tumour-sample cells) ---------------------------
  feat <- build_features(cn_ss, vaf_imp, vaf_scale, cfg)
  e1   <- embed_cluster(feat, cfg)
  cl   <- e1$cluster
  pal  <- make_palette(levels(cl))

  # ---- 6. normal cluster call and re-referencing -----------------------------
  nrm <- call_normal(cn_ss, cl, cfg, override = unname(cfg$normal_override[name]))
  if (!is.null(nrm$stats))
    write.table(nrm$stats, fp("CLUSTER_STATS", "tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (is.na(nrm$cluster)) {
    message("[", name, "] no normal cluster identified -> spike-in-only reference (see CLUSTER_STATS.tsv)")
    cn <- cn_ss; ref_txt <- "spikein_only"; ref_title <- "log2 ratio vs spike-in (no normal cluster found)"
  } else {
    mem_txt <- paste(nrm$members, collapse = "+")
    cn <- sweep(cn_ss, 1, nrm$ref, "-"); ref_txt <- paste0("cluster_", mem_txt, "_", nrm$method)
    ref_title <- paste0("log2 ratio vs normal cluster ", mem_txt)
    message("[", name, "] normal cluster(s) = ", mem_txt, " (", nrm$method, ")")
  }
  is_normal <- as.character(cl) %in% nrm$members
  ctype <- ifelse(as.character(cl) == "unassigned", "unassigned", ifelse(is_normal, "normal", "tumor"))
  names(ctype) <- names(cl)

  # putative doublet clusters (intermediate CN + VAF between normal and a tumour cluster, more reads)
  dbl <- NULL; dcl <- character(0)
  if (isTRUE(cfg$detect_doublets)) {
    reads_v <- setNames(qc_t[names(cl), "total_reads"], names(cl))
    dbl <- detect_doublet_clusters(cn, cl, vaf_imp, vaf_scale, reads_v, nrm$members, cfg)
    if (!is.null(dbl)) {
      write.table(dbl, fp("DOUBLET_STATS", "tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      dcl <- dbl$cluster[dbl$is_doublet]
      ctype[as.character(cl) %in% dcl] <- "putative_doublet"
      if (length(dcl)) message("[", name, "] putative doublet cluster(s): ", paste(dcl, collapse = ", "))
    }
  }
  tags <- c(if (length(nrm$members)) setNames(rep("normal", length(nrm$members)), nrm$members),
            if (length(dcl)) setNames(rep("doublet", length(dcl)), dcl))

  cell_tbl <- data.frame(cell = names(cl),
                         total_reads = qc_t[names(cl), "total_reads"],
                         frac_detected = qc_t[names(cl), "frac_detected"],
                         cluster = as.character(cl), cell_type = ctype[names(cl)],
                         e1$umap[names(cl), , drop = FALSE],
                         tumor_cluster = NA_character_,
                         tumor_UMAP1 = NA_real_, tumor_UMAP2 = NA_real_,
                         row.names = names(cl), stringsAsFactors = FALSE)

  # ---- 7. plots: full analysis -----------------------------------------------
  df_u <- data.frame(e1$umap, cluster = cl)
  ggs(plot_umap(df_u, pal), "UMAP_CLUSTERED", 9, 9)

  pct <- as.data.frame(table(cluster = cl), responseName = "Freq")
  pct$Percentage <- round(pct$Freq / sum(pct$Freq) * 100, 1)
  ggs(ggbarplot(pct, x = "cluster", y = "Percentage", fill = "cluster", xlab = "Cluster", ylim = c(0, 100)) +
        geom_text(aes(label = Percentage), vjust = -0.5, colour = "red", size = 5) +
        scale_fill_manual(name = "cluster", values = pal), "CLUSTER_PERCENT", 7, 7)

  rd <- data.frame(cluster = cl, total_reads = cell_tbl$total_reads)
  if (nlevels(droplevels(rd$cluster)) > 1) {
    kw  <- kruskal.test(total_reads ~ cluster, data = rd)
    med <- rd %>% group_by(cluster) %>% summarise(m = median(total_reads), .groups = "drop")
    ggs(ggboxplot(rd, x = "cluster", y = "total_reads", fill = "cluster", outlier.shape = NA) +
          geom_jitter(position = position_jitter(0.1), alpha = 0.1, size = 0.5) +
          geom_text(data = med, aes(x = cluster, y = m, label = round(m)), vjust = -0.5,
                    colour = "red", size = 3.5, inherit.aes = FALSE) +
          scale_fill_manual(name = "cluster", values = pal) + scale_y_log10() +
          labs(y = "total reads / cell", subtitle = sprintf("Kruskal-Wallis p = %.2g", kw$p.value)),
        "READS_PER_CLUSTER", 8, 8)
    pw <- tryCatch({
      p1 <- rstatix::pairwise_wilcox_test(rd %>% mutate(cluster = droplevels(cluster)),
                                          total_reads ~ cluster, p.adjust.method = "BH")
      ef <- tryCatch(rstatix::wilcox_effsize(rd %>% mutate(cluster = droplevels(cluster)), total_reads ~ cluster),
                     error = function(e) NULL)
      if (!is.null(ef)) p1 <- left_join(p1, ef[, c("group1", "group2", "effsize", "magnitude")],
                                        by = c("group1", "group2"))
      p1
    }, error = function(e) NULL)
    if (!is.null(pw)) write.table(pw, fp("READS_PAIRWISE_WILCOXON", "tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if (has_vaf) plot_vaf_panels(df_u, vaf_raw, pal, vaf_scale, fp, ggs)

  plot_heatmap(cn, cl, vaf_raw, pal, fp("PHEATMAP"), cfg, title = paste(name, "-", ref_title), tags = tags)

  # ---- 8. tumour-only re-clustering ------------------------------------------
  cl2 <- NULL
  tum_cells <- names(ctype)[ctype == "tumor" | (ctype == "unassigned" & !cfg$exclude_unassigned_from_tumor)]
  if (length(tum_cells) >= cfg$min_cells_tumor_recluster && length(tum_cells) < length(cells)) {
    feat_t <- build_features(cn[, tum_cells, drop = FALSE], vaf_imp, vaf_scale, cfg)
    e2   <- embed_cluster(feat_t, cfg)
    cl2  <- e2$cluster
    pal2 <- make_palette(levels(cl2))
    ggs(plot_umap(data.frame(e2$umap, cluster = cl2), pal2), "UMAP_TUMOR_ONLY", 9, 9)
    plot_heatmap(cn[, tum_cells, drop = FALSE], cl2, vaf_raw, pal2, fp("PHEATMAP_TUMOR_ONLY"), cfg,
                 title = paste(name, "- tumour cells only"))
    cell_tbl[names(cl2), "tumor_cluster"] <- as.character(cl2)
    cell_tbl[names(cl2), c("tumor_UMAP1", "tumor_UMAP2")] <- e2$umap[names(cl2), ]
  } else {
    message("[", name, "] tumour-only re-clustering skipped (", length(tum_cells), " cells)")
  }

  # ---- 9. CN significance + summary tables (full and tumour-only) -----------
  cn_calls <- NULL; cn_calls_t <- NULL
  normal_cells <- names(cl)[as.character(cl) %in% nrm$members]
  if (length(normal_cells)) {
    kept_amps <- cnv$amplicon_qc$amplicon[cnv$amplicon_qc$kept]
    cnt <- tum_raw[colnames(cn), kept_amps, drop = FALSE]
    skip <- function(e) { message("[", name, "] CN significance skipped: ", conditionMessage(e)); NULL }
    cn_calls <- tryCatch(call_cn_significance(cn, cl, ctype, normal_cells, cnv$n_amplicons, cfg,
                                              counts = cnt, arm_amps = cnv$arm_amplicons), error = skip)
    if (!is.null(cn_calls) && !is.null(cl2)) {                  # tumour sub-clusters vs. the same normal cells
      ctype_t <- setNames(rep("tumor", length(cl2)), names(cl2))
      cn_calls_t <- tryCatch(call_cn_significance(cn, cl2, ctype_t, normal_cells, cnv$n_amplicons, cfg,
                                                  counts = cnt, arm_amps = cnv$arm_amplicons,
                                                  null_arm = attr(cn_calls, "null_arm")), error = skip)
    }
  } else {
    message("[", name, "] CN significance calls not made (no normal reference)")
  }
  if (!is.null(cn_calls))
    write.table(cbind(Sample = name, cn_calls), fp("CN_CALLS", "tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(cn_calls_t))
    write.table(cbind(Sample = name, cn_calls_t), fp("CN_CALLS_TUMOR_ONLY", "tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  write_summary_table(fp("SUMMARY_TABLE", "tsv"), name, cell_tbl, "cluster", cn, vaf_raw, feats, ref_txt, cn_calls)
  if (!is.null(cl2))                                             # e.g. SAMPLE_A_SUMMARY_TABLE_TUMOR_ONLY.tsv
    write_summary_table(fp("SUMMARY_TABLE_TUMOR_ONLY", "tsv"), name, cell_tbl, "tumor_cluster", cn, vaf_raw,
                        feats, ref_txt, cn_calls_t)

  # ---- 10. per-cell table + RDS ----------------------------------------------
  wide <- cbind(cell_tbl,
                setNames(as.data.frame(t(cn[, cell_tbl$cell, drop = FALSE])), paste0(rownames(cn), "_log2ratio")),
                as.data.frame(vaf_raw[cell_tbl$cell, , drop = FALSE], check.names = FALSE))
  write.table(wide, fp("PER_CELL_TABLE", "tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(list(cfg = cfg, cell_table = cell_tbl, cn_log2 = cn, cn_spikein_ref = cn_ss,
               n_amplicons_per_arm = cnv$n_amplicons, variants = feats, vaf_raw = vaf_raw,
               normal = nrm, doublets = dbl, bias_applied = b_full,
               normal_profile_raw = if (!is.null(nrm$ref) && !grepl("low_confidence", nrm$method)) nrm$ref + b_full[names(nrm$ref)] else NULL,
               session = sessionInfo()),
          fp("RESULTS", "rds"))

  message("[", name, "] done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s")
  data.frame(sample = name, status = "OK", cells_in_file = nrow(tum_raw), cells_analysed = length(cells),
             n_clusters = nlevels(droplevels(cl)) - ("unassigned" %in% levels(droplevels(cl))),
             normal_reference = ref_txt, n_variants = nrow(feats),
             doublet_clusters = paste(dcl, collapse = ","), stringsAsFactors = FALSE)
}

# ---- summary table per cluster (used for the full AND the tumour-only analysis) ------
# cells_df: one row per cell with columns cell, total_reads, cell_type and the clustering column
# `cluster_col` (NA = cell not in this analysis). Columns: counts, median reads, median log2 ratio per
# arm vs. normal, CN calls (if provided), median VAF and % cells with VAF > 0 for oncogenic variants.
write_summary_table <- function(file, name, cells_df, cluster_col, cn, vaf_raw, feats, ref_txt,
                                cn_calls = NULL) {
  d <- cells_df
  d$cl <- as.character(d[[cluster_col]])
  d <- d[!is.na(d$cl), , drop = FALSE]
  summ <- d %>% group_by(cluster = cl) %>%
    summarise(Cells_in_cluster = n(), Cell_type = first(cell_type),
              Median_Reads = median(total_reads), .groups = "drop") %>%
    mutate(Sample = name, Total_cells = nrow(d), Normal_reference = ref_txt) %>%
    select(Sample, Total_cells, cluster, Cells_in_cluster, Cell_type, Median_Reads, Normal_reference)
  summ <- summ[order(suppressWarnings(as.numeric(summ$cluster)), summ$cluster), ]    # 1,2,...,10; unassigned last

  cn_df <- as.data.frame(t(cn[, d$cell, drop = FALSE]), check.names = FALSE)
  cn_df$cluster <- d$cl
  cn_med <- cn_df %>% group_by(cluster) %>% summarise(across(everything(), median), .groups = "drop") %>%
    rename_with(~ paste0(.x, "_log2ratio"), -cluster)
  summ <- left_join(summ, cn_med, by = "cluster")

  if (!is.null(cn_calls)) {
    wide_calls <- cn_calls %>% select(cluster, arm, call) %>%
      pivot_wider(names_from = arm, values_from = call) %>%
      rename_with(~ paste0(.x, "_call"), -cluster)
    summ <- left_join(summ, wide_calls, by = "cluster")
  }

  onc <- feats$label[feats$oncogenic]
  if (length(onc) > 0) {
    vd <- data.frame(cluster = d$cl, vaf_raw[d$cell, onc, drop = FALSE], check.names = FALSE)
    vmed <- vd %>% group_by(cluster) %>%
      summarise(across(all_of(onc), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
      rename_with(~ paste0(.x, "_VAF"), -cluster)
    vpct <- vd %>% group_by(cluster) %>%
      summarise(across(all_of(onc), ~ round(mean(.x > 0, na.rm = TRUE) * 100, 1)), .groups = "drop") %>%
      rename_with(~ paste0(.x, "_percent"), -cluster)
    summ <- summ %>% left_join(vmed, by = "cluster") %>% left_join(vpct, by = "cluster")
  }
  write.table(summ, file, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(summ)
}

# ---- how much "difference" does clustering alone create among NORMAL cells? ----
pseudocluster_null <- function(cn, nc, cfg) {
  arms <- rownames(cn)
  null_arm <- setNames(rep(NA_real_, length(arms)), arms)
  if (isTRUE(cfg$sig_pseudocluster_null) && length(nc) >= 2 * cfg$sig_min_group) {
    pc <- tryCatch({
      f0 <- matrix(0, nrow = length(nc), ncol = 0, dimnames = list(nc, NULL))
      embed_cluster(build_features(cn[, nc, drop = FALSE], f0, 1, cfg), cfg)$cluster
    }, error = function(e) NULL)
    if (!is.null(pc)) {
      g <- setdiff(levels(droplevels(pc)), "unassigned")
      g <- g[vapply(g, function(x) sum(pc == x), numeric(1)) >= cfg$sig_min_group]
      if (length(g) >= 2) {
        M <- sapply(g, function(x) apply(cn[, names(pc)[pc == x], drop = FALSE], 1, median))
        null_arm[] <- apply(M, 1, function(v) max(v) - min(v))
      }
    }
  }
  null_arm
}

# ---- significance of arm-level CN shifts vs. normal (tables only) --------------
# For every non-unassigned cluster and arm:
#   shift       = median(cluster) - median(normal cells)               (log2)
#   CI / p      = cell-level bootstrap CI and BH-adjusted Wilcoxon vs normal cells. These only
#                 reflect cell sampling; with thousands of cells nearly everything passes.
#   amplicon    = shift of every amplicon of the arm (log2 of cluster median / normal median of
#                 depth-scaled counts). Amplicons are the real replicates of an arm: one-sample
#                 t-test of these shifts against 0 (arm needs >= 3 amplicons) and the smallest
#                 leave-one-amplicon-out mean must keep the sign of the arm shift. This stops a
#                 single shifted amplicon (e.g. in a 3-amplicon arm) from producing a "loss".
#   threshold   = max(sig_min_effect, pseudo-cluster null of that arm). The null comes from
#                 clustering the normal cells alone with the same UMAP/HDBSCAN settings: the range
#                 of arm medians across the pseudo-clusters = difference manufactured by clustering.
call_cn_significance <- function(cn, cl, ctype, normal_cells, n_amp, cfg,
                                 counts = NULL, arm_amps = NULL, null_arm = NULL) {
  set.seed(cfg$seed)
  cap  <- function(x) if (length(x) > cfg$sig_max_cells) sample(x, cfg$sig_max_cells) else x
  nc   <- normal_cells                      # reference = ALL normal cells, whatever clusters are tested
  if (length(nc) < cfg$normal_min_cells) return(NULL)
  arms <- rownames(cn)
  nv   <- cn[, cap(nc), drop = FALSE]
  med_n <- apply(cn[, nc, drop = FALSE], 1, median)

  # (a) null from pseudo-clusters of the normal cells (can be passed in to avoid recomputing)
  if (is.null(null_arm)) null_arm <- pseudocluster_null(cn, nc, cfg)
  thr <- pmax(cfg$sig_min_effect, null_arm, na.rm = TRUE)

  # (b) bootstrap of the normal medians (reused for every cluster)
  B <- cfg$sig_n_boot
  boot_n <- t(replicate(B, apply(nv[, sample.int(ncol(nv), replace = TRUE), drop = FALSE], 1, median)))

  # (c) amplicon level: depth-scaled counts (pseudocount keeps zeros finite)
  use_amp <- !is.null(counts) && !is.null(arm_amps)
  if (use_amp) {
    lib <- rowSums(counts)
    cs  <- sweep(counts + cfg$pseudocount, 1, lib, "/") * median(lib)
    med_n_amp <- apply(cs[nc, , drop = FALSE], 2, median)
  }
  amp_cols <- c("n", "mean", "p", "conc", "loo")

  ks <- setdiff(levels(droplevels(cl)), "unassigned")
  res <- lapply(ks, function(k) {
    ck    <- names(cl)[cl == k]
    xk    <- cn[, cap(ck), drop = FALSE]
    delta <- apply(cn[, ck, drop = FALSE], 1, median) - med_n
    boot_k <- t(replicate(B, apply(xk[, sample.int(ncol(xk), replace = TRUE), drop = FALSE], 1, median)))
    ci   <- apply(boot_k - boot_n, 2, quantile, probs = c(0.025, 0.975))
    p    <- vapply(arms, function(a) tryCatch(suppressWarnings(wilcox.test(xk[a, ], nv[a, ])$p.value),
                                              error = function(e) NA_real_), numeric(1))
    padj <- p.adjust(p, "BH")
    sig_cells <- !is.na(padj) & padj < cfg$sig_alpha & (ci[1, ] > 0 | ci[2, ] < 0)

    if (use_amp) {
      s_k <- log2(apply(cs[ck, , drop = FALSE], 2, median) / med_n_amp)       # shift per amplicon
      at <- t(vapply(arms, function(a) {
        amp <- intersect(arm_amps[[a]], names(s_k)); n <- length(amp)
        if (n == 0) return(setNames(c(0, NA, NA, NA, NA), amp_cols))
        v <- s_k[amp]; m <- mean(v)
        pv <- if (n >= 3) tryCatch(suppressWarnings(t.test(v)$p.value), error = function(e) NA_real_) else NA_real_
        loo <- if (n >= 2) { lm <- (sum(v) - v) / (n - 1); lm[which.min(abs(lm))] } else m   # least extreme leave-one-out mean
        setNames(c(n, m, pv, mean(sign(v) == sign(m)), loo), amp_cols)
      }, numeric(5)))
      amp_n  <- at[, "n"]
      amp_ok <- (!is.na(at[, "p"]) & at[, "p"] < cfg$sig_amp_alpha & sign(at[, "loo"]) == sign(delta)) %in% TRUE
    } else {
      at <- matrix(NA_real_, nrow = length(arms), ncol = 5, dimnames = list(arms, amp_cols))
      amp_n <- rep(Inf, length(arms)); amp_ok <- rep(TRUE, length(arms))
    }

    scale_k <- median(sort(abs(delta), decreasing = TRUE)[seq_len(min(cfg$sig_scale_top_n, length(delta)))])
    thr_k   <- pmax(thr, cfg$sig_rel_effect * scale_k)           # cluster-adaptive threshold
    big <- abs(delta) >= thr_k
    cn_call <- dplyr::case_when(
      sig_cells & big & amp_n >= 3 & amp_ok  ~ ifelse(delta > 0, "gain", "loss"),
      sig_cells & big & amp_n < 3            ~ "uncertain_few_amplicons",
      sig_cells & big                        ~ "uncertain_amplicon_inconsistent",
      sig_cells & !big                       ~ "uncertain_small_effect",
      !sig_cells & big                       ~ "uncertain_not_significant",
      TRUE                                   ~ "neutral")
    data.frame(cluster = k, cell_type = unname(ctype[ck[1]]), arm = arms,
               n_cells = length(ck), n_normal_cells = length(nc), n_amplicons = unname(n_amp[arms]),
               shift_log2 = unname(delta), ci_low = unname(ci[1, ]), ci_high = unname(ci[2, ]),
               p_wilcoxon = unname(p), p_adj_BH = unname(padj),
               amp_mean_shift = unname(at[, "mean"]), amp_p_ttest = unname(at[, "p"]),
               amp_frac_concordant = unname(at[, "conc"]), amp_loo_shift = unname(at[, "loo"]),
               effect_threshold = unname(thr_k), pseudocluster_null = unname(null_arm),
               cluster_scale = scale_k,
               cluster_baseline_shift = median(delta),      # diagnostic only (closure / ploidy effects)
               call = unname(cn_call), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res)
  attr(out, "null_arm") <- null_arm          # reused for the tumour-only run
  out
}

# ---- learn systematic arm bias from the called normal clusters of all samples ----
# normal_profile_raw = median log2 profile of the normal cluster(s) relative to the spike-in
# BEFORE any bias subtraction. A per-arm median across samples is robust to a few wrong calls.
learn_arm_bias <- function(cfg, samples) {
  profs <- lapply(samples, function(nm) {
    f <- file.path(cfg$outdir, nm, paste0(nm, "_RESULTS.rds"))
    if (!file.exists(f)) return(NULL)
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(r)) NULL else r$normal_profile_raw
  })
  profs <- Filter(Negate(is.null), profs)
  if (length(profs) < cfg$min_bias_samples) {
    message("Arm bias NOT learned: only ", length(profs), " samples with a normal call (need ",
            cfg$min_bias_samples, ")")
    return(invisible(NULL))
  }
  arms <- unique(unlist(lapply(profs, names)))
  M    <- sapply(profs, function(p) unname(p[arms]))
  rownames(M) <- arms
  n_ok <- rowSums(!is.na(M))
  bias <- apply(M, 1, median, na.rm = TRUE)
  bias[n_ok < cfg$min_bias_samples] <- 0
  out <- data.frame(arm = names(bias), bias = round(unname(bias), 4), n_samples = unname(n_ok))
  f <- file.path(cfg$outdir, "ARM_BIAS.tsv")
  write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Arm bias learned from ", length(profs), " samples -> ", f,
          "\n  Set cfg$arm_bias_file to this path and re-run all samples.")
  invisible(bias)
}

# =============================================================================
# 3. DRIVER
# =============================================================================
dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
sh      <- load_shared(cfg)
samples <- unique(read.delim(cfg$samples_file, header = FALSE)[[1]])

safe_run <- function(nm) {
  tryCatch(run_sample(nm, cfg, sh), error = function(e) {
    call_txt <- paste(deparse(conditionCall(e)), collapse = " ")
    msg <- paste0(conditionMessage(e), " [in: ", substr(call_txt, 1, 150), "]")
    message("[", nm, "] ERROR: ", msg)
    data.frame(sample = nm, status = paste("ERROR:", msg), cells_in_file = NA, cells_analysed = NA,
               n_clusters = NA, normal_reference = NA, n_variants = NA, stringsAsFactors = FALSE)
  })
}

res <- parallel::mclapply(samples, safe_run, mc.cores = cfg$n_cores)
overview <- bind_rows(lapply(res, function(r) if (is.data.frame(r)) r else
  data.frame(sample = NA, status = paste("ERROR:", as.character(r)))))
write.table(overview, file.path(cfg$outdir, "RUN_OVERVIEW.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(cfg$outdir, "sessionInfo.txt"))
if (isTRUE(cfg$learn_arm_bias) && is.null(sh$arm_bias)) learn_arm_bias(cfg, samples)
print(overview)
