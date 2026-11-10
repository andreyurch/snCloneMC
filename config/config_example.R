# =============================================================================
# Example configuration for scripts/03_snclonemc.R
#
# Copy this file (e.g. to config/my_config.R), edit the paths and run:
#   Rscript scripts/03_snclonemc.R config/my_config.R
#
# Only `cfg_user` is read. Any parameter defined in section 0 of
# 03_snclonemc.R can be overridden here (QC thresholds, clustering, ...).
# =============================================================================
base <- "/path/to/project"

cfg_user <- list(
  # ---- required paths --------------------------------------------------------
  samples_file   = file.path(base, "data/SAMPLES_LIST.txt"),
  data_dir       = file.path(base, "data/parsed"),
  spikein_file   = file.path(base, "reference/normal_cells_read_counts.tsv"),
  arm_annot_file = file.path(base, "reference/panel_amplicon_arm_annotation.tsv"),
  idt_file       = file.path(base, "reference/known_variants.oncoKB.txt"),
  blacklist_file = file.path(base, "reference/variant_blacklist.txt"),
  outdir         = file.path(base, "results/run_01"),

  # ---- run settings ----------------------------------------------------------
  n_cores = 4,
  seed    = 12345

  # ---- second pass: subtract the learned spike-in/sample arm bias ------------
  # After a first run, uncomment the line below (and add a comma above):
  # arm_bias_file = file.path(base, "results/run_01/ARM_BIAS.tsv")
)
