# snCloneMC

**s**ingle-**n**ucleus **Clone**s from **M**utations and **C**opy number: clustering and clonal inference for **targeted single-nucleus DNA sequencing (snDNA-seq)** data based on mutations and copy-number alterations.

snCloneMC computes arm-level and focal copy-number profiles for every nucleus against a reference of normal cells, combines them with per-cell allele frequencies of selected somatic mutations, and clusters cells into normal, tumour and putative doublet populations to reconstruct the clonal architecture of each tumour. Each variant is treated as a distinct event, not collapsed to the gene level.

It was developed for, and applied to, 228 early-stage HR+/HER2– breast tumours (see [Citation](#citation)).

## Overview

```
 .h5 files ─► 1. 01_parse_h5.py ─────────► per-sample TSV files (read counts, AF, NGT, DP, GQ, variant metrics)
                                                │
                                                ▼
              2. 02_annotate_variants.sh ─► VCF ─► ANNOVAR (hg19) ─► MAF ─► OncoKB
                                                │
                                                ▼
              3. 03_snclonemc.R ─────────────► QC, copy-number profiles, clustering, normal/doublet calls,
                                              gain/loss calls, tables and figures per sample
```

Main steps of `03_snclonemc.R`, per sample:

1. **Cell QC**: total reads (absolute floor and MAD-based low/high filter), fraction of amplicons detected, VAF missingness.
2. **Amplicon QC** on the normal reference cells (detection rate and noise); failing amplicons are removed.
3. **Copy number**: raw counts are summed per chromosome arm or focal gene, converted to fractions, divided by the median fraction of the normal reference and log2-transformed.
4. **Mutation selection**: genotyping rate, allele frequency, heterozygous fraction and a germline VAF window, with rescue of OncoKB-oncogenic or known variants and removal of blacklisted artefacts.
5. **Clustering**: UMAP (Manhattan distance) and HDBSCAN on copy-number log2 ratios plus weighted VAFs.
6. **Normal cluster call** (flattest median profile), re-referencing of all cells to the called normal cells, and **putative doublet** detection (normal + tumour mixtures supported by copy number, VAF and read counts).
7. **Tumour-only re-clustering** and **gain/loss calls** per cluster (cell-level and amplicon-level tests with effect-size thresholds).
8. A systematic per-arm bias between the normal reference and the samples is **learned across samples** (`ARM_BIAS.tsv`) and subtracted in a second run.

## Requirements

| Step | Software |
|---|---|
| 1 | Python ≥ 3.8 with `pandas`, `numpy` (`requirements.txt`) and Mosaic (`missionbio.mosaic`), which reads the `.h5` files |
| 2 | bash, [ANNOVAR](https://annovar.openbioinformatics.org/) with hg19 databases `refGene`, `avsnp150`, `gnomad211_genome`, `clinvar_20221231`, `dbnsfp35a`; the [OncoKB annotator](https://github.com/oncokb/oncokb-annotator) and an OncoKB API token; optionally GNU parallel |
| 3 | R ≥ 4.1 with `dplyr`, `tidyr`, `ggplot2`, `ggpubr`, `uwot`, `dbscan`, `pheatmap`, `rstatix` |

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "ggpubr", "uwot", "dbscan", "pheatmap", "rstatix"))
```

Step 3 runs samples in parallel with `parallel::mclapply` (Linux and macOS; use `n_cores = 1` on Windows).

## Usage

### 1. Export the `.h5` files

```bash
python scripts/01_parse_h5.py --input-dir /path/to/h5_files --output-dir /path/to/data
```

For each `<sample>.dna.h5`, this writes `<sample>.dna.CNV.tsv`, `<sample>.dna.AF.tsv`, `<sample>.dna.var_stat.tsv` and the NGT, NGT_FILTERED, DP, GQ and META files. Errors are logged in `error_log.txt` and do not stop the other samples.

### 2. Annotate the variants

```bash
export ONCOKB_TOKEN=your-oncokb-token          # never write the token in a file of the repository
scripts/02_annotate_variants.sh -i /path/to/data \
    -a /path/to/annovar/table_annovar.pl -d /path/to/annovar/humandb -j 8 \
    -k /path/to/oncokb-annotator/MafAnnotator.py -t BRCA
```

The four sub-steps (`vcf`, `annovar`, `maf`, `oncokb`) can be run separately with `-s`, for example ANNOVAR on a computing cluster and the others locally. All files stay in the same directory, and the final file of each sample is:

```
<sample>.dna.var_stat.tsv.vcf.hg19_multianno.txt.maf.oncoKB.txt
```

The `maf` sub-step selects ANNOVAR columns by position and names them with `resources/maf_header.txt`; the column positions assume the ANNOVAR protocol above. Use `-t` to set the OncoTree tumour type (default `BRCA`). Run `scripts/02_annotate_variants.sh -h` for all options.

### 3. Clustering and clonal inference

```bash
cp config/config_example.R config/my_config.R   # edit the paths
Rscript scripts/03_snclonemc.R config/my_config.R
```

**Two-pass run (recommended).** The first run writes `<outdir>/ARM_BIAS.tsv`, the systematic per-arm offset between the normal reference and the samples, learned from the normal clusters of at least `min_bias_samples` samples. Set `arm_bias_file` to this file in your config and run again.

All parameters (QC thresholds, clustering, normal and doublet calls, gain/loss significance) are documented in section 0 of `scripts/03_snclonemc.R` and can be overridden in the config file. Use `normal_override` to set the normal cluster manually for a sample.

## Input files of step 3

**Per sample** (in `data_dir`):

| File | Content |
|---|---|
| `<sample>.dna.CNV.tsv` | Read counts, cells × amplicons (step 1) |
| `<sample>.dna.AF.tsv` | Allele frequencies, cells × variants (step 1) |
| `<sample>.dna.var_stat.tsv.vcf.hg19_multianno.txt.maf.oncoKB.txt` | Annotated variants (step 2). Required columns: `variant` (variant ID, as in the AF file), `AF`, `AF_popmax`, `Genotyped_cells`, `WT_cells`, `HET_cells`, `HOM_cells`, `VAF_by_cell_count`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `ONCOGENIC`, `ONCOKB_HUGO_SYMBOL`, `ONCOKB_PROTEIN_CHANGE` |

**Shared** (paths set in the config):

| Config entry | Content |
|---|---|
| `samples_file` | One sample name per line, no header (e.g. `SAMPLE_A` for `SAMPLE_A.dna.CNV.tsv`) |
| `spikein_file` | Read counts (cells × amplicons, cell IDs in the first column) of trusted normal cells; the first `n_spikein` rows form the normal reference |
| `arm_annot_file` | Amplicon annotation: amplicon ID in column 1 and chromosome arm (e.g. `8q`) or focal gene (e.g. `ERBB2_CN`) in column 8 |
| `idt_file` | Known variants annotated with OncoKB (columns `CHROMOSOME`, `START_POSITION`, `END_POSITION`, `REFERENCE_ALLELE`, `TUMOR_SEQ_ALLELE2`, `ONCOGENIC`, `ONCOKB_HUGO_SYMBOL`, `ONCOKB_PROTEIN_CHANGE`) |
| `blacklist_file` | Recurrent artefactual variants (columns `id` = `GENE_ProteinChange`, `fdr`, `log`); variants with `fdr < 0.1` and `log > 0` are removed |

The arms and focal genes analysed, and their order, are set by `arm_order`.

## Output

In `<outdir>/<sample>/`:

| File | Content |
|---|---|
| `_CELL_QC.tsv`, `_AMPLICON_QC.tsv` | Cell and amplicon QC |
| `_CLUSTER_STATS.tsv`, `_DOUBLET_STATS.tsv` | Normal cluster call and doublet evaluation |
| `_CN_CALLS.tsv`, `_CN_CALLS_TUMOR_ONLY.tsv` | Gain/loss calls per cluster and arm |
| `_SUMMARY_TABLE.tsv`, `_SUMMARY_TABLE_TUMOR_ONLY.tsv` | Per-cluster summary: cell counts, reads, median log2 ratios, gain/loss calls, VAFs of oncogenic variants |
| `_PER_CELL_TABLE.tsv` | One row per cell: cluster, cell type, log2 ratios, VAFs |
| `_READS_PAIRWISE_WILCOXON.tsv` | Read-count comparisons between clusters |
| `_RESULTS.rds` | All results and the configuration used, with `sessionInfo()` |
| `_UMAP_CLUSTERED.pdf`, `_UMAP_TUMOR_ONLY.pdf`, `_PHEATMAP.pdf`, `_PHEATMAP_TUMOR_ONLY.pdf`, `_VAF_ON_UMAP.pdf`, `_VAF_BOXPLOT.pdf`, `_CLUSTER_PERCENT.pdf`, `_READS_PER_CLUSTER.pdf` | Figures |

In `<outdir>/`: `RUN_OVERVIEW.tsv` (status of every sample), `ARM_BIAS.tsv` (first run) and `sessionInfo.txt`.

## Citation

If you use snCloneMC, please cite:

> Yurchenko AA, Verret B, Cerbone L, Tran-Dien A, Jelin R, Droin N, Scott V, Leite B, Francillette M, Aglave M, Scoazec J-Y, Ayadi C, Machavoine C, Camara-Clayette V, Lerebours F, Lemonnier J, André F, Tran DTN, Nikolaev S. Single-nucleus genomic profiling of HR+/HER2– breast tumors revealed limited intratumoral genetic heterogeneity in early-stage disease. **In preparation.**


## Licence

MIT, see [LICENSE](LICENSE).

## Contact

Andrey A. Yurchenko, Inserm U1361 Cancer Data Science, Gustave Roussy, Villejuif, France. Questions and bug reports: please open a GitHub issue.
