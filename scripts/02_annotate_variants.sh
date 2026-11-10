#!/usr/bin/env bash
# =============================================================================
# snCloneMC step 2: annotate the variants exported in step 1
#
#   <sample>.dna.var_stat.tsv
#     -> .vcf                                  (vcf)
#     -> .vcf.hg19_multianno.txt               (annovar: ANNOVAR, hg19)
#     -> .vcf.hg19_multianno.txt.maf           (maf: MAF-like table)
#     -> .vcf.hg19_multianno.txt.maf.oncoKB.txt (oncokb: OncoKB annotator)
#
# All files are read from and written to the same directory (-i), which is
# the `data_dir` used in step 3. Steps can be run separately, e.g. ANNOVAR on
# a cluster and the other steps locally.
#
# The OncoKB API token is read from the environment variable ONCOKB_TOKEN:
#   export ONCOKB_TOKEN=your-token
# Never write the token in this script or commit it to the repository.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- defaults ---------------------------------------------------------------
STEP="all"
JOBS=4
MAF_HEADER="${SCRIPT_DIR}/../resources/maf_header.txt"
TUMOR_TYPE="BRCA"          # OncoTree code passed to the OncoKB annotator
REF_GENOME="GRCh37"
ANNOVAR=""                 # path to table_annovar.pl
HUMANDB=""                 # ANNOVAR humandb directory
ONCOKB_ANNOTATOR=""        # path to oncokb-annotator MafAnnotator.py
IN_DIR=""

usage() {
  cat <<EOF
Usage: $(basename "$0") -i DIR [options]

  -i DIR    directory with the *.var_stat.tsv files from step 1 (required)
  -s STEP   vcf | annovar | maf | oncokb | all          [default: all]
  -a FILE   ANNOVAR table_annovar.pl                    (annovar step)
  -d DIR    ANNOVAR humandb directory (hg19)            (annovar step)
  -j N      parallel ANNOVAR jobs                       [default: ${JOBS}]
  -H FILE   MAF header file                             [default: resources/maf_header.txt]
  -k FILE   oncokb-annotator MafAnnotator.py            (oncokb step)
  -t CODE   OncoTree tumour type for OncoKB             [default: ${TUMOR_TYPE}]
  -r REF    reference genome for OncoKB                 [default: ${REF_GENOME}]
  -h        show this help

The OncoKB token must be set in the environment variable ONCOKB_TOKEN.
EOF
}

while getopts "i:s:a:d:j:H:k:t:r:h" opt; do
  case "$opt" in
    i) IN_DIR="$OPTARG" ;;
    s) STEP="$OPTARG" ;;
    a) ANNOVAR="$OPTARG" ;;
    d) HUMANDB="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    H) MAF_HEADER="$OPTARG" ;;
    k) ONCOKB_ANNOTATOR="$OPTARG" ;;
    t) TUMOR_TYPE="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

[[ -n "$IN_DIR" && -d "$IN_DIR" ]] || { echo "ERROR: -i DIR is required and must exist" >&2; usage; exit 1; }
MAF_HEADER="$(cd "$(dirname "$MAF_HEADER")" && pwd)/$(basename "$MAF_HEADER")"   # absolute, before cd
cd "$IN_DIR"
shopt -s nullglob

# ---- 1. var_stat.tsv -> VCF ---------------------------------------------------
# The variant ID (chr:pos:ref/alt) gives CHROM, POS, REF and ALT; the variant
# ID and the per-variant metrics are kept in the INFO field.
step_vcf() {
  local files=( *var_stat.tsv )
  (( ${#files[@]} )) || { echo "No *var_stat.tsv files in $IN_DIR" >&2; exit 1; }
  for f in "${files[@]}"; do
    sed 's/ /_/g' "$f" \
      | sed 's/[\/,:]/\t/g' \
      | grep -v "Genotyped_cells" \
      | awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,".",$3,$4,".","PASS",$1":"$2":"$3"/"$4";"$5";"$6";"$7";"$8";"$9";"$10";"$11";"$12";"$13";"$14}' \
      > "${f}.vcf"
  done
  echo "[vcf] ${#files[@]} files converted"
}

# ---- 2. ANNOVAR (hg19) -------------------------------------------------------
annovar_one() {
  "$ANNOVAR" "$1" "$HUMANDB" --vcfinput -buildver hg19 -out "$(basename "$1")" -remove \
    -protocol refGene,avsnp150,gnomad211_genome,clinvar_20221231,dbnsfp35a \
    -operation g,f,f,f,f -nastring . --thread 1 --maxgenethread 1
}

step_annovar() {
  [[ -x "$ANNOVAR" ]] || { echo "ERROR: -a must point to table_annovar.pl" >&2; exit 1; }
  [[ -d "$HUMANDB" ]] || { echo "ERROR: -d must point to the ANNOVAR humandb directory" >&2; exit 1; }
  local files=( *.vcf )
  (( ${#files[@]} )) || { echo "No *.vcf files in $IN_DIR" >&2; exit 1; }
  if command -v parallel >/dev/null 2>&1; then
    export ANNOVAR HUMANDB
    export -f annovar_one
    parallel --gnu --jobs "$JOBS" annovar_one {} ::: "${files[@]}"
  else
    echo "[annovar] GNU parallel not found: running sequentially"
    for f in "${files[@]}"; do annovar_one "$f"; done
  fi
  echo "[annovar] ${#files[@]} files annotated"
}

# ---- 3. ANNOVAR table -> MAF-like table ----------------------------------------
# Keeps the selected ANNOVAR columns and splits the INFO field (column 114:
# variant ID and per-variant metrics) into separate columns. Column names come
# from the MAF header file; the column numbers assume the ANNOVAR protocol of step 2.
step_maf() {
  [[ -f "$MAF_HEADER" ]] || { echo "ERROR: MAF header file not found: $MAF_HEADER" >&2; exit 1; }
  local files=( *_multianno.txt )
  (( ${#files[@]} )) || { echo "No *_multianno.txt files in $IN_DIR" >&2; exit 1; }
  for f in "${files[@]}"; do
    paste <(cut -f1-7,9-13,30-33,75,101 "$f") <(cut -f114 "$f" | sed 's/;/\t/g') \
      | grep -v "CADD_phred" \
      | cat "$MAF_HEADER" - > "${f}.maf"
  done
  echo "[maf] ${#files[@]} files converted"
}

# ---- 4. OncoKB annotation -------------------------------------------------------
step_oncokb() {
  [[ -f "$ONCOKB_ANNOTATOR" ]] || { echo "ERROR: -k must point to MafAnnotator.py" >&2; exit 1; }
  [[ -n "${ONCOKB_TOKEN:-}" ]] || { echo "ERROR: set the OncoKB token: export ONCOKB_TOKEN=..." >&2; exit 1; }
  local files=( *.maf )
  (( ${#files[@]} )) || { echo "No *.maf files in $IN_DIR" >&2; exit 1; }
  for f in "${files[@]}"; do
    python3 "$ONCOKB_ANNOTATOR" -q Genomic_Change -i "$f" -o "${f}.oncoKB.txt" \
      -t "$TUMOR_TYPE" -b "$ONCOKB_TOKEN" -r "$REF_GENOME"
  done
  echo "[oncokb] ${#files[@]} files annotated"
}

case "$STEP" in
  vcf)     step_vcf ;;
  annovar) step_annovar ;;
  maf)     step_maf ;;
  oncokb)  step_oncokb ;;
  all)     step_vcf; step_annovar; step_maf; step_oncokb ;;
  *)       echo "ERROR: unknown step '$STEP'" >&2; usage; exit 1 ;;
esac
