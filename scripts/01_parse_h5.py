#!/usr/bin/env python3
"""
snCloneMC step 1: export targeted snDNA-seq .h5 files to per-sample TSV files.

Input: .h5 files in the Mission Bio format, read with the Mosaic package.

For every <sample>.dna.h5 file in the input directory, writes to the output directory:
  <sample>.dna.var_stat.tsv      variant bulk metrics (input of step 2, 02_annotate_variants.sh)
  <sample>.dna.NGT_FILTERED.tsv  filtered genotype calls (cells x variants)
  <sample>.dna.NGT.tsv           genotype calls (cells x variants)
  <sample>.dna.AF.tsv            variant allele frequencies (cells x variants)
  <sample>.dna.DP.tsv            read depth (cells x variants)
  <sample>.dna.GQ.tsv            genotype quality (cells x variants)
  <sample>.dna.CNV.tsv           read counts (cells x amplicons)
  <sample>.dna.META.tsv          sample metadata (JSON)
Errors are logged per sample in <output_dir>/error_log.txt; other samples continue.

Requires Mosaic (missionbio.mosaic), pandas and numpy.

Usage:
  python scripts/01_parse_h5.py --input-dir /path/to/h5_files --output-dir /path/to/parsed
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import missionbio.mosaic as ms

DNA_LAYERS = ["NGT_FILTERED", "NGT", "AF", "DP", "GQ"]


def to_native(metadata):
    """Convert numpy values in the metadata dict to JSON-serialisable Python types."""
    for key, value in metadata.items():
        if isinstance(value, np.ndarray):
            metadata[key] = value.tolist()
        elif isinstance(value, np.generic):
            metadata[key] = value.item()
    return metadata


def parse_sample(sample_path, out_prefix):
    sample = ms.load(sample_path, filter_variants=True, filter_cells=True)

    # Keep variants passing Mosaic's default filters
    dna_vars = sample.dna.filter_variants()
    sample.dna = sample.dna[sample.dna.barcodes(), dna_vars]
    print(f"{os.path.basename(sample_path)} shape: {sample.dna.shape}")

    # Variant bulk metrics
    stat = sample.dna.get_variant_bulk_metrics()
    stat.to_csv(f"{out_prefix}.var_stat.tsv", sep="\t", index=True)

    # Per-cell DNA layers (cells x variants)
    for layer in DNA_LAYERS:
        df = pd.DataFrame(sample.dna.layers[layer], columns=dna_vars, index=sample.dna.barcodes())
        df.to_csv(f"{out_prefix}.{layer}.tsv", sep="\t", index=True)

    # Read counts per amplicon (cells x amplicons)
    cnv = pd.DataFrame(sample.cnv.layers["read_counts"],
                       index=sample.cnv.barcodes(), columns=sample.cnv.col_attrs["id"])
    cnv.to_csv(f"{out_prefix}.CNV.tsv", sep="\t", index=True)

    # Metadata
    with open(f"{out_prefix}.META.tsv", "w") as fh:
        json.dump(to_native(sample.dna.metadata), fh, indent=4)


def main():
    parser = argparse.ArgumentParser(description="snCloneMC step 1: export snDNA-seq .h5 files to per-sample TSV files.")
    parser.add_argument("--input-dir", required=True, help="Directory containing the .h5 files")
    parser.add_argument("--output-dir", required=True, help="Directory for the exported TSV files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    sample_files = sorted(f for f in os.listdir(args.input_dir) if f.endswith(".h5"))
    if not sample_files:
        raise SystemExit(f"No .h5 files found in {args.input_dir}")

    error_log_file = os.path.join(args.output_dir, "error_log.txt")
    n_ok = 0
    with open(error_log_file, "w") as error_log:
        for sample_file in sample_files:
            # e.g. SAMPLE.dna.h5 -> <output_dir>/SAMPLE.dna
            out_prefix = os.path.join(args.output_dir, sample_file[: -len(".h5")])
            try:
                parse_sample(os.path.join(args.input_dir, sample_file), out_prefix)
                n_ok += 1
            except Exception as e:
                msg = f"Error processing {sample_file}: {e}\n"
                print(msg)
                error_log.write(msg)

    print(f"Done: {n_ok}/{len(sample_files)} samples exported to {args.output_dir}")


if __name__ == "__main__":
    main()
