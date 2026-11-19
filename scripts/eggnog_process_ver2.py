#!/usr/bin/env python3
"""Run eggNOG-mapper once on all de-duplicated .faa files under a base directory, then split the result into one
file per genomic annotation file (<taxon>/eggnog/<genome>_genomic_emapper_results.tsv)."""
import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

TEMP_DIR_NAME = 'temp_annotations'
CONSOLIDATED_FAA = 'consolidated_all.faa'
OUTPUT_PREFIX = 'consolidated_results'
EMAPPER_SUFFIX = '.emapper.annotations'
EGGNOG_DATA_DIR = ('data', 'eggnog_data')
FAA_PATTERN = '*/faa/*.faa'
GENOMIC_ANNOTATION_PATTERN = '*/genomic_annotations/*.tsv'
EMAPPER_HEADER_LINES = 4
WRITE_WAIT_SECONDS = 5


def consolidate_faa(faa_files, output):
    """Writes every protein once (first occurrence wins) and returns the number of unique proteins."""
    seen = set()
    with open(output, 'w', encoding='utf-8') as out:
        for faa in faa_files:
            record, protein_id = [], None
            with open(faa, encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith('>'):
                        if protein_id and protein_id not in seen:
                            out.write(''.join(record))
                            seen.add(protein_id)
                        record, protein_id = [line], line.split()[0][1:]
                    else:
                        record.append(line)
            if protein_id and protein_id not in seen:
                out.write(''.join(record))
                seen.add(protein_id)
    return len(seen)


def run_emapper(faa, output_prefix, data_dir, cpu, tax_scope):
    command = ['emapper.py', '-i', str(faa), '-o', str(output_prefix), '--data_dir', str(data_dir),
               '--cpu', str(cpu), '--override', '--dbmem']
    if tax_scope:
        command += ['--tax_scope', tax_scope]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError:
        sys.exit("[ERROR] eggNOG-mapper failed.")
    finally:
        time.sleep(WRITE_WAIT_SECONDS)


def load_emapper(path):
    try:
        df = pd.read_csv(path, sep='\t', skiprows=EMAPPER_HEADER_LINES)
    except FileNotFoundError:
        sys.exit(f"Error: emapper file not found: {path}")
    except Exception as e:
        sys.exit(f"Error: failed to load {path}: {e}")
    df = df.rename(columns={'#query': 'Protein_ID'}).drop_duplicates(subset=['Protein_ID'])
    return df.set_index('Protein_ID')


def extract_for_annotation_file(emapper, annotation_path, output_path):
    print(f"-> {annotation_path}")
    try:
        annotation = pd.read_csv(annotation_path, sep='\t', header=0)
    except Exception as e:
        print(f"  !! Failed to read {annotation_path}: {e}")
        return
    if 'Protein_ID' not in annotation.columns:
        print("  >> Skipped: no 'Protein_ID' column.")
        return
    matched = emapper.loc[emapper.index.isin(annotation['Protein_ID'].drop_duplicates())]
    if matched.empty:
        print(f"  >> No matches for {annotation_path.name}.")
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    matched.to_csv(output_path, sep='\t')
    print(f"  >> {len(matched)} matches saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-dir', required=True, help="Directory containing the taxon directories.")
    parser.add_argument('--cpu', type=int, default=10, help="CPUs for eggNOG-mapper.")
    parser.add_argument('--tax_scope', default='Bacteria', help="Taxonomic scope for eggNOG-mapper.")
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    if not base_dir.is_dir():
        sys.exit(f"Error: base directory not found: {base_dir}")

    temp_dir = base_dir / TEMP_DIR_NAME
    temp_dir.mkdir(exist_ok=True)
    consolidated = temp_dir / CONSOLIDATED_FAA
    output_prefix = temp_dir / OUTPUT_PREFIX
    emapper_file = output_prefix.with_suffix(EMAPPER_SUFFIX)

    if emapper_file.exists():
        print(f"{emapper_file} exists; eggNOG-mapper not run.")
    else:
        faa_files = list(base_dir.rglob(FAA_PATTERN))
        if not faa_files:
            print("No .faa files found.", file=sys.stderr)
            sys.exit(0)
        print(f"Consolidating {len(faa_files)} files -> {consolidated}")
        print(f"{consolidate_faa(faa_files, consolidated)} unique proteins.")
        data_dir = base_dir.parent.joinpath(*EGGNOG_DATA_DIR)
        if not data_dir.is_dir():
            sys.exit(f"Error: eggNOG data directory not found: {data_dir}")
        run_emapper(consolidated, output_prefix, data_dir, args.cpu, args.tax_scope)

    emapper = load_emapper(emapper_file)
    print(f"Loaded {len(emapper)} unique emapper entries.")
    annotation_files = list(base_dir.glob(GENOMIC_ANNOTATION_PATTERN))
    if not annotation_files:
        print(f"No annotation files found for {base_dir / GENOMIC_ANNOTATION_PATTERN}")
        sys.exit(0)
    for path in annotation_files:
        output = path.parent.parent / 'eggnog' / f"{path.stem}_emapper_results.tsv"
        extract_for_annotation_file(emapper, path, output)
    print(f"Done. {consolidated.name} was kept.")


if __name__ == '__main__':
    main()
