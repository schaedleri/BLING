#!/usr/bin/env python3
"""
External (community-wide) BLASTP screen.

Queries every protein that passed the internal (within-genus) screen against
the unified BLAST+ database (all reference genomes merged via
blastdb_aliastool), restricted to taxids outside the target taxonomic group.
Proteins with no hit are putatively taxon-restricted.
"""
import sys
import os
import re
import argparse
from pathlib import Path
import shlex
from multiprocessing import Pool, cpu_count


def read_taxids_from_tsv(tsv_path):
    """Reads taxids from a TSV file's 4th column."""
    taxids = set()
    if not tsv_path.is_file():
        return taxids
    try:
        with open(tsv_path, 'r', encoding='utf-8') as f:
            try:
                next(f)
            except StopIteration:
                return taxids
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= 4 and cols[3]:
                    taxids.add(cols[3].strip())
    except IOError as e:
        print(f"Warning: Could not read {tsv_path}: {e}", file=sys.stderr)
    return taxids


def load_ids_from_txt(file_path):
    """Loads IDs from a simple text file, one ID per line."""
    ids = set()
    if not file_path.is_file():
        return ids
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            ids.update(line.strip() for line in f if line.strip())
    except IOError as e:
        print(f"Warning: Could not read {file_path}: {e}", file=sys.stderr)
    return ids


def parse_clstr_file(clstr_path):
    """Parses a CD-HIT .clstr file into a list of {star_id, ids} clusters."""
    clusters = []
    if not clstr_path.is_file():
        return clusters
    try:
        with open(clstr_path, 'r', encoding='utf-8') as f:
            current_members, star_id = [], None
            for line in f:
                if line.startswith('>'):
                    if current_members:
                        clusters.append({'star_id': star_id, 'ids': current_members})
                    current_members, star_id = [], None
                else:
                    match = re.search(r'>\s*(WP_\S+)\.\.\.', line)
                    if match:
                        member_id = match.group(1)
                        current_members.append(member_id)
                        if line.strip().endswith('*'):
                            star_id = member_id
            if current_members:
                clusters.append({'star_id': star_id, 'ids': current_members})
    except IOError as e:
        print(f"Warning: Could not read cluster file {clstr_path}: {e}", file=sys.stderr)
    return clusters


def parse_hit_query_pairs(pairs_path):
    """Reads hit_query_pairs.tsv (from the internal BLASTP step): representative -> members."""
    pairs = {}
    if not pairs_path.is_file():
        return pairs
    try:
        with open(pairs_path, 'r', encoding='utf-8') as f:
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= 2:
                    pairs.setdefault(cols[0], []).append(cols[1])
    except IOError as e:
        print(f"Warning: Could not read {pairs_path}: {e}", file=sys.stderr)
    return pairs


def blast_task_worker(args_tuple):
    """Runs a single BLASTP search against the unified reference database."""
    fasta_path, unified_db_path, taxid_str, genus_out_dir, species_out_dir, mode, cdhit_ids = args_tuple

    protein_id = fasta_path.stem
    if protein_id not in cdhit_ids:
        return

    output_dir = genus_out_dir if mode == 'genus' else species_out_dir
    blast_out_file = output_dir / f"{protein_id}_Unified.tsv"

    q_query = shlex.quote(str(fasta_path))
    q_db = shlex.quote(str(unified_db_path))
    q_out = shlex.quote(str(blast_out_file))
    outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq'
    q_outfmt = shlex.quote(outfmt)
    taxids_arg = f"-taxids '{taxid_str}'" if taxid_str else ""

    # -num_threads 1 avoids CPU thrashing across many parallel worker processes;
    # -seg yes masks low-complexity regions to avoid spurious hits.
    cmd_str = (
        f"blastp -query {q_query} "
        f"-db {q_db} "
        f"-task blastp-fast "
        f"-num_threads 1 "
        f"-max_target_seqs 10 "
        f"-evalue 1e-10 "
        f"-seg yes "
        f"-outfmt {q_outfmt} "
        f"-out {q_out} "
        f"{taxids_arg}"
    )

    try:
        exit_code = os.system(cmd_str)
        if exit_code != 0:
            print(f"Warning: blastp failed for {protein_id} with exit code {exit_code}.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred for {protein_id}: {e}", file=sys.stderr)

    return f"Finished BLAST for {protein_id}"


def main():
    parser = argparse.ArgumentParser(description="Run the external (community-wide) BLASTP screen against the unified database.")
    parser.add_argument('--base-dir', required=True, help="Project root directory.")
    parser.add_argument('--target-dir', required=True, help="Sample directory name (e.g., 'level-7').")
    parser.add_argument('taxon_dir', help="The specific Genus_* or Species_* directory to process.")
    args = parser.parse_args()

    project_root = Path(args.base_dir).resolve()
    sample_dir = project_root / args.target_dir
    taxon_path = sample_dir / args.taxon_dir
    mode = 'genus' if args.taxon_dir.startswith('Genus_') else 'species'

    print(f"--- Running external BLASTP for: {args.taxon_dir} (mode: {mode}) ---")

    print("Step 1: Calculating taxid exclusion list...")
    all_taxids = set()
    for pattern in ["Species_*", "Genus_*", "Family_*", "Order_*"]:
        for tsv_file in sample_dir.glob(f"{pattern}/*.tsv"):
            all_taxids.update(read_taxids_from_tsv(tsv_file))

    selected_taxids = read_taxids_from_tsv(taxon_path / f"{args.taxon_dir}.tsv")
    taxid_diff = all_taxids - selected_taxids
    taxid_str = ",".join(sorted(taxid_diff)) if taxid_diff else ""
    print(f"Found {len(all_taxids)} total taxids, {len(selected_taxids)} for target, {len(taxid_diff)} for exclusion.")
    if not taxid_str:
        print("Warning: taxid exclusion list is empty.", file=sys.stderr)

    cdhit_dir = taxon_path / "CDhit"
    fasta_root = taxon_path / "fasta"
    output_dir = taxon_path / "BLAST_specify"
    output_dir.mkdir(exist_ok=True)

    cdhit_ids = load_ids_from_txt(cdhit_dir / f"{mode}.txt")
    if not cdhit_ids:
        print(f"Warning: No IDs found in {cdhit_dir / f'{mode}.txt'}. Nothing to do.", file=sys.stderr)
        sys.exit(0)

    print(f"Step 2: Preparing to run BLAST for {len(cdhit_ids)} protein IDs...")
    unified_db_path = project_root / "DB" / "bacteria_strain_taxid_DB" / "Unified_Bacteria_DB"
    if not unified_db_path.with_suffix(".pal").is_file():
        print(f"[ERROR] Unified database not found at: {unified_db_path}.pal", file=sys.stderr)
        print("Build it with blastdb_aliastool before running this script.", file=sys.stderr)
        sys.exit(1)

    all_fasta_files = sorted(fasta_root.rglob("WP_*.fasta"))
    print(f"Found {len(all_fasta_files)} total FASTA files to check.")

    tasks = []
    queued_ids = set()
    for fasta_file in all_fasta_files:
        protein_id = fasta_file.stem
        if protein_id not in cdhit_ids or protein_id in queued_ids:
            continue
        queued_ids.add(protein_id)

        species_dir_name = fasta_file.parent.name
        base_out = output_dir / species_dir_name
        genus_out = base_out / "genus"
        species_out = base_out / "species"
        genus_out.mkdir(parents=True, exist_ok=True)
        species_out.mkdir(parents=True, exist_ok=True)

        tasks.append((fasta_file, unified_db_path, taxid_str, genus_out, species_out, mode, cdhit_ids))

    print(f"Reduced {len(all_fasta_files)} FASTA files to {len(tasks)} unique query tasks.")

    num_cpus = min(cpu_count(), 20)
    print(f"Running {len(tasks)} BLAST tasks in parallel using {num_cpus} CPUs...")
    with Pool(processes=num_cpus) as pool:
        for i, result in enumerate(pool.imap_unordered(blast_task_worker, tasks), 1):
            if result:
                sys.stdout.write(f"\rProgress: {i}/{len(tasks)} tasks completed.")
                sys.stdout.flush()
    print("\nAll BLAST tasks finished.")

    print("Step 3: Aggregating BLAST results...")
    result_ids = {'genus': set(), 'species': set()}
    for species_subdir in output_dir.glob("*/"):
        if not species_subdir.is_dir():
            continue
        for id_subdir_mode in ['genus', 'species']:
            for f in (species_subdir / id_subdir_mode).glob("*_Unified.tsv"):
                # An empty result file means no hit was found -> taxon-restricted.
                if f.stat().st_size == 0:
                    result_ids[id_subdir_mode].add(f.stem.replace("_Unified", ""))

    print("Step 4: Expanding results with cluster and pair information...")
    cdhit_clusters = []
    clstr_files = list(fasta_root.glob("*.clstr"))
    if clstr_files:
        cdhit_clusters = parse_clstr_file(clstr_files[0])

    blast_pairs = parse_hit_query_pairs(cdhit_dir / 'hit_query_pairs.tsv')
    print(f"  -> Loaded {len(blast_pairs)} BLAST pairs from the internal screen.")

    for m in ['genus', 'species']:
        if not result_ids[m]:
            continue
        current_ids = result_ids[m].copy()
        ids_to_add = set()
        for c in cdhit_clusters:
            if c.get('star_id') in current_ids:
                ids_to_add.update(c.get('ids', []))
        for rep_id in current_ids:
            if rep_id in blast_pairs:
                ids_to_add.update(blast_pairs[rep_id])
        result_ids[m].update(ids_to_add)

    print("Step 5: Writing final result files...")
    for m, ids in result_ids.items():
        result_path = output_dir / f"result_{m}.txt"
        with open(result_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(sorted(ids)) + '\n')
        print(f"  -> Created {result_path} with {len(ids)} IDs.")

    print(f"\nSuccessfully finished external BLAST for {args.taxon_dir}.")


if __name__ == "__main__":
    main()
