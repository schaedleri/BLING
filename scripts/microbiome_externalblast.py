#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse
from pathlib import Path
import subprocess
import shlex
from multiprocessing import Pool, cpu_count

# --- Helper Functions ---

def read_taxids_from_tsv(tsv_path):
    """Reads taxids from a TSV file's 4th column, safely handling empty files."""
    taxids = set()
    if not tsv_path.is_file():
        return taxids
    try:
        with open(tsv_path, 'r', encoding='utf-8') as f:
            try:
                next(f)  # Try to skip the header
            except StopIteration:
                return taxids # If file is empty, return the empty set immediately
            
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
    """Parses a CD-HIT .clstr file."""
    clusters = []
    if not clstr_path.is_file():
        return clusters
    try:
        with open(clstr_path, 'r', encoding='utf-8') as f:
            current_members, star_id = [], None
            for line in f:
                if line.startswith('>'):
                    if current_members: clusters.append({'star_id': star_id, 'ids': current_members})
                    current_members, star_id = [], None
                else:
                    if match := re.search(r'>\s*(WP_\S+)\.\.\.', line):
                        member_id = match.group(1)
                        current_members.append(member_id)
                        if line.strip().endswith('*'): star_id = member_id
            if current_members: clusters.append({'star_id': star_id, 'ids': current_members})
    except IOError as e:
        print(f"Warning: Could not read cluster file {clstr_path}: {e}", file=sys.stderr)
    return clusters

def blast_task_worker(args_tuple):
    """
    Worker function for running a set of BLAST searches for a single FASTA file using os.system.
    This ensures the command string is passed directly to the shell.
    """
    fasta_path, dbs_config, taxid_str, genus_out_dir, species_out_dir, mode, cdhit_ids = args_tuple
    
    protein_id = fasta_path.stem
    if protein_id not in cdhit_ids:
        return # Skip if this ID is not in the target list

    output_dir = genus_out_dir if mode == 'genus' else species_out_dir
    
    for db_info in dbs_config:
        db_path = db_info["path"]
        db_number = db_info["number"]
        blast_out_file = output_dir / f"{protein_id}_{db_number}.tsv"
        
        # Step 1: Build the command string for the shell.
        # Use shlex.quote() to make sure file paths and arguments with spaces are safe.
        q_query = shlex.quote(str(fasta_path))
        q_db = shlex.quote(str(db_path))
        q_out = shlex.quote(str(blast_out_file))
        
        # The -outfmt argument string contains spaces, so it must be quoted as well.
        outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq'
        q_outfmt = shlex.quote(outfmt)

        # Step 2: Build the -taxids argument part.
        # This creates the '-taxids '...' string part only if taxid_str exists.
        taxids_arg = f"-taxids '{taxid_str}'" if taxid_str else ""
        
        # Step 3: Assemble all parts into the final command string.
        # The f-string builds the exact command that will be passed to the shell.
        final_cmd_str = (
            f"blastp -query {q_query} "
            f"-db {q_db} "
            f"-max_target_seqs 10 "
            f"-evalue 1e-10 "
            f"-outfmt {q_outfmt} "
            f"-out {q_out} "
            f"{taxids_arg}"
        )

        
        try:
            # Step 4: Execute the command using os.system.
            # os.system passes the string directly to the system's shell.
            exit_code = os.system(final_cmd_str)
            
            # Check the exit code to see if the command was successful.
            if exit_code != 0:
                print(f"Warning: blastp failed for {protein_id} against {db_number} with exit code {exit_code}.", file=sys.stderr)
                
        except Exception as e:
            print(f"An unexpected error occurred while running os.system for {protein_id}: {e}", file=sys.stderr)
            
    return f"Finished BLAST for {protein_id}"

# --- Main Execution Block ---
def main():
    parser = argparse.ArgumentParser(description="Run external BLAST step.")
    parser.add_argument('--base-dir', required=True, help="Project root directory.")
    parser.add_argument('--target-dir', required=True, help="Sample directory name (e.g., 'depression').")
    parser.add_argument('taxon_dir', help="The specific Genus_* or Species_* directory to process.")
    args = parser.parse_args()

    project_root = Path(args.base_dir).resolve()
    sample_dir = project_root / args.target_dir
    taxon_path = sample_dir / args.taxon_dir
    mode = 'genus' if args.taxon_dir.startswith('Genus_') else 'species'
    
    print(f"--- Running External BLAST for: {args.taxon_dir} (Mode: {mode}) ---")

    print("Step 1: Calculating TaxID exclusion list...")
    all_taxids = set()
    for tsv_file in sample_dir.glob("Genus_*/*.tsv"):
        all_taxids.update(read_taxids_from_tsv(tsv_file))
    for tsv_file in sample_dir.glob("Species_*/*.tsv"):
        all_taxids.update(read_taxids_from_tsv(tsv_file))

    selected_taxon_tsv = taxon_path / f"{args.taxon_dir}.tsv"
    selected_taxids = read_taxids_from_tsv(selected_taxon_tsv)
    
    taxid_diff = all_taxids - selected_taxids
    taxid_str = ",".join(sorted(list(taxid_diff))) if taxid_diff else ""
    print(f"Found {len(all_taxids)} total taxids in sample, {len(selected_taxids)} for target, {len(taxid_diff)} for exclusion.")
    if not taxid_str:
        print("Warning: TaxID exclusion list is empty. BLAST may not run as expected.", file=sys.stderr)

    cdhit_dir = taxon_path / "CDhit"
    fasta_root = taxon_path / "fasta"
    output_dir = taxon_path / "BLAST_specify"
    output_dir.mkdir(exist_ok=True)
    
    cdhit_ids = load_ids_from_txt(cdhit_dir / f"{mode}.txt")
    if not cdhit_ids:
        print(f"Warning: No IDs found in {cdhit_dir / f'{mode}.txt'}. Nothing to do.", file=sys.stderr)
        sys.exit(0)
    
    print(f"Step 2: Preparing to run BLAST for {len(cdhit_ids)} protein IDs...")
    
    dbs_config = []
    for i in range(1, 11):
        db_num = 0 if i == 10 else i
        db_name = f"BSDB{db_num}"
        dbs_config.append({
            "path": project_root / "DB" / "bacteria_strain_taxid_DB" / f"DB{db_num}" / db_name,
            "number": f"DB{i}"
        })

    tasks = []
    all_fasta_files = list(fasta_root.rglob("WP_*.fasta"))
    print(f"Found {len(all_fasta_files)} total FASTA files to check.")

    for fasta_file in all_fasta_files:
        species_dir_name = fasta_file.parent.name
        base_out = output_dir / species_dir_name
        genus_out = base_out / "genus"
        species_out = base_out / "species"
        genus_out.mkdir(parents=True, exist_ok=True)
        species_out.mkdir(parents=True, exist_ok=True)
        
        tasks.append((fasta_file, dbs_config, taxid_str, genus_out, species_out, mode, cdhit_ids))

    num_cpus = min(cpu_count(), 20)
    print(f"Running {len(tasks)} BLAST tasks in parallel using {num_cpus} CPUs...")
    with Pool(processes=num_cpus) as pool:
        for i, result in enumerate(pool.imap_unordered(blast_task_worker, tasks), 1):
            if result:
                sys.stdout.write(f"\rProgress: {i}/{len(tasks)} tasks completed. Last: {result}")
                sys.stdout.flush()
    print("\nAll BLAST tasks finished.")

    print("Step 3: Aggregating BLAST results...")
    result_ids = {'genus': set(), 'species': set()}
    
    for species_subdir in output_dir.glob("*/"):
        if not species_subdir.is_dir(): continue
        for id_subdir_mode in ['genus', 'species']:
            id_list = {f.stem.replace("_DB1", "") for f in (species_subdir / id_subdir_mode).glob("*_DB1.tsv")}
            for protein_id in id_list:
                if all(not (f.exists() and f.stat().st_size > 0) for i in range(1, 11) if (f := species_subdir / id_subdir_mode / f"{protein_id}_DB{i}.tsv")):
                    result_ids[id_subdir_mode].add(protein_id)
    
    print("Step 4: Expanding results with cluster information...")
    if clstr_files := list(fasta_root.glob("*.clstr")):
        clusters = parse_clstr_file(clstr_files[0])
        for m in ['genus', 'species']:
            if result_ids[m]:
                ids_to_add = {pid for c in clusters if c.get('star_id') in result_ids[m] for pid in c.get('ids', [])}
                result_ids[m].update(ids_to_add)

    print("Step 5: Writing final result files...")
    for m, ids in result_ids.items():
        result_path = output_dir / f"result_{m}.txt"
        with open(result_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(sorted(list(ids))) + '\n')
        print(f"  -> Created {result_path} with {len(ids)} IDs.")

    print(f"\nSuccessfully finished externalblast for {args.taxon_dir}.")

if __name__ == "__main__":
    main()