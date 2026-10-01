#!/usr/bin/env python3
"""
Internal (within-genus) BLASTP screen.

For each genus directory, queries every protein against the ten partitioned
reference databases (BSDB0-BSDB9), restricted to other taxids within the same
genus, and records proteins with no significant hit as candidates for the
external (community-wide) screen.
"""
import os
import sys
import re
import argparse
import fcntl
import glob
from multiprocessing import Pool
from pathlib import Path
import shlex

DBS_CONFIG = []


def setup_dbs_config(project_root):
    """Initializes the global DBS_CONFIG list with the ten partitioned databases."""
    global DBS_CONFIG
    db_map = {f"DB{i}": f"BSDB{i}" for i in range(1, 10)}
    db_map["DB10"] = "BSDB0"

    for i in range(1, 11):
        db_num_str = f"DB{i}"
        db_name = db_map[db_num_str]
        dir_num = 0 if i == 10 else i
        db_path = os.path.join(project_root, "DB", "bacteria_strain_taxid_DB", f"DB{dir_num}", db_name)
        DBS_CONFIG.append({"db": db_path, "number": db_num_str})


def load_ids(filepath):
    """Loads protein IDs from a file into a set for quick lookups."""
    ids = set()
    if not os.path.exists(filepath):
        return ids
    try:
        with open(filepath, 'r') as f:
            for line in f:
                match = re.search(r'WP_\d+\.\d+', line)
                if match:
                    ids.add(match.group(0))
    except IOError as e:
        print(f"[ERROR] Could not read file {filepath}: {e}", file=sys.stderr)
    return ids


def append_ids_locked(filepath, ids_to_add):
    """Appends a list of IDs to a file with an exclusive lock."""
    if not ids_to_add:
        return
    try:
        with open(filepath, 'a') as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            for item in ids_to_add:
                f.write(f"{item}\n")
            fcntl.flock(f, fcntl.LOCK_UN)
    except IOError as e:
        print(f"[ERROR] Could not write to file {filepath}: {e}", file=sys.stderr)


def process_species_dir(genus_dir):
    """
    Handles directories prefixed with 'Species_': extracts protein IDs from the
    CD-HIT FASTA directly, without running BLAST (no within-species screen needed).
    """
    genus_name = os.path.basename(genus_dir)
    multi_fasta = os.path.join(genus_dir, "fasta", f"{genus_name}_all_sequences_cdhit")
    cdhit_dir = os.path.join(genus_dir, "CDhit")
    species_txt = os.path.join(cdhit_dir, "species.txt")

    os.makedirs(cdhit_dir, exist_ok=True)

    if not os.path.exists(multi_fasta):
        print(f"[ERROR] Fasta file not found: {multi_fasta}", file=sys.stderr)
        return

    count = 0
    try:
        with open(multi_fasta, 'r') as infile, open(species_txt, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    match = re.search(r'(WP_\d+\.\d+)', line)
                    if match:
                        outfile.write(f"{match.group(1)}\n")
                        count += 1
        print(f"[INFO] Wrote {count} species IDs to {species_txt}")
    except IOError as e:
        print(f"[ERROR] File operation failed for {genus_name}: {e}", file=sys.stderr)


def process_record(args):
    """Worker function: runs BLASTP for a single protein record against all ten partitions."""
    header, seq, org_label, output_dir, species2taxid, genus_taxids = args

    id_match = re.search(r'(WP_\d+\.\d+)', header)
    if not id_match:
        return
    protein_id = id_match.group(1)

    org_out_dir = os.path.join(output_dir, org_label)
    os.makedirs(org_out_dir, exist_ok=True)

    tmp_fasta = os.path.join(org_out_dir, f"{protein_id}.fasta")
    try:
        with open(tmp_fasta, 'w') as tf:
            tf.write(f">{header}\n{seq}\n")
    except IOError as e:
        print(f"[ERROR] Could not write temporary fasta {tmp_fasta}: {e}", file=sys.stderr)
        return

    new_hits = []
    new_queries = []

    for db_info in DBS_CONFIG:
        db_path = db_info["db"]
        db_number = db_info["number"]
        blast_out_file = os.path.join(org_out_dir, f"{protein_id}_{db_number}.tsv")

        genus_name, species_name_part = (org_label.split('_', 1) + [None])[:2]
        self_taxid = species2taxid.get(f"{genus_name}_{species_name_part}")
        other_taxids = [taxid for taxid in genus_taxids.get(genus_name, []) if taxid != self_taxid]
        taxid_str = ",".join(other_taxids)

        if not taxid_str:
            genus_txt_path = os.path.join(output_dir, "genus.txt")
            append_ids_locked(genus_txt_path, [protein_id])
            continue

        outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq'
        # seg yes masks low-complexity regions in the query, preventing a short
        # compositionally biased stretch from producing a spurious significant hit.
        cmd = [
            'blastp',
            '-query', tmp_fasta,
            '-db', db_path,
            '-task', 'blastp-fast',
            '-max_target_seqs', '10',
            '-evalue', '1e-10',
            '-seg', 'yes',
            '-taxids', taxid_str,
            '-outfmt', outfmt,
            '-out', blast_out_file
        ]

        cmd_str = ' '.join(shlex.quote(s) for s in cmd)
        exit_code = os.system(cmd_str)
        if exit_code != 0:
            print(f"[ERROR] blastp failed for {protein_id} DB {db_number}", file=sys.stderr)

        if os.path.exists(blast_out_file) and os.path.getsize(blast_out_file) > 0:
            try:
                with open(blast_out_file, 'r') as f:
                    for line in f:
                        cols = line.strip().split('\t')
                        if len(cols) < 16:
                            continue
                        try:
                            qcovs = float(cols[14])
                            pident = float(cols[15])
                            qlen = int(cols[4])
                            slen = int(cols[9])
                        except (ValueError, IndexError):
                            continue

                        is_good = False
                        if pident >= 95 and qcovs >= 95 and slen > 0:
                            ratio = qlen / slen
                            if 0.95 <= ratio <= 1.05:
                                is_good = True

                        if is_good:
                            s_match = re.search(r'(WP_\d+\.\d+)', cols[5])
                            if s_match:
                                new_hits.append(s_match.group(1))
                                new_queries.append(cols[0])
            except IOError as e:
                print(f"[ERROR] Could not read blast output {blast_out_file}: {e}", file=sys.stderr)

    os.remove(tmp_fasta)

    used_hits_file = os.path.join(output_dir, "used_hits.txt")
    used_queries_file = os.path.join(output_dir, "used_queries.txt")
    append_ids_locked(used_hits_file, new_hits)
    append_ids_locked(used_queries_file, new_queries)


def write_hit_query_pairs(cdhit_dir):
    """Scans all BLAST output TSV files and writes a summary of hit-query pairs."""
    tsv_files = glob.glob(os.path.join(cdhit_dir, '*', '*.tsv'))
    pairs = []
    for file in tsv_files:
        try:
            with open(file, 'r') as f:
                for line in f:
                    cols = line.strip().split('\t')
                    if len(cols) < 6:
                        continue
                    subject_match = re.search(r'(WP_\d+\.\d+)', cols[5])
                    if subject_match:
                        pairs.append((subject_match.group(1), cols[0]))
        except IOError:
            continue

    out_file = os.path.join(cdhit_dir, 'hit_query_pairs.tsv')
    try:
        with open(out_file, 'w') as out:
            for pair in pairs:
                out.write(f"{pair[0]}\t{pair[1]}\n")
        print(f"[INFO] Wrote {len(pairs)} pairs to {out_file}")
    except IOError as e:
        print(f"[ERROR] Cannot open {out_file}: {e}", file=sys.stderr)


def generate_genus_txt(cdhit_root):
    """
    Updates genus.txt: an ID is added if at least one of its ten per-partition
    BLAST result files is non-empty (i.e. a hit was found somewhere).
    """
    if not os.path.isdir(cdhit_root):
        return

    genus_txt_path = os.path.join(cdhit_root, "genus.txt")
    already_present_ids = load_ids(genus_txt_path)

    potential_records = []
    for subdir in glob.glob(os.path.join(cdhit_root, "*")):
        if not os.path.isdir(subdir):
            continue
        for file in glob.glob(os.path.join(subdir, "*_DB1.tsv")):
            match = re.search(r'([^\/\\]+)_DB1\.tsv$', os.path.basename(file))
            if match:
                protein_id = match.group(1)
                if protein_id not in already_present_ids:
                    potential_records.append({'id': protein_id, 'dir': subdir})

    new_ids_to_add = []
    for rec in potential_records:
        protein_id = rec['id']
        has_nonzero_hit = any(
            os.path.exists(f_path) and os.path.getsize(f_path) > 0
            for f_path in (os.path.join(rec['dir'], f"{protein_id}_DB{i}.tsv") for i in range(1, 11))
        )
        if has_nonzero_hit:
            new_ids_to_add.append(protein_id)

    all_ids = sorted(already_present_ids.union(new_ids_to_add))
    try:
        with open(genus_txt_path, 'w') as f:
            for protein_id in all_ids:
                f.write(f"{protein_id}\n")
        print(f"[INFO] Updated {genus_txt_path} with {len(new_ids_to_add)} new IDs.")
    except IOError as e:
        print(f"[ERROR] Could not write to {genus_txt_path}: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Run the internal (within-genus) BLASTP screen in parallel.",
        usage="%(prog)s --base-dir <sample_dir> [--cpu <n>] <Genus_or_Species_dir> ..."
    )
    parser.add_argument('--base-dir', required=True, help='Path to the sample directory')
    parser.add_argument('--cpu', type=int, default=10, help='Number of CPUs to use')
    parser.add_argument('target_dirs', nargs='+', help='One or more target directories (e.g., Genus_Xxx)')
    args, _ = parser.parse_known_args()

    sample_dir = Path(args.base_dir).resolve()
    project_root = sample_dir.parent

    if not sample_dir.is_dir():
        print(f"Error: Sample directory '{sample_dir}' is not a valid directory.", file=sys.stderr)
        sys.exit(1)
    if not (project_root / "DB").is_dir():
        print(f"Error: DB directory not found under project root '{project_root}'", file=sys.stderr)
        sys.exit(1)

    setup_dbs_config(project_root)

    for genus_name in args.target_dirs:
        print(f"\n--- Processing target: {genus_name} ---")
        genus_path = sample_dir / genus_name

        if genus_name.lower().startswith('species_'):
            process_species_dir(str(genus_path))
            continue

        fasta_file = genus_path / "fasta" / f"{genus_name}_all_sequences_cdhit"
        tsv_file = genus_path / f"{genus_name}.tsv"
        output_dir = genus_path / "CDhit"
        output_dir.mkdir(exist_ok=True)

        used_hits_file = output_dir / "used_hits.txt"
        used_queries_file = output_dir / "used_queries.txt"
        used_hits_file.touch()
        used_queries_file.touch()
        used_hits = load_ids(str(used_hits_file))

        species2taxid = {}
        genus_taxids = {}
        if not tsv_file.is_file():
            print(f"[ERROR] Taxid TSV file not found: {tsv_file}", file=sys.stderr)
            continue

        with open(tsv_file, 'r') as f:
            next(f)
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 4:
                    continue
                org_name, taxid = cols[1], cols[3]
                if not (org_name and taxid):
                    continue
                parts = org_name.split()
                if len(parts) >= 2:
                    genus, species = parts[0], parts[1]
                    species2taxid[f"{genus}_{species}"] = taxid
                    genus_taxids.setdefault(genus, []).append(taxid)

        if not fasta_file.is_file():
            print(f"[ERROR] Fasta file not found: {fasta_file}", file=sys.stderr)
            continue

        tasks = []
        with open(fasta_file, 'r') as mf:
            header, seq = '', ''
            for line in mf:
                if line.startswith('>'):
                    if header and seq:
                        id_match = re.search(r'(WP_\d+\.\d+)', header)
                        if id_match and id_match.group(1) not in used_hits:
                            org_match = re.search(r'\[([^\]]+)\]', header)
                            org_label = org_match.group(1).replace(' ', '_') if org_match else 'other'
                            tasks.append((header.strip(), seq, org_label, str(output_dir), species2taxid, genus_taxids))
                    header = line[1:]
                    seq = ''
                else:
                    seq += line.strip()
            if header and seq:
                id_match = re.search(r'(WP_\d+\.\d+)', header)
                if id_match and id_match.group(1) not in used_hits:
                    org_match = re.search(r'\[([^\]]+)\]', header)
                    org_label = org_match.group(1).replace(' ', '_') if org_match else 'other'
                    tasks.append((header.strip(), seq, org_label, str(output_dir), species2taxid, genus_taxids))

        print(f"[INFO] Found {len(tasks)} new records to process.")
        if tasks:
            with Pool(processes=args.cpu) as pool:
                pool.map(process_record, tasks)

        write_hit_query_pairs(str(output_dir))

    for genus_name in args.target_dirs:
        if not genus_name.lower().startswith('species_'):
            generate_genus_txt(str(sample_dir / genus_name / "CDhit"))


if __name__ == '__main__':
    main()
