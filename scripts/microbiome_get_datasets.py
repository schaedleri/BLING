#!/usr/bin/env python3
import sys
import os
import re
import subprocess
from pathlib import Path
from multiprocessing import Pool, cpu_count

# --- Global Configuration ---
# Set the number of parallel processes.
# The original Perl script used 1. cpu_count() uses all available cores.
NUM_PROCESSES = 1
MAX_RETRIES = 3

def find_tsv_files(root_dir):
    """
    Recursively finds all files ending with .tsv in a given directory.
    """
    tsv_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower().endswith('.tsv'):
                tsv_files.append(os.path.join(dirpath, filename))
    return tsv_files

def download_and_unzip_worker(job):
    """
    A worker function for a single download/unzip job.
    Designed to be run in a multiprocessing Pool.
    It returns a status tuple: (is_success, accession, message).
    """
    accession, organism, rawgbff_dir = job
    
    # Define paths for the zip file and the final unzipped directory.
    zipfile_path = os.path.join(rawgbff_dir, f"{organism}_{accession}.zip")
    unzip_dir_path = os.path.join(rawgbff_dir, f"{organism}_{accession}")

    # --- Skip if already processed ---
    if os.path.isdir(unzip_dir_path):
        # Using print directly in worker is okay for simple progress, but can be messy.
        # print(f"Skip existing: {unzip_dir_path}")
        return (True, accession, "Skipped, already exists")

    # --- Download Command with Retries ---
    download_cmd = [
        "datasets", "download", "genome", "accession", accession,
        "--assembly-source", "RefSeq", "--reference",
        "--filename", zipfile_path, "--include", "gbff"
    ]
    
    download_success = False
    for attempt in range(1, MAX_RETRIES + 1):
        if attempt > 1:
            print(f"Warning: Attempt {attempt} for {accession} ({organism})", file=sys.stderr)
        
        result = subprocess.run(download_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            download_success = True
            break
        else:
            # Keep the last error message
            last_error = result.stderr.strip()
            
    if not download_success:
        message = f"Failed to download after {MAX_RETRIES} attempts. Last error: {last_error}"
        return (False, accession, message)

    # --- Unzip Command ---
    os.makedirs(unzip_dir_path, exist_ok=True)
    unzip_cmd = ["unzip", "-q", zipfile_path, "-d", unzip_dir_path]
    unzip_result = subprocess.run(unzip_cmd, capture_output=True, text=True)

    if unzip_result.returncode == 0:
        print(f"Unzipped: {zipfile_path} -> {unzip_dir_path}")
        try:
            os.remove(zipfile_path) # Cleanup the zip file
        except OSError as e:
            # This is not a fatal error, just warn.
            print(f"Warning: Failed to remove zip file {zipfile_path}: {e}", file=sys.stderr)
        return (True, accession, "Success")
    else:
        message = f"Failed to unzip {zipfile_path}. Error: {unzip_result.stderr.strip()}"
        return (False, accession, message)


def main():
    """
    Main function to find TSVs and orchestrate the download process.
    """
    # === Argument Handling ===
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <sample_output_directory>", file=sys.stderr)
        sys.exit(1)

    sample_output_dir = sys.argv[1]
    if not os.path.isdir(sample_output_dir):
        print(f"Error: Sample output directory not found: {sample_output_dir}", file=sys.stderr)
        sys.exit(1)

    # === Path Configuration ===
    try:
        # Determine project root assuming this script is in a 'scripts' subdirectory.
        project_root = Path(__file__).resolve().parent.parent
    except NameError:
        # Fallback for interactive environments like Jupyter.
        project_root = Path(os.getcwd())

    # Define shared and sample-specific paths.
    rawgbff_dir = project_root / "data" / "rawgbff"
    failed_file_path = os.path.join(sample_output_dir, "failed_downloads.txt")

    print(f"? Genome data will be stored in shared directory: {rawgbff_dir}")
    os.makedirs(rawgbff_dir, exist_ok=True)

    # === Find Target TSV Files ===
    target_tsvs = find_tsv_files(sample_output_dir)
    if not target_tsvs:
        print(f"Error: No .tsv files found in {sample_output_dir}. Did step 1 (assembly) run correctly?", file=sys.stderr)
        sys.exit(1)

    # === Prepare Download Jobs ===
    jobs = []
    seen_accessions = set()
    for tsv_file in target_tsvs:
        print(f"Processing accession file: {tsv_file}")
        try:
            with open(tsv_file, 'r', encoding='utf-8') as f:
                next(f) # Skip header
                for line in f:
                    if not line.strip(): continue
                    cols = line.strip().split('\t')
                    
                    accession = cols[0]
                    # Validate accession format and check for duplicates.
                    if re.match(r'^GC[AF]_\d+\.\d+$', accession) and accession not in seen_accessions:
                        organism = cols[1] if len(cols) > 1 else 'unknown'
                        # Sanitize organism name for use in filenames.
                        safe_organism = re.sub(r'[^\w.-]', '_', organism)
                        
                        jobs.append((accession, safe_organism, str(rawgbff_dir)))
                        seen_accessions.add(accession)
        except Exception as e:
            print(f"Warning: Could not process file {tsv_file}: {e}", file=sys.stderr)
            
    if not jobs:
        print("No valid accessions found to download.")
        return

    # === Run Jobs in Parallel ===
    print(f"Starting download of {len(jobs)} accessions using {NUM_PROCESSES} processes...")
    with Pool(processes=NUM_PROCESSES) as pool:
        results = pool.map(download_and_unzip_worker, jobs)

    # === Log Failures ===
    failures = [res for res in results if not res[0]]
    if failures:
        print(f"\nEncountered {len(failures)} failed jobs. Logging to {failed_file_path}")
        try:
            with open(failed_file_path, 'a', encoding='utf-8') as f_fail:
                for _, accession, message in failures:
                    f_fail.write(f"Accession: {accession} - {message}\n")
        except IOError as e:
            print(f"Error: Could not write to failed log {failed_file_path}: {e}", file=sys.stderr)
            
    print("\nDownload & unzip complete.")
    print(f"Genomic data is stored in: {rawgbff_dir}")
    if failures:
        print(f"Failed list is stored in: {failed_file_path}")
    else:
        print("All jobs completed successfully.")


if __name__ == "__main__":
    main()