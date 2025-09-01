#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import argparse
import subprocess
import shutil
from pathlib import Path

def main():
    """
    Main function to concatenate FASTA files and run CD-HIT.
    """
    parser = argparse.ArgumentParser(description="Concatenate FASTA files and run CD-HIT.")
    
    # --- ▼▼▼ 修正箇所 ▼▼▼ ---
    # --sample-dir を --base-dir に変更
    parser.add_argument('--base-dir', type=str, required=True,
                        help="The path to the sample's base directory (e.g., 'depression').")
    parser.add_argument('target_dirs', nargs='+',
                        help="One or more target directories to process (e.g., 'Genus_Xxx').")
    
    # Jupyterなどからの未知の引数を無視するために parse_known_args() を使用
    args, unknown = parser.parse_known_args()
    
    base_dir = Path(args.base_dir).resolve()
    # --- ▲▲▲ 修正箇所 ▲▲▲ ---

    for group_name in args.target_dirs:
        print(f"\n--- Processing group: {group_name} ---")
        group_dir = base_dir / group_name
        fasta_dir = group_dir / "fasta"

        if not fasta_dir.is_dir():
            print(f"Warning: Fasta directory not found, skipping: {fasta_dir}", file=sys.stderr)
            continue

        print(f"Searching for .fasta files in {fasta_dir}...")
        fasta_files = list(fasta_dir.rglob('*.fasta'))
        
        if not fasta_files:
            print(f"Warning: No .fasta files found under {fasta_dir}, skipping.", file=sys.stderr)
            continue
            
        print(f"Found {len(fasta_files)} FASTA files.")

        output_fasta_path = fasta_dir / f"{group_name}_all_sequences.fasta"
        cdhit_output_path = fasta_dir / f"{group_name}_all_sequences_cdhit"

        print(f"Concatenating files into {output_fasta_path}...")
        try:
            with open(output_fasta_path, 'wb') as f_out:
                for file_path in fasta_files:
                    with open(file_path, 'rb') as f_in:
                        shutil.copyfileobj(f_in, f_out)
        except IOError as e:
            print(f"Error: Cannot write to {output_fasta_path}: {e}", file=sys.stderr)
            continue

        print(f"Running CD-HIT...")
        cmd = ["cd-hit", "-i", str(output_fasta_path), "-o", str(cdhit_output_path), "-c", "0.9"]
        
        try:
            # 実行結果を直接表示するように変更
            subprocess.run(cmd, check=True, text=True)
            print(f"? CD-HIT successfully completed for {group_name}")
        except FileNotFoundError:
            print("Error: 'cd-hit' command not found. Is it in your PATH?", file=sys.stderr)
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            print(f"Error: Failed to run cd-hit for {group_name}", file=sys.stderr)
            # Stderrはsubprocess.runが自動で表示するので、ここでは表示しない
            sys.exit(1)

if __name__ == "__main__":
    main()