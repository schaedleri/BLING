#!/usr/bin/env python3
"""Run InterProScan on the .faa files of every taxon directory (<taxon>/faa/*.faa -> <taxon>/interpro/*.tsv).
A result already produced for the same file name in another taxon directory is copied instead of recomputed."""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

TAXON_LEVELS = ('Class', 'Order', 'Family', 'Genus', 'Species')
TOOLS_DIR = 'tools'
INTERPROSCAN_DIR_PATTERN = 'interproscan-*'
INTERPROSCAN_SCRIPT = 'interproscan.sh'
DEFAULT_APPLICATIONS = 'NCBIfam,SUPERFAMILY,PANTHER,ProSitePatterns,Pfam'


def find_interproscan(base_dir, given_dir):
    if given_dir:
        directory = Path(given_dir).resolve()
    else:
        tools_dir = base_dir.parent / TOOLS_DIR
        candidates = sorted(p for p in tools_dir.glob(INTERPROSCAN_DIR_PATTERN) if p.is_dir())
        if not candidates:
            sys.exit(f"Error: no uncompressed InterProScan directory under {tools_dir}.")
        directory = candidates[0]
    script = directory / INTERPROSCAN_SCRIPT
    if not (script.is_file() and os.access(script, os.X_OK)):
        sys.exit(f"Error: '{script}' does not exist or is not executable.")
    return script


def run_interproscan(script, faa, output, cpu, applications):
    command = [str(script), '-i', str(faa), '-f', 'tsv', '-o', str(output), '--goterms',
               '-cpu', str(cpu), '-appl', applications]
    try:
        result = subprocess.run(command, capture_output=True, text=True)
    except Exception as e:
        print(f"[!] Could not run InterProScan for {faa.name}: {e}", file=sys.stderr)
        return False
    if result.returncode != 0:
        print(f"[!] InterProScan failed for {faa.name}", file=sys.stderr)
        if result.stderr:
            print(f"[!] Stderr: {result.stderr.strip()}", file=sys.stderr)
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--base-dir', required=True, help="Sample directory containing the taxon directories.")
    parser.add_argument('--interproscan-dir', help=f"InterProScan directory; default: first '{INTERPROSCAN_DIR_PATTERN}' "
                                                   f"directory under <project root>/{TOOLS_DIR}.")
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--appl', default=DEFAULT_APPLICATIONS, help="Comma-separated InterProScan applications.")
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    if not base_dir.is_dir():
        sys.exit(f"Error: --base-dir '{args.base_dir}' not found or not a directory.")
    script = find_interproscan(base_dir, args.interproscan_dir)
    print(f"Using InterProScan at: {script}")

    done = {}
    for level in TAXON_LEVELS:
        for faa in base_dir.glob(f"{level}_*/faa/*.faa"):
            out_dir = faa.parent.parent / 'interpro'
            out_dir.mkdir(exist_ok=True)
            output = out_dir / f"{faa.stem}.tsv"
            if output.exists():
                print(f"[?] Existing output for {faa.stem}; skipped.")
                done.setdefault(faa.stem, output)
            elif faa.stem in done:
                print(f"[->] Copying {done[faa.stem]} -> {output}")
                try:
                    shutil.copy(done[faa.stem], output)
                except Exception as e:
                    print(f"[!] Copy failed: {e}", file=sys.stderr)
            else:
                print(f"[*] InterProScan: {faa.stem}")
                if run_interproscan(script, faa, output, args.cpu, args.appl):
                    done[faa.stem] = output
    print("InterProScan processing finished.")


if __name__ == '__main__':
    main()
