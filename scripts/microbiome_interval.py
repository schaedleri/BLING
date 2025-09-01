#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import re
import glob
import argparse
from pathlib import Path
from collections import defaultdict

# (このスクリプトの他のヘルパー関数部分は変更ありません。main関数のみを修正・掲載します)
def load_go_classifications(obo_file):
    """Parses a go.obo file to get namespaces, needed for filtering InterPro GO terms."""
    classifications = {}
    if not Path(obo_file).is_file(): return classifications
    try:
        with open(obo_file, 'r', encoding='utf-8') as f:
            content = f.read()
            term_blocks = content.split('[Term]\n')[1:]
            for block in term_blocks:
                goid_match = re.search(r'^id: (GO:\d+)', block, re.MULTILINE)
                ns_match = re.search(r'^namespace: (\S+)', block, re.MULTILINE)
                if goid_match:
                    goid = goid_match.group(1)
                    namespace = ns_match.group(1) if ns_match else 'unknown'
                    classifications[goid] = namespace
    except IOError:
        print(f"Warning: Could not read go.obo file at {obo_file}", file=sys.stderr)
    return classifications

def read_go_data(file_path):
    """Reads a _go_annotations.tsv file."""
    go_data = defaultdict(list)
    if not file_path.is_file(): return go_data
    with open(file_path, 'r', encoding='utf-8') as f:
        next(f) # Skip header
        for line in f:
            if not line.strip(): continue
            try:
                acc, org, go_id, go_type = line.strip().split('\t')
                if go_type in ['biological_process', 'molecular_function', 'cellular_component']:
                     go_data[acc].append({'id': go_id, 'type': go_type})
            except ValueError: continue
    return go_data

def read_ip_data(file_path, go_classifications):
    """Reads an InterProScan TSV and returns a set of protein IDs with valid IP hits."""
    ip_data = set()
    if not file_path.is_file(): return ip_data
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 14: continue
            prot_id, n_col, go_field = cols[0], cols[12], cols[13]
            go_terms = re.findall(r'(GO:\d+)', go_field)
            ns_count = defaultdict(int)
            for go_id in go_terms:
                ns = go_classifications.get(go_id, "unknown")
                ns_count[ns] += 1
            if list(ns_count.keys()) == ['cellular_component']: continue
            if n_col != '-': ip_data.add(prot_id)
    return ip_data

def read_ip_go_data(file_path, go_classifications):
    """Reads InterProScan TSV and extracts functional GO terms."""
    ip_go_data = defaultdict(list)
    if not file_path.is_file(): return ip_go_data
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 14: continue
            prot_id, go_field = cols[0], cols[13]
            go_terms = re.findall(r'(GO:\d+)', go_field)
            for go_id in go_terms:
                ns = go_classifications.get(go_id, "")
                if ns == 'molecular_function': ip_go_data[prot_id].append({'id': go_id, 'type': 'function'})
                elif ns == 'biological_process': ip_go_data[prot_id].append({'id': go_id, 'type': 'process'})
    return ip_go_data

def read_genomic_data(file_path):
    """Reads a _genomic.tsv file."""
    genomic_data = defaultdict(list)
    with open(file_path, 'r', encoding='utf-8') as f:
        header = next(f).strip().split('\t')
        for line in f:
            if not line.strip(): continue
            values = line.strip().split('\t')
            gene_info = dict(zip(header, values))
            try:
                gene_info['begin'] = int(gene_info['Begin'])
                gene_info['end'] = int(gene_info['End'])
            except (ValueError, KeyError): continue
            genomic_data[gene_info['Accession']].append(gene_info)
    return genomic_data

def generate_clusters(genomic_data, interval=100):
    """Groups genes into clusters."""
    clusters = []
    for accession in genomic_data:
        genes = sorted(genomic_data[accession], key=lambda x: x['begin'])
        if not genes: continue
        current_cluster = [genes[0]]
        for i in range(1, len(genes)):
            gap = genes[i]['begin'] - genes[i-1]['end']
            same_strand = genes[i]['Strand'] == genes[i-1]['Strand']
            if gap <= interval and same_strand:
                current_cluster.append(genes[i])
            else:
                if len(current_cluster) > 2: clusters.append(current_cluster)
                current_cluster = [genes[i]]
        if len(current_cluster) > 2: clusters.append(current_cluster)
    return clusters

def annotate_cluster(cluster, go_data, ip_data):
    """Annotates each gene in a cluster with a status: 'GO', 'IP', or 'nonGO'."""
    for gene in cluster:
        pid = gene['Protein_ID']
        if pid in go_data: gene['status'] = 'GO'
        elif pid in ip_data: gene['status'] = 'IP'
        else: gene['status'] = 'nonGO'

def is_mixed_cluster(cluster):
    """Checks if a cluster contains both 'nonGO' and other genes."""
    statuses = {gene['status'] for gene in cluster}
    return 'nonGO' in statuses and ('GO' in statuses or 'IP' in statuses)

def print_cluster(file_handle, cluster):
    """Prints all genes in a cluster to a file handle."""
    header = ['Accession', 'Organism', 'Begin', 'End', 'Strand', 'Product', 'Gene', 'Locus_Tag', 'Protein_Length', 'Protein_ID', 'status']
    for gene in cluster:
        row = [str(gene.get(h, '')) for h in header]
        file_handle.write("\t".join(row) + "\n")
    file_handle.write("\n")

def find_nearest_neighbors(cluster, cluster_id, file_handle, max_dist=300):
    """Finds and prints pairs of nonGO genes and their nearby GO/IP neighbors."""
    for gene in cluster:
        if gene['status'] != 'nonGO': continue
        for other in cluster:
            if gene is other or other['status'] == 'nonGO': continue
            dist = 0
            if gene['begin'] > other['end']: dist = gene['begin'] - other['end']
            elif gene['end'] < other['begin']: dist = other['begin'] - gene['end']
            if dist <= max_dist:
                file_handle.write(f"{cluster_id}\t{gene['Protein_ID']}\tneighbor_{other['status']}\t{dist}\t{other['Protein_ID']}\n")

def output_go_pairs(cluster, file_handle, go_data, ip_go_data, max_dist=300):
    """Finds nonGO genes and outputs pairs with the GO terms of their neighbors."""
    for gene in cluster:
        if gene['status'] != 'nonGO': continue
        for other in cluster:
            if gene is other or other['status'] == 'nonGO': continue
            dist = 0
            if gene['begin'] > other['end']: dist = gene['begin'] - other['end']
            elif gene['end'] < other['begin']: dist = other['begin'] - gene['end']
            if dist > max_dist: continue
            neighbor_pid = other['Protein_ID']
            go_terms_to_print = []
            if other['status'] == 'GO' and neighbor_pid in go_data:
                go_terms_to_print.extend(go_data[neighbor_pid])
            if other['status'] == 'IP' and neighbor_pid in ip_go_data:
                go_terms_to_print.extend(ip_go_data[neighbor_pid])
            for go in go_terms_to_print:
                 if go.get('type') in ['biological_process', 'molecular_function', 'function', 'process']:
                    go_type = go['type']
                    if go_type == 'molecular_function': go_type = 'function'
                    if go_type == 'biological_process': go_type = 'process'
                    file_handle.write(f"{gene['Protein_ID']}\t{go['id']}\t{go_type}\n")

def output_all_ip_go_pairs(ip_file, go_classifications, file_handle):
    """Reads an InterProScan file and prints all Protein-GO pairs."""
    if not ip_file.is_file(): return
    with open(ip_file, 'r', encoding='utf-8') as f_in:
        for line in f_in:
            cols = line.strip().split('\t')
            if len(cols) < 14: continue
            prot_id, go_field = cols[0], cols[13]
            go_terms = re.findall(r'(GO:\d+)', go_field)
            for go_id in go_terms:
                ns = go_classifications.get(go_id, "unknown")
                file_handle.write(f"{prot_id}\t{go_id}\t{ns}\n")

def main():
    parser = argparse.ArgumentParser(description="Analyze gene intervals and proximity to GO-annotated genes.")
    parser.add_argument('--base-dir', type=str, required=True, help="The base directory for a specific sample.")
    args, unknown = parser.parse_known_args() # For Jupyter compatibility

    sample_dir = Path(args.base_dir).resolve()
    project_root = sample_dir.parent
    go_obo_file = project_root / "data" / "go.obo"

    print("Loading GO classifications (for filtering)...")
    go_classifications = load_go_classifications(go_obo_file)

    # === ★★★ この部分を修正 ★★★ ===
    # Find all classification-level directories (e.g., Species_..., Family_...).
    # The previous glob pattern was incorrect for Python.
    class_dirs = []
    for item in sample_dir.iterdir():
        if item.is_dir() and re.match(r'^(Class|Family|Genus|Order|Species)_', item.name):
            class_dirs.append(item)
    # === ★★★★★★★★★★★★★★★ ===
    
    if not class_dirs:
        print(f"Warning: No classification directories found in {sample_dir}", file=sys.stderr)
        return

    for class_dir in class_dirs:
        genomic_dir = class_dir / 'genomic_annotations'
        if not genomic_dir.is_dir(): continue
        
        for genomic_file in genomic_dir.glob('*.tsv'):
            basename = genomic_file.stem.replace('_genomic', '')
            if not basename: continue
            
            print(f"Processing {basename} in {class_dir.name}...")
            
            output_dir = class_dir / "interval_test" / basename
            output_dir.mkdir(parents=True, exist_ok=True)
            
            go_file = class_dir / 'go_annotations' / f"{basename}_go_annotations.tsv"
            ip_file = class_dir / 'interpro' / f"{basename}.tsv"

            go_data = read_go_data(go_file)
            ip_data = read_ip_data(ip_file, go_classifications)
            ip_go_data = read_ip_go_data(ip_file, go_classifications)
            genomic_data = read_genomic_data(genomic_file)

            clusters = generate_clusters(genomic_data)

            if not clusters:
                print(f"  -> No valid gene clusters found for {basename}. Skipping output file generation.")
                continue

            try:
                with open(output_dir / "interval_afterpsi100.tsv", 'w') as fh_all, \
                     open(output_dir / "filtered_interval_afterpsi100.tsv", 'w') as fh_filtered, \
                     open(output_dir / "nearest_results100.tsv", 'w') as fh_nearest, \
                     open(output_dir / "nonGO_GO_pairs100.tsv", 'w') as fh_pairs, \
                     open(output_dir / "IP_GO_pairs.tsv", 'w') as fh_ipgo:
                    
                    fh_ipgo.write("IP_Protein_ID\tGO_ID\tGO_Type\n")
                    fh_all.write("Accession\tOrganism\tBegin\tEnd\tStrand\tProduct\tGene\tLocus_Tag\tProtein_Length\tProtein_ID\tstatus\n")
                    fh_filtered.write("Accession\tOrganism\tBegin\tEnd\tStrand\tProduct\tGene\tLocus_Tag\tProtein_Length\tProtein_ID\tstatus\n")
                    fh_nearest.write("Cluster_ID\tNonGO_Protein_ID\tNeighbor_Status\tDistance\tNeighbor_Protein_ID\n")
                    fh_pairs.write("NonGO_Protein_ID\tGO_Protein_ID\tGO_Type\n")
                    
                    output_all_ip_go_pairs(ip_file, go_classifications, fh_ipgo)
                    
                    cluster_id = 1
                    for cluster in clusters:
                        annotate_cluster(cluster, go_data, ip_data)
                        print_cluster(fh_all, cluster)
                        if is_mixed_cluster(cluster): print_cluster(fh_filtered, cluster)
                        find_nearest_neighbors(cluster, cluster_id, fh_nearest)
                        output_go_pairs(cluster, fh_pairs, go_data, ip_go_data)
                        cluster_id += 1
            except IOError as e:
                print(f"Error writing output files for {basename}: {e}", file=sys.stderr)
                continue
            print(f"Finished processing {basename} in {class_dir.name}.")

if __name__ == "__main__":
    main()