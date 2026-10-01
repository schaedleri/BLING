#!/usr/bin/env python3
"""
Fixed-threshold gene-neighborhood operon prediction and GO propagation.

Clusters adjacent genes on the same strand separated by <= 100 bp, then
propagates Biological Process GO terms from annotated genes to unannotated
("nonGO") neighbors within 300 bp, matching the GMM-based method's scope
(see predict_and_analyze_operons.py).
"""
import sys
import os
import re
import argparse
from pathlib import Path
from collections import defaultdict


def load_go_classifications(obo_file):
    """Parses a go.obo file to get GO-term namespaces."""
    classifications = {}
    if not Path(obo_file).is_file():
        return classifications
    try:
        with open(obo_file, 'r', encoding='utf-8') as f:
            content = f.read()
            for block in content.split('[Term]\n')[1:]:
                goid_match = re.search(r'^id: (GO:\d+)', block, re.MULTILINE)
                ns_match = re.search(r'^namespace: (\S+)', block, re.MULTILINE)
                if goid_match:
                    classifications[goid_match.group(1)] = ns_match.group(1) if ns_match else 'unknown'
    except IOError:
        print(f"Warning: Could not read go.obo file at {obo_file}", file=sys.stderr)
    return classifications


def read_go_data(file_path):
    """Reads a _go_annotations.tsv file (direct GBFF-derived GO annotations)."""
    go_data = defaultdict(list)
    if not file_path.is_file():
        return go_data
    with open(file_path, 'r', encoding='utf-8') as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            try:
                acc, org, go_id, go_type = line.strip().split('\t')
                if go_type in ['biological_process', 'molecular_function', 'cellular_component']:
                    go_data[acc].append({'id': go_id, 'type': go_type})
            except ValueError:
                continue
    return go_data


def read_ip_data(file_path, go_classifications):
    """Returns the set of protein IDs with a valid InterPro domain hit (excluding CC-only hits)."""
    ip_data = set()
    if not file_path.is_file():
        return ip_data
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 14:
                continue
            prot_id, n_col, go_field = cols[0], cols[12], cols[13]
            go_terms = re.findall(r'(GO:\d+)', go_field)
            ns_count = defaultdict(int)
            for go_id in go_terms:
                ns_count[go_classifications.get(go_id, "unknown")] += 1
            if list(ns_count.keys()) == ['cellular_component']:
                continue
            if n_col != '-':
                ip_data.add(prot_id)
    return ip_data


def read_ip_go_data(file_path, go_classifications):
    """Reads InterProScan output and extracts BP/MF GO terms per protein."""
    ip_go_data = defaultdict(list)
    if not file_path.is_file():
        return ip_go_data
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 14:
                continue
            prot_id, go_field = cols[0], cols[13]
            for go_id in re.findall(r'(GO:\d+)', go_field):
                ns = go_classifications.get(go_id, "")
                if ns == 'molecular_function':
                    ip_go_data[prot_id].append({'id': go_id, 'type': 'function'})
                elif ns == 'biological_process':
                    ip_go_data[prot_id].append({'id': go_id, 'type': 'process'})
    return ip_go_data


def read_eggnog_data(file_path, go_classifications):
    """Reads a pre-processed eggNOG-mapper result file (*_emapper_results.tsv)."""
    eggnog_go_data = defaultdict(list)
    if not file_path.is_file():
        return eggnog_go_data
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            try:
                prot_id_idx = header.index('Protein_ID')
                go_idx = header.index('GOs')
            except ValueError:
                print(f"Warning: 'Protein_ID' or 'GOs' column not found in {file_path}", file=sys.stderr)
                return eggnog_go_data

            for line in f:
                cols = line.strip().split('\t')
                if len(cols) <= max(prot_id_idx, go_idx):
                    continue
                prot_id, go_field = cols[prot_id_idx], cols[go_idx]
                if not go_field or go_field == '-':
                    continue
                for go_id in go_field.split(','):
                    go_id = go_id.strip()
                    if not go_id:
                        continue
                    ns = go_classifications.get(go_id, "unknown")
                    if ns in ['biological_process', 'molecular_function']:
                        eggnog_go_data[prot_id].append({'id': go_id, 'type': ns})
    except IOError as e:
        print(f"Warning: Could not read pre-processed eggNOG file at {file_path}: {e}", file=sys.stderr)
    return eggnog_go_data


def read_genomic_data(file_path):
    """Reads a _genomic.tsv file (per-gene coordinates and metadata)."""
    genomic_data = defaultdict(list)
    with open(file_path, 'r', encoding='utf-8') as f:
        header = next(f).strip().split('\t')
        for line in f:
            if not line.strip():
                continue
            gene_info = dict(zip(header, line.strip().split('\t')))
            try:
                gene_info['begin'] = int(gene_info['Begin'])
                gene_info['end'] = int(gene_info['End'])
            except (ValueError, KeyError):
                continue
            genomic_data[gene_info['Accession']].append(gene_info)
    return genomic_data


def generate_clusters(genomic_data, interval=100):
    """Groups adjacent, same-strand genes separated by <= `interval` bp into clusters of >2 genes."""
    clusters = []
    for accession in genomic_data:
        genes = sorted(genomic_data[accession], key=lambda x: x['begin'])
        if not genes:
            continue
        current_cluster = [genes[0]]
        for i in range(1, len(genes)):
            gap = genes[i]['begin'] - genes[i - 1]['end']
            same_strand = genes[i]['Strand'] == genes[i - 1]['Strand']
            if gap <= interval and same_strand:
                current_cluster.append(genes[i])
            else:
                if len(current_cluster) > 2:
                    clusters.append(current_cluster)
                current_cluster = [genes[i]]
        if len(current_cluster) > 2:
            clusters.append(current_cluster)
    return clusters


def annotate_cluster(cluster, go_data, ip_data):
    """Labels each gene in a cluster as 'GO', 'IP', or 'nonGO'."""
    for gene in cluster:
        pid = gene['Protein_ID']
        gene['status'] = 'GO' if pid in go_data else ('IP' if pid in ip_data else 'nonGO')


def is_mixed_cluster(cluster):
    """True if a cluster contains both annotated and unannotated ('nonGO') genes."""
    statuses = {gene['status'] for gene in cluster}
    return 'nonGO' in statuses and ('GO' in statuses or 'IP' in statuses)


def print_cluster(file_handle, cluster):
    """Writes every gene in a cluster to a file handle."""
    header = ['Accession', 'Organism', 'Begin', 'End', 'Strand', 'Product', 'Gene',
              'Locus_Tag', 'Protein_Length', 'Protein_ID', 'status']
    for gene in cluster:
        file_handle.write("\t".join(str(gene.get(h, '')) for h in header) + "\n")
    file_handle.write("\n")


def find_nearest_neighbors(cluster, cluster_id, file_handle, max_dist=300):
    """Writes pairs of unannotated genes and their nearby annotated neighbors."""
    for gene in cluster:
        if gene['status'] != 'nonGO':
            continue
        for other in cluster:
            if gene is other or other['status'] == 'nonGO':
                continue
            dist = max(gene['begin'] - other['end'], other['begin'] - gene['end'], 0)
            if dist <= max_dist:
                file_handle.write(f"{cluster_id}\t{gene['Protein_ID']}\tneighbor_{other['status']}\t{dist}\t{other['Protein_ID']}\n")


def output_go_pairs(cluster, cluster_id, genome_id, file_handle, go_data, ip_go_data, max_dist=300):
    """
    Propagates Biological Process GO terms from annotated genes to unannotated
    ('nonGO') neighbors within `max_dist` bp in the same cluster.

    The entire cluster is treated as a single evidence unit (one
    Evidence_Group_ID per cluster), matching the GMM-side design so that
    downstream enrichment analysis does not double-count annotations copied
    from a single source. Evidence_Group_ID includes genome_id to stay
    globally unique across the merged, community-wide annotation file.

    Only Biological Process terms are propagated, matching the scope of the
    GMM-based method (predict_and_analyze_operons.py); Molecular Function and
    Cellular Component annotations are not propagated across clusters.
    """
    evidence_group_id = f"{genome_id}_cluster{cluster_id}"
    for gene in cluster:
        if gene['status'] != 'nonGO':
            continue
        for other in cluster:
            if gene is other or other['status'] == 'nonGO':
                continue
            dist = max(gene['begin'] - other['end'], other['begin'] - gene['end'], 0)
            if dist > max_dist:
                continue
            neighbor_pid = other['Protein_ID']
            go_terms_to_print = []
            if other['status'] == 'GO' and neighbor_pid in go_data:
                go_terms_to_print.extend(go_data[neighbor_pid])
            if other['status'] == 'IP' and neighbor_pid in ip_go_data:
                go_terms_to_print.extend(ip_go_data[neighbor_pid])
            for go in go_terms_to_print:
                if go.get('type') in ['biological_process', 'process']:
                    file_handle.write(
                        f"{gene['Protein_ID']}\t{go['id']}\tprocess\t{cluster_id}\t{evidence_group_id}\n"
                    )


def output_all_ip_go_pairs(ip_file, go_classifications, file_handle):
    """Writes every Protein-GO pair found in an InterProScan output file."""
    if not ip_file.is_file():
        return
    with open(ip_file, 'r', encoding='utf-8') as f_in:
        for line in f_in:
            cols = line.strip().split('\t')
            if len(cols) < 14:
                continue
            prot_id, go_field = cols[0], cols[13]
            for go_id in re.findall(r'(GO:\d+)', go_field):
                ns = go_classifications.get(go_id, "unknown")
                file_handle.write(f"{prot_id}\t{go_id}\t{ns}\n")


def main():
    parser = argparse.ArgumentParser(description="Fixed-threshold gene-neighborhood operon prediction and GO propagation.")
    parser.add_argument('--base-dir', type=str, required=True, help="The base directory for a specific sample.")
    args, _ = parser.parse_known_args()

    sample_dir = Path(args.base_dir).resolve()
    project_root = sample_dir.parent
    go_obo_file = project_root / "data" / "go.obo"

    print("Loading GO classifications...")
    go_classifications = load_go_classifications(go_obo_file)

    class_dirs = [item for item in sample_dir.iterdir()
                  if item.is_dir() and re.match(r'^(Class|Family|Genus|Order|Species)_', item.name)]
    if not class_dirs:
        print(f"Warning: No classification directories found in {sample_dir}", file=sys.stderr)
        return

    for class_dir in class_dirs:
        genomic_dir = class_dir / 'genomic_annotations'
        if not genomic_dir.is_dir():
            continue

        for genomic_file in genomic_dir.glob('*.tsv'):
            basename = genomic_file.stem.replace('_genomic', '')
            if not basename:
                continue

            print(f"Processing {basename} in {class_dir.name}...")
            output_dir = class_dir / "interval_test" / basename
            output_dir.mkdir(parents=True, exist_ok=True)

            go_file = class_dir / 'go_annotations' / f"{basename}_go_annotations.tsv"
            ip_file = class_dir / 'interpro' / f"{basename}.tsv"
            eggnog_file = class_dir / 'eggnog' / f"{basename}_genomic_emapper_results.tsv"

            go_data = read_go_data(go_file)
            eggnog_go_data = read_eggnog_data(eggnog_file, go_classifications)
            for prot_id, go_list in eggnog_go_data.items():
                existing_go_ids = {go['id'] for go in go_data.get(prot_id, [])}
                for go_item in go_list:
                    if go_item['id'] not in existing_go_ids:
                        go_data[prot_id].append(go_item)
                        existing_go_ids.add(go_item['id'])

            ip_data = read_ip_data(ip_file, go_classifications)
            ip_go_data = read_ip_go_data(ip_file, go_classifications)
            genomic_data = read_genomic_data(genomic_file)
            clusters = generate_clusters(genomic_data)

            if not clusters:
                print(f"  -> No valid gene clusters found for {basename}. Skipping.")
                continue

            try:
                with open(output_dir / "interval_afterpsi100.tsv", 'w') as fh_all, \
                     open(output_dir / "filtered_interval_afterpsi100.tsv", 'w') as fh_filtered, \
                     open(output_dir / "nearest_results100.tsv", 'w') as fh_nearest, \
                     open(output_dir / "nonGO_GO_pairs100.tsv", 'w') as fh_pairs, \
                     open(output_dir / "IP_GO_pairs.tsv", 'w') as fh_ipgo:

                    fh_ipgo.write("IP_Protein_ID\tGO_ID\tGO_Type\n")
                    header = "Accession\tOrganism\tBegin\tEnd\tStrand\tProduct\tGene\tLocus_Tag\tProtein_Length\tProtein_ID\tstatus\n"
                    fh_all.write(header)
                    fh_filtered.write(header)
                    fh_nearest.write("Cluster_ID\tNonGO_Protein_ID\tNeighbor_Status\tDistance\tNeighbor_Protein_ID\n")
                    fh_pairs.write("NonGO_Protein_ID\tGO_ID\tGO_Type\tCluster_ID\tEvidence_Group_ID\n")

                    output_all_ip_go_pairs(ip_file, go_classifications, fh_ipgo)

                    for cluster_id, cluster in enumerate(clusters, start=1):
                        annotate_cluster(cluster, go_data, ip_data)
                        print_cluster(fh_all, cluster)
                        if is_mixed_cluster(cluster):
                            print_cluster(fh_filtered, cluster)
                        find_nearest_neighbors(cluster, cluster_id, fh_nearest)
                        output_go_pairs(cluster, cluster_id, basename, fh_pairs, go_data, ip_go_data)
            except IOError as e:
                print(f"Error writing output files for {basename}: {e}", file=sys.stderr)
                continue
            print(f"Finished processing {basename} in {class_dir.name}.")


if __name__ == "__main__":
    main()
