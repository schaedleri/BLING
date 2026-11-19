#!/usr/bin/env python3
"""Fixed-threshold gene-neighborhood clusters and BP GO propagation.

Adjacent same-strand genes separated by <= --interval bp form a cluster (at least --min-cluster-size genes). Within a
cluster, Biological Process terms of annotated genes (GO or InterPro) are propagated to unannotated ('nonGO') genes
within --max-dist bp. One Evidence_Group_ID per cluster, prefixed with the genome id."""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

TAXON_DIR_PATTERN = r'^(Class|Family|Genus|Order|Species)_'
GO_OBO = ('data', 'go.obo')
GO_NAMESPACES = ('biological_process', 'molecular_function', 'cellular_component')
BP_TYPES = ('biological_process', 'process')
INTERPRO_MIN_COLUMNS = 14
CLUSTER_COLUMNS = ['Accession', 'Organism', 'Begin', 'End', 'Strand', 'Product', 'Gene', 'Locus_Tag',
                   'Protein_Length', 'Protein_ID', 'status']
OUTPUT_DIR = 'interval_test'
OUTPUT_FILES = {'all': 'interval_afterpsi100.tsv', 'mixed': 'filtered_interval_afterpsi100.tsv',
                'nearest': 'nearest_results100.tsv', 'pairs': 'nonGO_GO_pairs100.tsv', 'ip_go': 'IP_GO_pairs.tsv'}


def load_go_classifications(obo_path):
    classifications = {}
    if not obo_path.is_file():
        return classifications
    try:
        blocks = obo_path.read_text(encoding='utf-8').split('[Term]\n')[1:]
    except OSError:
        print(f"Warning: could not read {obo_path}", file=sys.stderr)
        return classifications
    for block in blocks:
        go_id = re.search(r'^id: (GO:\d+)', block, re.MULTILINE)
        if go_id:
            namespace = re.search(r'^namespace: (\S+)', block, re.MULTILINE)
            classifications[go_id.group(1)] = namespace.group(1) if namespace else 'unknown'
    return classifications


def read_go_data(path):
    go_data = defaultdict(list)
    if not path.is_file():
        return go_data
    with open(path, encoding='utf-8') as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            try:
                accession, _, go_id, go_type = line.strip().split('\t')
            except ValueError:
                continue
            if go_type in GO_NAMESPACES:
                go_data[accession].append({'id': go_id, 'type': go_type})
    return go_data


def read_interpro_rows(path):
    if not path.is_file():
        return
    with open(path, encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= INTERPRO_MIN_COLUMNS:
                yield cols[0], cols[12], re.findall(r'(GO:\d+)', cols[13])


def read_ip_data(path, go_classifications):
    """Proteins with an InterPro hit, excluding those whose GO terms are all cellular component."""
    proteins = set()
    for protein, pathway, gos in read_interpro_rows(path):
        namespaces = {go_classifications.get(g, 'unknown') for g in gos}
        if namespaces != {'cellular_component'} and pathway != '-':
            proteins.add(protein)
    return proteins


def read_ip_go_data(path, go_classifications):
    data = defaultdict(list)
    for protein, _, gos in read_interpro_rows(path):
        for go_id in gos:
            namespace = go_classifications.get(go_id, '')
            if namespace == 'molecular_function':
                data[protein].append({'id': go_id, 'type': 'function'})
            elif namespace == 'biological_process':
                data[protein].append({'id': go_id, 'type': 'process'})
    return data


def read_eggnog_data(path, go_classifications):
    data = defaultdict(list)
    if not path.is_file():
        return data
    try:
        with open(path, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            if 'Protein_ID' not in header or 'GOs' not in header:
                print(f"Warning: 'Protein_ID' or 'GOs' column not found in {path}", file=sys.stderr)
                return data
            protein_i, go_i = header.index('Protein_ID'), header.index('GOs')
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) <= max(protein_i, go_i) or cols[go_i] in ('', '-'):
                    continue
                for go_id in (g.strip() for g in cols[go_i].split(',')):
                    namespace = go_classifications.get(go_id, 'unknown')
                    if go_id and namespace in ('biological_process', 'molecular_function'):
                        data[cols[protein_i]].append({'id': go_id, 'type': namespace})
    except IOError as e:
        print(f"Warning: could not read {path}: {e}", file=sys.stderr)
    return data


def read_genomic_data(path):
    data = defaultdict(list)
    with open(path, encoding='utf-8') as f:
        header = next(f).strip().split('\t')
        for line in f:
            if not line.strip():
                continue
            gene = dict(zip(header, line.strip().split('\t')))
            try:
                gene['begin'], gene['end'] = int(gene['Begin']), int(gene['End'])
            except (ValueError, KeyError):
                continue
            data[gene['Accession']].append(gene)
    return data


def generate_clusters(genomic_data, interval, min_size):
    clusters = []
    for genes in genomic_data.values():
        genes = sorted(genes, key=lambda g: g['begin'])
        if not genes:
            continue
        current = [genes[0]]
        for prev, cur in zip(genes, genes[1:]):
            if cur['begin'] - prev['end'] <= interval and cur['Strand'] == prev['Strand']:
                current.append(cur)
                continue
            if len(current) >= min_size:
                clusters.append(current)
            current = [cur]
        if len(current) >= min_size:
            clusters.append(current)
    return clusters


def annotate_cluster(cluster, go_data, ip_data):
    for gene in cluster:
        pid = gene['Protein_ID']
        gene['status'] = 'GO' if pid in go_data else 'IP' if pid in ip_data else 'nonGO'


def is_mixed_cluster(cluster):
    statuses = {g['status'] for g in cluster}
    return 'nonGO' in statuses and bool(statuses & {'GO', 'IP'})


def print_cluster(handle, cluster):
    for gene in cluster:
        handle.write('\t'.join(str(gene.get(h, '')) for h in CLUSTER_COLUMNS) + '\n')
    handle.write('\n')


def neighbors_of_nongo(cluster, max_dist):
    """Yields (nonGO gene, annotated gene, distance) for annotated genes within max_dist bp."""
    for gene in cluster:
        if gene['status'] != 'nonGO':
            continue
        for other in cluster:
            if other is gene or other['status'] == 'nonGO':
                continue
            dist = max(gene['begin'] - other['end'], other['begin'] - gene['end'], 0)
            if dist <= max_dist:
                yield gene, other, dist


def write_nearest(cluster, cluster_id, handle, max_dist):
    for gene, other, dist in neighbors_of_nongo(cluster, max_dist):
        handle.write(f"{cluster_id}\t{gene['Protein_ID']}\tneighbor_{other['status']}\t{dist}\t{other['Protein_ID']}\n")


def write_go_pairs(cluster, cluster_id, genome_id, handle, go_data, ip_go_data, max_dist):
    group_id = f"{genome_id}_cluster{cluster_id}"
    sources = {'GO': go_data, 'IP': ip_go_data}
    for gene, other, _ in neighbors_of_nongo(cluster, max_dist):
        for go in sources[other['status']].get(other['Protein_ID'], []):
            if go.get('type') in BP_TYPES:
                handle.write(f"{gene['Protein_ID']}\t{go['id']}\tprocess\t{cluster_id}\t{group_id}\n")


def write_all_ip_go_pairs(ip_file, go_classifications, handle):
    for protein, _, gos in read_interpro_rows(ip_file):
        for go_id in gos:
            handle.write(f"{protein}\t{go_id}\t{go_classifications.get(go_id, 'unknown')}\n")


def process_genome(genomic_file, taxon_dir, go_classifications, args):
    name = genomic_file.stem.replace('_genomic', '')
    if not name:
        return
    print(f"Processing {name} in {taxon_dir.name}...")
    out_dir = taxon_dir / OUTPUT_DIR / name
    out_dir.mkdir(parents=True, exist_ok=True)

    ip_file = taxon_dir / 'interpro' / f"{name}.tsv"
    go_data = read_go_data(taxon_dir / 'go_annotations' / f"{name}_go_annotations.tsv")
    eggnog = read_eggnog_data(taxon_dir / 'eggnog' / f"{name}_genomic_emapper_results.tsv", go_classifications)
    for protein, terms in eggnog.items():
        known = {go['id'] for go in go_data.get(protein, [])}
        for term in terms:
            if term['id'] not in known:
                go_data[protein].append(term)
                known.add(term['id'])
    ip_data = read_ip_data(ip_file, go_classifications)
    ip_go_data = read_ip_go_data(ip_file, go_classifications)
    clusters = generate_clusters(read_genomic_data(genomic_file), args.interval, args.min_cluster_size)
    if not clusters:
        print(f"  -> No clusters for {name}; skipped.")
        return

    header = '\t'.join(CLUSTER_COLUMNS) + '\n'
    try:
        with open(out_dir / OUTPUT_FILES['all'], 'w') as fh_all, \
                open(out_dir / OUTPUT_FILES['mixed'], 'w') as fh_mixed, \
                open(out_dir / OUTPUT_FILES['nearest'], 'w') as fh_nearest, \
                open(out_dir / OUTPUT_FILES['pairs'], 'w') as fh_pairs, \
                open(out_dir / OUTPUT_FILES['ip_go'], 'w') as fh_ip_go:
            fh_ip_go.write("IP_Protein_ID\tGO_ID\tGO_Type\n")
            fh_all.write(header)
            fh_mixed.write(header)
            fh_nearest.write("Cluster_ID\tNonGO_Protein_ID\tNeighbor_Status\tDistance\tNeighbor_Protein_ID\n")
            fh_pairs.write("NonGO_Protein_ID\tGO_ID\tGO_Type\tCluster_ID\tEvidence_Group_ID\n")
            write_all_ip_go_pairs(ip_file, go_classifications, fh_ip_go)
            for cluster_id, cluster in enumerate(clusters, 1):
                annotate_cluster(cluster, go_data, ip_data)
                print_cluster(fh_all, cluster)
                if is_mixed_cluster(cluster):
                    print_cluster(fh_mixed, cluster)
                write_nearest(cluster, cluster_id, fh_nearest, args.max_dist)
                write_go_pairs(cluster, cluster_id, name, fh_pairs, go_data, ip_go_data, args.max_dist)
    except IOError as e:
        print(f"Error writing output for {name}: {e}", file=sys.stderr)
        return
    print(f"Finished {name} in {taxon_dir.name}.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-dir', required=True, help="Sample directory containing the taxon directories.")
    parser.add_argument('--interval', type=int, default=100, help="Maximum gap between adjacent genes of a cluster (bp).")
    parser.add_argument('--max-dist', type=int, default=300,
                        help="Maximum distance between an annotated and an unannotated gene for propagation (bp).")
    parser.add_argument('--min-cluster-size', type=int, default=3, help="Minimum number of genes in a cluster.")
    args, _ = parser.parse_known_args()

    sample_dir = Path(args.base_dir).resolve()
    go_classifications = load_go_classifications(sample_dir.parent.joinpath(*GO_OBO))
    taxon_dirs = [d for d in sample_dir.iterdir() if d.is_dir() and re.match(TAXON_DIR_PATTERN, d.name)]
    if not taxon_dirs:
        print(f"Warning: no taxon directories found in {sample_dir}", file=sys.stderr)
        return
    for taxon_dir in taxon_dirs:
        genomic_dir = taxon_dir / 'genomic_annotations'
        if genomic_dir.is_dir():
            for genomic_file in genomic_dir.glob('*.tsv'):
                process_genome(genomic_file, taxon_dir, go_classifications, args)


if __name__ == '__main__':
    main()
