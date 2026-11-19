#!/usr/bin/env python3
"""GMM-based operon-like cluster prediction and BP GO propagation for every genome of every taxon directory."""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from warnings import filterwarnings

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from sklearn.mixture import GaussianMixture

filterwarnings('ignore', category=UserWarning, module='sklearn')

TAXON_DIR_PATTERN = r'^(Class|Family|Genus|Order|Species)_'
GO_NAMESPACES = ('biological_process', 'molecular_function', 'cellular_component')
IP_GO_NAMESPACES = ('molecular_function', 'biological_process')
WILDCARD_COGS = {'S', 'R'}
ALLOWED_MISMATCH_PAIRS = {frozenset(p) for p in (
    'KE', 'KG', 'KF', 'KH', 'KI', 'KP', 'KC', 'TM', 'TN', 'TU', 'TD', 'TV', 'GC', 'EH',
    'PC', 'UM', 'JO', 'MO', 'LD', 'MD', 'CP', 'FI')}
CLUSTER_COLUMNS = ['Accession', 'Organism', 'Begin', 'End', 'Strand', 'Product', 'Gene', 'Locus_Tag',
                   'Protein_Length', 'Protein_ID', 'status', 'Predicted_Name', 'COG_Category', 'Description']
SUMMARY_COLUMNS = ['Genome', 'Threshold_bp', 'Optimal_K', 'N_Intergenic_Distances', 'GMM_Suggestion_bp',
                   'Threshold_Basis']
K_SELECTION_N_INIT = 5
FINAL_FIT_N_INIT = 10
RANDOM_STATE = 0
HISTOGRAM_BINS = 80
PLOT_DPI = 150


def load_go_classifications(obo_file):
    classifications = {}
    if not obo_file.is_file():
        return classifications
    try:
        blocks = obo_file.read_text(encoding='utf-8').split('[Term]\n')[1:]
    except OSError:
        print(f"Warning: could not read {obo_file}", file=sys.stderr)
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
            go_type = 'biological_process' if go_type == 'process' else go_type
            if go_type in GO_NAMESPACES:
                go_data[accession].append({'id': go_id, 'type': go_type})
    return go_data


def load_eggnog_annotations(path, go_classifications):
    required = ['Protein_ID', 'GOs', 'COG_category', 'Preferred_name', 'Description']
    data = {}
    if not path.is_file():
        return data
    try:
        with open(path, encoding='utf-8') as f:
            col = {name: i for i, name in enumerate(f.readline().strip().split('\t'))}
            if not all(c in col for c in required):
                return data

            def field(fields, name):
                i = col[name]
                return fields[i] if i < len(fields) and fields[i] != '-' else ''

            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                cog, gos = field(fields, 'COG_category'), field(fields, 'GOs')
                data[fields[col['Protein_ID']]] = {
                    'cog_categories': list(cog),
                    'go_terms': [{'id': g, 'type': go_classifications.get(g, 'unknown')}
                                 for g in gos.split(',') if gos],
                    'predicted_name': field(fields, 'Preferred_name'),
                    'description': field(fields, 'Description'),
                }
    except Exception as e:
        print(f"Error reading eggNOG file {path.name}: {e}", file=sys.stderr)
    return data


def read_interpro_rows(path):
    if not path.is_file():
        return
    with open(path, encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 14:
                yield cols[0], cols[12], re.findall(r'(GO:\d+)', cols[13])


def read_ip_data(path, go_classifications):
    return {protein for protein, pathway, gos in read_interpro_rows(path)
            if pathway != '-' or any(go_classifications.get(g, 'unknown') != 'cellular_component' for g in gos)}


def read_ip_go_data(path, go_classifications):
    data = defaultdict(list)
    for protein, _, gos in read_interpro_rows(path):
        for go_id in gos:
            namespace = go_classifications.get(go_id, '')
            if namespace in IP_GO_NAMESPACES:
                data[protein].append({'id': go_id, 'type': namespace})
    return data


def read_genomic_data(path):
    data = defaultdict(list)
    if not path.is_file():
        return data
    try:
        with open(path, encoding='utf-8', errors='ignore') as f:
            header = {name: i for i, name in enumerate(next(f).strip().split('\t'))}
            for line in f:
                values = line.strip().split('\t')
                if len(values) != len(header):
                    continue
                gene = {k: values[i] for k, i in header.items()}
                try:
                    gene['begin'], gene['end'] = int(gene['Begin']), int(gene['End'])
                except (ValueError, KeyError):
                    continue
                data[gene['Accession']].append(gene)
    except Exception as e:
        print(f"Error reading {path}: {e}", file=sys.stderr)
    return data


def ordered_genes(genes):
    """Genes of one replicon by start position, all strands together. Adjacent pairs are consecutive genes on
    the same strand, so a gene on the opposite strand in between breaks the adjacency."""
    return sorted(genes, key=lambda g: (g['begin'], g['end']))


def get_distances_and_offset(genomic_data, outlier_threshold, min_distances):
    distances = []
    for genes in genomic_data.values():
        genes = ordered_genes(genes)
        for prev, cur in zip(genes, genes[1:]):
            if prev['Protein_ID'] != cur['Protein_ID'] and prev['Strand'] == cur['Strand']:
                distances.append(cur['begin'] - prev['end'])
    distances = np.array(distances)
    distances = distances[distances < outlier_threshold]
    if len(distances) < min_distances:
        return None, None, None
    min_gap = distances.min()
    offset = abs(min_gap) + 1 if min_gap <= 0 else 1
    return distances, np.log(distances + offset).reshape(-1, 1), offset


def find_optimal_k(log_distances, max_k):
    ks = range(1, max_k + 1)
    fits = [GaussianMixture(n_components=k, random_state=RANDOM_STATE, n_init=K_SELECTION_N_INIT).fit(log_distances)
            for k in ks]
    bics = [g.bic(log_distances) for g in fits]
    aics = [g.aic(log_distances) for g in fits]
    return ks[int(np.argmin(bics))], ks, bics, aics


def plot_k_selection(ks, bics, aics, optimal_k, path):
    plt.figure(figsize=(10, 6))
    plt.plot(ks, bics, marker='o', label='BIC')
    plt.plot(ks, aics, marker='o', linestyle='--', label='AIC', alpha=0.7)
    plt.axvline(optimal_k, color='r', linestyle=':', label=f'Optimal K = {optimal_k}')
    plt.title(f'GMM Component Selection ({path.stem})')
    plt.xlabel('Number of Components (K)')
    plt.ylabel('Score')
    plt.xticks(ks)
    plt.legend()
    plt.grid(True, alpha=0.6)
    plt.savefig(path, dpi=PLOT_DPI)
    plt.close()


def select_operon_component(gmm):
    """Among components whose mean does not exceed the median of the component means, the one with the largest weight."""
    means, weights = gmm.means_.flatten(), gmm.weights_.flatten()
    if len(means) == 1:
        return 0
    candidates = np.where(means <= np.median(means))[0]
    return candidates[np.argmax(weights[candidates])]


def log_transform(x, offset):
    with np.errstate(invalid='ignore'):
        return np.nan_to_num(np.log(x + offset), nan=0.0, posinf=0.0, neginf=0.0)


def find_threshold_and_plot(distances, log_distances, offset, optimal_k, plot_path, outlier_threshold,
                            min_threshold, scan_start):
    gmm = GaussianMixture(n_components=optimal_k, random_state=RANDOM_STATE, n_init=FINAL_FIT_N_INIT).fit(log_distances)
    operon_idx = select_operon_component(gmm)
    print(f"     -> Operon component: index {operon_idx} (log-mean {gmm.means_[operon_idx][0]:.2f}, "
          f"weight {gmm.weights_[operon_idx]:.2%})")

    scan = np.arange(scan_start, outlier_threshold)
    in_operon = scan[gmm.predict(log_transform(scan, offset).reshape(-1, 1)) == operon_idx]
    gmm_threshold = in_operon.max() if len(in_operon) else -1
    threshold = max(gmm_threshold, min_threshold)

    x_min = distances.min() + scan_start
    x = np.linspace(x_min, outlier_threshold, 1000).reshape(-1, 1)
    log_x = log_transform(x, offset)
    plt.figure(figsize=(12, 7))
    plt.hist(distances, bins=HISTOGRAM_BINS, density=True, alpha=0.6, label=f'Distances (< {outlier_threshold} bp)')
    for i, (w, m, c) in enumerate(zip(gmm.weights_, gmm.means_, gmm.covariances_)):
        pdf = np.nan_to_num(w * norm.pdf(log_x, m[0], np.sqrt(c[0][0])) / (x + offset),
                            nan=0.0, posinf=0.0, neginf=0.0)
        plt.plot(x, pdf, lw=2, label='Operon Comp.' if i == operon_idx else f'Comp. {i + 1}')
    plt.axvline(threshold, color='r', linestyle='--', label=f'Final Threshold: {threshold} bp')
    if threshold != gmm_threshold and gmm_threshold >= scan_start:
        plt.axvline(gmm_threshold, color='gray', linestyle=':', label=f'GMM Suggestion: {gmm_threshold} bp')
    plt.title(f'Final GMM Fit (K={optimal_k}) for {plot_path.stem}')
    plt.xlabel('Distance (bp)')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, alpha=0.6)
    plt.xlim(x_min, outlier_threshold)
    plt.ylim(bottom=0)
    plt.savefig(plot_path, dpi=PLOT_DPI)
    plt.close()
    return threshold, gmm_threshold


def cogs_compatible(cogs1, cogs2):
    return bool(cogs1 & cogs2 or cogs1 & WILDCARD_COGS or cogs2 & WILDCARD_COGS
                or any(frozenset({a, b}) in ALLOWED_MISMATCH_PAIRS for a in cogs1 for b in cogs2))


def generate_clusters(genomic_data, threshold, eggnog_data):
    def cogs(gene):
        return set(eggnog_data.get(gene['Protein_ID'], {}).get('cog_categories', []))

    clusters = []
    for genes in genomic_data.values():
        genes = ordered_genes(genes)
        if not genes:
            continue
        current = [genes[0]]
        for prev, cur in zip(genes, genes[1:]):
            if prev['Protein_ID'] == cur['Protein_ID']:
                continue
            gap = cur['begin'] - prev['end']
            joined = (prev['Strand'] == cur['Strand'] and gap <= threshold
                      and (gap <= 0 or cogs_compatible(cogs(prev), cogs(cur))))
            if joined:
                current.append(cur)
                continue
            if len(current) >= 2:
                clusters.append(current)
            current = [cur]
        if len(current) >= 2:
            clusters.append(current)
    return clusters


def annotate_cluster(cluster, go_data, eggnog_data, ip_data):
    for gene in cluster:
        pid = gene['Protein_ID']
        info = eggnog_data.get(pid, {})
        cogs = info.get('cog_categories', ['S'])
        gene.update(cog_categories=cogs, predicted_name=info.get('predicted_name', ''),
                    description=info.get('description', ''))
        if go_data.get(pid):
            gene['status'] = 'GO'
        elif pid in ip_data:
            gene['status'] = 'IP'
        elif cogs and not all(c == 'S' for c in cogs):
            gene['status'] = 'COG'
        else:
            gene['status'] = 'nonGO'


def is_valid_cog_operon(cluster):
    annotated = [g for g in cluster if g['status'] in ('GO', 'COG', 'IP')]
    if len(annotated) < 2:
        return True
    distinct = sorted(set.union(*[set(g['cog_categories']) for g in annotated]) - WILDCARD_COGS)
    return len(distinct) <= 1 or (len(distinct) == 2 and frozenset(distinct) in ALLOWED_MISMATCH_PAIRS)


def is_mixed_cluster(cluster):
    statuses = {g['status'] for g in cluster}
    return 'nonGO' in statuses and bool(statuses & {'GO', 'COG', 'IP'})


def print_cluster(handle, cluster):
    for gene in cluster:
        gene['COG_Category'] = ''.join(gene.get('cog_categories', []))
        handle.write('\t'.join(str(gene.get(h.lower(), gene.get(h, ''))) for h in CLUSTER_COLUMNS) + '\n')
    handle.write('\n')


def output_go_pairs(cluster, cluster_id, genome_id, handle, go_data):
    """Every BP term observed in the cluster is assigned to all its genes; one Evidence_Group_ID per cluster,
    prefixed with the genome id because cluster ids restart in every genome."""
    bp_ids = {go['id'] for gene in cluster for go in go_data.get(gene['Protein_ID'], [])
              if go['type'] == 'biological_process'}
    group_id = f"{genome_id}_cluster{cluster_id}"
    for gene in cluster:
        for go_id in sorted(bp_ids):
            handle.write(f"{gene['Protein_ID']}\t{go_id}\tprocess\t{cluster_id}\t{group_id}\n")


def threshold_basis(distances, optimal_k, gmm_suggestion, min_threshold):
    if distances is None:
        return 'default_too_few_distances'
    if optimal_k == 1:
        return 'single_component_K1'
    return 'floor_applied' if gmm_suggestion < min_threshold else 'GMM'


def merge_go_sources(*sources):
    merged = defaultdict(list)
    for source in sources:
        for pid, terms in source.items():
            seen = {g['id'] for g in merged[pid]}
            merged[pid].extend(g for g in terms if g['id'] not in seen)
    return merged


def process_genome(genomic_file, taxon_dir, out_dir, go_classifications, args):
    name = genomic_file.stem.replace('_genomic', '')
    print(f"\n  -> Genome: {name}")
    eggnog_data = load_eggnog_annotations(taxon_dir / 'eggnog' / f"{name}_genomic_emapper_results.tsv",
                                          go_classifications)
    if not eggnog_data:
        print(f"     -> Warning: no eggNOG data for {name}; skipped.")
        return None
    genomic_data = read_genomic_data(genomic_file)
    distances, log_distances, offset = get_distances_and_offset(genomic_data, args.outlier_threshold,
                                                                args.min_distances)
    optimal_k, gmm_suggestion = None, None
    if distances is None:
        print(f"     -> Too few intergenic distances for {name}; using the minimum threshold.")
        threshold = args.min_threshold
    else:
        optimal_k, ks, bics, aics = find_optimal_k(log_distances, args.max_k)
        print(f"     -> Optimal K: {optimal_k} (minimum BIC)")
        plot_k_selection(ks, bics, aics, optimal_k, out_dir / f"{name}_k_selection.png")
        threshold, gmm_suggestion = find_threshold_and_plot(
            distances, log_distances, offset, optimal_k, out_dir / f"{name}_gmm_fit.png",
            args.outlier_threshold, args.min_threshold, args.scan_start)
    record = {'genome': name, 'threshold_bp': threshold, 'optimal_k': optimal_k,
              'n_distances': 0 if distances is None else len(distances),
              'gmm_suggestion_bp': 'NA' if gmm_suggestion is None else gmm_suggestion,
              'basis': threshold_basis(distances, optimal_k, gmm_suggestion, args.min_threshold)}
    if optimal_k == 1:
        print("     -> K = 1: threshold not estimable; genome excluded from clustering and propagation.")
        return record

    print(f"     -> Threshold used for clustering: {threshold:.1f} bp")
    clusters = generate_clusters(genomic_data, threshold, eggnog_data)
    if not clusters:
        print("        -> No clusters predicted.")
        return record
    print(f"     -> {len(clusters)} clusters predicted.")

    go_data = merge_go_sources(
        read_go_data(taxon_dir / 'go_annotations' / f"{name}_go_annotations.tsv"),
        read_ip_go_data(taxon_dir / 'interpro' / f"{name}.tsv", go_classifications),
        {pid: d['go_terms'] for pid, d in eggnog_data.items()})
    ip_data = read_ip_data(taxon_dir / 'interpro' / f"{name}.tsv", go_classifications)

    genome_dir = out_dir / name
    genome_dir.mkdir(parents=True, exist_ok=True)
    header = '\t'.join(CLUSTER_COLUMNS) + '\n'
    with open(genome_dir / 'all_predicted_operons.tsv', 'w') as fh_all, \
            open(genome_dir / 'filtered_cog_coherent_operons.tsv', 'w') as fh_cog, \
            open(genome_dir / 'filtered_mixed_operons.tsv', 'w') as fh_mixed, \
            open(genome_dir / 'propagated_GO_pairs_BP.tsv', 'w') as fh_pairs:
        for handle in (fh_all, fh_cog, fh_mixed):
            handle.write(header)
        fh_pairs.write('Protein_ID\tGO_ID\tGO_Type\tCluster_ID\tEvidence_Group_ID\n')
        coherent = []
        for cluster in clusters:
            annotate_cluster(cluster, go_data, eggnog_data, ip_data)
            print_cluster(fh_all, cluster)
            if is_mixed_cluster(cluster):
                print_cluster(fh_mixed, cluster)
            if is_valid_cog_operon(cluster):
                coherent.append(cluster)
        for cluster_id, cluster in enumerate(coherent, 1):
            print_cluster(fh_cog, cluster)
            output_go_pairs(cluster, cluster_id, name, fh_pairs, go_data)
    return record


def write_threshold_summary(records, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(SUMMARY_COLUMNS) + '\n')
        for r in sorted(records, key=lambda r: r['genome']):
            f.write(f"{r['genome']}\t{r['threshold_bp']}\t{r['optimal_k']}\t{r['n_distances']}\t"
                    f"{r['gmm_suggestion_bp']}\t{r['basis']}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-dir', required=True, help="Sample directory containing the taxon directories.")
    parser.add_argument('--max-k', type=int, default=8)
    parser.add_argument('--outlier-threshold', type=int, default=1000,
                        help="Intergenic distances >= this value (bp) are excluded from the GMM fit.")
    parser.add_argument('--min-threshold', type=int, default=20, help="Lower bound of the distance threshold (bp).")
    parser.add_argument('--min-distances', type=int, default=20,
                        help="Minimum number of intergenic distances required to fit a GMM.")
    parser.add_argument('--scan-start', type=int, default=-50,
                        help="Smallest distance (bp) scanned when locating the threshold.")
    args = parser.parse_args()

    sample_dir = Path(args.base_dir).resolve()
    go_classifications = load_go_classifications(sample_dir.parent / 'data' / 'go.obo')
    taxon_dirs = [d for d in sample_dir.iterdir() if d.is_dir() and re.match(TAXON_DIR_PATTERN, d.name)]
    if not taxon_dirs:
        sys.exit(f"No taxon directories found in {sample_dir}.")

    for taxon_dir in taxon_dirs:
        print(f"\n--- {taxon_dir.name} ---")
        out_dir = taxon_dir / 'operon_analysis'
        out_dir.mkdir(parents=True, exist_ok=True)
        genomic_dir = taxon_dir / 'genomic_annotations'
        if not genomic_dir.exists():
            print(f"  -> Warning: {genomic_dir} not found; skipped.")
            continue
        records = [r for f in genomic_dir.glob('*_genomic.tsv')
                   if (r := process_genome(f, taxon_dir, out_dir, go_classifications, args))]
        if records:
            summary = out_dir / 'gmm_threshold_summary.tsv'
            write_threshold_summary(records, summary)
            print(f"  -> GMM threshold summary: {summary}")
    print("\nDone.")


if __name__ == '__main__':
    main()
