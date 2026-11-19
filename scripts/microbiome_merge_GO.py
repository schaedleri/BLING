#!/usr/bin/env python3
"""Aggregate GO and InterPro annotations, expand GO terms to their ancestors and write three annotation tiers.

  ncbi_go_gene.tsv / ncbi_background_list.txt      GenBank annotations only
  direct_go_gene.tsv / direct_background_list.txt  GenBank + InterProScan + eggNOG (no propagation)
  fixed_go_gene.tsv / fixed_background_list.txt    direct evidence + fixed-threshold propagation
  fixed_ipr_gene.tsv                               InterPro entries (used by the direct and fixed tiers)

Every pair carries an Evidence_Group_ID: direct sources give one group per gene and term, propagated pairs keep the
group read from nonGO_GO_pairs*.tsv, and ancestors inherit the group of the term they were expanded from."""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

PARENT_RELATIONSHIPS = {'is_a', 'part_of', 'regulates', 'positively_regulates', 'negatively_regulates',
                        'occurs_in', 'enables', 'contributes_to'}
GO_OBO = ('data', 'go.obo')
INTERPRO_GLOB = '*/interpro/*.tsv'
INTERPRO_IPR_COLUMN = 11
# source -> (glob pattern, in NCBI tier, in direct tier)
SOURCES = {
    'gbff': ('*/go_annotations/*.tsv', True, True),
    'interpro': ('*/interval_test/*/IP_GO_pairs.tsv', False, True),
    'propagated': ('*/interval_test/*/nonGO_GO_pairs*.tsv', False, False),
    'eggnog': ('*/eggnog/*_emapper_results.tsv', False, True),
}
PROPAGATED_COLUMNS = ('nongo_protein_id', 'go_id', 'evidence_group_id')
EGGNOG_COLUMNS = ('protein_id', 'gos')
GO_HEADER = "GO_ID\tGene_ID\tEvidence_Group_ID\n"
IPR_HEADER = "IPR_ID\tGene_ID\n"


def load_go_graph(obo_path):
    """Returns (child -> parents, obsolete term -> replacement/suggested terms)."""
    graph, obsolete_map = defaultdict(list), defaultdict(list)
    try:
        f = open(obo_path, encoding='utf-8')
    except FileNotFoundError:
        sys.exit(f"Error: GO OBO file not found: {obo_path}")
    with f:
        current, obsolete = None, False
        for line in (l.strip() for l in f):
            if not line:
                continue
            if line == '[Term]':
                current, obsolete = None, False
            elif line.startswith('id: '):
                current = line.split(' ')[1]
            elif line.startswith('is_obsolete: true'):
                obsolete = True
            elif current and not obsolete:
                parent = None
                if line.startswith('relationship: '):
                    parts = line.split(' ')
                    if parts[1] in PARENT_RELATIONSHIPS:
                        parent = parts[2]
                elif ':' in line:
                    tag, value = line.split(':', 1)
                    if tag in PARENT_RELATIONSHIPS:
                        parent = value.strip().split(' ')[0]
                if parent:
                    graph[current].append(parent)
            elif current and (line.startswith('replaced_by: ') or line.startswith('consider: ')):
                obsolete_map[current].append(line.split(' ')[1].split('!')[0].strip())
    return graph, obsolete_map


def get_all_ancestors(go_id, graph, memo):
    if go_id not in memo:
        parents = graph.get(go_id, [])
        ancestors = set(parents)
        for parent in parents:
            ancestors |= get_all_ancestors(parent, graph, memo)
        memo[go_id] = ancestors
    return memo[go_id]


def read_rows(path):
    with open(path, encoding='utf-8') as f:
        for line in f:
            yield line.strip().split('\t')


def pairs_from_file(source, path):
    """Yields (gene, go, evidence_group_id) from one input file."""
    table = read_rows(path)
    if source in ('propagated', 'eggnog'):
        wanted = PROPAGATED_COLUMNS if source == 'propagated' else EGGNOG_COLUMNS
        header = [c.lower() for c in next(table, [''])]
        if not all(c in header for c in wanted):
            print(f"Warning: expected columns {wanted} not found in {path}; skipped.", file=sys.stderr)
            return
        idx = [header.index(c) for c in wanted]
    for cols in table:
        if source == 'gbff' and len(cols) >= 3 and cols[2].startswith('GO:'):
            yield cols[0], cols[2], f"direct_gbff_{cols[0]}_{cols[2]}"
        elif source == 'interpro' and len(cols) >= 2 and cols[1].startswith('GO:'):
            yield cols[0], cols[1], f"direct_interpro_{cols[0]}_{cols[1]}"
        elif source == 'propagated' and len(cols) > max(idx) and cols[idx[1]].startswith('GO:'):
            yield cols[idx[0]], cols[idx[1]], cols[idx[2]]
        elif source == 'eggnog' and len(cols) > max(idx) and cols[idx[0]] and cols[idx[1]] != '-':
            for go in (g.strip() for g in cols[idx[1]].split(',')):
                if go.startswith('GO:'):
                    yield cols[idx[0]], go, f"direct_eggnog_{cols[idx[0]]}_{go}"


def collect_go_triples(sample_dir, graph, obsolete_map):
    """Returns the (GO_ID, Gene_ID, Evidence_Group_ID) sets of the NCBI, direct and full tiers."""
    tiers = {'ncbi': set(), 'direct': set(), 'all': set()}
    memo = {}
    for source, (pattern, in_ncbi, in_direct) in SOURCES.items():
        for path in sample_dir.glob(pattern):
            try:
                for gene, go, group in pairs_from_file(source, path):
                    for term in obsolete_map.get(go, {go}):
                        expanded = {(t, gene, group) for t in {term} | get_all_ancestors(term, graph, memo)}
                        tiers['all'] |= expanded
                        if in_direct:
                            tiers['direct'] |= expanded
                        if in_ncbi:
                            tiers['ncbi'] |= expanded
            except Exception as e:
                print(f"Warning: could not process {path}: {e}", file=sys.stderr)
    return tiers


def collect_ipr_pairs(sample_dir):
    pairs = set()
    for path in sample_dir.glob(INTERPRO_GLOB):
        try:
            for cols in read_rows(path):
                if len(cols) > INTERPRO_IPR_COLUMN and cols[INTERPRO_IPR_COLUMN] != '-':
                    pairs.update((ipr, cols[0]) for ipr in re.findall(r'IPR\d+', cols[INTERPRO_IPR_COLUMN]))
        except Exception as e:
            print(f"Warning: could not process {path}: {e}", file=sys.stderr)
    return pairs


def write_lines(path, header, lines):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(header + '\n'.join(sorted(lines)) + '\n')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--base-dir', required=True, help="Sample directory containing the taxon directories.")
    args = parser.parse_args()

    sample_dir = Path(args.base_dir).resolve()
    obo_path = sample_dir.parent.joinpath(*GO_OBO)
    print(f"Parsing GO hierarchy from {obo_path}...")
    graph, obsolete_map = load_go_graph(obo_path)
    print(f"{len(graph)} GO terms with parents; {len(obsolete_map)} obsolete terms.")

    tiers = collect_go_triples(sample_dir, graph, obsolete_map)
    ipr_pairs = collect_ipr_pairs(sample_dir)
    ipr_genes = {gene for _, gene in ipr_pairs}

    written = [write_lines(sample_dir / 'fixed_ipr_gene.tsv', IPR_HEADER, (f"{i}\t{g}" for i, g in ipr_pairs))]
    for tier, prefix, with_ipr in (('ncbi', 'ncbi', False), ('direct', 'direct', True), ('all', 'fixed', True)):
        genes = {gene for _, gene, _ in tiers[tier]} | (ipr_genes if with_ipr else set())
        written.append(write_lines(sample_dir / f"{prefix}_go_gene.tsv", GO_HEADER,
                                   (f"{go}\t{gene}\t{g}" for go, gene, g in tiers[tier])))
        written.append(write_lines(sample_dir / f"{prefix}_background_list.txt", '', genes))
    print("Written:\n" + '\n'.join(f" - {p}" for p in written))


if __name__ == '__main__':
    main()
