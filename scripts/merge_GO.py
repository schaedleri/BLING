#!/usr/bin/env python3
"""Aggregate GO and InterPro annotations of all taxa, expand GO terms to their ancestors and write
<prefix>_go_gene.tsv, <prefix>_ipr_gene.tsv and <prefix>_background_list.txt.

Every GO-gene pair carries an Evidence_Group_ID. Direct sources (GenBank, InterProScan, eggNOG) give one group per
gene and term; propagated pairs keep the group read from propagated_GO_pairs_BP.tsv; ancestors inherit the group of
the term they were expanded from."""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

PARENT_RELATIONSHIPS = {'is_a', 'part_of', 'regulates', 'positively_regulates', 'negatively_regulates',
                        'occurs_in', 'enables', 'contributes_to'}
GO_OBO = ('data', 'go.obo')
GO_ANNOTATION_GLOB = '*/go_annotations/*.tsv'
INTERPRO_GLOB = '*/interpro/*.tsv'
PROPAGATED_GLOB = '*/operon_analysis/*/propagated_GO_pairs_BP.tsv'
EGGNOG_GLOB = '*/eggnog/*_emapper_results.tsv'
PROPAGATED_COLUMNS = ('Protein_ID', 'GO_ID', 'Evidence_Group_ID')
INTERPRO_MIN_COLUMNS = 14
INTERPRO_GO_COLUMN = 13
INTERPRO_IPR_COLUMN = 11


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


def rows(path):
    with open(path, encoding='utf-8') as f:
        for line in f:
            yield line.strip().split('\t')


def collect_direct_and_propagated(sample_dir):
    triples = set()
    for path in sample_dir.glob(GO_ANNOTATION_GLOB):
        for cols in list(rows(path))[1:]:
            if len(cols) >= 4 and cols[2].startswith('GO:'):
                triples.add((cols[0], cols[2], f"direct_gbff_{cols[0]}_{cols[2]}"))

    for path in sample_dir.glob(INTERPRO_GLOB):
        for cols in rows(path):
            if len(cols) >= INTERPRO_MIN_COLUMNS:
                for go in re.findall(r'GO:\d+', cols[INTERPRO_GO_COLUMN]):
                    triples.add((cols[0], go, f"direct_interpro_{cols[0]}_{go}"))

    for path in sample_dir.glob(PROPAGATED_GLOB):
        table = rows(path)
        header = next(table, [''])
        if not all(c in header for c in PROPAGATED_COLUMNS):
            print(f"Warning: expected columns not found in {path}; skipped.", file=sys.stderr)
            continue
        gene_i, go_i, group_i = (header.index(c) for c in PROPAGATED_COLUMNS)
        for cols in table:
            if len(cols) > max(gene_i, go_i, group_i) and cols[go_i].startswith('GO:'):
                triples.add((cols[gene_i], cols[go_i], cols[group_i]))

    for path in sample_dir.glob(EGGNOG_GLOB):
        try:
            table = rows(path)
            header = {name: i for i, name in enumerate(next(table, ['']))}
            if 'Protein_ID' not in header or 'GOs' not in header:
                continue
            gene_i, go_i = header['Protein_ID'], header['GOs']
            for cols in table:
                if len(cols) > max(gene_i, go_i) and cols[go_i] != '-':
                    for go in (g.strip() for g in cols[go_i].split(',')):
                        if go:
                            triples.add((cols[gene_i], go, f"direct_eggnog_{cols[gene_i]}_{go}"))
        except Exception as e:
            print(f"Warning: could not process {path}: {e}", file=sys.stderr)
    return triples


def expand_to_ancestors(triples, graph, obsolete_map):
    expanded, memo = set(), {}
    for gene, go, group in triples:
        for term in obsolete_map.get(go, {go}):
            expanded.add((term, gene, group))
            expanded.update((a, gene, group) for a in get_all_ancestors(term, graph, memo))
    return expanded


def collect_ipr_pairs(sample_dir):
    pairs = set()
    for path in sample_dir.glob(INTERPRO_GLOB):
        try:
            for cols in rows(path):
                if len(cols) > INTERPRO_IPR_COLUMN and cols[INTERPRO_IPR_COLUMN] != '-':
                    pairs.update((ipr, cols[0]) for ipr in re.findall(r'IPR\d+', cols[INTERPRO_IPR_COLUMN]))
        except Exception as e:
            print(f"Warning: could not process {path}: {e}", file=sys.stderr)
    return pairs


def write_lines(path, header, lines):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(header)
        f.write('\n'.join(sorted(lines)) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--base-dir', required=True, help="Sample directory containing the taxon directories.")
    parser.add_argument('--prefix', default='GMM', help="Prefix of the output files.")
    args = parser.parse_args()

    sample_dir = Path(args.base_dir).resolve()
    obo_path = sample_dir.parent.joinpath(*GO_OBO)
    print(f"Parsing GO hierarchy from {obo_path}...")
    graph, obsolete_map = load_go_graph(obo_path)
    print(f"{len(graph)} GO terms with parents; {len(obsolete_map)} obsolete terms.")

    triples = collect_direct_and_propagated(sample_dir)
    print(f"Expanding {len(triples)} GO triples to ancestors...")
    expanded = expand_to_ancestors(triples, graph, obsolete_map)
    ipr_pairs = collect_ipr_pairs(sample_dir)

    go_path = sample_dir / f"{args.prefix}_go_gene.tsv"
    ipr_path = sample_dir / f"{args.prefix}_ipr_gene.tsv"
    background_path = sample_dir / f"{args.prefix}_background_list.txt"
    write_lines(go_path, "GO_ID\tGene_ID\tEvidence_Group_ID\n", (f"{go}\t{gene}\t{g}" for go, gene, g in expanded))
    write_lines(ipr_path, "IPR_ID\tGene_ID\n", (f"{ipr}\t{gene}" for ipr, gene in ipr_pairs))
    write_lines(background_path, '', {gene for _, gene, _ in expanded} | {gene for _, gene in ipr_pairs})
    print(f"Written:\n - {go_path}\n - {ipr_path}\n - {background_path}")


if __name__ == '__main__':
    main()
