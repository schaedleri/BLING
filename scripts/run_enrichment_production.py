#!/usr/bin/env python3
"""GO / InterPro enrichment with evidence units (a propagation cluster counts as one unit).

Per GO namespace: M = background units with an annotation, n = those carrying the term, N = units with a study
gene, k = study units carrying the term. p = hypergeom.sf(k - 1, M, n, N); fold = (k / N) / (n / M); BH over all
terms with k > 0 and min_term_size <= n <= max_term_fraction * M.
"""
import argparse
import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from goatools.obo_parser import GODag
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

CLUSTER_TAG = '_cluster'
ANNOTATION_FILES = {  # annotation file -> (background list, output suffix)
    'ncbi_go_gene.tsv': ('ncbi_background_list.txt', '_from_ncbi'),
    'direct_go_gene.tsv': ('direct_background_list.txt', '_from_direct'),
    'fixed_go_gene.tsv': ('fixed_background_list.txt', '_from_fixed'),
    'GMM_go_gene.tsv': ('GMM_background_list.txt', '_from_GMM'),
}
NAMESPACES = (('biological_process', 'BP'), ('molecular_function', 'MF'), ('cellular_component', 'CC'))
GO_OBO = ('data', 'go.obo')
IPR_ENTRY_LIST = ('data', 'entry.list.txt')
IPR_ANNOTATION_FILE = 'ipr_gene.tsv'
GENE_LIST_DIR = 'BLAST_specify'
OUTPUT_DIR = 'Enrichment_Analysis'
TAXON_PREFIXES = ('Genus_', 'Species_')
REQUIRED_GO_COLUMNS = {'GO_ID', 'Gene_ID', 'Evidence_Group_ID'}
BOOTSTRAP_CI = (2.5, 97.5)


def is_cluster_row(df):
    return df['Evidence_Group_ID'].astype(str).str.contains(CLUSTER_TAG, na=False)


def build_units(df, mode):
    """Unit id per row. 'exclusive': a gene with cluster rows is represented by its cluster (the smallest cluster
    id if identical proteins sit in clusters of several genomes), other genes by themselves. 'legacy': cluster rows
    belong to the cluster, other rows to the gene. Also returns the genes assigned to more than one cluster."""
    is_cl = is_cluster_row(df).values
    if mode == 'legacy':
        return np.where(is_cl, df['Evidence_Group_ID'].values, df['Gene'].values), set()
    clusters = df.loc[is_cl, ['Gene', 'Evidence_Group_ID']].drop_duplicates()
    n_clusters = clusters.groupby('Gene')['Evidence_Group_ID'].nunique()
    unit = df['Gene'].map(clusters.groupby('Gene')['Evidence_Group_ID'].min())
    return unit.where(unit.notna(), df['Gene']).values, set(n_clusters[n_clusters > 1].index)


def run_go_enrichment(study_genes, background_genes, annotations, go_dag, unit_mode, log, args):
    study, background = set(study_genes), set(background_genes)
    d = annotations[annotations['Gene'].isin(background)].copy()
    if d.empty:
        return pd.DataFrame()
    d['Unit'], multi = build_units(d, unit_mode)
    if multi:
        log(f"       [NOTE] {len(multi)} background genes carry more than one cluster id; the smallest was used. "
            f"Study genes affected: {len(multi & study)}.")
    n_bg_genes = d['Gene'].nunique()
    study_in = study & set(d['Gene'])
    M = d['Unit'].nunique()
    study_units = set(d.loc[d['Gene'].isin(study_in), 'Unit'])
    N = len(study_units)
    log(f"       - Background genes: {n_bg_genes}; study genes with annotation: {len(study_in)}")
    log(f"       - M = {M}; N = {N}  [unit mode: {unit_mode}]")
    if N == 0:
        return pd.DataFrame()

    unit_terms = d[['GO_ID', 'Unit']].drop_duplicates()
    n_by_term = unit_terms.groupby('GO_ID').size()
    max_size = int(M * args.max_term_fraction)
    valid = set(n_by_term[(n_by_term >= args.min_term_size) & (n_by_term <= max_size)].index)
    log(f"       - Term-size filter (units): {args.min_term_size} <= n <= {max_size}; terms before/after: "
        f"{len(n_by_term)}/{len(valid)}")
    study_unit_terms = unit_terms[unit_terms['Unit'].isin(study_units)]
    k_by_term = study_unit_terms.groupby('GO_ID').size()
    terms = [t for t in k_by_term.index if t in valid]
    if not terms:
        return pd.DataFrame()

    gene_terms = d[['GO_ID', 'Gene']].drop_duplicates()
    n_genes = gene_terms.groupby('GO_ID').size()
    study_gene_terms = gene_terms[gene_terms['Gene'].isin(study_in)]
    k_genes = study_gene_terms.groupby('GO_ID').size()
    gene_ids = study_gene_terms.groupby('GO_ID')['Gene'].apply(lambda s: '/'.join(sorted(s)))

    k = k_by_term.loc[terms].values.astype(int)
    n = n_by_term.loc[terms].values.astype(int)
    pvalues = hypergeom.sf(k - 1, M, n, N)
    bg_ratio = n / M
    fold = (k / N) / bg_ratio

    position = {u: i for i, u in enumerate(sorted(study_units))}
    counts = np.random.default_rng(args.random_seed).multinomial(N, np.full(N, 1.0 / N), size=args.n_bootstrap)
    unit_idx = study_unit_terms.assign(_i=study_unit_terms['Unit'].map(position)).groupby('GO_ID')['_i'].apply(
        lambda s: s.values)
    ci = np.array([np.percentile((counts[:, unit_idx[t]].sum(axis=1) / N) / r, BOOTSTRAP_CI)
                   for t, r in zip(terms, bg_ratio)])

    df = pd.DataFrame({
        'ID': terms,
        'Description': [go_dag[t].name if t in go_dag else 'Unknown' for t in terms],
        'GeneRatio': [f"{int(k_genes.get(t, 0)):,}/{len(study_in):,}" for t in terms],
        'BgRatio': [f"{int(n_genes.get(t, 0)):,}/{n_bg_genes:,}" for t in terms],
        'EvidenceUnitRatio': [f"{a:,}/{N:,}" for a in k],
        'BgEvidenceUnitRatio': [f"{a:,}/{M:,}" for a in n],
        'FoldEnrichment': fold,
        'FoldEnrichment_CI95_low': ci[:, 0],
        'FoldEnrichment_CI95_high': ci[:, 1],
        'pvalue': pvalues,
        'geneID': [gene_ids.get(t, '') for t in terms],
        'Count': [int(k_genes.get(t, 0)) for t in terms],
    })
    df['p.adjust'] = multipletests(df['pvalue'], alpha=args.fdr_alpha, method='fdr_bh')[1]
    return df[df['p.adjust'] < args.fdr_alpha].sort_values('pvalue')


def filter_terminal_go(enrichment, go_dag):
    """Drops terms that are ancestors of another significant term."""
    if enrichment.empty:
        return enrichment
    ids = set(enrichment['ID'])
    ancestors = {g: go_dag[g].get_all_parents() for g in ids if g in go_dag}
    redundant = {g for g in ids for other in ids - {g} if other in ancestors and g in ancestors[other]}
    return enrichment[~enrichment['ID'].isin(redundant)].copy()


def run_ipr_enrichment(study_genes, background_genes, annotations, args):
    universe = set(background_genes)
    study = set(study_genes) & universe
    term_genes = annotations[annotations['Gene'].isin(universe)].groupby('IPR_ID')['Gene'].apply(set)
    M, N = len(universe), len(study)
    rows = []
    for term_id, genes in term_genes.items():
        hits = study & genes
        if hits:
            n, k = len(genes), len(hits)
            rows.append({'ID': term_id, 'Description': 'Unknown', 'GeneRatio': f"{k:,}/{N:,}",
                         'BgRatio': f"{n:,}/{M:,}", 'pvalue': hypergeom.sf(k - 1, M, n, N),
                         'geneID': '/'.join(sorted(hits)), 'Count': k})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    if args.ipr_legacy_p:
        return df[df['pvalue'] < args.fdr_alpha]
    df['p.adjust'] = multipletests(df['pvalue'], alpha=args.fdr_alpha, method='fdr_bh')[1]
    return df[df['p.adjust'] < args.fdr_alpha].sort_values('pvalue')


def worked_example(go_id, study_genes, background_genes, bp_annotations, go_dag, log):
    study = set(study_genes)
    d = bp_annotations[bp_annotations['Gene'].isin(set(background_genes))].copy()
    d['is_cl'] = is_cluster_row(d).values
    log(f"\n    ===== Worked example: {go_id} ({go_dag[go_id].name if go_id in go_dag else 'Unknown'}), "
        f"Biological Process =====")
    for mode in ('legacy', 'exclusive'):
        d['Unit'], _ = build_units(d, mode)
        M = d['Unit'].nunique()
        study_units = set(d.loc[d['Gene'].isin(study), 'Unit'])
        N = len(study_units)
        carriers = set(d.loc[d['GO_ID'] == go_id, 'Unit'])
        n, k = len(carriers), len(carriers & study_units)
        p = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
        log(f"    [{mode:9}] M = {M}, n = {n}, N = {N}, k = {k}, expected k = {N * n / M:.3f}, p = {p:.3g}")
    sub = d[(d['GO_ID'] == go_id) & d['Unit'].isin(study_units)]
    log(f"    Study units carrying {go_id} (exclusive): {sub['Unit'].nunique()}")
    for unit, g in sub.groupby('Unit'):
        members = sorted(set(d.loc[(d['Unit'] == unit) & d['Gene'].isin(study), 'Gene']))
        direct = sorted(set(g.loc[~g['is_cl'], 'Gene']))
        propagated = bool(g['is_cl'].any())
        kind = 'cluster' if CLUSTER_TAG in str(unit) else 'single gene'
        how = 'direct and propagated' if direct and propagated else 'direct only' if direct else 'propagated only'
        log(f"      unit {unit} ({kind}); study genes: {', '.join(members)}; evidence: {how}"
            + (f"; genes with direct evidence: {', '.join(direct)}" if direct else ''))
    both = set(sub.loc[~sub['is_cl'], 'Gene']) & set(sub.loc[sub['is_cl'], 'Gene']) & study
    log(f"    Study genes with direct and propagated rows for this term (counted twice under 'legacy'): {len(both)}")


def read_ids(path):
    return pd.read_csv(path, header=None)[0].astype(str).tolist()


def load_go_annotations(path):
    df = pd.read_csv(path, sep='\t')
    if not REQUIRED_GO_COLUMNS.issubset(df.columns):
        raise ValueError(f"{path.name} lacks columns {REQUIRED_GO_COLUMNS}; found {list(df.columns)}")
    return df.rename(columns={'Gene_ID': 'Gene'})


def drop_propagation(go_df, prefixes, log):
    group = go_df['Evidence_Group_ID'].astype(str)
    drop = np.zeros(len(go_df), dtype=bool)
    for prefix in prefixes:
        drop |= (is_cluster_row(go_df) & group.str.startswith(prefix)).values
    log(f"    [INFO] Dropping {int(drop.sum())} propagated rows from genomes {prefixes}.")
    if not drop.any():
        log(f"    [WARNING] No row matched; example cluster ids: {group[is_cluster_row(go_df)].head(3).tolist()}")
    return go_df[~drop]


def run_tag(args):
    tag = f"_u{'X' if args.unit_mode == 'exclusive' else 'L'}_s{args.propagation_scope}" \
          f"_b{'foc' if args.background_mode == 'focal' else 'com'}"
    return tag + ('_dropK1' if args.drop_propagation_from_genomes else '') + ('_excl' if args.exclude_genes else '')


def process_one(sample_dir, project_root, taxon_dir, annotation_file, background_file, args, log):
    suffix = ANNOTATION_FILES[annotation_file][1]
    tag = run_tag(args)
    mode = 'genus' if taxon_dir.name.startswith('Genus_') else 'species'
    out_dir = taxon_dir / OUTPUT_DIR
    out_dir.mkdir(exist_ok=True)
    excel_path = out_dir / f"{taxon_dir.name}_enrichment{suffix}{tag}.xlsx"
    record = {'taxon': taxon_dir.name, 'annotation_file': annotation_file, 'background_file': background_file,
              'tag': tag, 'status': None, 'output_excel': str(excel_path), 'error': None}
    if args.skip_existing and excel_path.is_file():
        record['status'] = 'skipped_existing'
        return record
    try:
        background_path = Path(args.focal_background_file) if args.background_mode == 'focal' \
            else sample_dir / background_file
        paths = {'gene_list': taxon_dir / GENE_LIST_DIR / f"result_{mode}.txt", 'background': background_path,
                 'go': sample_dir / annotation_file, 'ipr': sample_dir / IPR_ANNOTATION_FILE,
                 'obo': project_root.joinpath(*GO_OBO), 'ipr_entries': project_root.joinpath(*IPR_ENTRY_LIST)}
        missing = [str(p) for p in paths.values() if not p.is_file()]
        if missing:
            raise FileNotFoundError(f"Missing file(s): {missing}")

        study = read_ids(paths['gene_list'])
        if args.exclude_genes:
            excluded = set(read_ids(Path(args.exclude_genes)))
            log(f"    [INFO] Excluding {len(excluded & set(study))} genes from the foreground (listed: {len(excluded)}).")
            study = [g for g in study if g not in excluded]
        background = read_ids(paths['background'])
        log(f"    [INFO] Foreground genes: {len(study)}; background: {len(background)} ({args.background_mode}); "
            f"foreground genes missing from the background: {len(set(study) - set(background))}")

        go_df = load_go_annotations(paths['go'])
        if args.drop_propagation_from_genomes:
            prefixes = [x.strip() for x in args.drop_propagation_from_genomes.split(',') if x.strip()]
            go_df = drop_propagation(go_df, prefixes, log)
        go_dag = GODag(str(paths['obo']), load_obsolete=False, prt=None)
        namespace = {g: (go_dag[g].namespace if g in go_dag else None) for g in go_df['GO_ID'].unique()}
        go_df['Namespace'] = go_df['GO_ID'].map(namespace)
        go_df = go_df.dropna(subset=['Namespace'])
        if args.propagation_scope == 'BP':
            drop = is_cluster_row(go_df) & (go_df['Namespace'] != 'biological_process')
            log(f"    [INFO] Propagation scope BP: dropping {int(drop.sum())} propagated MF/CC rows.")
            go_df = go_df[~drop]
        ipr_df = pd.read_csv(paths['ipr'], sep='\t', names=['IPR_ID', 'Gene'], skiprows=1)

        raw, filtered = {}, {}
        for ns, label in NAMESPACES:
            log(f"    [INFO] GO {label} enrichment...")
            sub = go_df[go_df['Namespace'] == ns]
            raw[label] = run_go_enrichment(study, background, sub, go_dag, args.unit_mode, log, args)
            filtered[label] = filter_terminal_go(raw[label], go_dag)
            for table in (raw[label], filtered[label]):
                if not table.empty:
                    table['Category'] = label
            if label == 'BP' and args.worked_example:
                worked_example(args.worked_example, study, background, sub, go_dag, log)
        log("    [INFO] InterPro enrichment...")
        ipr = run_ipr_enrichment(study, background, ipr_df, args)
        if not ipr.empty:
            entries = pd.read_csv(paths['ipr_entries'], sep='\t', header=None,
                                  names=['IPR_ID', 'ENTRY_TYPE', 'ENTRY_NAME'])
            ipr['Description'] = ipr['ID'].map(entries.set_index('IPR_ID')['ENTRY_NAME']).fillna('Unknown')

        labels = [label for _, label in NAMESPACES]
        tables = {f"GO_{l}_filtered": filtered[l] for l in labels}
        tables.update({f"GO_{l}_raw": raw[l] for l in labels})
        tables['ALL_GO_filtered'] = pd.concat([filtered[l] for l in labels], ignore_index=True)
        tables['ALL_GO_raw'] = pd.concat([raw[l] for l in labels], ignore_index=True)
        tables['IPR'] = ipr
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            non_empty = [name for name in sorted(tables) if not tables[name].empty]
            for name in non_empty:
                tables[name].to_excel(writer, sheet_name=name, index=False)
            if not non_empty:
                pd.DataFrame({'note': ['no significant terms']}).to_excel(writer, sheet_name='EMPTY', index=False)
        for name in (n for n in tables if not tables[n].empty):
            tables[name].to_csv(out_dir / f"{name}{suffix if 'GO' in name else ''}{tag}.tsv", sep='\t', index=False)
        record.update(status='success', n_sig_bp=len(filtered['BP']), n_sig_mf=len(filtered['MF']),
                      n_sig_cc=len(filtered['CC']), n_sig_ipr=len(ipr))
        log(f"    [INFO] Significant terms: BP {len(filtered['BP'])}, MF {len(filtered['MF'])}, "
            f"CC {len(filtered['CC'])}, InterPro {len(ipr)}")
    except Exception as e:
        log(f"    [ERROR] {taxon_dir.name} / {annotation_file}: {e}")
        log(traceback.format_exc())
        record.update(status='failed', error=str(e))
    return record


def background_for(args, annotation_file):
    if args.annotation_file != 'all' and args.background_file:
        return args.background_file
    return ANNOTATION_FILES[annotation_file][0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--base-dir', required=True)
    ap.add_argument('--sample-dir', required=True)
    ap.add_argument('--taxa', required=True, help="Comma-separated taxon directories, or 'all'.")
    ap.add_argument('--annotation-file', required=True, choices=list(ANNOTATION_FILES) + ['all'])
    ap.add_argument('--background-file', default=None)
    ap.add_argument('--unit-mode', choices=['exclusive', 'legacy'], default='exclusive')
    ap.add_argument('--propagation-scope', choices=['BP', 'all'], default='BP',
                    help="BP drops propagated MF/CC rows (Evidence_Group_ID containing '_cluster').")
    ap.add_argument('--background-mode', choices=['community', 'focal'], default='community')
    ap.add_argument('--focal-background-file', default=None, help="Protein ids of the focal genome, one per line.")
    ap.add_argument('--exclude-genes', default=None, help="Protein ids to remove from the foreground.")
    ap.add_argument('--drop-propagation-from-genomes', default=None,
                    help="Comma-separated Evidence_Group_ID prefixes whose propagated rows are removed.")
    ap.add_argument('--worked-example', default=None, help="BP GO id for which the counts are printed in detail.")
    ap.add_argument('--ipr-legacy-p', action='store_true', help="Raw p < alpha for InterPro instead of BH.")
    ap.add_argument('--min-term-size', type=int, default=5)
    ap.add_argument('--max-term-fraction', type=float, default=0.1)
    ap.add_argument('--fdr-alpha', type=float, default=0.05)
    ap.add_argument('--n-bootstrap', type=int, default=2000)
    ap.add_argument('--random-seed', type=int, default=0)
    ap.add_argument('--skip-existing', action='store_true')
    ap.add_argument('--manifest-path', default=None)
    args = ap.parse_args()
    if args.background_mode == 'focal' and not args.focal_background_file:
        sys.exit("Error: --focal-background-file is required with --background-mode focal.")

    project_root = Path(args.base_dir).resolve()
    sample_dir = project_root / args.sample_dir
    if not sample_dir.is_dir():
        sys.exit(f"Error: sample directory not found: {sample_dir}")
    if args.taxa.strip().lower() == 'all':
        taxon_dirs = sorted(d for d in sample_dir.iterdir() if d.is_dir() and d.name.startswith(TAXON_PREFIXES))
    else:
        taxon_dirs = [sample_dir / t.strip() for t in args.taxa.split(',') if t.strip()]
        for d in taxon_dirs:
            if not d.is_dir():
                sys.exit(f"Error: taxon directory not found: {d}")
    combos = list(ANNOTATION_FILES) if args.annotation_file == 'all' else [args.annotation_file]
    for annotation_file in combos:
        background_file = background_for(args, annotation_file)
        if not (sample_dir / annotation_file).is_file():
            sys.exit(f"Error: annotation file not found: {sample_dir / annotation_file}")
        if args.background_mode == 'community' and not (sample_dir / background_file).is_file():
            sys.exit(f"Error: background file not found: {sample_dir / background_file}")

    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    manifest_path = Path(args.manifest_path) if args.manifest_path else sample_dir / f"enrichment_v2_manifest_{stamp}.json"
    log_lines = []

    def log(msg):
        print(msg)
        log_lines.append(msg)

    log(f"[INFO] Run started {stamp}; options: {vars(args)}")
    results = []
    for annotation_file in combos:
        background_file = background_for(args, annotation_file)
        for taxon_dir in taxon_dirs:
            log(f"\n{'=' * 20} {taxon_dir.name} [{annotation_file} / {background_file}] {'=' * 20}")
            results.append(process_one(sample_dir, project_root, taxon_dir, annotation_file, background_file,
                                       args, log))
    manifest = {'timestamp_utc': stamp, 'options': vars(args), 'results': results, 'log': log_lines}
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False, default=str)
    log(f"[INFO] Manifest: {manifest_path}")
    if any(r['status'] == 'failed' for r in results):
        sys.exit(1)


if __name__ == '__main__':
    main()
