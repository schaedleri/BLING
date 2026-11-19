#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GO / InterPro enrichment runner, version 2 (non-interactive).

Changes with respect to run_enrichment_production.py
  1. Evidence units are an EXCLUSIVE partition of the genes (--unit-mode exclusive,
     default). A gene that belongs to a propagation cluster is represented ONLY by its
     cluster; every other gene represents itself. A unit carries a GO term if any of its
     members carries it, whether directly or by propagation. A gene with direct evidence
     for a term that is also propagated within its cluster is therefore counted once.
     (--unit-mode legacy reproduces the earlier definition, in which such a gene was
     counted both as a gene unit and as a cluster unit.)
  2. Propagated annotations can be restricted to Biological Process
     (--propagation-scope BP, default): propagated rows (Evidence_Group_ID containing
     '_cluster') of the MF and CC namespaces, which only arise from ancestor expansion of
     propagated BP terms, are dropped. --propagation-scope all keeps them (earlier behaviour).
  3. Background: --background-mode community (default; the background list of the
     annotation tier) or focal (all proteins of the focal genome, --focal-background-file).
  4. Bootstrap confidence intervals resample the evidence units with multiplicity
     (the earlier version collapsed resampled units into a set, which lowered N).
  5. InterPro enrichment uses Benjamini-Hochberg FDR over all tested entries
     (--ipr-legacy-p restores the raw p < 0.05 rule).
  6. --exclude-genes removes genes from the foreground; --drop-propagation-from-genomes
     removes the propagated rows generated from the listed genomes (sensitivity analysis,
     e.g. for genomes with a degenerate single-component GMM fit).
  7. --worked-example GO:xxxxxxx prints the counts M, n, N, k under both unit definitions
     and the units contributing to k, for one GO term (BP).

Statistics (per namespace)
  M : evidence units with at least one annotation in the namespace among the background genes
  n : units among them carrying the term
  N : units containing at least one study gene
  k : study units carrying the term
  p = hypergeom.sf(k - 1, M, n, N); fold = (k / N) / (n / M); BH over all terms with k > 0 that
  pass the term-size filter 5 <= n <= 0.1 * M (in evidence units).
"""
import sys
import os
import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from goatools.obo_parser import GODag
from statsmodels.stats.multitest import multipletests

ANNOTATION_BACKGROUND_MAP = {
    "ncbi_go_gene.tsv": "ncbi_background_list.txt",
    "direct_go_gene.tsv": "direct_background_list.txt",
    "fixed_go_gene.tsv": "fixed_background_list.txt",
    "GMM_go_gene.tsv": "GMM_background_list.txt",
}
OUTPUT_SUFFIX_MAP = {
    "ncbi_go_gene.tsv": "_from_ncbi",
    "direct_go_gene.tsv": "_from_direct",
    "fixed_go_gene.tsv": "_from_fixed",
    "GMM_go_gene.tsv": "_from_GMM",
}


def is_cluster_series(df):
    return df['Evidence_Group_ID'].astype(str).str.contains('_cluster', na=False)


def build_units(df, mode):
    """Returns (array of unit ids aligned with df, set of genes assigned to more than one cluster).

    Exclusive mode: a gene that carries cluster rows is represented by its cluster. Identical RefSeq
    proteins (same WP_ accession) occurring in several genomes can belong to clusters of more than one
    genome; such a gene is assigned to the lexicographically smallest cluster id, so that the result
    does not depend on the row order of the input file."""
    is_cl = is_cluster_series(df).values
    if mode == 'legacy':
        unit = np.where(is_cl, df['Evidence_Group_ID'].values, df['Gene'].values)
        return unit, set()
    cl = df.loc[is_cl, ['Gene', 'Evidence_Group_ID']].drop_duplicates()
    nun = cl.groupby('Gene')['Evidence_Group_ID'].nunique()
    multi = set(nun[nun > 1].index)
    g2c = cl.groupby('Gene')['Evidence_Group_ID'].min()
    unit = df['Gene'].map(g2c)
    unit = unit.where(unit.notna(), df['Gene'])
    return unit.values, multi


def run_go_enrichment(gene_list, background_list, term2gene_df, go_dag, unit_mode='exclusive',
                      log=print, n_bootstrap=2000, random_seed=0):
    bg, st = set(background_list), set(gene_list)
    d = term2gene_df[term2gene_df['Gene'].isin(bg)].copy()
    if d.empty:
        return pd.DataFrame()
    d['Unit'], multi = build_units(d, unit_mode)
    if multi:
        log(f"       [NOTE] {len(multi)} background genes carry more than one cluster id (identical proteins in "
            f"several genomes); the smallest cluster id was used. Study genes affected: {len(multi & st)}.")
    eff_bg_genes = d['Gene'].nunique()
    study_genes_in = st & set(d['Gene'])
    M = d['Unit'].nunique()
    study_units = set(d.loc[d['Gene'].isin(study_genes_in), 'Unit'])
    N = len(study_units)
    log(f"       - Background genes: {eff_bg_genes}; study genes with annotation: {len(study_genes_in)}")
    log(f"       - M (background evidence units) = {M}; N (study evidence units) = {N}  [unit mode: {unit_mode}]")
    if N == 0:
        return pd.DataFrame()

    ut = d[['GO_ID', 'Unit']].drop_duplicates()
    n_by = ut.groupby('GO_ID').size()
    min_gs, max_gs = 5, int(M * 0.1)
    valid = set(n_by[(n_by >= min_gs) & (n_by <= max_gs)].index)
    log(f"       - Term-size filter (evidence units): {min_gs} <= n <= {max_gs}; terms before/after: {len(n_by)}/{len(valid)}")
    ut_s = ut[ut['Unit'].isin(study_units)]
    k_by = ut_s.groupby('GO_ID').size()
    terms = [t for t in k_by.index if t in valid]
    if not terms:
        return pd.DataFrame()

    gt = d[['GO_ID', 'Gene']].drop_duplicates()
    n_genes_by = gt.groupby('GO_ID').size()
    gt_s = gt[gt['Gene'].isin(study_genes_in)]
    k_genes_by = gt_s.groupby('GO_ID').size()
    genes_by = gt_s.groupby('GO_ID')['Gene'].apply(lambda s: '/'.join(sorted(s)))

    k = k_by.loc[terms].values.astype(int)
    n = n_by.loc[terms].values.astype(int)
    pvals = hypergeom.sf(k - 1, M, n, N)
    bg_ratio = n / M
    fold = (k / N) / bg_ratio

    # bootstrap of the study units WITH multiplicity
    units_list = sorted(study_units)
    pos = {u: i for i, u in enumerate(units_list)}
    rng = np.random.default_rng(random_seed)
    counts = rng.multinomial(N, np.full(N, 1.0 / N), size=n_bootstrap)
    idx_by = ut_s.assign(_i=ut_s['Unit'].map(pos)).groupby('GO_ID')['_i'].apply(lambda s: s.values)
    lo, hi = [], []
    for t, br in zip(terms, bg_ratio):
        kb = counts[:, idx_by[t]].sum(axis=1)
        fb = (kb / N) / br
        a, b = np.percentile(fb, [2.5, 97.5])
        lo.append(a)
        hi.append(b)

    df = pd.DataFrame({
        "ID": terms,
        "Description": [go_dag[t].name if t in go_dag else "Unknown" for t in terms],
        "GeneRatio": [f"{int(k_genes_by.get(t, 0)):,}/{len(study_genes_in):,}" for t in terms],
        "BgRatio": [f"{int(n_genes_by.get(t, 0)):,}/{eff_bg_genes:,}" for t in terms],
        "EvidenceUnitRatio": [f"{int(a):,}/{N:,}" for a in k],
        "BgEvidenceUnitRatio": [f"{int(a):,}/{M:,}" for a in n],
        "FoldEnrichment": fold,
        "FoldEnrichment_CI95_low": lo,
        "FoldEnrichment_CI95_high": hi,
        "pvalue": pvals,
        "geneID": [genes_by.get(t, '') for t in terms],
        "Count": [int(k_genes_by.get(t, 0)) for t in terms],
    })
    _, padj, _, _ = multipletests(df['pvalue'], alpha=0.05, method='fdr_bh')
    df['p.adjust'] = padj
    return df[df['p.adjust'] < 0.05].sort_values(by='pvalue')


def filter_terminal_go(enrich_df, go_dag):
    if enrich_df.empty:
        return enrich_df
    go_ids = set(enrich_df['ID'])
    anc = {g: go_dag[g].get_all_parents() for g in go_ids if g in go_dag}
    redundant = {g for g in go_ids for o in (go_ids - {g}) if o in anc and g in anc[o]}
    return enrich_df[~enrich_df['ID'].isin(redundant)]


def run_ipr_enrichment(gene_list, background_list, term2gene_df, legacy_p=False):
    universe, study = set(background_list), set(gene_list) & set(background_list)
    t2g = term2gene_df[term2gene_df['Gene'].isin(universe)]
    term_to_genes = t2g.groupby('IPR_ID')['Gene'].apply(set)
    M, N = len(universe), len(study)
    rows = []
    for term_id, genes in term_to_genes.items():
        k = len(study & genes)
        if k == 0:
            continue
        n = len(genes)
        rows.append({"ID": term_id, "Description": "Unknown", "GeneRatio": f"{k:,}/{N:,}",
                     "BgRatio": f"{n:,}/{M:,}", "pvalue": hypergeom.sf(k - 1, M, n, N),
                     "geneID": '/'.join(sorted(study & genes)), "Count": k})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    if legacy_p:
        return df[df['pvalue'] < 0.05]
    _, padj, _, _ = multipletests(df['pvalue'], alpha=0.05, method='fdr_bh')
    df['p.adjust'] = padj
    return df[df['p.adjust'] < 0.05].sort_values('pvalue')


def worked_example(go_id, gene_list, background_list, df_bp, go_dag, log):
    bg, st = set(background_list), set(gene_list)
    d = df_bp[df_bp['Gene'].isin(bg)].copy()
    name = go_dag[go_id].name if go_id in go_dag else 'Unknown'
    log(f"\n    ===== Worked example: {go_id} ({name}), Biological Process =====")
    d['is_cl'] = is_cluster_series(d).values
    for mode in ('legacy', 'exclusive'):
        d['Unit'], _ = build_units(d, mode)
        M = d['Unit'].nunique()
        su = set(d.loc[d['Gene'].isin(st), 'Unit'])
        N = len(su)
        has = set(d.loc[d['GO_ID'] == go_id, 'Unit'])
        n, k = len(has), len(has & su)
        p = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0
        log(f"    [{mode:9}] M = {M}, n = {n}, N = {N}, k = {k}, expected k = {N*n/M:.3f}, p = {p:.3g}")
    d['Unit'], _ = build_units(d, 'exclusive')
    su = set(d.loc[d['Gene'].isin(st), 'Unit'])
    sub = d[(d['GO_ID'] == go_id) & (d['Unit'].isin(su))]
    log(f"    Study units carrying {go_id} (exclusive definition): {sub['Unit'].nunique()}")
    for u, g in sub.groupby('Unit'):
        members = sorted(set(d.loc[(d['Unit'] == u) & d['Gene'].isin(st), 'Gene']))
        direct = sorted(set(g.loc[~g['is_cl'], 'Gene']))
        prop = bool(g['is_cl'].any())
        kind = 'cluster' if str(u).find('_cluster') >= 0 else 'single gene'
        how = ('direct and propagated' if direct and prop else 'direct only' if direct else 'propagated only')
        log(f"      unit {u} ({kind}); study genes: {', '.join(members)}; evidence: {how}"
            + (f"; genes with direct evidence: {', '.join(direct)}" if direct else ''))
    both = set(sub.loc[~sub['is_cl'], 'Gene']) & set(sub.loc[sub['is_cl'], 'Gene']) & st
    log(f"    Study genes with both direct and propagated rows for this term (counted twice under the legacy definition): {len(both)}")


def load_go_annotation_file(path, log=print):
    df = pd.read_csv(path, sep='\t')
    need = {"GO_ID", "Gene_ID", "Evidence_Group_ID"}
    if not need.issubset(df.columns):
        raise ValueError(f"{path.name} lacks columns {need}; found {list(df.columns)}")
    return df.rename(columns={"Gene_ID": "Gene"})


def read_ids(path):
    return pd.read_csv(path, header=None)[0].astype(str).tolist()


def process_one(sample_dir, project_root, taxon_dir, annotation_file, background_file, args, run_log):
    suffix = OUTPUT_SUFFIX_MAP.get(annotation_file, "_from_unknown")
    tag = f"_u{'X' if args.unit_mode == 'exclusive' else 'L'}_s{args.propagation_scope}_b{'foc' if args.background_mode == 'focal' else 'com'}"
    if args.drop_propagation_from_genomes:
        tag += "_dropK1"
    if args.exclude_genes:
        tag += "_excl"
    mode = 'genus' if taxon_dir.name.startswith('Genus_') else 'species'
    out_dir = taxon_dir / "Enrichment_Analysis"
    out_dir.mkdir(exist_ok=True)
    excel_path = out_dir / f"{taxon_dir.name}_enrichment{suffix}{tag}.xlsx"
    rec = {"taxon": taxon_dir.name, "annotation_file": annotation_file, "background_file": background_file,
           "tag": tag, "status": None, "output_excel": str(excel_path), "error": None}
    if args.skip_existing and excel_path.is_file():
        rec["status"] = "skipped_existing"
        return rec
    try:
        bg_path = Path(args.focal_background_file) if args.background_mode == 'focal' else sample_dir / background_file
        paths = {"gene_list": taxon_dir / "BLAST_specify" / f"result_{mode}.txt", "background_list": bg_path,
                 "go_annot": sample_dir / annotation_file, "ipr_annot": sample_dir / "ipr_gene.tsv",
                 "go_obo": project_root / "data" / "go.obo", "ipr_entry": project_root / "data" / "entry.list.txt"}
        missing = [str(p) for p in paths.values() if not p.is_file()]
        if missing:
            raise FileNotFoundError(f"Missing file(s): {missing}")
        gene_list = read_ids(paths["gene_list"])
        if args.exclude_genes:
            ex = set(read_ids(Path(args.exclude_genes)))
            run_log(f"    [INFO] Excluding {len(ex & set(gene_list))} genes from the foreground (listed: {len(ex)}).")
            gene_list = [g for g in gene_list if g not in ex]
        background_list = read_ids(paths["background_list"])
        run_log(f"    [INFO] Foreground genes: {len(gene_list)}; background list: {len(background_list)} "
                f"({args.background_mode}); foreground genes missing from the background: "
                f"{len(set(gene_list) - set(background_list))}")

        go_df = load_go_annotation_file(paths["go_annot"], run_log)
        is_cl = is_cluster_series(go_df)
        if args.drop_propagation_from_genomes:
            prefixes = [x.strip() for x in args.drop_propagation_from_genomes.split(',') if x.strip()]
            m = np.zeros(len(go_df), dtype=bool)
            egid = go_df['Evidence_Group_ID'].astype(str)
            for pfx in prefixes:
                m |= (is_cl & egid.str.startswith(pfx)).values
            run_log(f"    [INFO] Dropping {int(m.sum())} propagated rows from genomes {prefixes}.")
            if m.sum() == 0:
                run_log(f"    [WARNING] No row matched; example cluster ids: {egid[is_cl].head(3).tolist()}")
            go_df = go_df[~m]
        go_dag = GODag(str(paths["go_obo"]), load_obsolete=False, prt=None)
        ns_map = {g: (go_dag[g].namespace if g in go_dag else None) for g in go_df['GO_ID'].unique()}
        go_df['Namespace'] = go_df['GO_ID'].map(ns_map)
        go_df = go_df.dropna(subset=['Namespace'])
        if args.propagation_scope == 'BP':
            drop = is_cluster_series(go_df) & (go_df['Namespace'] != 'biological_process')
            run_log(f"    [INFO] Propagation scope BP: dropping {int(drop.sum())} propagated MF/CC rows.")
            go_df = go_df[~drop]
        ipr_df = pd.read_csv(paths["ipr_annot"], sep='\t', names=["IPR_ID", "Gene"], skiprows=1)

        res, raw = {}, {}
        for ns, label in (('biological_process', 'BP'), ('molecular_function', 'MF'), ('cellular_component', 'CC')):
            run_log(f"    [INFO] Running GO {label} enrichment...")
            sub = go_df[go_df['Namespace'] == ns]
            r = run_go_enrichment(gene_list, background_list, sub, go_dag, args.unit_mode, run_log,
                                  random_seed=args.random_seed)
            raw[label] = r.copy()
            res[label] = filter_terminal_go(r.copy(), go_dag)
            for dd in (raw[label], res[label]):
                if not dd.empty:
                    dd.loc[:, 'Category'] = label
            if label == 'BP' and args.worked_example:
                worked_example(args.worked_example, gene_list, background_list, sub, go_dag, run_log)
        run_log("    [INFO] Running InterPro enrichment...")
        ipr_res = run_ipr_enrichment(gene_list, background_list, ipr_df, args.ipr_legacy_p)
        if not ipr_res.empty:
            ent = pd.read_csv(paths["ipr_entry"], sep='\t', header=None, names=["IPR_ID", "ENTRY_TYPE", "ENTRY_NAME"])
            ipr_res['Description'] = ipr_res['ID'].map(ent.set_index('IPR_ID')['ENTRY_NAME']).fillna("Unknown")

        out = {"GO_BP_filtered": res['BP'], "GO_MF_filtered": res['MF'], "GO_CC_filtered": res['CC'],
               "ALL_GO_filtered": pd.concat([res['BP'], res['MF'], res['CC']], ignore_index=True), "IPR": ipr_res,
               "GO_BP_raw": raw['BP'], "GO_MF_raw": raw['MF'], "GO_CC_raw": raw['CC'],
               "ALL_GO_raw": pd.concat([raw['BP'], raw['MF'], raw['CC']], ignore_index=True)}
        with pd.ExcelWriter(excel_path, engine='openpyxl') as w:
            wrote = False
            for name in sorted(out):
                if not out[name].empty:
                    out[name].to_excel(w, sheet_name=name, index=False)
                    wrote = True
            if not wrote:
                pd.DataFrame({"note": ["no significant terms"]}).to_excel(w, sheet_name="EMPTY", index=False)
        for name, dd in out.items():
            if not dd.empty:
                dd.to_csv(out_dir / f"{name}{suffix if 'GO' in name else ''}{tag}.tsv", sep='\t', index=False)
        rec.update(status="success", n_sig_bp=len(res['BP']), n_sig_mf=len(res['MF']), n_sig_cc=len(res['CC']),
                   n_sig_ipr=len(ipr_res))
        run_log(f"    [INFO] Significant terms: BP {len(res['BP'])}, MF {len(res['MF'])}, CC {len(res['CC'])}, "
                f"InterPro {len(ipr_res)}")
    except Exception as e:
        import traceback
        run_log(f"    [ERROR] {taxon_dir.name} / {annotation_file}: {e}")
        run_log(traceback.format_exc())
        rec.update(status="failed", error=str(e))
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--base-dir', required=True)
    ap.add_argument('--sample-dir', required=True)
    ap.add_argument('--taxa', required=True, help="Comma-separated taxon directories, or 'all'.")
    ap.add_argument('--annotation-file', required=True, choices=list(ANNOTATION_BACKGROUND_MAP) + ["all"])
    ap.add_argument('--background-file', default=None)
    ap.add_argument('--unit-mode', choices=['exclusive', 'legacy'], default='exclusive')
    ap.add_argument('--propagation-scope', choices=['BP', 'all'], default='BP')
    ap.add_argument('--background-mode', choices=['community', 'focal'], default='community')
    ap.add_argument('--focal-background-file', default=None,
                    help="Text file with all protein ids of the focal genome (one per line).")
    ap.add_argument('--exclude-genes', default=None, help="Text file with protein ids to remove from the foreground.")
    ap.add_argument('--drop-propagation-from-genomes', default=None,
                    help="Comma-separated genome ids (prefix of Evidence_Group_ID) whose propagated rows are removed.")
    ap.add_argument('--worked-example', default=None, help="GO id (BP) for which counts are printed in detail.")
    ap.add_argument('--ipr-legacy-p', action='store_true')
    ap.add_argument('--random-seed', type=int, default=0)
    ap.add_argument('--skip-existing', action='store_true')
    ap.add_argument('--cpu', type=int, default=1, help="Accepted for compatibility; not used.")
    ap.add_argument('--manifest-path', default=None)
    args = ap.parse_args()
    if args.background_mode == 'focal' and not args.focal_background_file:
        sys.exit("Error: --focal-background-file is required with --background-mode focal.")

    project_root = Path(args.base_dir).resolve()
    sample_dir = project_root / args.sample_dir
    if not sample_dir.is_dir():
        sys.exit(f"Error: sample directory not found: {sample_dir}")
    if args.taxa.strip().lower() == 'all':
        taxon_dirs = sorted(d for d in sample_dir.iterdir()
                            if d.is_dir() and d.name.startswith(('Genus_', 'Species_')))
    else:
        taxon_dirs = [sample_dir / t.strip() for t in args.taxa.split(',') if t.strip()]
        for d in taxon_dirs:
            if not d.is_dir():
                sys.exit(f"Error: taxon directory not found: {d}")
    combos = list(ANNOTATION_BACKGROUND_MAP.items()) if args.annotation_file == "all" else \
        [(args.annotation_file, args.background_file or ANNOTATION_BACKGROUND_MAP[args.annotation_file])]
    for a, b in combos:
        if not (sample_dir / a).is_file():
            sys.exit(f"Error: annotation file not found: {sample_dir / a}")
        if args.background_mode == 'community' and not (sample_dir / b).is_file():
            sys.exit(f"Error: background file not found: {sample_dir / b}")

    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = Path(args.manifest_path) if args.manifest_path else sample_dir / f"enrichment_v2_manifest_{ts}.json"
    log_lines = []

    def run_log(msg):
        print(msg)
        log_lines.append(msg)

    run_log(f"[INFO] Run started {ts}; options: {vars(args)}")
    results = []
    for a, b in combos:
        for td in taxon_dirs:
            run_log(f"\n{'='*20} {td.name} [{a} / {b}] {'='*20}")
            results.append(process_one(sample_dir, project_root, td, a, b, args, run_log))
    nf = sum(r['status'] == 'failed' for r in results)
    manifest = {"timestamp_utc": ts, "options": vars(args), "results": results, "log": log_lines,
                "package_versions": {"python": platform.python_version(), "pandas": pd.__version__,
                                     "numpy": np.__version__}}
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False, default=str)
    run_log(f"[INFO] Manifest: {manifest_path}")
    if nf:
        sys.exit(1)


if __name__ == "__main__":
    main()
