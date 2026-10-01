#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Production GO/InterPro enrichment runner (non-interactive).

Every run parameter is given as an explicit CLI argument, and each run
writes a run_manifest.json into the sample directory recording the exact
annotation/background combination(s), random seed, package versions, and a
per-taxon success/failure log, so that any result can be traced back to the
configuration that produced it. Pass --annotation-file all to run every
annotation/background combination in one invocation, for cross-tier
sensitivity comparisons.

Statistical design: genes within one operon-like cluster are not
independent observations, since a GO term propagated within a cluster is a
single piece of evidence copied across its members. The cluster is
therefore treated as a single evidence unit: predict_and_analyze_operons.py
and microbiome_interval.py write one Evidence_Group_ID per cluster for
propagated annotations, and one evidence group per (gene, GO term) pair for
direct evidence (GBFF/InterPro/eggNOG). A classical hypergeometric test
computed over the number of distinct evidence units (not genes) accounts
for this non-independence without requiring a permutation test.
GeneRatio/BgRatio are reported as gene-level counts for readability; the
p-value computation itself uses evidence-unit counts.
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
from collections import defaultdict
from multiprocessing import get_context
from scipy.stats import hypergeom
from goatools.obo_parser import GODag
from statsmodels.stats.multitest import multipletests

# Worker globals for parallelizing the bootstrap CI computation across GO
# terms (the hypergeometric test itself is closed-form and fast; the
# bootstrap resampling loop, run independently per term, is what benefits
# from parallelization). Populated via the Pool initializer and relied upon
# through fork's copy-on-write semantics on Linux.
_BW_STUDY_EG_LIST = None
_BW_GO2EG_BG = None


def _init_bootstrap_worker(study_evidence_units_list, go2evidence_units_bg):
    global _BW_STUDY_EG_LIST, _BW_GO2EG_BG
    _BW_STUDY_EG_LIST = study_evidence_units_list
    _BW_GO2EG_BG = go2evidence_units_bg


def _bootstrap_worker(args):
    """
    Computes a bootstrap 95% CI on fold enrichment for one GO term.
    Resamples (with replacement) the evidence units touching the study set
    (genes for direct evidence, clusters for propagated evidence), using an
    independent SeedSequence child so results are reproducible regardless
    of how work is split across worker processes.
    """
    go_id, seed_seq, n_bootstrap, n_study_units, bg_ratio = args
    rng = np.random.default_rng(seed_seq)
    bg_units_with_go = _BW_GO2EG_BG[go_id]

    boot_folds = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        sampled_idx = rng.integers(0, n_study_units, size=n_study_units)
        boot_units = {_BW_STUDY_EG_LIST[idx] for idx in sampled_idx}
        N_b = len(boot_units)
        k_b = len(boot_units & bg_units_with_go)
        boot_folds[b] = (k_b / N_b) / bg_ratio if N_b > 0 else np.nan

    valid_boot = boot_folds[~np.isnan(boot_folds)]
    if len(valid_boot) > 0:
        ci_low, ci_high = np.percentile(valid_boot, [2.5, 97.5])
    else:
        ci_low, ci_high = float('nan'), float('nan')
    return go_id, ci_low, ci_high


# Four nested annotation tiers, each with its own matching background:
#   ncbi_go_gene.tsv   -> ncbi_background_list.txt    (GBFF/NCBI only)
#   direct_go_gene.tsv -> direct_background_list.txt  (+ InterPro + eggNOG, no propagation)
#   fixed_go_gene.tsv  -> fixed_background_list.txt   (+ fixed-threshold propagation)
#   GMM_go_gene.tsv    -> GMM_background_list.txt     (+ GMM-based propagation)
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


def run_go_enrichment(gene_list, background_list, term2gene_df, go_dag, log=print,
                       n_bootstrap=2000, random_seed=0, n_jobs=1):
    """
    Classical hypergeometric enrichment test, using evidence units as the
    statistical unit rather than genes or raw Evidence_Group_ID instances.

    Evidence_Group_ID is a provenance record: for direct evidence
    (GBFF/InterPro/eggNOG) it is per (gene, GO term) pair; for propagated
    evidence it is per operon/cluster (shared across every GO term that
    cluster produces). Using Evidence_Group_ID directly as the
    population-size unit would over-count, since a single gene with direct
    evidence for five different GO terms would contribute five distinct
    Evidence_Group_IDs, inflating the population as if it were five
    different individuals.

    Evidence_Unit_ID instead identifies the actual sampling individual,
    term-independently:
      - direct evidence rows  -> the gene itself (Gene column)
      - propagated rows       -> the cluster (Evidence_Group_ID, already
                                  term-independent for clusters)
    M, N, n, k are all computed over Evidence_Unit_ID. This also correctly
    collapses the case where a single gene has independent direct evidence
    for the same GO term from two different sources (e.g. GBFF and eggNOG
    both separately support GO:X): Evidence_Group_ID would count that as
    two evidence instances, but Evidence_Unit_ID correctly treats it as one
    individual (one gene).

    GeneRatio/BgRatio in the output remain gene-level counts (for
    readability). Fold enrichment and its bootstrap 95% CI are computed
    consistently in Evidence_Unit_ID terms. The bootstrap loop is
    parallelized across GO terms via multiprocessing when n_jobs > 1.
    """
    background_set = set(background_list)
    study_set = set(gene_list)

    term2gene_bg = term2gene_df[term2gene_df['Gene'].isin(background_set)]
    effective_background_size = term2gene_bg['Gene'].nunique()
    log(f"       - Effective background (genes): {effective_background_size}")

    has_evidence_groups = 'Evidence_Group_ID' in term2gene_bg.columns
    if not has_evidence_groups:
        log("       [WARNING] No Evidence_Group_ID column found; falling back to "
            "treating each (GO_ID, Gene) pair as its own evidence group.")
        term2gene_bg = term2gene_bg.copy()
        term2gene_bg['Evidence_Group_ID'] = term2gene_bg['Gene'] + '_' + term2gene_bg['GO_ID']
    else:
        term2gene_bg = term2gene_bg.copy()

    # Derive the term-agnostic sampling individual for every row.
    # Cluster-level Evidence_Group_IDs (propagated evidence) already
    # identify a term-independent operon/cluster; direct-evidence IDs are
    # per (gene, term) and must be collapsed back to the gene itself.
    term2gene_bg['Evidence_Unit_ID'] = np.where(
        term2gene_bg['Evidence_Group_ID'].str.contains('_cluster'),
        term2gene_bg['Evidence_Group_ID'],
        term2gene_bg['Gene'],
    )

    term2gene_study = term2gene_df[term2gene_df['Gene'].isin(study_set)]
    effective_study_size = term2gene_study['Gene'].nunique()
    log(f"       - Effective study set (for GeneRatio, genes): {effective_study_size}")

    if term2gene_bg.empty:
        return pd.DataFrame()

    go2evidence_units_bg = defaultdict(set)
    go2genes_bg = defaultdict(set)
    for row in term2gene_bg.itertuples(index=False):
        go2evidence_units_bg[row.GO_ID].add(row.Evidence_Unit_ID)
        go2genes_bg[row.GO_ID].add(row.Gene)

    M = term2gene_bg['Evidence_Unit_ID'].nunique()
    log(f"       - Effective background (M, evidence units): {M}")

    study_genes_in_category = study_set.intersection(term2gene_bg['Gene'].unique())
    study_rows = term2gene_bg[term2gene_bg['Gene'].isin(study_genes_in_category)]
    study_evidence_units = set(study_rows['Evidence_Unit_ID'])
    N = len(study_evidence_units)
    log(f"       - Study set size (N, evidence units): {N}")

    min_gs_size = 5
    max_gs_size = int(M * 0.1)
    log(f"       - Filtering criteria (evidence-unit counts): minGSSize >= {min_gs_size}, maxGSSize <= {max_gs_size}")

    valid_go_ids = {
        go_id for go_id, units in go2evidence_units_bg.items()
        if min_gs_size <= len(units) <= max_gs_size
    }
    log(f"       - Term count before filtering: {len(go2evidence_units_bg)}, after filtering: {len(valid_go_ids)}")

    if not valid_go_ids:
        return pd.DataFrame()

    study_evidence_units_list = list(study_evidence_units)
    n_study_units = len(study_evidence_units_list)

    # Pass 1: closed-form hypergeometric test for every term, which also
    # determines which terms have k > 0 so the slower bootstrap pass below
    # only runs for terms that will actually be kept.
    term_stats = {}
    for go_id in valid_go_ids:
        bg_units_with_go = go2evidence_units_bg[go_id]
        study_units_with_go = bg_units_with_go & study_evidence_units

        k = len(study_units_with_go)
        if k == 0:
            continue
        n = len(bg_units_with_go)

        p_val = hypergeom.sf(k - 1, M, n, N)

        study_genes_with_go = study_genes_in_category & go2genes_bg[go_id]
        k_genes = len(study_genes_with_go)
        n_genes = len(go2genes_bg[go_id])

        bg_ratio = n / M if M > 0 else float('nan')
        fold_enrichment = (k / N) / bg_ratio if N > 0 and bg_ratio > 0 else float('nan')

        term_stats[go_id] = {
            "k": k, "n": n, "p_val": p_val, "k_genes": k_genes, "n_genes": n_genes,
            "fold_enrichment": fold_enrichment, "bg_ratio": bg_ratio,
            "study_genes_with_go": study_genes_with_go,
        }

    if not term_stats:
        return pd.DataFrame()

    # Pass 2: bootstrap 95% CI on fold enrichment, per kept term.
    # Parallelized across terms via multiprocessing when n_jobs > 1;
    # independent SeedSequence children guarantee reproducible,
    # non-overlapping random draws regardless of how work is distributed.
    n_jobs = max(1, int(n_jobs))
    base_seed_seq = np.random.SeedSequence(random_seed)
    go_ids_for_bootstrap = list(term_stats.keys())
    child_seed_seqs = base_seed_seq.spawn(len(go_ids_for_bootstrap))
    worker_args = [
        (go_id, seed_seq, n_bootstrap, n_study_units, term_stats[go_id]["bg_ratio"])
        for go_id, seed_seq in zip(go_ids_for_bootstrap, child_seed_seqs)
    ]

    log(f"       - Computing bootstrap 95% CIs for {len(worker_args)} terms (n_jobs={n_jobs})...")
    if n_jobs > 1 and len(worker_args) > 1:
        ctx = get_context("fork")
        with ctx.Pool(
            processes=n_jobs,
            initializer=_init_bootstrap_worker,
            initargs=(study_evidence_units_list, go2evidence_units_bg),
        ) as pool:
            bootstrap_results = pool.map(_bootstrap_worker, worker_args)
    else:
        _init_bootstrap_worker(study_evidence_units_list, go2evidence_units_bg)
        bootstrap_results = [_bootstrap_worker(a) for a in worker_args]

    for go_id, ci_low, ci_high in bootstrap_results:
        term_stats[go_id]["fold_ci_low"] = ci_low
        term_stats[go_id]["fold_ci_high"] = ci_high

    results = []
    for go_id, stats in term_stats.items():
        results.append({
            "ID": go_id,
            "Description": go_dag[go_id].name if go_id in go_dag else "Unknown",
            "GeneRatio": f"{stats['k_genes']:,}/{effective_study_size:,}",
            "BgRatio": f"{stats['n_genes']:,}/{effective_background_size:,}",
            "EvidenceUnitRatio": f"{stats['k']:,}/{N:,}",
            "BgEvidenceUnitRatio": f"{stats['n']:,}/{M:,}",
            "FoldEnrichment": stats["fold_enrichment"],
            "FoldEnrichment_CI95_low": stats["fold_ci_low"],
            "FoldEnrichment_CI95_high": stats["fold_ci_high"],
            "pvalue": stats["p_val"],
            "geneID": '/'.join(sorted(stats["study_genes_with_go"])),
            "Count": stats["k_genes"],
        })

    df = pd.DataFrame(results)
    reject_bh, pvals_corr_bh, _, _ = multipletests(df['pvalue'], alpha=0.05, method='fdr_bh')
    df['p.adjust'] = pvals_corr_bh

    significant_df = df[df['p.adjust'] < 0.05].copy()
    return significant_df.sort_values(by='pvalue')


def filter_terminal_go(enrich_df, go_dag):
    """Removes a term if a more specific (child) significant term is also present."""
    if enrich_df.empty:
        return enrich_df
    go_ids = set(enrich_df['ID'])
    all_ancestors = {go_id: go_dag[go_id].get_all_parents() for go_id in go_ids if go_id in go_dag}
    redundant_terms = {
        go_id for go_id in go_ids for other_id in (go_ids - {go_id})
        if other_id in all_ancestors and go_id in all_ancestors[other_id]
    }
    return enrich_df[~enrich_df['ID'].isin(redundant_terms)]


def run_ipr_enrichment(gene_list, background_list, term2gene_df):
    """
    Classical hypergeometric test for InterPro domain enrichment. IPR
    annotations are not produced via operon/cluster propagation anywhere in
    the pipeline -- they come directly from each protein's own InterProScan
    hit -- so no evidence-unit correction is needed here.
    """
    universe_genes, study_genes = set(background_list), set(gene_list)
    term2gene_df = term2gene_df[term2gene_df['Gene'].isin(universe_genes)]
    term_to_genes = term2gene_df.groupby('IPR_ID')['Gene'].apply(set)
    M = len(universe_genes)
    N = len(study_genes)
    results = []
    for term_id, term_genes in term_to_genes.items():
        k = len(study_genes.intersection(term_genes))
        if k == 0:
            continue
        n = len(term_genes)
        p_val = hypergeom.sf(k - 1, M, n, N)
        if p_val < 0.05:
            results.append({
                "ID": term_id,
                "Description": "Unknown",
                "GeneRatio": f"{k:,}/{N:,}",
                "BgRatio": f"{n:,}/{M:,}",
                "pvalue": p_val,
                "geneID": '/'.join(sorted(study_genes.intersection(term_genes))),
                "Count": k,
            })
    return pd.DataFrame(results)


def load_go_annotation_file(go_annot_path, log=print):
    """
    Loads a *_go_gene.tsv file. Expects the (GO_ID, Gene_ID, Evidence_Group_ID)
    format; falls back to an older two-column (GO_ID, Gene) format with a
    warning, treating each (GO_ID, Gene) pair as its own evidence group in
    that case. Returns a DataFrame with a 'Gene' column (renamed from
    'Gene_ID') and a boolean flag indicating whether Evidence_Group_ID was
    present.
    """
    go_annot_df = pd.read_csv(go_annot_path, sep='\t')
    expected_cols = {"GO_ID", "Gene_ID", "Evidence_Group_ID"}
    if expected_cols.issubset(set(go_annot_df.columns)):
        go_annot_df = go_annot_df.rename(columns={"Gene_ID": "Gene"})
        return go_annot_df, True
    else:
        log(f"[WARNING] '{go_annot_path.name}' does not contain the expected "
            f"(GO_ID, Gene_ID, Evidence_Group_ID) columns (found: {list(go_annot_df.columns)}). "
            f"Falling back to treating each (GO_ID, Gene) pair as an independent "
            f"evidence group, which may overstate significance for propagated annotations.")
        go_annot_df = pd.read_csv(go_annot_path, sep='\t', names=["GO_ID", "Gene"], skiprows=1)
        go_annot_df['Evidence_Group_ID'] = go_annot_df['Gene'] + '_' + go_annot_df['GO_ID']
        return go_annot_df, False


def process_one_combination(sample_dir, project_root, taxon_dir, annotation_file, background_file,
                             random_seed, skip_existing, run_log, n_jobs=1):
    """
    Runs GO (BP/MF/CC) and InterPro enrichment for one taxon, using one
    specific (annotation_file, background_file) combination. Writes results
    into taxon_dir/Enrichment_Analysis and returns a dict summarizing the
    outcome (for the overall run manifest).
    """
    output_suffix = OUTPUT_SUFFIX_MAP.get(annotation_file, "_from_unknown")
    mode = 'genus' if taxon_dir.name.startswith('Genus_') else 'species'

    output_dir = taxon_dir / "Enrichment_Analysis"
    output_dir.mkdir(exist_ok=True)
    excel_path = output_dir / f"{taxon_dir.name}_enrichment{output_suffix}.xlsx"

    combo_record = {
        "taxon": taxon_dir.name,
        "annotation_file": annotation_file,
        "background_file": background_file,
        "random_seed": random_seed,
        "status": None,
        "output_excel": str(excel_path),
        "had_evidence_group_id": None,
        "error": None,
    }

    if skip_existing and excel_path.is_file():
        run_log(f"    -> [SKIP] {excel_path.name} already exists.")
        combo_record["status"] = "skipped_existing"
        return combo_record

    try:
        paths = {
            "gene_list": taxon_dir / "BLAST_specify" / f"result_{mode}.txt",
            "background_list": sample_dir / background_file,
            "go_annot": sample_dir / annotation_file,
            "ipr_annot": sample_dir / "ipr_gene.tsv",
            "go_obo": project_root / "data" / "go.obo",
            "ipr_entry": project_root / "data" / "entry.list.txt",
        }
        missing = [str(p) for p in paths.values() if not p.is_file()]
        if missing:
            raise FileNotFoundError(f"Missing required file(s): {missing}")

        gene_list = pd.read_csv(paths["gene_list"], header=None)[0].tolist()
        background_list = pd.read_csv(paths["background_list"], header=None)[0].tolist()

        go_annot_df, had_eg = load_go_annotation_file(paths["go_annot"], log=run_log)
        combo_record["had_evidence_group_id"] = had_eg

        ipr_annot_df = pd.read_csv(paths["ipr_annot"], sep='\t', names=["IPR_ID", "Gene"], skiprows=1)

        go_dag = GODag(str(paths["go_obo"]))
        go_annot_df['Namespace'] = go_annot_df['GO_ID'].map(
            lambda go_id: go_dag[go_id].namespace if go_id in go_dag else None)
        go_annot_df.dropna(subset=['Namespace'], inplace=True)

        run_log("    [INFO] Running GO BP enrichment...")
        go_bp = go_annot_df[go_annot_df['Namespace'] == 'biological_process']
        res_bp_raw = run_go_enrichment(gene_list, background_list, go_bp, go_dag,
                                        random_seed=random_seed, log=run_log, n_jobs=n_jobs)
        res_bp_filtered = filter_terminal_go(res_bp_raw.copy(), go_dag)
        if not res_bp_raw.empty:
            res_bp_raw.loc[:, 'Category'] = 'BP'
        if not res_bp_filtered.empty:
            res_bp_filtered.loc[:, 'Category'] = 'BP'

        run_log("    [INFO] Running GO MF enrichment...")
        go_mf = go_annot_df[go_annot_df['Namespace'] == 'molecular_function']
        res_mf_raw = run_go_enrichment(gene_list, background_list, go_mf, go_dag,
                                        random_seed=random_seed, log=run_log, n_jobs=n_jobs)
        res_mf_filtered = filter_terminal_go(res_mf_raw.copy(), go_dag)
        if not res_mf_raw.empty:
            res_mf_raw.loc[:, 'Category'] = 'MF'
        if not res_mf_filtered.empty:
            res_mf_filtered.loc[:, 'Category'] = 'MF'

        run_log("    [INFO] Running GO CC enrichment...")
        go_cc = go_annot_df[go_annot_df['Namespace'] == 'cellular_component']
        res_cc_raw = run_go_enrichment(gene_list, background_list, go_cc, go_dag,
                                        random_seed=random_seed, log=run_log, n_jobs=n_jobs)
        res_cc_filtered = filter_terminal_go(res_cc_raw.copy(), go_dag)
        if not res_cc_raw.empty:
            res_cc_raw.loc[:, 'Category'] = 'CC'
        if not res_cc_filtered.empty:
            res_cc_filtered.loc[:, 'Category'] = 'CC'

        run_log("    [INFO] Running InterPro enrichment...")
        res_ipr_df = run_ipr_enrichment(gene_list, background_list, ipr_annot_df)
        if not res_ipr_df.empty:
            ipr_entry_df = pd.read_csv(paths["ipr_entry"], sep='\t', header=None,
                                        names=["IPR_ID", "ENTRY_TYPE", "ENTRY_NAME"])
            ipr_id_to_name = ipr_entry_df.set_index('IPR_ID')['ENTRY_NAME']
            res_ipr_df['Description'] = res_ipr_df['ID'].map(ipr_id_to_name).fillna("Unknown")

        results_to_save = {
            "GO_BP_filtered": res_bp_filtered, "GO_MF_filtered": res_mf_filtered, "GO_CC_filtered": res_cc_filtered,
            "ALL_GO_filtered": pd.concat([res_bp_filtered, res_mf_filtered, res_cc_filtered], ignore_index=True),
            "IPR": res_ipr_df,
            "GO_BP_raw": res_bp_raw, "GO_MF_raw": res_mf_raw, "GO_CC_raw": res_cc_raw,
            "ALL_GO_raw": pd.concat([res_bp_raw, res_mf_raw, res_cc_raw], ignore_index=True),
        }

        run_log(f"    [INFO] Writing results to {output_dir}...")
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            for sheet_name in sorted(results_to_save.keys()):
                df = results_to_save[sheet_name]
                if not df.empty:
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
        run_log(f"      - Saved Excel file: {excel_path.name}")

        for name, df in results_to_save.items():
            if not df.empty:
                current_suffix = output_suffix if "GO" in name else ""
                tsv_path = output_dir / f"{name}{current_suffix}.tsv"
                df.to_csv(tsv_path, sep='\t', index=False)
                run_log(f"      - Saved TSV file: {tsv_path.name}")

        combo_record["status"] = "success"
        combo_record["n_significant_bp"] = int(len(res_bp_filtered))
        combo_record["n_significant_mf"] = int(len(res_mf_filtered))
        combo_record["n_significant_cc"] = int(len(res_cc_filtered))
        combo_record["n_significant_ipr"] = int(len(res_ipr_df))

    except Exception as e:
        import traceback
        run_log(f"    [ERROR] {taxon_dir.name} / {annotation_file}: {e}")
        run_log(traceback.format_exc())
        combo_record["status"] = "failed"
        combo_record["error"] = str(e)

    return combo_record


def main():
    parser = argparse.ArgumentParser(
        description="Run GO and InterPro enrichment analysis (production / non-interactive version).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--base-dir', type=str, required=True,
                        help="Project root directory (contains 'data/' and the sample directory).")
    parser.add_argument('--sample-dir', type=str, required=True,
                        help="Sample directory name under --base-dir (e.g. 'level-7').")
    parser.add_argument('--taxa', type=str, required=True,
                        help="Comma-separated list of taxon directory names (e.g. 'Species_schaedleri,Genus_Bacteroides'), "
                             "or 'all' to process every Genus_*/Species_* directory found.")
    parser.add_argument('--annotation-file', type=str, required=True,
                        choices=list(ANNOTATION_BACKGROUND_MAP.keys()) + ["all"],
                        help="Which *_go_gene.tsv file to use, or 'all' to run every "
                             "annotation/background combination (for sensitivity analysis).")
    parser.add_argument('--background-file', type=str, default=None,
                        help="Which *background_list.txt file to use. If omitted, the natural "
                             "match for --annotation-file is used automatically "
                             "(e.g. ncbi_go_gene.tsv -> ncbi_background_list.txt, "
                             "GMM_go_gene.tsv -> GMM_background_list.txt). Ignored if "
                             "--annotation-file all is given.")
    parser.add_argument('--random-seed', type=int, default=0,
                        help="Random seed for the fold-enrichment bootstrap CI (for reproducibility).")
    parser.add_argument('--skip-existing', action='store_true',
                        help="Skip a (taxon, annotation_file) combination if its output Excel "
                             "file already exists.")
    parser.add_argument('--manifest-path', type=str, default=None,
                        help="Where to write the run manifest JSON. Defaults to "
                             "<sample-dir>/enrichment_run_manifest_<timestamp>.json")
    parser.add_argument('--cpu', type=int, default=os.cpu_count() or 1,
                        help="Number of CPU cores to use for the fold-enrichment bootstrap CI "
                             "computation (parallelized across GO terms via multiprocessing, "
                             "Linux fork-based). The hypergeometric test itself is fast and "
                             "unaffected. Defaults to all detected cores. Use 1 to disable "
                             "parallelization.")

    args = parser.parse_args()

    project_root = Path(args.base_dir).resolve()
    sample_dir = project_root / args.sample_dir
    if not sample_dir.is_dir():
        sys.exit(f"Error: sample directory not found: {sample_dir}")

    if args.taxa.strip().lower() == 'all':
        taxon_dirs = sorted([d for d in sample_dir.iterdir()
                              if d.is_dir() and (d.name.startswith('Genus_') or d.name.startswith('Species_'))])
    else:
        requested = [t.strip() for t in args.taxa.split(',') if t.strip()]
        taxon_dirs = []
        for t in requested:
            d = sample_dir / t
            if not d.is_dir():
                sys.exit(f"Error: requested taxon directory not found: {d}")
            taxon_dirs.append(d)

    if not taxon_dirs:
        sys.exit("Error: no taxon directories to process.")

    if args.annotation_file == "all":
        combinations = list(ANNOTATION_BACKGROUND_MAP.items())
    else:
        bg_file = args.background_file or ANNOTATION_BACKGROUND_MAP[args.annotation_file]
        combinations = [(args.annotation_file, bg_file)]

    for annotation_file, background_file in combinations:
        if not (sample_dir / annotation_file).is_file():
            sys.exit(f"Error: annotation file not found: {sample_dir / annotation_file}")
        if not (sample_dir / background_file).is_file():
            sys.exit(f"Error: background file not found: {sample_dir / background_file}")

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = Path(args.manifest_path) if args.manifest_path else \
        sample_dir / f"enrichment_run_manifest_{timestamp}.json"

    log_lines = []

    def run_log(msg):
        print(msg)
        log_lines.append(msg)

    run_log(f"[INFO] Run started at {timestamp} (UTC)")
    run_log(f"[INFO] Sample directory: {sample_dir}")
    run_log(f"[INFO] Taxa to process ({len(taxon_dirs)}): {[d.name for d in taxon_dirs]}")
    run_log(f"[INFO] Annotation/background combinations: {combinations}")
    run_log(f"[INFO] random_seed={args.random_seed}, skip_existing={args.skip_existing}, cpu={args.cpu}")

    all_results = []
    for annotation_file, background_file in combinations:
        for taxon_dir in taxon_dirs:
            run_log(f"\n{'='*20} Processing: {taxon_dir.name} "
                    f"[{annotation_file} / {background_file}] {'='*20}")
            combo_record = process_one_combination(
                sample_dir=sample_dir,
                project_root=project_root,
                taxon_dir=taxon_dir,
                annotation_file=annotation_file,
                background_file=background_file,
                random_seed=args.random_seed,
                skip_existing=args.skip_existing,
                run_log=run_log,
                n_jobs=args.cpu,
            )
            all_results.append(combo_record)

    n_success = sum(1 for r in all_results if r["status"] == "success")
    n_failed = sum(1 for r in all_results if r["status"] == "failed")
    n_skipped = sum(1 for r in all_results if r["status"] == "skipped_existing")
    run_log(f"\n[INFO] Run complete: {n_success} succeeded, {n_failed} failed, {n_skipped} skipped "
            f"(of {len(all_results)} total combinations).")

    manifest = {
        "timestamp_utc": timestamp,
        "base_dir": str(project_root),
        "sample_dir": str(sample_dir),
        "taxa_requested": args.taxa,
        "taxa_processed": [d.name for d in taxon_dirs],
        "annotation_background_combinations": combinations,
        "random_seed": args.random_seed,
        "cpu": args.cpu,
        "skip_existing": args.skip_existing,
        "package_versions": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "summary": {"success": n_success, "failed": n_failed, "skipped": n_skipped},
        "results": all_results,
        "log": log_lines,
    }
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    run_log(f"[INFO] Run manifest written to: {manifest_path}")

    if n_failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
