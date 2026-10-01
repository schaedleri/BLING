# BLING

**B**acterial **L**ineage-**I**nformed **N**eighborhood **G**enomics — an integrated computational
framework for predicting the functions of putatively taxon-restricted bacterial proteins by
combining 16S rRNA amplicon-based taxonomic profiling with reference-genome retrieval,
homology-based annotation, and genomic-context (operon-based) functional propagation.

This repository contains the pipeline used in the accompanying manuscript to (1) identify
putatively taxon-restricted proteins within a microbial community, (2) expand functional
annotation coverage for those proteins using sequence homology and genomic context, and
(3) test for functional enrichment among the resulting annotation sets.

## Overview of the pipeline

The pipeline runs in the following stages. Each stage's script(s) are listed with their
required inputs and produced outputs.

```
1. 16S rRNA amplicon processing (QIIME2 / DADA2)
        |
        v
2. Reference genome/proteome retrieval + unified BLAST+ database construction
        |
        v
3. Protein clustering (CD-HIT, -c 0.95 -aS 0.95 -g 1)
        |
        v
4. Taxon-restriction screen (two-step BLASTP, both steps against the same unified database)
   4a. Internal BLASTP (within-genus screen)      -> cdhit_blastscreen.py
   4b. External BLASTP (community-wide screen)    -> microbiome_externalblast_ver2.py
        |
        v
5. Functional annotation of all community proteins
   5a. InterProScan (domain-based annotation)
   5b. eggNOG-mapper (orthology-based annotation)  -> eggnog_process_ver2.py
        |
        v
6. Operon-like cluster prediction + genomic-context GO propagation (Biological Process only)
   6a. GMM-based (statistically derived threshold)      -> predict_and_analyze_operons.py
   6b. Fixed-threshold (100 bp, exploratory alternative) -> microbiome_interval.py
        |
        v
7. Aggregation into annotation + background files
   7a. GMM-side aggregation     -> merge_GO.py
   7b. Fixed-side aggregation (also produces the NCBI-only and direct-evidence-only tiers)
                                -> microbiome_merge_GO.py
        |
        v
8. GO / InterPro enrichment analysis
        -> run_enrichment_production.py (scriptable, CLI-driven)
        -> run_enrichment_interactive.py (interactive, for exploratory use)
```

### Reporting utilities (built on top of the pipeline output)

These are not part of the core annotation pipeline above, but operate on its output to
produce the summary figures and statistics reported in the manuscript.

- `aggregate_gmm_thresholds.py` — collects the per-genome `gmm_threshold_summary.tsv` files
  written under every taxonomic group's `operon_analysis/` directory into one deduplicated,
  community-wide table (one row per unique genome), flagging genomes with a degenerate GMM
  fit (`Optimal_K == 1`) rather than excluding them. Used to produce the genome-wide GMM
  threshold distribution figure.
- `microbiome_summary_ver2.py` / `microbiome_summary_fixed.py` — compute the community-wide
  annotation-coverage statistics (baseline GO coverage, newly annotated proteins contributed
  by InterPro, eggNOG, and operon-based propagation) for the GMM and fixed-threshold methods,
  respectively.

## Directory structure (per sample)

```
<sample_dir>/
├── Genus_<name>/ or Species_<name>/         # one directory per taxonomic group of interest
│   ├── <name>.tsv                           # taxonomy ID list for this group
│   ├── fasta/                               # per-genome FASTA files + concatenated/CD-HIT output,
│   │                                          and the per-protein query FASTA files used for BLASTP
│   ├── CDhit/                               # CD-HIT clustering output, internal BLASTP results
│   │                                          (*_Unified.tsv, one file per query against the
│   │                                          unified database)
│   ├── BLAST_specify/                       # external BLASTP results, taxon-restriction screen output
│   ├── genomic_annotations/                 # per-genome *_genomic.tsv (Protein_ID, Begin, End, Strand, ...)
│   ├── go_annotations/                      # GBFF-derived direct GO annotations
│   ├── interpro/                            # InterProScan output
│   ├── eggnog/                              # eggNOG-mapper output (per-genome extracts)
│   ├── operon_analysis/                     # GMM operon prediction output (per genome), incl.
│   │                                          gmm_threshold_summary.tsv (fitted threshold per genome)
│   └── interval_test/                       # fixed-threshold (100 bp) operon prediction output,
│                                              incl. nonGO_GO_pairs100.tsv (propagated BP terms)
├── ncbi_go_gene.tsv / ncbi_background_list.txt        # Tier 1: NCBI/GBFF annotations only
├── direct_go_gene.tsv / direct_background_list.txt    # Tier 2: + InterPro + eggNOG, no propagation
├── fixed_go_gene.tsv / fixed_background_list.txt       # Tier 3: + fixed-threshold propagation
├── fixed_ipr_gene.tsv
├── GMM_go_gene.tsv / GMM_background_list.txt          # Tier 4: + GMM-based propagation
└── */Enrichment_Analysis/                             # per-taxon enrichment results (Excel/TSV)
```

## Requirements

- Python 3.10+
- `scikit-learn`, `numpy`, `scipy`, `pandas`, `statsmodels`, `goatools`
- BLAST+ (`blastp`, `makeblastdb`, `blastdb_aliastool`)
- CD-HIT
- InterProScan
- eggNOG-mapper (with a local eggNOG database)
- QIIME2 (for the 16S amplicon processing stage)

See `env.yml` for a pinned dependency list.

## Usage

Each stage is run independently via its own script, using `--base-dir` (project root) and
`--target-dir` / taxon-directory arguments to specify the sample and taxonomic group being
processed. Example, for a single taxonomic group:

```bash
# 1. Cluster sequences with CD-HIT
python3 microbiome_cd-hit.py --base-dir /path/to/project/level-7 Species_example

# 2. Internal (within-genus) BLASTP screen, against the unified database
python3 cdhit_blastscreen.py --base-dir /path/to/project/level-7 --cpu 10 Genus_example

# 3. External (community-wide) BLASTP screen, against the same unified database
python3 microbiome_externalblast_ver2.py --base-dir /path/to/project --target-dir level-7 Species_example

# 4. eggNOG-mapper annotation (consolidated across all taxa in the sample)
python3 eggnog_process_ver2.py --base-dir /path/to/project/level-7 --cpu 10

# 5. GMM-based operon prediction and GO propagation (Biological Process only)
python3 predict_and_analyze_operons.py --base-dir /path/to/project/level-7 [options]

# 6. Fixed-threshold (100 bp) operon prediction and GO propagation (Biological Process only)
python3 microbiome_interval.py --base-dir /path/to/project/level-7 [options]

# 7. Aggregate annotation/background files
python3 merge_GO.py --base-dir /path/to/project/level-7            # GMM side
python3 microbiome_merge_GO.py --base-dir /path/to/project/level-7  # fixed side + NCBI/direct tiers

# 8. GO/InterPro enrichment analysis
python3 run_enrichment_production.py \
    --base-dir /path/to/project --sample-dir level-7 \
    --taxa Species_example --annotation-file all --cpu 8

# Reporting: aggregate GMM thresholds across all genomes
python3 aggregate_gmm_thresholds.py \
    --sample-dir /path/to/project/level-7 \
    --output all_genomes_gmm_threshold_summary.tsv

# Reporting: community-wide annotation coverage statistics
python3 microbiome_summary_ver2.py --base-dir /path/to/project/level-7 \
    --background-file /path/to/project/level-7/GMM_background_list.txt
python3 microbiome_summary_fixed.py --base-dir /path/to/project/level-7 \
    --background-file /path/to/project/level-7/fixed_background_list.txt
```

Run any script with `--help` for the full list of options.

## Reproducibility notes

- Both BLASTP steps (internal and external) query the same unified BLAST+ database, built by
  merging the ten reference-genome partitions (BSDB0–BSDB9) into one alias with
  `blastdb_aliastool`. Querying the partitions separately would understate their effective
  e-value (since e-value scales with the size of the searched database), so both steps are
  deliberately kept consistent.
- `seg yes` (low-complexity masking) is applied in both BLASTP steps.
- BLASTP searches use a fixed, lexicographically sorted processing order to ensure
  deterministic output regardless of filesystem enumeration order.
- Genomic-context propagation (both the GMM-based and fixed-threshold methods) is restricted
  to Biological Process (BP) GO terms. Molecular Function and Cellular Component annotations
  from direct sequence evidence (GBFF, InterProScan, eggNOG-mapper) are retained and used in
  enrichment analysis, but are not propagated across operon-like clusters.
- The GMM operon-prediction threshold is fit independently per genome (see
  `gmm_threshold_summary.tsv` under each taxon's `operon_analysis/` directory, or the
  community-wide aggregation produced by `aggregate_gmm_thresholds.py`, for the distribution
  of fitted thresholds across all genomes analyzed). Genomes for which the fit is degenerate
  (a single Gaussian component selected) are flagged rather than silently excluded.
- The permitted COG-category mismatch pairs used during operon cluster refinement are
  hard-coded as the `ALLOWED_MISMATCH_PAIRS` constant in `predict_and_analyze_operons.py`.
- Every propagated GO annotation carries a provenance identifier (`Evidence_Group_ID`) that
  traces it back to either its originating direct-evidence source (GBFF / InterPro / eggNOG)
  or the specific operon-like cluster it was propagated within. Enrichment analysis treats
  each such group as a single evidence unit, so co-propagated annotations are not counted as
  independent observations.

## License

This project is released under the MIT License. See `LICENSE` for details.

## Citation

If you use this pipeline, please cite the accompanying manuscript (citation to be added
upon publication).
