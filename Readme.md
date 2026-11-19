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
   5a. InterProScan (domain-based annotation)      -> microbiome_run_interpro.py
   5b. eggNOG-mapper (orthology-based annotation)  -> eggnog_process_ver2.py
        |
        v
6. Operon-like cluster prediction + genomic-context GO propagation (Biological Process only)
   6a. GMM-based (statistically derived threshold)      -> predict_and_analyze_operons.py
   6b. Fixed-threshold (100 bp, exploratory alternative) -> microbiome_interval.py
        |
        v
7. Aggregation into annotation + background files (GO terms expanded to their ancestors)
   7a. GMM-side aggregation     -> merge_GO.py
   7b. Fixed-side aggregation (also produces the NCBI-only and direct-evidence-only tiers)
                                -> microbiome_merge_GO.py
        |
        v
8. GO / InterPro enrichment analysis
        -> run_enrichment_v2.py
```

### Reporting utilities (built on top of the pipeline output)

These are not part of the core annotation pipeline above, but operate on its output to
produce the summary statistics reported in the manuscript.

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
├── temp_annotations/                                  # consolidated protein FASTA + eggNOG-mapper output
├── ncbi_go_gene.tsv / ncbi_background_list.txt        # Tier 1: NCBI/GBFF annotations only
├── direct_go_gene.tsv / direct_background_list.txt    # Tier 2: + InterPro + eggNOG, no propagation
├── fixed_go_gene.tsv / fixed_background_list.txt       # Tier 3: + fixed-threshold propagation
├── fixed_ipr_gene.tsv
├── GMM_go_gene.tsv / GMM_background_list.txt          # Tier 4: + GMM-based propagation
├── GMM_ipr_gene.tsv
├── enrichment_v2_manifest_<timestamp>.json            # options, log and results of each enrichment run
└── */Enrichment_Analysis/                             # per-taxon enrichment results (Excel/TSV)
```

## Requirements

- Python 3.10+
- `scikit-learn`, `numpy`, `scipy`, `pandas`, `matplotlib`, `statsmodels`, `goatools`, `openpyxl`
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

# 4. InterProScan annotation (all taxa of the sample; identical file names are computed once)
python3 microbiome_run_interpro.py --base-dir /path/to/project/level-7 --cpu 8

# 5. eggNOG-mapper annotation (consolidated across all taxa in the sample)
python3 eggnog_process_ver2.py --base-dir /path/to/project/level-7 --cpu 10

# 6. GMM-based operon prediction and GO propagation (Biological Process only)
python3 predict_and_analyze_operons.py --base-dir /path/to/project/level-7 [options]

# 7. Fixed-threshold (100 bp) operon prediction and GO propagation (Biological Process only)
python3 microbiome_interval.py --base-dir /path/to/project/level-7 [options]

# 8. Aggregate annotation/background files
python3 merge_GO.py --base-dir /path/to/project/level-7            # GMM side (--prefix GMM by default)
python3 microbiome_merge_GO.py --base-dir /path/to/project/level-7  # fixed side + NCBI/direct tiers

# 9. GO/InterPro enrichment analysis
python3 run_enrichment_v2.py \
    --base-dir /path/to/project --sample-dir level-7 \
    --taxa Species_example --annotation-file all


# Reporting: community-wide annotation coverage statistics
python3 microbiome_summary_ver2.py --base-dir /path/to/project/level-7 \
    --background-file /path/to/project/level-7/GMM_background_list.txt
python3 microbiome_summary_fixed.py --base-dir /path/to/project/level-7 \
    --background-file /path/to/project/level-7/fixed_background_list.txt
```

Run any script with `--help` for the full list of options. The main options of the last three stages are:

| Script | Option | Default | Meaning |
|---|---|---|---|
| `predict_and_analyze_operons.py` | `--max-k` | 8 | Largest number of GMM components tested (BIC) |
| | `--outlier-threshold` | 1000 | Intergenic distances >= this value (bp) are excluded from the fit |
| | `--min-threshold` | 20 | Lower bound of the distance threshold (bp) |
| | `--min-distances` | 20 | Minimum number of intergenic distances needed to fit a genome |
| | `--scan-start` | -50 | Smallest distance (bp) scanned when locating the threshold |
| `microbiome_interval.py` | `--interval` | 100 | Maximum gap between adjacent genes of a cluster (bp) |
| | `--max-dist` | 300 | Maximum distance for propagation to an unannotated gene (bp) |
| | `--min-cluster-size` | 3 | Minimum number of genes in a cluster |
| `run_enrichment_v2.py` | `--unit-mode` | `exclusive` | `exclusive`: a gene in a cluster is represented by its cluster only; `legacy`: counted as gene and as cluster |
| | `--propagation-scope` | `BP` | `BP` drops propagated MF/CC rows; `all` keeps them |
| | `--background-mode` | `community` | `community`: background list of the annotation tier; `focal`: all proteins of the focal genome (`--focal-background-file`) |
| | `--min-term-size`, `--max-term-fraction` | 5, 0.1 | Term-size filter in evidence units (n >= 5 and n <= 10% of M) |
| | `--fdr-alpha`, `--n-bootstrap` | 0.05, 2000 | Benjamini-Hochberg threshold; bootstrap replicates for the fold-enrichment interval |
| | `--exclude-genes`, `--drop-propagation-from-genomes` | none | Sensitivity analyses: remove genes from the foreground, or remove the propagated rows of listed genomes |
| | `--worked-example GO:xxxxxxx` | none | Print M, n, N, k and the contributing units for one BP term under both unit definitions |

The output file names of `run_enrichment_v2.py` carry the settings used (e.g. `_uX_sBP_bcom` = exclusive units, BP-only
propagation, community background).

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
- The GMM operon-prediction threshold is fit independently per genome. A one-dimensional GMM
  is fitted to log(d + offset) of the distances between adjacent same-strand genes (distances
  below `--outlier-threshold`); the number of components is chosen by BIC. The intra-operon
  component is the one with the largest weight among the components whose mean does not
  exceed the median of the component means, and the threshold is the largest distance assigned
  to it (never below `--min-threshold`). Adjacency is strict: a gene on the opposite strand
  between two genes breaks their adjacency. `gmm_threshold_summary.tsv` under each taxon's
  `operon_analysis/` directory lists the threshold, K and the basis of the threshold
  (`GMM`, `floor_applied`, `single_component_K1`, `default_too_few_distances`) for every
  genome. Genomes for which a single Gaussian component is selected (K = 1) are flagged in
  this file and excluded from clustering and propagation.
- The permitted COG-category mismatch pairs used during operon cluster refinement are
  hard-coded as the `ALLOWED_MISMATCH_PAIRS` constant in `predict_and_analyze_operons.py`.
- Every GO annotation carries a provenance identifier (`Evidence_Group_ID`): one per gene and
  term for direct evidence (GBFF / InterPro / eggNOG), and one per cluster (`<genome>_cluster<n>`)
  for propagated annotations. Ancestor terms inherit the identifier of the term they were
  expanded from. In the enrichment analysis the statistical unit is the evidence unit: a gene
  that belongs to a propagation cluster is represented by its cluster, every other gene by
  itself, and a unit carries a term if any member carries it. Co-propagated annotations are
  therefore not counted as independent observations, and a gene with both direct and
  propagated evidence for a term is counted once.
- Enrichment is a one-sided hypergeometric test per GO namespace (and for InterPro entries),
  with Benjamini-Hochberg correction within each namespace; the confidence interval of the fold
  enrichment comes from a bootstrap over evidence units.

## License

This project is released under the MIT License. See `LICENSE` for details.

## Citation

If you use this pipeline, please cite the accompanying manuscript (citation to be added
upon publication).
