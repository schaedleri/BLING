#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import re
import glob
import argparse
from pathlib import Path
from collections import defaultdict

def load_go_graph(obo_file_path):
    """
    Parses a go.obo file and builds the following two dictionaries:
    1. go_graph: A Directed Acyclic Graph (DAG) representing parent-child relationships in GO.
    2. obsolete_map: A dictionary mapping obsolete GO IDs to their replacement or suggested terms.
    """
    go_graph = defaultdict(list)
    obsolete_map = defaultdict(list)
    
    # Define the relationship types to be traced as parent-child relationships.
    PARENT_RELATIONSHIPS = {
        'is_a',
        'part_of',
        'regulates',
        'positively_regulates',
        'negatively_regulates',
        'occurs_in',
        'enables',
        'contributes_to'
    }
    
    try:
        with open(obo_file_path, 'r', encoding='utf-8') as f:
            current_go_id = None
            is_obsolete = False
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line == '[Term]':
                    current_go_id = None
                    is_obsolete = False
                elif line.startswith('id: '):
                    current_go_id = line.split(' ')[1]
                elif line.startswith('is_obsolete: true'):
                    is_obsolete = True
                
                # Process parent-child relationships for non-obsolete terms.
                elif current_go_id and not is_obsolete:
                    parent_go_id = None
                    # Handle the generic format like 'relationship: part_of GO:...'.
                    if line.startswith('relationship: '):
                        parts = line.split(' ')
                        relationship_type = parts[1]
                        if relationship_type in PARENT_RELATIONSHIPS:
                            parent_go_id = parts[2]
                    # Handle the direct format like 'is_a: GO:...'.
                    elif ':' in line:
                        tag, value = line.split(':', 1)
                        if tag in PARENT_RELATIONSHIPS:
                            parent_go_id = value.strip().split(' ')[0]
                    
                    if parent_go_id:
                        go_graph[current_go_id].append(parent_go_id)

                # Process replacement information for obsolete terms.
                elif current_go_id and is_obsolete:
                    if line.startswith('replaced_by: ') or line.startswith('consider: '):
                        replacement_id = line.split(' ')[1].split('!')[0].strip()
                        obsolete_map[current_go_id].append(replacement_id)

    except FileNotFoundError:
        print(f"Error: GO OBO file not found at {obo_file_path}", file=sys.stderr)
        sys.exit(1)
    
    return go_graph, obsolete_map

def get_all_ancestors(go_id, go_graph, memo=None):
    """Recursively finds all ancestor GO terms for a given GO ID."""
    if memo is None: memo = {}
    if go_id in memo: return memo[go_id]
    
    parents = go_graph.get(go_id, [])
    ancestors = set(parents)
    for parent in parents:
        ancestors.update(get_all_ancestors(parent, go_graph, memo))
        
    memo[go_id] = ancestors
    return ancestors

def main():
    """
    Main function to aggregate GO/IPR annotations, expand GO terms to their
    ancestors, and create final summary files for enrichment analysis.
    """
    parser = argparse.ArgumentParser(
        description="Aggregate GO/IPR annotations and create summary files."
    )
    parser.add_argument('--base-dir', type=str, required=True,
                        help="The base directory for a specific sample. Required.")
    args = parser.parse_args()

    sample_dir = Path(args.base_dir).resolve()
    project_root = sample_dir.parent
    go_obo_file = project_root / "data" / "go.obo"
    
    print(f"Parsing GO hierarchy from {go_obo_file}...")
    go_graph, obsolete_map = load_go_graph(go_obo_file)
    print(f"Found relationships for {len(go_graph)} GO terms.")
    print(f"Found {len(obsolete_map)} obsolete GO terms to be replaced.")

    go_gene_pairs = set()
    ipr_gene_pairs = set()

    go_file_patterns = [
        str(sample_dir / "*/go_annotations/*.tsv"),
        str(sample_dir / "*/interval_test/*/IP_GO_pairs.tsv"),
        str(sample_dir / "*/interval_test/*/nonGO_GO_pairs*.tsv")
    ]
    
    print("Collecting GO-gene pairs...")
    ancestor_cache = {}
    for pattern in go_file_patterns:
        for file_path in glob.glob(pattern):
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    # Skip header if it exists
                    first_line = next(f, '').lower()
                    if 'go_id' in first_line: pass
                    else: f.seek(0)
                    
                    for line in f:
                        cols = line.strip().split('\t')
                        gene, go = None, None
                        
                        if '/go_annotations/' in file_path:
                            if len(cols) >= 3 and cols[2].startswith('GO:'): gene, go = cols[0], cols[2]
                        elif file_path.endswith('IP_GO_pairs.tsv'):
                            if len(cols) >= 2 and cols[1].startswith('GO:'): gene, go = cols[0], cols[1]
                        elif 'nonGO_GO_pairs' in file_path:
                            if len(cols) >= 2 and cols[1].startswith('GO:'): gene, go = cols[0], cols[1]
                        
                        if gene and go:
                            # Resolve obsolete GO terms
                            go_terms_to_process = set()
                            if go in obsolete_map:
                                go_terms_to_process.update(obsolete_map[go])
                            else:
                                go_terms_to_process.add(go)

                            # Process each valid/replaced GO term and its ancestors
                            for term in go_terms_to_process:
                                go_gene_pairs.add(f"{term}\t{gene}")
                                ancestors = get_all_ancestors(term, go_graph, ancestor_cache)
                                for ancestor in ancestors:
                                    go_gene_pairs.add(f"{ancestor}\t{gene}")
            except Exception as e:
                print(f"Warning: Could not process file {file_path}: {e}", file=sys.stderr)

    print("Collecting IPR-gene pairs...")
    interpro_files = glob.glob(str(sample_dir / "*/interpro/*.tsv"))
    for file_path in interpro_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    cols = line.strip().split('\t')
                    if len(cols) > 11:
                        gene, ipr_field = cols[0], cols[11]
                        if ipr_field and ipr_field != '-':
                            ipr_terms = re.findall(r'IPR\d+', ipr_field)
                            for ipr in ipr_terms:
                                ipr_gene_pairs.add(f"{ipr}\t{gene}")
        except Exception as e:
            print(f"Warning: Could not process file {file_path}: {e}", file=sys.stderr)
            
    # --- Write Output Files ---
    print("Writing output summary files...")
    go_out_path = sample_dir / "go_gene.tsv"
    with open(go_out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(sorted(list(go_gene_pairs))) + '\n')

    ipr_out_path = sample_dir / "ipr_gene.tsv"
    with open(ipr_out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(sorted(list(ipr_gene_pairs))) + '\n')
        
    background_genes = set()
    for pair in go_gene_pairs: background_genes.add(pair.split('\t')[1])
    for pair in ipr_gene_pairs: background_genes.add(pair.split('\t')[1])
        
    bg_out_path = sample_dir / "background_list.txt"
    with open(bg_out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(sorted(list(background_genes))) + '\n')
        
    print(f"\nSuccessfully created summary files in {sample_dir}")
    print(f" - {go_out_path}\n - {ipr_out_path}\n - {bg_out_path}")

if __name__ == "__main__":
    main()