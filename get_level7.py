import pandas as pd
from biom import load_table
from ete3 import NCBITaxa
import re

# --- 設定項目 ---
BIOM_FILE = 'exported-feature-table/feature-table.biom'
TAXONOMY_FILE = 'exported-taxonomy/taxonomy.tsv'
OUTPUT_FILE = 'final_aggregated_summary_table_final.tsv'
# ----------------

# --- 関数（get_taxid_hierarchicallyを大幅に改良）---
def get_taxid_hierarchically(taxon_string, ncbi_db):
    """
    QIIMEの分類文字列からNCBI TaxIDを取得する。
    種レベルの検索では、属名と種小名を結合して完全な学名を作成する。
    """
    skip_keywords = ['uncultured', 'unidentified', 'metagenome', 'bacterium', 'sp.']
    tax_levels = [part.strip() for part in taxon_string.split(';')]

    # 1. ★★★ 種レベルの検索を特別に処理 ★★★
    genus_name, species_epithet = None, None
    for level in tax_levels:
        if level.startswith('g__'):
            genus_name = level[3:].strip()
        elif level.startswith('s__'):
            species_epithet = level[3:].strip()
            
    # 属名と種小名の両方が存在する場合、完全な学名で検索
    if genus_name and species_epithet and all(kw not in species_epithet.lower() for kw in skip_keywords):
        full_species_name = f"{genus_name.replace('_', ' ')} {species_epithet.replace('_', ' ')}"
        name2taxid = ncbi_db.get_name_translator([full_species_name])
        if full_species_name in name2taxid:
            return name2taxid[full_species_name][0]

    # 2. 種で見つからない場合、属以上の階層を検索（これまでのロジック）
    ranks = ['g__', 'f__', 'o__', 'c__', 'p__']
    for rank_prefix in ranks:
        for level in reversed(tax_levels):
            if level.startswith(rank_prefix):
                name = level[len(rank_prefix):].strip()
                if not name or any(keyword in name.lower() for keyword in skip_keywords): continue
                search_terms = set()
                search_terms.add(name.replace('_', ' ')); search_terms.add(name)
                match = re.search(r'\[(.+?)\]_([^_]+)', name)
                if match:
                    bracket_part, second_part = match.group(1), match.group(2)
                    search_terms.add(f"{bracket_part} {second_part}"); search_terms.add(f"{bracket_part}_{second_part}")
                    search_terms.add(f"[{bracket_part}] {second_part}"); search_terms.add(f"[{bracket_part}]_{second_part}")
                if '_' in name:
                    base_name = name.split('_', 1)[0].replace('[', '').replace(']', '')
                    if base_name: search_terms.add(base_name)
                for term in sorted(list(search_terms), key=len, reverse=True):
                    term_cleaned = term.replace(' group', '').replace(' family', '').strip()
                    if not term_cleaned or len(term_cleaned) < 3: continue
                    name2taxid = ncbi_db.get_name_translator([term_cleaned])
                    if term_cleaned in name2taxid: return name2taxid[term_cleaned][0]
                break
    return None

def get_new_ncbi_lineage(taxid, ncbi_db):
    ranks_to_fill = ['Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    if pd.isna(taxid): return {rank: pd.NA for rank in ranks_to_fill}
    try:
        lineage = ncbi_db.get_lineage(taxid)
        ranks, names = ncbi_db.get_rank(lineage), ncbi_db.get_taxid_translator(lineage)
        lineage_dict = {rank: pd.NA for rank in ranks_to_fill}
        for tax_id in lineage:
            rank, name = ranks.get(tax_id), names.get(tax_id)
            if rank == 'superkingdom': lineage_dict['Domain'] = name
            elif rank == 'kingdom': lineage_dict['Kingdom'] = name
            elif rank == 'phylum': lineage_dict['Phylum'] = name
            elif rank == 'class': lineage_dict['Class'] = name
            elif rank == 'order': lineage_dict['Order'] = name
            elif rank == 'family': lineage_dict['Family'] = name
            elif rank == 'genus': lineage_dict['Genus'] = name
            elif rank == 'species': lineage_dict['Species'] = name
        return lineage_dict
    except: return {rank: pd.NA for rank in ranks_to_fill}

# --- メイン処理 ---
# （メイン処理部分は変更がないため、簡潔に表示します。実際には前回のコードをそのままお使いください）
print("--- 最終処理開始 (種名検索 精度向上版) ---")
# 1. データの読み込み
print("1/5: Loading data...")
table = load_table(BIOM_FILE)
feature_df = table.to_dataframe(dense=True).astype(int)
taxonomy_df = pd.read_csv(TAXONOMY_FILE, sep='\t', index_col=0, header=0, names=['Feature ID', 'Taxon', 'Confidence'])
# 2. リード数の集計
print("2/5: Aggregating counts by taxon string...")
merged_df = feature_df.join(taxonomy_df['Taxon'])
merged_df.dropna(subset=['Taxon'], inplace=True)
aggregated_df = merged_df.groupby('Taxon').sum()
# 3. TaxIDの紐付けと最新分類階層の再構築
print("3/5: Getting TaxIDs and reconstructing lineage...")
ncbi_db = NCBITaxa()
unique_taxa = aggregated_df.index.to_series()
aggregated_df['taxID'] = unique_taxa.apply(lambda x: get_taxid_hierarchically(x, ncbi_db))
new_lineages = pd.DataFrame(aggregated_df['taxID'].apply(lambda x: get_new_ncbi_lineage(x, ncbi_db)).tolist(), index=aggregated_df.index)
# 4. 全データの結合
print("4/5: Creating pre-final table...")
pre_final_df = new_lineages.join(aggregated_df)
pre_final_df['ID'] = pre_final_df.index
# 5. 最終集計ステップ
print("5/5: Final aggregation...")
grouping_columns = ['Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'taxID']
sample_columns = feature_df.columns.tolist()
pre_final_df[grouping_columns] = pre_final_df[grouping_columns].fillna('Unassigned')
agg_rules = {col: 'sum' for col in sample_columns}
agg_rules['ID'] = lambda x: ' ; '.join(x.unique())
final_df = pre_final_df.groupby(grouping_columns).agg(agg_rules).reset_index()
# 列の並び替えと保存
output_columns = ['ID', 'Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'] + sample_columns + ['taxID']
final_df = final_df[output_columns]
final_df.to_csv(OUTPUT_FILE, sep='\t', index=False, na_rep='')

print(f"--- 全処理完了 ---\n最終的な集計・統合テーブルを '{OUTPUT_FILE}' に保存しました。")
