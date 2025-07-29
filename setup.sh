#!/bin/bash
#
# BLINGパイプラインの環境をセットアップするスクリプト
#

# --- 設定 (この部分をあなたの環境に合わせて変更してください) ---
# あなたのGitHubユーザー名/リポジトリ名を指定
GH_USER="schaedleri"
GH_REPO="BLING"
RELEASE_TAG="v1.0.0" # DBファイルなどを置いているGitHub Releaseのタグ名
# ---------------------------------------------------------

# GitHub RawのベースURL
BASE_URL="https://raw.githubusercontent.com/${GH_USER}/${GH_REPO}/master"

# ダウンロードするファイルリスト
ROOT_FILES=(
    "run_route_blastcluster.pl"
    "run_route_interpro.pl"
    "run_route_preprocess.pl"
    "get_level7.py"
)
SCRIPT_FILES=(
    "microbiome_assembly.pl"
    "microbiome_cd-hit.pl"
    "microbiome_cdhit_blastscreening.pl"
    "microbiome_externalblast.pl"
    "microbiome_extract.pl"
    "microbiome_get_datasets.pl"
    "microbiome_interval.pl"
    "microbiome_merge_GO.pl"
    "microbiome_run_interpro.pl"
)
DATA_FILES=(
    "go.obo"
    "taxonomy.tsv"
)
DONE_FILES=(
    "setup.done"
)

# 外部ツールの情報
IPRSCAN_FILE="interproscan-5.75-106.0-64-bit.tar.gz"
TOOLS_DIR="tools"

# BLAST DBの情報
BLAST_DB_FILE="bacteria_strain_taxid_DB.tar.gz"
BLAST_DB_DIR="DB"


# --- セットアップ開始 ---
echo "セットアップを開始します..."
echo "リポジトリ: https://github.com/${GH_USER}/${GH_REPO}"
echo ""

# 1. ディレクトリの作成
echo "ディレクトリを作成します: data, scripts, done, tools, DB"
mkdir -p data scripts done "${TOOLS_DIR}" "${BLAST_DB_DIR}"
echo ""

# 2. スクリプト等のダウンロード
echo "スクリプトとデータをダウンロードします..."

# ルートファイルのダウンロード
for file in "${ROOT_FILES[@]}"; do
    echo "  - ${file}"
    curl -fsSL -o "${file}" "${BASE_URL}/${file}"
    if [ $? -ne 0 ]; then
        echo "エラー: ${file} のダウンロードに失敗しました。"
        exit 1
    fi
done

# scriptファイルのダウンロード
for file in "${SCRIPT_FILES[@]}"; do
    echo "  - scripts/${file}"
    curl -fsSL -o "scripts/${file}" "${BASE_URL}/scripts/${file}"
    if [ $? -ne 0 ]; then
        echo "エラー: scripts/${file} のダウンロードに失敗しました。"
        exit 1
    fi
done

# dataファイルのダウンロード
for file in "${DATA_FILES[@]}"; do
    echo "  - data/${file}"
    curl -fsSL -o "data/${file}" "${BASE_URL}/data/${file}"
    if [ $? -ne 0 ]; then
        echo "エラー: data/${file} のダウンロードに失敗しました。"
        exit 1
    fi
done

# doneファイルのダウンロード
for file in "${DONE_FILES[@]}"; do
    echo "  - done/${file}"
    curl -fsSL -o "done/${file}" "${BASE_URL}/done/${file}"
    if [ $? -ne 0 ]; then
        echo "エラー: done/${file} のダウンロードに失敗しました。"
        exit 1
    fi
done
echo ""

# 3. ツール(InterProScan)のダウンロードと展開
echo "ツール(InterProScan)をダウンロードします... (サイズが大きいため時間がかかります)"
IPRSCAN_URL="https://github.com/${GH_USER}/${GH_REPO}/releases/download/${RELEASE_TAG}/${IPRSCAN_FILE}"
curl -L "${IPRSCAN_URL}" -o "${TOOLS_DIR}/${IPRSCAN_FILE}"
if [ $? -ne 0 ]; then
    echo "エラー: InterProScan のダウンロードに失敗しました。"
    exit 1
fi

echo "InterProScanを展開します..."
tar -xzvf "${TOOLS_DIR}/${IPRSCAN_FILE}" -C "${TOOLS_DIR}/"
if [ $? -ne 0 ]; then
    echo "エラー: InterProScan の展開に失敗しました。"
    exit 1
fi

echo "ダウンロードした圧縮ファイルを削除します..."
rm "${TOOLS_DIR}/${IPRSCAN_FILE}"
echo ""

# 4. BLAST DBのダウンロードと展開
echo "BLAST DBをダウンロードします..."
DB_URL="https://github.com/${GH_USER}/${GH_REPO}/releases/download/${RELEASE_TAG}/${BLAST_DB_FILE}"
curl -L "${DB_URL}" -o "${BLAST_DB_DIR}/${BLAST_DB_FILE}"
if [ $? -ne 0 ]; then
    echo "エラー: BLAST DB のダウンロードに失敗しました。"
    exit 1
fi

echo "BLAST DBを展開します..."
tar -xzvf "${BLAST_DB_DIR}/${BLAST_DB_FILE}" -C "${BLAST_DB_DIR}/"
if [ $? -ne 0 ]; then
    echo "エラー: BLAST DB の展開に失敗しました。"
    exit 1
fi

echo "ダウンロードした圧縮ファイルを削除します..."
rm "${BLAST_DB_DIR}/${BLAST_DB_FILE}"
echo ""

# 5. 実行権限の付与
echo "スクリプトに実行権限を付与します..."
chmod +x *.pl *.py
chmod +x scripts/*.pl
echo ""

# --- 完了 ---
echo "✅ セットアップが完了しました。"
