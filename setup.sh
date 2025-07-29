#!/bin/bash
#
# BLINGパイプラインの環境をセットアップするスクリプト
#

# --- 設定 (この部分をあなたの環境に合わせて変更してください) ---
# あなたのGitHubユーザー名/リポジトリ名を指定
GH_USER="schaedleri"
GH_REPO="BLING"
# ---------------------------------------------------------

# GitHub RawのベースURL
BASE_URL="https://raw.githubusercontent.com/${GH_USER}/${GH_REPO}/main"

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
    "go.obo"
)

# InterProScanの情報
IPRSCAN_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz"
IPRSCAN_FILE="interproscan-5.75-106.0-64-bit.tar.gz"
TOOLS_DIR="tools"


# --- セットアップ開始 ---
echo "セットアップを開始します..."
echo "リポジトリ: https://github.com/${GH_USER}/${GH_REPO}"
echo ""

# 1. ディレクトリの作成
echo "ディレクトリを作成します: data, scripts, done, tools"
mkdir -p data scripts done "${TOOLS_DIR}"
echo ""

# 2. ファイルのダウンロード
echo "ファイルをダウンロードします..."

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

# 4. 実行権限の付与
echo "スクリプトに実行権限を付与します..."
chmod +x *.pl *.py
chmod +x scripts/*.pl
echo ""

# --- 完了 ---
echo "✅ セットアップが完了しました。"
