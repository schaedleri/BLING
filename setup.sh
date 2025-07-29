#!/bin/bash
#
# BLINGパイプラインの環境をセットアップするスクリプト
#

# --- 設定 (この部分をあなたの環境に合わせて変更してください) ---
# あなたのGitHubユーザー名/リポジトリ名を指定
GH_USER="your-username"
GH_REPO="your-repo-name"
# ---------------------------------------------------------

# GitHub RawのベースURL
BASE_URL="https://raw.githubusercontent.com/${GH_USER}/${GH_REPO}/main"

# ダウンロードするファイルリスト
ROOT_FILES=(
    "run_route_blastcluster.pl"
    "run_route_interpro.pl"
    "run_route_preprocess.pl"
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
)

# --- セットアップ開始 ---
echo "セットアップを開始します..."
echo "リポジトリ: https://github.com/${GH_USER}/${GH_REPO}"
echo ""

# 1. ディレクトリの作成
echo "ディレクトリを作成します: data, scripts"
mkdir -p data
mkdir -p scripts
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
echo ""

# 3. 実行権限の付与
echo "Perlスクリプトに実行権限を付与します..."
chmod +x *.pl
chmod +x scripts/*.pl
echo ""

# --- 完了 ---
echo "✅ セットアップが完了しました。"
