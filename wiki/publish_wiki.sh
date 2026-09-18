#!/bin/bash
set -e

WIKI_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WIKI_REPO_URL="https://github.com/Riperedo/AnalyticalStructureFactors.jl.wiki.git"
TEMP_DIR="$(mktemp -d)"

echo "========================================================================"
echo "Publishing AnalyticalStructureFactors.jl Documentation to GitHub Wiki"
echo "========================================================================"

echo "Cloning wiki repository from: $WIKI_REPO_URL..."
if git clone "$WIKI_REPO_URL" "$TEMP_DIR" 2>/dev/null; then
    echo "Copying wiki markdown pages..."
    cp "$WIKI_DIR"/*.md "$TEMP_DIR/"
    cd "$TEMP_DIR"
    git add *.md
    if git diff --staged --quiet; then
        echo "No changes to commit. GitHub Wiki is already up to date!"
    else
        git commit -m "docs: Update GitHub Wiki documentation"
        git push origin master 2>/dev/null || git push origin main
        echo "GitHub Wiki successfully published!"
    fi
    rm -rf "$TEMP_DIR"
else
    echo "------------------------------------------------------------------------"
    echo "Note: The GitHub Wiki repository is created on-demand by GitHub."
    echo "To initialize the Wiki on GitHub:"
    echo "  1. Go to https://github.com/Riperedo/AnalyticalStructureFactors.jl/wiki"
    echo "  2. Click 'Create the first page' and click 'Save Page'."
    echo "  3. Re-run this script: bash wiki/publish_wiki.sh"
    echo "------------------------------------------------------------------------"
fi
