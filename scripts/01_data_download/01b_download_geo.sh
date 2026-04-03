#!/bin/bash
# Download GEO scRNA-seq datasets for thyroid cancer
set -e

BASE="/data2/projects/Carcinoma_of_thyroid/data/geo_scrna"

echo "=== Downloading GSE184362 ==="
cd "$BASE/GSE184362"

# Try GEO supplementary files
wget -q "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184362/suppl/" -O index.html 2>/dev/null || true

# Download using GEO FTP
for f in $(grep -oP 'GSE184362_[^"<]+' index.html 2>/dev/null | sort -u | head -10); do
    echo "Downloading $f..."
    wget -q "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184362/suppl/$f" -O "$f" 2>/dev/null || true
done
rm -f index.html

# If h5ad or h5 files, decompress if needed
for f in *.gz; do
    [ -f "$f" ] && gunzip -f "$f" 2>/dev/null || true
done

echo "GSE184362 files:"
ls -lh "$BASE/GSE184362/" 2>/dev/null || echo "(empty)"

echo ""
echo "=== Downloading GSE191288 ==="
cd "$BASE/GSE191288"
wget -q "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE191nnn/GSE191288/suppl/" -O index.html 2>/dev/null || true
for f in $(grep -oP 'GSE191288_[^"<]+' index.html 2>/dev/null | sort -u | head -10); do
    echo "Downloading $f..."
    wget -q "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE191nnn/GSE191288/suppl/$f" -O "$f" 2>/dev/null || true
done
rm -f index.html

for f in *.gz; do
    [ -f "$f" ] && gunzip -f "$f" 2>/dev/null || true
done

echo "GSE191288 files:"
ls -lh "$BASE/GSE191288/" 2>/dev/null || echo "(empty)"

echo ""
echo "=== GEO downloads complete ==="
