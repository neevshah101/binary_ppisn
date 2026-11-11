#!/usr/bin/env bash
# Make one .mp4 per pgstar/pgbinary png series in ./png (e.g. png/grid1_*.png -> grid1.mp4).
# Usage: ./make_movies.sh [png_dir]     (run from the run folder)
set -euo pipefail
cd "$(dirname "$0")"

png_dir=${1:-png}
command -v images_to_movie >/dev/null || { echo "images_to_movie not found (load the MESA SDK)"; exit 1; }
[ -d "$png_dir" ] || { echo "No $png_dir/ directory"; exit 1; }

# series prefixes: filename with trailing digits and .png removed, e.g. grid1_000010.png -> grid1_
prefixes=$(find "$png_dir" -maxdepth 1 -name '*.png' -printf '%f\n' | sed -E 's/[0-9]+\.png$//' | sort -u)
[ -n "$prefixes" ] || { echo "No png files in $png_dir/"; exit 1; }

for p in $prefixes; do
    out="${p%_}.mp4"
    echo "$png_dir/${p}*.png -> $out"
    images_to_movie "$png_dir/${p}*.png" "$out"
done
