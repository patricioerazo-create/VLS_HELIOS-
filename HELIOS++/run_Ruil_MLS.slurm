#!/bin/bash
#SBATCH --job-name=Ruil_MLS
#SBATCH --output=output/slurm_%j.log
#SBATCH --error=output/slurm_%j.log
#SBATCH --time=48:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=16
#SBATCH --partition=main

# Ajustar rutas de Conda/binarios según tu cluster (desde HELIOS++/).
set -euo pipefail
HELIOS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HELIOS_ROOT"
OUT="$HELIOS_ROOT/output"
mkdir -p "$OUT"

module load Miniconda3/4.9.2
HELIOS=/home/perazo/.conda/envs/helios/lib/python3.13/site-packages/pyhelios/bin/helios++
LASMERGE=/home/perazo/.conda/envs/helios/bin/lasmerge
LAS2LAS=/home/perazo/.conda/envs/helios/bin/las2las
echo "=== Ruil_MLS | vox 3.5cm PAD 1.0 | z=1.3m | walk 0.5 m/s $(date) ==="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Node: $(hostname)"

$HELIOS data/surveys/Ruil_MLS.xml \
  --lasOutput \
  --assets "$HELIOS_ROOT" \
  --output "$OUT" \
  --rebuildScene \
  --njobs 16

RUN_DIR="$(find "$OUT/MLS" -mindepth 1 -maxdepth 1 -type d -name '20*' 2>/dev/null | sort | tail -1)"
if [[ -n "$RUN_DIR" ]]; then
  mapfile -t LEGS < <(find "$RUN_DIR" -maxdepth 1 -name 'leg*_points.las' | sort)
  if [[ ${#LEGS[@]} -gt 0 ]]; then
    MERGED="${OUT}/Ruil_MLS.las"
    echo "=== Merge ${#LEGS[@]} legs -> ${MERGED} ==="
    "$LASMERGE" -i "${LEGS[@]}" -o "$MERGED" -verbose
    ls -lh "$MERGED"

    echo "=== Compress merged -> ${OUT}/Ruil_MLS.laz ==="
    "$LAS2LAS" -i "$MERGED" -o "${OUT}/Ruil_MLS.laz" -set_version 1.4
    ls -lh "${OUT}/Ruil_MLS.laz"

    rm -f "$MERGED"
  fi
fi

echo "=== Remove HELIOS run subdirs ==="
rm -rf "$OUT/MLS"

echo "=== Done $(date) ==="
