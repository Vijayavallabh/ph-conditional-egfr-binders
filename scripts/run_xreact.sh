#!/bin/bash
# Cross-reactivity for final candidates: mutate target to mouse DIII + human-repack control, af2ig both.
# Usage: run_xreact.sh <cand_dir> <tag> <gpus>
set -uo pipefail
REPO=.; cd "$REPO"
export TMPDIR="$REPO/runs/tmp"; mkdir -p "$TMPDIR"; export PATH="$REPO/envs/pyrosetta/bin:$PATH"
CAND=$1; TAG=${2:-fin}; GPUS=${3:-1-7}
PR="$REPO/envs/pyrosetta/bin/python"; MD="$REPO/envs/md/bin/python"
echo "=== [$TAG] mouse-DIII mutate+repack ==="
"$PR" scripts/mouse_mutate.py --indir "$CAND" --outdir "runs/xreact_${TAG}_mouse" --workers 24
echo "=== [$TAG] human-repack control ==="
"$PR" scripts/mouse_mutate.py --indir "$CAND" --outdir "runs/xreact_${TAG}_humanctrl" --workers 24 --control
echo "=== af2ig mouse + human control ==="
"$MD" scripts/run_af2ig.py --pdbdir "runs/xreact_${TAG}_mouse" --out "runs/af2ig_xreact_${TAG}_mouse" --gpus "$GPUS" --mode fast
"$MD" scripts/run_af2ig.py --pdbdir "runs/xreact_${TAG}_humanctrl" --out "runs/af2ig_xreact_${TAG}_human" --gpus "$GPUS" --mode fast
echo "XREACT_${TAG}_DONE"
