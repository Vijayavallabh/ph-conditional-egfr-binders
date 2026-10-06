#!/bin/bash
# Fold-validate designed binder-target complexes with AF2 initial-guess on ONE GPU.
# Baker-lab standard de-novo binder filter: pae_interaction (<10 good), plddt_binder (>80),
# binder_aligned_rmsd (<2A = the design was realised). pH-AGNOSTIC: confirms the binder, NOT
# the pH switch (that is PROPKA linkage, scripts/obj1_switch.py).
#
# usage: fold_af2ig.sh GPU PDBDIR OUTDIR [MODE]
set -uo pipefail
R=.
if [ $# -lt 3 ]; then echo "usage: fold_af2ig.sh GPU PDBDIR OUTDIR [MODE]" >&2; exit 2; fi
GPU=$1; PDBDIR=$(readlink -f "$2"); OUT=$(mkdir -p "$3" && readlink -f "$3"); MODE=${4:-exact}  # absolute: run.sh cd's into the kit dir
unset VIRTUAL_ENV 2>/dev/null || true
# shellcheck disable=SC1091
source "$R/envs/uplifting-af2ig/bin/activate"
set -a; source "$R/envs/uplifting.env"; set +a
export AF2IG_DIR="$R/third_party/uplifting-biomolecular-modeling/af2ig/dl_binder_design/af2_initial_guess"
export CUDA_VISIBLE_DEVICES=$GPU
export MODEL_OPT_JIT_ROOT="$R/.jit"
mkdir -p "$OUT"
cd "$R/third_party/uplifting-biomolecular-modeling/af2ig" || exit 1
echo "[af2ig] gpu=$GPU mode=$MODE pdbdir=$PDBDIR out=$OUT start=$(date -u +%FT%TZ)"
bash run.sh pred --config h100 --mode "$MODE" --pdbdir "$PDBDIR" --out "$OUT"
rc=$?
echo "[af2ig] DONE rc=$rc end=$(date -u +%FT%TZ)"
exit $rc
