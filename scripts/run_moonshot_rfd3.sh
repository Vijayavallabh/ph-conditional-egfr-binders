#!/bin/bash
# WAVE 22 MOONSHOT: 8-GPU RFd3 generation of de-novo binders carrying a binder-His that
# salt-bridges a CONSERVED (human=mouse) antigen Asp on EGFR Domain III (inverted pH switch;
# removes the His433-pKa dependence). Two hotspots: Asp347 (A323, exposed) + Asp460 (A436, central).
# Each GPU runs the combined spec with a unique seed -> distinct designs; outputs per-GPU to avoid
# name collisions. Binder=chain A, target=chain B in the output CIFs.
#   Usage: run_moonshot_rfd3.sh [inputs.json] [n_batches] [diffusion_batch_size]
set -uo pipefail
REPO=.
INPUTS=${1:-$REPO/specs/rfd3/egfr_hisasp_combined.json}
NBATCH=${2:-6}
BSIZE=${3:-16}
GEN=${4:-$REPO/runs/moonshot/gen}          # out_dir (default = WAVE 22)
GPUS=${5:-"0 1 2 3 4 5 6 7"}               # GPU list (default = all 8)
mkdir -p "$GEN"
echo "=== MOONSHOT RFd3 gen: inputs=$INPUTS n_batches=$NBATCH bsize=$BSIZE out=$GEN gpus=[$GPUS]  $(date -u +%FT%TZ) ==="
pids=()
for i in $GPUS; do
  CUDA_VISIBLE_DEVICES=$i \
  PATH=$REPO/envs/uplifting-rfd3-kit/bin:$PATH \
  RFDIFFUSION3_STOCK_PYTHON=$REPO/envs/uplifting-rfd3-stock/bin/python \
  RFD3_CKPT=$REPO/weights/rfd3/rfd3_latest.ckpt \
  MODEL_OPT_JIT_ROOT=$REPO/.jit \
  bash $REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion3/run.sh design --config h100 --mode fast \
    inputs="$INPUTS" out_dir="$GEN/gpu$i" seed=$((100*(i+1))) \
    diffusion_batch_size=$BSIZE n_batches=$NBATCH > "$GEN/gpu$i.log" 2>&1 &
  pids+=($!)
  echo "launched GPU $i pid ${pids[-1]} -> $GEN/gpu$i"
done
echo "all 8 launched; waiting..."
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
n=$(find "$GEN" -name '*.cif.gz' 2>/dev/null | wc -l)
echo "MOONSHOT_RFD3_DONE $(date -u +%FT%TZ)  structures=$n  failed_shards=$fail"