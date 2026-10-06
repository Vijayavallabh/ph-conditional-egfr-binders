#!/bin/bash
# WAVE 22c: RFd1 PARTIAL DIFFUSION of moon01 (egfr_phsw_inv02 parent) to drop its novelty (oracle TM
# 0.804 -> <0.80 with margin) while preserving binding + the inverted binder-His -> conserved-Asp460
# switch. Shards a partial_T sweep across 8 GPUs (endgame cap lifted). Each GPU gets its own out_dir +
# schedule dir (RFd1 races on a shared schedule). Downstream re-locks the His by PROXIMITY to Asp460
# (RFd1 is backbone-only -> no His survives diffusion), via moonshot_prep_partial.py.
# Usage: run_pd_inv02.sh [ndes_per_gpu]
set -uo pipefail
REPO=.; cd "$REPO"
INP="$REPO/runs/moonshot/pd_inv02/moon01_complex.pdb"
NDES=${1:-96}
[ -s "$INP" ] || { echo "run_pd_inv02: input complex missing: $INP" >&2; exit 2; }
# GPU -> partial_T map: 0,1 gentle(12); 2,3,4 moderate(18); 5,6,7 aggressive(24)
declare -A T=( [0]=12 [1]=12 [2]=18 [3]=18 [4]=18 [5]=24 [6]=24 [7]=24 )
pids=()
for g in 0 1 2 3 4 5 6 7; do
  pt=${T[$g]}
  out="pd_inv02/gpu${g}_T${pt}"
  sched="$REPO/runs/rfd/pd_inv02/sched_gpu${g}"
  bash scripts/run_rfd_partial.sh "$out" "$INP" "$NDES" "$g" "$pt" "$sched" \
      > "$REPO/runs/rfd/pd_inv02_gpu${g}.log" 2>&1 &
  pids+=($!)
  sleep 1
done
echo "[pd_inv02] launched ${#pids[@]} partial-diffusion shards, ${NDES}/gpu, T={12,18,24} $(date -u +%FT%TZ)"
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
tot=$(find "$REPO/runs/rfd/pd_inv02" -name "des_*.pdb" 2>/dev/null | wc -l)
echo "PD_INV02_DONE rc=$rc backbones=$tot $(date -u +%FT%TZ)"
