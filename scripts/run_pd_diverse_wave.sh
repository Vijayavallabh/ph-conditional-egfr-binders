#!/bin/bash
# Diverse partial-diffusion wave: amplify the NON-A_s baseline winners (C_s, C_l, A_l, H_s, conserved)
# to rebalance the submission's backbone diversity away from A_s dominance. One worker per GPU, each
# processing its round-robin share of the 14 seeds sequentially. partial_T=10 keeps the binding mode.
# Usage: run_pd_diverse_wave.sh [ndes_per_seed] [partial_T] "<gpu list>"
set -uo pipefail
REPO=.; cd "$REPO"
NDES=${1:-150}; PT=${2:-10}; GPUS=${3:-"0 1 2 3 4 5 6 7"}
mkdir -p logs/rfd
mapfile -t SEEDS < runs/pd_diverse_seeds.txt
ga=($GPUS); ng=${#ga[@]}
echo "PD-diverse: ${#SEEDS[@]} seeds x $NDES designs, partial_T=$PT, gpus=(${ga[*]}) $(date -u +%T)"
for gi in "${!ga[@]}"; do
  g=${ga[$gi]}
  (
    for si in "${!SEEDS[@]}"; do
      [ $((si % ng)) -eq $gi ] || continue
      line="${SEEDS[$si]}"; name=$(echo "$line"|cut -f1); pdb=$(echo "$line"|cut -f2)
      out="pd2_${name}"
      echo "[gpu$g] PD $name -> runs/rfd/$out $(date -u +%T)"
      bash scripts/run_rfd_partial.sh "$out" "$pdb" "$NDES" "$g" "$PT" "runs/rfd/sched_$out" \
        >> "logs/rfd/pd2_g${g}.log" 2>&1 || echo "[gpu$g] FAILED $name"
    done
  ) &
done
wait
echo "PD_DIVERSE_WAVE2_DONE $(date -u +%TZ)"
echo "backbones generated: $(find runs/rfd -path '*pd2_*' -name 'des_*.pdb' 2>/dev/null | wc -l)"