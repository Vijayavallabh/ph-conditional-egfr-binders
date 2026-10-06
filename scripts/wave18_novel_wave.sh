#!/bin/bash
# WAVE 18 novelty-recovery generation: generalized partial-diffusion wave.
# Generalizes run_pd_diverse_wave.sh (which hardcoded the seeds path + "pd2_" prefix) to take a
# seeds file, ndesigns, partial_T, output-prefix, and GPU list. Used to diversify the validated
# complexes AGGRESSIVELY (partial_T=20 "more diverse") to push designs off common folds and raise
# the Adaptyv STRUCTURAL-novelty score (the confirmed blocker: 12/20 fail on fold-commonness, not
# sequence). One worker per GPU, round-robin seeds. Prints the marker the wave17 downstream expects.
# Usage: wave18_novel_wave.sh <seeds_file> <ndes> <partial_T> <out_prefix> "<gpu list>"
set -uo pipefail
REPO=.; cd "$REPO"
SEEDS=${1:-runs/pd_diverse_seeds.txt}; NDES=${2:-400}; PT=${3:-20}; PREFIX=${4:-pdnov}; GPUS=${5:-"0 1 2 3 4 5 6 7"}
mkdir -p logs/rfd
mapfile -t SEEDLINES < "$SEEDS"
ga=($GPUS); ng=${#ga[@]}
echo "WAVE18 [$PREFIX]: ${#SEEDLINES[@]} seeds x $NDES designs, partial_T=$PT, gpus=(${ga[*]}) $(date -u +%FT%TZ)"
for gi in "${!ga[@]}"; do
  g=${ga[$gi]}
  (
    for si in "${!SEEDLINES[@]}"; do
      [ $((si % ng)) -eq "$gi" ] || continue
      line="${SEEDLINES[$si]}"; name=$(echo "$line"|cut -f1); pdb=$(echo "$line"|cut -f2)
      out="${PREFIX}_${name}"
      echo "[gpu$g] $out (T=$PT) $(date -u +%T)"
      bash scripts/run_rfd_partial.sh "$out" "$pdb" "$NDES" "$g" "$PT" "runs/rfd/sched_$out" \
        >> "logs/rfd/${PREFIX}_g${g}.log" 2>&1 || echo "[gpu$g] FAILED $name"
    done
  ) &
done
wait
echo "backbones: $(find runs/rfd -path "*${PREFIX}_*" -name 'des_*.pdb' 2>/dev/null | grep -v _traj | wc -l)"
echo "WAVE18_${PREFIX}_DONE $(date -u +%FT%TZ)"
echo "PD_DIVERSE_WAVE2_DONE $(date -u +%FT%TZ)"   # marker the wave17 downstream driver waits for
