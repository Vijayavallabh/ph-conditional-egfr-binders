#!/bin/bash
# WAVE 22 moonshot stage 2-4: motif-(His)-locked soluble ProteinMPNN -> full-atom-target af2ig inputs
# -> AF2-ig fold filter. The binder His that salt-bridges the conserved antigen Asp is FIXED during
# MPNN (fixed_residues_multi from moonshot_prep_mpnn.py); Cys is forbidden (0-Cys expression gate).
#   Usage: moonshot_mpnn_af2ig.sh [staged_dir] [gpus lo-hi] [nseq]
set -uo pipefail
REPO=.
cd "$REPO"
STAGED=${1:-runs/moonshot/staged}
GPUS=${2:-0-7}
NSEQ=${3:-2}
# output dirs sit beside the staged dir: runs/moonshot/staged -> runs/moonshot/{mpnn,af2in,af2ig}
# (original behavior), runs/moonshot/pd_inv02/staged -> runs/moonshot/pd_inv02/{mpnn,af2in,af2ig}.
BASE=$(dirname "$STAGED")
WORK=$BASE/mpnn
AF2IN=$BASE/af2in
AF2IG=$BASE/af2ig
L="$REPO/third_party/LigandMPNN"
PY="$L/.venv/bin/python"
CKPT="$L/model_params/solublempnn_v_48_020.pt"
MDPY="$REPO/envs/md/bin/python"
IFS=- read -r lo hi <<< "$GPUS"; GPULIST=$(seq "$lo" "$hi"); NG=$(echo $GPULIST | wc -w)

echo "=== [moonshot] His-locked soluble MPNN  $(date -u +%FT%TZ) ==="
mkdir -p "$WORK"
"$MDPY" - "$STAGED/pdbs.json" "$WORK" "$NG" <<'PY'
import json,sys,os
pdbs=json.load(open(sys.argv[1])); work=sys.argv[2]; ng=int(sys.argv[3])
for i in range(ng):
    os.makedirs(f"{work}/gpu{i}",exist_ok=True)
    json.dump(pdbs[i::ng], open(f"{work}/gpu{i}/shard.json","w"))
print(f"sharded {len(pdbs)} pdbs across {ng} gpus")
PY
pids=(); gi=0
for g in $GPULIST; do
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py \
      --model_type soluble_mpnn --checkpoint_soluble_mpnn "$CKPT" \
      --pdb_path_multi "$REPO/$WORK/gpu$gi/shard.json" \
      --fixed_residues_multi "$REPO/$STAGED/fixed_residues_multi.json" \
      --chains_to_design A --bias_AA "C:-100.0" \
      --out_folder "$REPO/$WORK/gpu$gi" \
      --batch_size "$NSEQ" --number_of_batches 1 --temperature 0.2 ) > "$WORK/gpu$gi.log" 2>&1 &
  pids+=($!); gi=$((gi+1))
done
for p in "${pids[@]}"; do wait "$p" || true; done
NBB=$(find "$WORK" -path '*/backbones/*.pdb' 2>/dev/null | wc -l)
echo "MPNN done: $NBB designed complexes ($(date -u +%FT%TZ))"
[ "$NBB" -gt 0 ] || { echo "FATAL: no MPNN backbones"; exit 1; }

echo "=== [moonshot] full-atom-target af2ig inputs  $(date -u +%FT%TZ) ==="
"$MDPY" scripts/build_fullatom_af2in.py --mpnn "$WORK/gpu*/backbones/*.pdb" \
  --target data/targets/6ARU_domainIII.pdb --target-chain A --out "$AF2IN"
NIN=$(find "$AF2IN" -name '*.pdb' 2>/dev/null | wc -l)
echo "af2in built: $NIN ($(date -u +%FT%TZ))"
[ "$NIN" -gt 0 ] || { echo "FATAL: no af2in inputs"; exit 1; }

echo "=== [moonshot] AF2-ig fold filter  $(date -u +%FT%TZ) ==="
"$MDPY" scripts/run_af2ig.py --pdbdir "$AF2IN" --out "$AF2IG" --gpus "$GPUS" --mode fast \
  || { echo "FATAL: run_af2ig.py exited non-zero"; exit 1; }
NPRED=$(find "$AF2IG" -name '*_af2pred.pdb' 2>/dev/null | wc -l)
[ "$NPRED" -gt 0 ] || { echo "FATAL: af2ig produced no *_af2pred.pdb"; exit 1; }
echo "MOONSHOT_MPNN_AF2IG_DONE ($NPRED af2pred) $(date -u +%FT%TZ)"