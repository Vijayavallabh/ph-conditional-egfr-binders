#!/bin/bash
# RFdiffusion PARTIAL DIFFUSION of a binder:target complex: diversify a validated binder while
# keeping the binding mode (Vazquez-Torres). Input: binder = chain A (len Lb), target = chain B (1..191).
# Contig [Lb-Lb/0 B1-191] diffuses the binder, keeps the target. diffuser.partial_T controls noise
# (T=50 default; ~10 gentle/retains binding, ~20 more diverse).
# Usage: run_rfd_partial.sh <out_subdir> <input_complex.pdb> <num_designs> <gpu> [partial_T] [sched_dir]
set -euo pipefail
REPO=.
OUTSUB=$1; INP=$(readlink -f "$2"); NDES=$3; GPU=$4; PT=${5:-10}; SCHED=${6:-}
source "$REPO/envs/uplifting-rfdiffusion1/bin/activate"
export RFD_ROOT="$REPO/third_party/RFdiffusion" WEIGHTS="$REPO/weights/rfdiffusion1"
export MODEL_OPT_JIT_ROOT="$REPO/.jit" CUDA_VISIBLE_DEVICES=$GPU
# Cap CPU threads per process: torch/OMP default to ALL cores, so 8 concurrent waves (one per GPU)
# oversubscribe a 192-core box to load ~400 and thrash (0% GPU util during SE3 prep, esp. on larger
# binders). ~24 threads x 8 procs = 192, keeping every GPU fed. Override by exporting RFD_CPU_THREADS.
export OMP_NUM_THREADS=${RFD_CPU_THREADS:-24} MKL_NUM_THREADS=${RFD_CPU_THREADS:-24} OPENBLAS_NUM_THREADS=${RFD_CPU_THREADS:-24}
Lb=$(awk 'substr($0,1,4)=="ATOM" && substr($0,22,1)=="A" && substr($0,13,4)~/CA/{n++}END{print n}' "$INP")
Lt=$(awk 'substr($0,1,4)=="ATOM" && substr($0,22,1)=="B" && substr($0,13,4)~/CA/{n++}END{print n}' "$INP")
OUTDIR="$REPO/runs/rfd/$OUTSUB"; mkdir -p "$OUTDIR"
EXTRA=(); [ -n "$SCHED" ] && { mkdir -p "$SCHED"; SCHED=$(readlink -f "$SCHED"); EXTRA+=("inference.schedule_directory_path=$SCHED"); }
echo "[partial] $OUTSUB Lb=$Lb Lt=$Lt partial_T=$PT ndes=$NDES gpu=$GPU"
cd "$REPO/third_party/uplifting-biomolecular-modeling/rfdiffusion1"
exec bash run.sh design --config h100 --mode fast \
  inference.input_pdb="$INP" \
  "contigmap.contigs=[${Lb}-${Lb}/0 B1-${Lt}]" \
  diffuser.partial_T="$PT" \
  inference.num_designs="$NDES" \
  "${EXTRA[@]}" \
  inference.output_prefix="$OUTDIR/des"
