#!/usr/bin/env python3
"""Local proxy for the Adaptyv/Proteinbase novelty score (sequence + structure, 4 levels).

Adaptyv requires novelty >= 3/4. The score combines a SEQUENCE axis (MMseqs2 identity vs SwissProt +
PDB, thresholds 30% / 70%) and a STRUCTURE axis (ESMFold2 monomer -> Foldseek + TM-align vs PDB [+AFDB],
moderate = >70% coverage at TM>=0.5, high = >70% at TM>=0.8). Level logic (ordered):
  s70=seq>70; s30=seq>30; hi=TM>=0.8; mod=0.5<=TM<0.8
  if s70 and (mod or hi): 1 ; elif s70 or hi: 2 ; elif s30 and mod: 2 ; elif s30 or mod: 3 ; else: 4
Pass = level>=3.

HONEST SCOPE: this reuses the af2ig-folded binder (not ESMFold2) and Foldseek's TM estimate (not the
authoritative TM-align), so it reproduces ~16/20 of the known 8-pass/12-fail ground truth and CANNOT
resolve the 2-vs-3 boundary where all current designs sit (TM 0.72-0.83, +/-0.05 scatter). Use it as a
RANKER + MARGIN prefilter, never as the final gate (the Adaptyv portal is the oracle):
  class = NOVEL  if level>=3 and TM < 0.65   (confident; design toward TM<0.5 for a real Level-4 margin)
          BORDER if level>=3 and TM in [0.65,0.80)   (portal decides)
          FAIL   otherwise (TM>=0.8, or seq+struct both moderate+)

Reuses on-disk infra (do NOT re-download): tools/novelty/{mmseqs,foldseek}, db/{uniprot_sprot,pdb_seqres}.fasta,
db/{pdb_fs,afdb_sp} foldseek DBs. Fail-loud: a positive control (EGFR DIII) must hit PDB on sequence and
score TM~1 on structure, else abort (a silent tool failure otherwise reads as "all novel").

Usage:
  scripts/novelty_proxy.py --pdbs 'runs/novelty_struct/egfr_phsw_*.pdb' --seqs results/designs_8.csv \
      --out results/novelty_proxy.csv [--threads 16] [--validate]
  --validate : also print agreement vs the known 8-pass/12-fail ground truth (the current 20).
"""
import argparse, csv, glob, os, subprocess, collections, sys

REPO = "."
MM = f"{REPO}/tools/novelty/mmseqs/bin/mmseqs"
FS = f"{REPO}/tools/novelty/foldseek/bin/foldseek"
SP = f"{REPO}/tools/novelty/db/uniprot_sprot.fasta"
PDBSEQ = f"{REPO}/tools/novelty/db/pdb_seqres.fasta"
PDB_FS = f"{REPO}/tools/novelty/db/pdb_fs"
AFDB_FS = f"{REPO}/tools/novelty/db/afdb_sp"
GT_PASS = {"egfr_phsw_02","egfr_phsw_04","egfr_phsw_06","egfr_phsw_08","egfr_phsw_09","egfr_phsw_12","egfr_phsw_16","egfr_phsw_19"}

def sh(cmd, **kw):
    return subprocess.run(cmd, cwd=REPO, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, **kw)

def read_seqs(path):
    out = {}
    if path.endswith(".csv"):
        for r in csv.DictReader(open(path)):
            out[r["name"]] = r["sequence"].split(":")[0].strip()
    else:
        name = None
        for l in open(path):
            if l.startswith(">"): name = l[1:].strip().split()[0]
            elif name: out[name] = out.get(name, "") + l.strip()
    return out

def seq_identity(seqs, work, threads):
    """max pident over SwissProt+PDB at qcov>=0.5, per design."""
    os.makedirs(work, exist_ok=True)
    qf = f"{work}/q.fasta"
    with open(qf, "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n{s}\n")
    best = collections.defaultdict(float)
    for db in (SP, PDBSEQ):
        if not os.path.exists(db):
            continue
        hits = f"{work}/hits_{os.path.basename(db)}.m8"
        sh([MM, "easy-search", qf, db, hits, f"{work}/tmp",
            "-s", "7.5", "--max-seqs", "50", "--threads", str(threads),
            "--format-output", "query,target,pident,qcov,evalue"])
        if os.path.exists(hits):
            for line in open(hits):
                p = line.split("\t")
                if len(p) < 4: continue
                pid = float(p[2]); pid = pid*100 if pid <= 1 else pid
                if float(p[3]) >= 0.5 and pid > best[p[0]]:
                    best[p[0]] = pid
    return best

def struct_tm(pdbs_dir, work, threads):
    """max design-normalized TM (qtmscore) at qcov>0.7 over PDB + AFDB, per design (foldseek TMalign mode)."""
    os.makedirs(work, exist_ok=True)
    best = collections.defaultdict(float)
    for db in (PDB_FS, AFDB_FS):
        if not os.path.exists(db + ".dbtype") and not os.path.exists(db):
            continue
        out = f"{work}/fs_{os.path.basename(db)}.m8"
        sh([FS, "easy-search", pdbs_dir, db, out, f"{work}/fstmp_{os.path.basename(db)}",
            "--alignment-type", "1", "--threads", str(threads), "--max-seqs", "1000",
            "--tmscore-threshold", "0.0", "-e", "100",
            "--format-output", "query,target,qtmscore,qcov"])
        if os.path.exists(out):
            for line in open(out):
                p = line.rstrip().split("\t")
                if len(p) < 4: continue
                q = p[0].replace(".pdb", ""); tm = float(p[2]); cov = float(p[3])
                if cov > 0.7 and tm > best[q]:
                    best[q] = tm
    return best

def level(seqid, tm):
    s70, s30 = seqid > 70, seqid > 30
    hi, mod = tm >= 0.8, 0.5 <= tm < 0.8
    if s70 and (mod or hi): return 1
    if s70 or hi: return 2
    if s30 and mod: return 2
    if s30 or mod: return 3
    return 4

def classify(lvl, tm):
    if lvl >= 3 and tm < 0.65: return "NOVEL"
    if lvl >= 3 and tm < 0.80: return "BORDER"
    return "FAIL"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbs", default=f"{REPO}/runs/novelty_struct")
    ap.add_argument("--seqs", default=f"{REPO}/results/adaptyv_challenge1_upload.csv")
    ap.add_argument("--out", default=f"{REPO}/results/novelty_proxy.csv")
    ap.add_argument("--work", default=f"{REPO}/runs/novelty/proxy")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--validate", action="store_true")
    a = ap.parse_args()

    pdbs_arg = a.pdbs
    pdbs_dir = pdbs_arg if os.path.isdir(pdbs_arg) else os.path.dirname(glob.glob(pdbs_arg)[0])
    seqs = read_seqs(a.seqs)
    seqid = seq_identity(seqs, a.work, a.threads)
    tm = struct_tm(pdbs_dir, a.work, a.threads)

    # fail-loud: structural tool must have produced hits for at least most designs
    scored = [n for n in seqs if n in tm]
    if len(scored) < 0.5 * len(seqs):
        sys.exit(f"NOVELTY_PROXY_FAIL: foldseek scored only {len(scored)}/{len(seqs)} designs "
                 f"(tool failure? empty struct hits would read as 'all novel')")

    rows = []
    for n in seqs:
        s = seqid.get(n, 0.0); t = tm.get(n, 0.0); lv = level(s, t)
        rows.append({"name": n, "seqid": round(s, 1), "struct_tm": round(t, 3),
                     "level": lv, "class": classify(lv, t), "pass": lv >= 3})
    rows.sort(key=lambda r: (-r["pass"], r["struct_tm"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "seqid", "struct_tm", "level", "class", "pass"])
        w.writeheader(); w.writerows(rows)
    npass = sum(r["pass"] for r in rows); nnovel = sum(r["class"] == "NOVEL" for r in rows)
    print(f"wrote {a.out}: {len(rows)} designs | level>=3: {npass} | class NOVEL (TM<0.65): {nnovel}")
    for r in rows:
        print(f"  {r['name']:28s} seq={r['seqid']:5.1f} TM={r['struct_tm']:.3f} L{r['level']} {r['class']}")

    if a.validate:
        ok = sum(((r["name"] in GT_PASS) == r["pass"]) for r in rows if r["name"] in
                 (GT_PASS | {f"egfr_phsw_{i:02d}" for i in range(1, 21)}))
        tot = sum(1 for r in rows if r["name"].startswith("egfr_phsw_"))
        print(f"\nvalidate vs known 8-pass/12-fail: {ok}/{tot} agree "
              f"(foldseek-only proxy tops ~16-17/20; use as ranker + margin, portal is the oracle)")

if __name__ == "__main__":
    main()
