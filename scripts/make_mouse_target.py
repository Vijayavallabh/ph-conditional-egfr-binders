#!/usr/bin/env python3
"""Mouse EGFR Domain III counterpart of the human 6ARU crop (for cross-reactivity, objective 2).

Aligns the human crop sequence (from 6ARU_domainIII.pdb) to mouse EGFR (UniProt Q01279) and
writes the mouse Domain III sequence aligned 1:1 to the crop positions, plus a per-position
human/mouse map. Verifies the conserved anchors match (human==mouse) and the divergent
cetuximab contacts differ. Run with: uv run python scripts/make_mouse_target.py
"""
from __future__ import annotations

import json
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

REPO = Path(__file__).resolve().parents[1]
CROP = REPO / "data" / "targets" / "6ARU_domainIII.pdb"
MOUSE_FA = REPO / "data" / "targets" / "EGFR_mouse_Q01279.fasta"
CROP_MAP = REPO / "data" / "targets" / "6ARU_domainIII_map.json"
OUT_FA = REPO / "data" / "targets" / "EGFR_mouse_domainIII.fasta"
OUT_MAP = REPO / "data" / "targets" / "mouse_domainIII_map.json"

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
}


def crop_residues():
    out, seen = [], set()
    for line in CROP.read_text().splitlines():
        if line[:4] != "ATOM":
            continue
        if line[12:16].strip() != "CA":
            continue
        author = int(line[22:26])
        if author in seen:
            continue
        seen.add(author)
        out.append((author, THREE2ONE.get(line[17:20].strip(), "X")))
    return out


def main() -> None:
    res = crop_residues()
    human_seq = "".join(aa for _, aa in res)
    mouse_full = str(next(SeqIO.parse(str(MOUSE_FA), "fasta")).seq)

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    aln = aligner.align(human_seq, mouse_full)[0]

    human_to_mouse: dict[int, str] = {}
    hblocks, mblocks = aln.aligned
    for (h0, h1), (m0, m1) in zip(hblocks, mblocks):
        for k in range(h1 - h0):
            human_to_mouse[h0 + k] = mouse_full[m0 + k]

    cm = json.loads(CROP_MAP.read_text())
    conserved = {int(v["ordinal"]) for v in cm["conserved_binding"].values()}
    divergent = {int(v["ordinal"]) for v in cm["divergent_not_binding"].values()}

    rows, mouse_chars = [], []
    n_ident = n_aligned = 0
    for i, (author, haa) in enumerate(res):
        maa = human_to_mouse.get(i, "-")
        mouse_chars.append(maa if maa != "-" else "X")
        if maa != "-":
            n_aligned += 1
            n_ident += (maa == haa)
        rows.append({
            "ordinal": i + 1, "author_6aru": author, "p00533": author + 24,
            "human": haa, "mouse": maa, "same": maa == haa,
        })
    mouse_seq = "".join(mouse_chars)

    def check(ordset):
        return {o: {"human": rows[o - 1]["human"], "mouse": rows[o - 1]["mouse"],
                    "same": rows[o - 1]["same"]} for o in sorted(ordset)}

    OUT_FA.write_text(
        f">EGFR_mouse_Q01279_domainIII_alignedto_6ARUcrop len={len(mouse_seq)}\n{mouse_seq}\n"
    )
    result = {
        "note": "mouse residues aligned 1:1 to human 6ARU DomainIII crop (ordinal = crop sequential position).",
        "crop_len": len(res),
        "identity_pct_over_crop": round(100 * n_ident / max(1, n_aligned), 2),
        "mouse_seq": mouse_seq,
        "conserved_anchor_check": check(conserved),
        "divergent_check": check(divergent),
        "positions": rows,
    }
    OUT_MAP.write_text(json.dumps(result, indent=2))

    print(f"crop {len(res)} res | mouse identity over crop = {result['identity_pct_over_crop']}%")
    print("conserved anchors (expect SAME):")
    for o, d in result["conserved_anchor_check"].items():
        print(f"  ord{o}: human {d['human']} mouse {d['mouse']} same={d['same']}")
    print("divergent (expect DIFFER):")
    for o, d in result["divergent_check"].items():
        print(f"  ord{o}: human {d['human']} mouse {d['mouse']} same={d['same']}")
    print(f"wrote {OUT_FA}")
    print(f"wrote {OUT_MAP}")


if __name__ == "__main__":
    main()
