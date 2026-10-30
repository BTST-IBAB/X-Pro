#!/usr/bin/env python3
"""
xpro_summarise.py
═════════════════
Reads the canonical interaction JSON produced by X-Pro's
interaction_comparison.py and automatically generates a
research-paper-style 2-paragraph structural interpretation
using a local Ollama model (default: phi4-mini-local).

X-Pro calls this automatically at the end of each job if
--auto is passed, or you can run it manually at any time.

Usage
─────
  # Manual run after X-Pro completes
  python3 xpro_summarise.py <canonical_json> [options]

  # Auto mode (called by X-Pro at job end)
  python3 xpro_summarise.py <canonical_json> --auto

Options
───────
  --model      Ollama model        (default: phi4-mini-local)
  --ollama     Ollama URL          (default: http://localhost:11434)
  --protein    Protein name        (default: read from JSON)
  --disease    Disease context     (optional)
  --out        Save output to file (default: summary.txt beside the JSON)
  --auto       Suppress progress output, write file only

Example
───────
  python3 xpro_summarise.py \
    /var/www/html/fileUpload/XPRO-20261006-.../canonical.json \
    --model phi4-mini-local \
    --protein "GALNS" \
    --disease "Mucopolysaccharidosis IVA (Morquio A syndrome)"
"""

import argparse
import json
import sys
import requests
from pathlib import Path


# ── CLI ───────────────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument("canonical_json",  help="Path to X-Pro canonical interaction JSON")
ap.add_argument("--model",   default="phi4-mini-local")
ap.add_argument("--ollama",  default="http://localhost:11434")
ap.add_argument("--protein", default=None)
ap.add_argument("--disease", default=None)
ap.add_argument("--out",     default=None)
ap.add_argument("--auto",      action="store_true", help="Suppress progress output")
ap.add_argument("--user-hint", default=None,
                help="User's belief about the verdict (overrides auto-detect hint)")
args = ap.parse_args()


def log(msg):
    if not args.auto:
        print(msg)


# ── Load JSON ─────────────────────────────────────────────────────────────────
json_path = Path(args.canonical_json)
if not json_path.exists():
    print(f"[ERROR] File not found: {json_path}", file=sys.stderr)
    sys.exit(1)

try:
    dataset = json.loads(json_path.read_text(encoding="utf-8"))
except Exception as e:
    print(f"[ERROR] Could not parse JSON: {e}", file=sys.stderr)
    sys.exit(1)

# Default output path beside the JSON
out_path = Path(args.out) if args.out else json_path.parent / "summary.txt"


# ── Extract mutation info ─────────────────────────────────────────────────────
mutation  = dataset.get("mutation", {})
wt_res    = mutation.get("residue_wt",     "?")
mut_res   = mutation.get("residue_mutant", "?")
position  = mutation.get("residue_number", "?")
chain     = mutation.get("chain",          "?")
protein   = args.protein or dataset.get("protein_name", "the protein")
disease   = args.disease or ""
mut_label = f"{wt_res}{position}{mut_res}"

log(f"[INFO] Mutation  : {mut_label}  (chain {chain})")
log(f"[INFO] Protein   : {protein}")
log(f"[INFO] Model     : {args.model}")
log(f"[INFO] Records   : {len(dataset.get('records', []))}")


# ── Parse all records from canonical JSON ────────────────────────────────────
records = dataset.get("records", [])

# Categorise every record
hbonds_retained   = []   # (partner_label, wt_dist, mut_dist)
hbonds_lost       = []   # (partner_label, wt_dist)
hbonds_gained     = []   # (partner_label, mut_dist)

contacts_retained_changed = []  # contact count changed (retained but Δ ≠ 0)
contacts_retained_same    = []  # retained, count unchanged
contacts_lost             = []  # completely lost
contacts_gained           = []  # newly gained

for rec in records:
    result      = rec.get("result", "")
    iclass      = rec.get("interaction_class", "")
    partner     = rec.get("partner_residue_name") or "UNK"
    partner_num = rec.get("partner_residue", "")
    p_label     = f"{partner}{partner_num}"

    wt_dist   = rec.get("wt_distance")
    mut_dist  = rec.get("mutant_distance")
    wt_cnt    = rec.get("wt_contact_count",    0)
    mut_cnt   = rec.get("mutant_contact_count", 0)
    delta     = rec.get("delta_contact_count",  0)

    if iclass == "hydrogen_bond":
        if result == "Retained":
            hbonds_retained.append((p_label, wt_dist, mut_dist))
        elif result == "Lost":
            hbonds_lost.append((p_label, wt_dist))
        elif result == "Gained":
            hbonds_gained.append((p_label, mut_dist))

    elif iclass == "non_bonded_contact":
        if result == "Lost":
            contacts_lost.append((p_label, wt_cnt))
        elif result == "Gained":
            contacts_gained.append((p_label, mut_cnt))
        elif result == "Retained":
            if delta != 0:
                contacts_retained_changed.append((p_label, wt_cnt, mut_cnt, delta))
            else:
                contacts_retained_same.append((p_label, wt_cnt))


# ── Format data summary for prompt ───────────────────────────────────────────
lines = []
lines.append(f"Protein: {protein}")
lines.append(f"Mutation: {mut_label} (chain {chain})")
if disease:
    lines.append(f"Disease relevance: {disease}")
lines.append("")

# WT H-bonds
if hbonds_retained or hbonds_lost:
    wt_hb_parts = []
    for label, wt_d, _ in hbonds_retained:
        wt_hb_parts.append(f"{label} ({wt_d:.2f} Å)" if wt_d else label)
    for label, wt_d in hbonds_lost:
        wt_hb_parts.append(f"{label} ({wt_d:.2f} Å)" if wt_d else label)
    lines.append(f"Wild type H-bonds: {', '.join(wt_hb_parts)}")

# WT non-bonded contacts
all_wt_contacts = (
    [(l, wt) for l, wt, _, _ in contacts_retained_changed] +
    [(l, wt) for l, wt in contacts_retained_same] +
    [(l, wt) for l, wt in contacts_lost]
)
if all_wt_contacts:
    wt_c_parts = [f"{l} ({wt} contacts)" for l, wt in all_wt_contacts]
    if contacts_gained:
        wt_c_parts += [f"{l} (not present in WT)" for l, _ in contacts_gained]
    lines.append(f"Wild type non-bonded contacts: {', '.join(wt_c_parts)}")

lines.append("")

# Mutant H-bonds
if hbonds_retained:
    ret_parts = []
    for label, wt_d, mut_d in hbonds_retained:
        if mut_d:
            ret_parts.append(f"{label} ({mut_d:.2f} Å)")
        else:
            ret_parts.append(label)
    lines.append(f"Retained H-bonds in mutant: {', '.join(ret_parts)}")

if hbonds_lost:
    lines.append(f"Lost H-bonds: {', '.join(l for l, _ in hbonds_lost)}")

if hbonds_gained:
    gained_parts = []
    for label, mut_d in hbonds_gained:
        gained_parts.append(f"{label} ({mut_d:.2f} Å)" if mut_d else label)
    lines.append(f"Newly gained H-bonds in mutant: {', '.join(gained_parts)}")

# Contact changes
if contacts_retained_changed:
    for label, wt_cnt, mut_cnt, delta in contacts_retained_changed:
        direction = "rises" if delta > 0 else "falls"
        lines.append(
            f"Non-bonded contact change — {label}: {direction} from "
            f"{wt_cnt} to {mut_cnt} atom pairs (Δ {delta:+d})"
        )

if contacts_lost:
    for label, wt_cnt in contacts_lost:
        lines.append(f"Non-bonded contact lost — {label}: {wt_cnt} contact(s) → 0 (completely lost)")

if contacts_gained:
    for label, mut_cnt in contacts_gained:
        lines.append(f"Non-bonded contact gained — {label}: 0 → {mut_cnt} new contact(s)")

if contacts_retained_same:
    same_parts = [f"{l} ({wt} contacts, unchanged)" for l, wt in contacts_retained_same]
    lines.append(f"Unchanged contacts: {', '.join(same_parts)}")

# Summary flags for net conclusion guidance
n_hb_lost    = len(hbonds_lost)
n_hb_gained  = len(hbonds_gained)
n_hb_ret     = len(hbonds_retained)
n_ct_lost    = len(contacts_lost)
n_ct_gained  = len(contacts_gained)
total_gained_contacts = sum(mut for _, mut in contacts_gained)
total_lost_contacts   = sum(wt for _, wt in contacts_lost)
gains_in_retained     = sum(d for _, _, _, d in contacts_retained_changed if d > 0)
losses_in_retained    = sum(abs(d) for _, _, _, d in contacts_retained_changed if d < 0)

# Net statistics are kept SEPARATE from data_summary so the model
# does not copy them verbatim into the output paragraphs.
# They are only used internally to build the _verdict_hint below.
data_summary = "\n".join(lines)

log("\n[INFO] Extracted data summary:")
log(data_summary)
log("")


# ── Build prompt ──────────────────────────────────────────────────────────────
disease_line = f"The mutation is associated with {disease}.\n" if disease else ""

# Determine verdict guidance from net statistics to help the model reason correctly
_hb_lost    = len(hbonds_lost)
_hb_gained  = len(hbonds_gained)
_ct_lost    = len(contacts_lost)
_ct_gained  = len(contacts_gained)
_delta_sum  = sum(d for _, _, _, d in contacts_retained_changed)

# Build a reasoning hint so the model can make a consistent verdict
if _hb_lost == 0 and _ct_lost == 0 and _delta_sum >= 0:
    _verdict_hint = (
        "All hydrogen bonds are retained and no contacts are lost. "
        "Any additional contacts are a gain. "
        "Based on this, the mutation is likely BENIGN or mildly destabilising."
    )
elif _hb_lost == 0 and _ct_lost == 0 and _delta_sum < 0:
    _verdict_hint = (
        "All hydrogen bonds are retained but some non-bonded contacts are reduced. "
        "Based on this, the mutation may be mildly DELETERIOUS (subtle packing defect)."
    )
elif _hb_lost > 0 and _ct_lost == 0:
    _verdict_hint = (
        f"{_hb_lost} hydrogen bond(s) are lost. "
        "Loss of hydrogen bonds typically destabilises protein structure. "
        "Based on this, the mutation is likely DELETERIOUS."
    )
elif _hb_lost == 0 and _ct_lost > 0:
    _verdict_hint = (
        f"{_ct_lost} non-bonded contact pair(s) are completely lost. "
        "Based on this, the mutation is likely DELETERIOUS due to reduced packing."
    )
else:
    _verdict_hint = (
        f"{_hb_lost} hydrogen bond(s) lost and {_ct_lost} contact pair(s) lost. "
        "Multiple interaction losses indicate the mutation is likely DELETERIOUS."
    )

# If the user provided a hint via --user-hint, use it instead of the auto-detected one.
# This lets the user guide the verdict when they already know the expected outcome.
if args.user_hint:
    _verdict_hint = args.user_hint

prompt = f"""/no_think You are a structural biologist writing a supplementary section for a research paper.

Write exactly 4 paragraphs interpreting the structural effect of the following mutation.
{disease_line}
Follow this EXACT writing style and structure:
- Paragraph 1: Classify the amino acid change. State the physicochemical nature of both the wild-type and mutant residues (choose from: aliphatic/hydrophobic, aromatic, polar/hydrophilic, acidic/negatively charged, basic/positively charged). Then state whether the mutant is larger, smaller, or similar in size. End with one sentence on what structural consequence this physicochemical and size difference may introduce.
- Paragraph 2: Describe what the wild-type residue does — its position in the protein, all contacts it makes, and its hydrogen bonds with specific Å distances
- Paragraph 3: Describe what changes in the mutant — which H-bonds are retained (with new distances), which are lost, which are gained, and how non-bonded contact counts changed with exact atom pair numbers
- Paragraph 4: State clearly whether the mutation is BENIGN or DELETERIOUS and explain why in 2-3 sentences based on the structural evidence above

Strict rules:
1. Use ONLY the data provided below — do not invent any contacts, residues or distances
2. Write in concise, precise, third-person research paper style
3. No bullet points, no headers, no markdown — plain prose paragraphs only
4. Always include specific Å distances for every hydrogen bond mentioned
5. For non-bonded contacts always state the exact atom pair counts (e.g. "rises from 3 to 8 atom pairs")
6. Paragraph 1 MUST classify both residues using standard biochemistry terminology and state the size change
7. Paragraph 4 MUST start with either "This mutation is classified as BENIGN." or "This mutation is classified as DELETERIOUS." — no other opening is acceptable
8. The final sentence of paragraph 4 must summarise the net structural mechanism
9. NEVER reproduce any labels, headers or statistics from the Data block — write only prose paragraphs

Amino acid classification reference (for paragraph 1 only):
- Aliphatic/hydrophobic: Gly, Ala, Val, Leu, Ile, Met, Pro
- Aromatic: Phe, Tyr, Trp, His
- Polar/hydrophilic: Ser, Thr, Cys, Asn, Gln
- Acidic (negatively charged): Asp, Glu
- Basic (positively charged): Lys, Arg, His
- Size order (smallest to largest): Gly < Ala < Ser/Cys < Val/Thr < Leu/Ile/Asn/Asp < Gln/Glu/Lys < Met < Phe/Arg < Tyr < Trp

Reasoning hint for verdict (guide paragraph 4 only — do not copy this text):
{_verdict_hint}

Data (use values only — do not copy labels or headers):
{data_summary}
"""


# ── Call Ollama ───────────────────────────────────────────────────────────────
log(f"[INFO] Calling Ollama ({args.model})...")

try:
    response = requests.post(
        f"{args.ollama}/api/generate",
        json={
            "model":  args.model,
            "prompt": prompt,
            "stream": False,
            "think":  False,
            "options": {"temperature": 0.2},
        },
        proxies={"http": "", "https": ""},   # bypass proxy for localhost
        timeout=300,
    )
    response.raise_for_status()
    paragraph = response.json().get("response", "").strip()

except Exception as e:
    print(f"[ERROR] Ollama call failed: {e}", file=sys.stderr)
    sys.exit(1)


# ── Output ────────────────────────────────────────────────────────────────────
if not args.auto:
    print("\n" + "═" * 60)
    print("GENERATED INTERPRETATION")
    print("═" * 60)
    print(paragraph)
    print("═" * 60 + "\n")

out_path.write_text(paragraph, encoding="utf-8")
log(f"[INFO] Summary saved to: {out_path}")
