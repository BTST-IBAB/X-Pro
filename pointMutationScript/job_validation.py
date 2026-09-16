#!/usr/bin/env python3
"""Scientific input/output validation shared by the X-Pro job worker."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
RESIDUE_ID = re.compile(r"^(-?\d+)([A-Za-z]?)$")


class ValidationError(RuntimeError):
    pass


def sha256_file(path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_residue_id(value: str) -> tuple[int, str]:
    match = RESIDUE_ID.fullmatch(str(value).strip())
    if not match:
        raise ValidationError("Residue must be an integer with an optional insertion code.")
    return int(match.group(1)), match.group(2).upper()


def parse_pdb(path) -> dict:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValidationError("A generated structure file is missing or empty.")
    residues = {}
    atom_count = 0
    model_index = 0
    active = True
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                model_index += 1
                active = model_index == 1
                continue
            if record == "ENDMDL":
                if active:
                    active = False
                continue
            if not active or record != "ATOM" or len(line) < 54:
                continue
            altloc = line[16:17].strip()
            if altloc not in {"", "A"}:
                continue
            chain = line[21:22].strip()
            try:
                number = int(line[22:26])
            except ValueError:
                continue
            insertion = line[26:27].strip().upper()
            residue_name = line[17:20].strip().upper()
            if residue_name not in AA3_TO_1:
                continue
            key = (chain, number, insertion)
            residue = residues.setdefault(key, {"name": residue_name, "atoms": set()})
            if residue["name"] != residue_name:
                raise ValidationError(f"Ambiguous residue identity at {chain}:{number}{insertion}.")
            residue["atoms"].add(line[12:16].strip())
            atom_count += 1
    if atom_count == 0 or not residues:
        raise ValidationError("The structure contains no valid standard protein ATOM records.")
    chains = {}
    for (chain, number, insertion), residue in residues.items():
        chains.setdefault(chain, []).append({
            "number": number, "insertion_code": insertion,
            "name": residue["name"], "atoms": residue["atoms"],
        })
    return {"atom_count": atom_count, "residues": residues, "chains": chains}


def validate_structure_mutation(path, chain: str, residue_id: str,
                                expected_wt: str, mutant: str) -> dict:
    structure = parse_pdb(path)
    number, insertion = parse_residue_id(residue_id)
    key = (str(chain).strip(), number, insertion)
    if key not in structure["residues"]:
        raise ValidationError(f"Mutation site {chain}:{residue_id} is absent from the selected structure model.")
    actual = structure["residues"][key]["name"]
    expected_wt = str(expected_wt).upper()
    mutant = str(mutant).upper()
    if expected_wt not in AA3_TO_1 or mutant not in AA3_TO_1:
        raise ValidationError("Only standard amino-acid substitutions are supported.")
    if actual != expected_wt:
        raise ValidationError(
            f"Mutation validation failed. Requested {expected_wt}{residue_id} → {mutant}; "
            f"the structure contains {actual}{residue_id}."
        )
    if actual == mutant:
        raise ValidationError("The mutant residue must differ from the WT residue.")
    missing = {"N", "CA", "C", "O"} - structure["residues"][key]["atoms"]
    if missing:
        raise ValidationError(
            f"Mutation site {chain}:{residue_id} lacks required backbone atoms: "
            + ", ".join(sorted(missing))
        )
    return {"actual_wt": actual, "atom_count": structure["atom_count"]}


def validate_canonical_pair(canonical_path, wt_pdb, mutant_pdb) -> dict:
    try:
        dataset = json.loads(Path(canonical_path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValidationError("The canonical interaction dataset is missing or invalid.") from exc
    pair = dataset.get("model_pair") or {}
    if pair.get("wt", {}).get("sha256") != sha256_file(wt_pdb):
        raise ValidationError("Canonical WT model hash does not match the analyzed WT PDB.")
    if pair.get("mutant", {}).get("sha256") != sha256_file(mutant_pdb):
        raise ValidationError("Canonical mutant model hash does not match the analyzed mutant PDB.")
    if not isinstance(dataset.get("records"), list):
        raise ValidationError("Canonical interaction records are invalid.")
    return dataset
