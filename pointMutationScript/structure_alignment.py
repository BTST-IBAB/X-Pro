"""Create display-only WT/mutant structural alignment for the X-Pro viewer."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


def _first_model_ca(path):
    atoms = {}
    saw_model = False
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record != "ATOM" or line[12:16].strip() != "CA" or len(line) < 54:
                continue
            altloc = line[16:17]
            if altloc not in {" ", "", "A"}:
                continue
            try:
                key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip().upper())
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            if key not in atoms or altloc.strip() == "":
                atoms[key] = xyz
    return atoms


def _transform_pdb(source, destination, rotation, moving_centroid, fixed_centroid):
    """Transform a copy for display; the downloaded source stays byte-for-byte unchanged."""
    lines = []
    model_seen = False
    with Path(source).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if model_seen:
                    break
                model_seen = True
                continue
            if record == "ENDMDL" and model_seen:
                break
            if record in {"ATOM", "HETATM"} and len(line) >= 54:
                try:
                    point = np.array([
                        float(line[30:38]), float(line[38:46]), float(line[46:54])
                    ])
                except ValueError:
                    continue
                transformed = (point - moving_centroid) @ rotation + fixed_centroid
                padded = line.rstrip("\n").ljust(80)
                line = (
                    padded[:30]
                    + f"{transformed[0]:8.3f}{transformed[1]:8.3f}{transformed[2]:8.3f}"
                    + padded[54:]
                    + "\n"
                )
            lines.append(line)
    if not any(line.startswith("END") for line in lines):
        lines.append("END\n")
    Path(destination).write_text("".join(lines), encoding="utf-8")


def create_display_alignment(wt_pdb, mutant_pdb, output_dir, mutation_chain,
                             mutation_residue, mutation_insertion_code=""):
    """Align mutant onto WT using corresponding CA atoms except the mutation site."""
    wt_pdb = Path(wt_pdb)
    mutant_pdb = Path(mutant_pdb)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    display_wt = output_dir / "display_wt.pdb"
    display_mutant = output_dir / "display_mutant_aligned.pdb"
    display_wt.write_bytes(wt_pdb.read_bytes())

    wt_atoms = _first_model_ca(wt_pdb)
    mutant_atoms = _first_model_ca(mutant_pdb)
    excluded = (str(mutation_chain).strip(), int(mutation_residue),
                str(mutation_insertion_code).strip().upper())
    shared = sorted((set(wt_atoms) & set(mutant_atoms)) - {excluded})
    metadata = {
        "method": "Kabsch least-squares superposition of corresponding CA atoms",
        "mutation_site_excluded": True,
        "correspondence_identity": "chain + residue number + insertion code",
        "display_only": True,
        "download_coordinates_unchanged": True,
        "matched_ca_count": len(shared),
        "wt_display_file": display_wt.name,
        "mutant_display_file": display_mutant.name,
        "warnings": [],
    }
    if len(shared) < 3:
        display_mutant.write_bytes(mutant_pdb.read_bytes())
        metadata.update({"aligned": False, "rmsd_angstrom": None})
        metadata["warnings"].append(
            "Fewer than three corresponding CA atoms were available; the viewer uses unaligned coordinates."
        )
    else:
        fixed = np.vstack([wt_atoms[key] for key in shared])
        moving = np.vstack([mutant_atoms[key] for key in shared])
        fixed_centroid = fixed.mean(axis=0)
        moving_centroid = moving.mean(axis=0)
        covariance = (moving - moving_centroid).T @ (fixed - fixed_centroid)
        u, _, vt = np.linalg.svd(covariance)
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            u[:, -1] *= -1
            rotation = u @ vt
        fitted = (moving - moving_centroid) @ rotation + fixed_centroid
        rmsd = math.sqrt(float(np.mean(np.sum((fitted - fixed) ** 2, axis=1))))
        _transform_pdb(mutant_pdb, display_mutant, rotation, moving_centroid, fixed_centroid)
        metadata.update({
            "aligned": True,
            "rmsd_angstrom": round(rmsd, 4),
            "rotation_matrix": rotation.round(8).tolist(),
            "moving_centroid": moving_centroid.round(8).tolist(),
            "fixed_centroid": fixed_centroid.round(8).tolist(),
        })
    (output_dir / "display_alignment.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return metadata
