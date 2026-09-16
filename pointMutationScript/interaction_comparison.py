#!/usr/bin/env python3
"""Canonical WT-vs-mutant interaction comparison for X-Pro.

Scientific interaction classes are not inferred here. Hydrogen bonds come from
the LigPlot ``.hhb`` file produced from HBPLUS output, and non-bonded contacts
come from the LigPlot ``.nnb`` file. Covalent ``.bonds`` records are
intentionally excluded, so peptide connectivity is never classified as a
gained/lost non-bonded interaction. This module only normalizes, compares,
exports, and visualizes those engine-defined interactions.

Residues are identified structurally by ``chain + residue number + insertion
code``.  Residue names are deliberately excluded from comparison identity so
that (for example) WT LEU A:87 and mutant TRP A:87 are the same position.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path


SCHEMA_VERSION = "1.0"
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})
INTERACTION_TYPES = {
    "hydrogen_bond": {
        "label": "Hydrogen bond",
        "source": "HBPLUS via the LigPlot .hhb output",
    },
    "non_bonded_contact": {
        "label": "Non-bonded contact",
        "source": "LigPlot/HBPLUS via the LigPlot .nnb output",
    },
}


def _sha256_file(filename):
    digest = hashlib.sha256()
    with Path(filename).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_residue_identifier(value):
    """Return ``(PDB residue number, insertion code)`` from values like 87A."""
    match = re.fullmatch(r"\s*(-?\d+)\s*([A-Za-z]?)\s*", str(value))
    if not match:
        raise ValueError(f"Invalid PDB residue identifier: {value!r}")
    return int(match.group(1)), match.group(2).upper()


def _position(chain, residue_number, insertion_code=""):
    return {
        "chain": (chain or "").strip(),
        "residue_number": int(residue_number),
        "insertion_code": (insertion_code or "").strip().upper(),
    }


def _position_key(position):
    return (
        position["chain"],
        int(position["residue_number"]),
        position.get("insertion_code", ""),
    )


def _position_from_key(key):
    return _position(key[0], key[1], key[2])


def _position_text(position):
    chain = position["chain"] or "(blank chain)"
    return f"{chain}:{position['residue_number']}{position.get('insertion_code', '')}"


def _residue_label(position, residue_name):
    name = residue_name or "UNK"
    chain = position["chain"] or "_"
    return f"{chain}:{name}{position['residue_number']}{position.get('insertion_code', '')}"


def _atom_role(atom_name):
    return "backbone" if atom_name.strip().upper() in BACKBONE_ATOMS else "side-chain"


def parse_pdb_structure(filename):
    """Read coordinates from only the first PDB model.

    Blank alternate locations are preferred over altloc A.  Other alternate
    locations are ignored so distances never silently mix conformers/models.
    Both ATOM and HETATM are read for diagnostics, although selected amino-acid
    mutation sites are expected to remain ATOM records.
    """
    atoms = {}
    residues = {}
    path = Path(filename)
    if not path.exists():
        return {"atoms": atoms, "residues": residues, "warnings": [f"PDB file missing: {path}"]}

    warnings = []
    saw_model = False
    with path.open("r", errors="replace") as handle:
        for line in handle:
            record = line[0:6].strip()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            altloc = line[16:17]
            if altloc not in {" ", "", "A"}:
                continue
            try:
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip().upper()
                chain = line[21:22].strip()
                residue_number = int(line[22:26])
                insertion_code = line[26:27].strip().upper()
                coordinates = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                warnings.append(f"Skipped malformed coordinate record in {path.name}: {line.rstrip()}")
                continue

            pos_key = (chain, residue_number, insertion_code)
            atom_key = pos_key + (atom_name,)
            atom_data = {
                "coordinates": coordinates,
                "residue_name": residue_name,
                "record_type": record,
                "altloc": altloc.strip(),
            }
            # A blank altloc is canonical; do not replace it with altloc A.
            if atom_key not in atoms or (not atom_data["altloc"] and atoms[atom_key]["altloc"]):
                atoms[atom_key] = atom_data

            residue = residues.setdefault(
                pos_key,
                {"residue_name": residue_name, "record_types": set(), "atoms": set()},
            )
            residue["record_types"].add(record)
            residue["atoms"].add(atom_name)

    return {"atoms": atoms, "residues": residues, "warnings": warnings}


def _parse_contact_line(line):
    """Parse LigPlot's fixed-width human-readable .hhb/.nnb contact row."""
    if len(line.rstrip("\n")) < 40:
        return None
    try:
        residue_name_1 = line[0:3].strip().upper()
        chain_1 = line[4:5].strip()
        residue_token_1 = line[5:10].strip()
        atom_1 = line[12:16].strip()
        residue_name_2 = line[21:24].strip().upper()
        chain_2 = line[25:26].strip()
        residue_token_2 = line[26:31].strip()
        atom_2 = line[33:37].strip()
        reported_distance = float(line[39:].strip().split()[0])
        residue_number_1, insertion_code_1 = split_residue_identifier(residue_token_1)
        residue_number_2, insertion_code_2 = split_residue_identifier(residue_token_2)
    except (ValueError, IndexError):
        return None
    if not all((residue_name_1, atom_1, residue_name_2, atom_2)):
        return None
    return {
        "residue_1": {
            **_position(chain_1, residue_number_1, insertion_code_1),
            "residue_name": residue_name_1,
        },
        "atom_1": atom_1,
        "residue_2": {
            **_position(chain_2, residue_number_2, insertion_code_2),
            "residue_name": residue_name_2,
        },
        "atom_2": atom_2,
        "engine_distance": reported_distance,
    }


def _coordinate_for(structure, residue, atom_name):
    key = _position_key(residue) + (atom_name,)
    atom = structure["atoms"].get(key)
    return atom["coordinates"] if atom else None


def parse_ligplot_contact_file(filename, interaction_class, structure):
    """Parse one LigPlot contact file and calculate coordinate distances."""
    if interaction_class not in INTERACTION_TYPES:
        raise ValueError(f"Unsupported interaction class: {interaction_class}")
    path = Path(filename)
    if not path.exists():
        return [], [f"Interaction file missing: {path}"]

    contacts = []
    warnings = []
    with path.open("r", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            parsed = _parse_contact_line(line)
            if parsed is None:
                continue
            coordinate_1 = _coordinate_for(structure, parsed["residue_1"], parsed["atom_1"])
            coordinate_2 = _coordinate_for(structure, parsed["residue_2"], parsed["atom_2"])
            if coordinate_1 is not None and coordinate_2 is not None:
                distance = math.dist(coordinate_1, coordinate_2)
                distance_source = "PDB coordinates (Euclidean)"
            else:
                # The engine distance is measured, not inferred.  Preserve it
                # when a coordinate lookup is impossible and label its source.
                distance = parsed["engine_distance"]
                distance_source = "LigPlot/HBPLUS reported distance"
                warnings.append(
                    f"Used engine-reported distance for {path.name}:{line_number}; "
                    "one or both PDB atom coordinates were unavailable."
                )

            contact = {
                **parsed,
                "interaction_class": interaction_class,
                "interaction_type": INTERACTION_TYPES[interaction_class]["label"],
                "interaction_source": INTERACTION_TYPES[interaction_class]["source"],
                "distance": round(distance, 4),
                "distance_source": distance_source,
                "atom_1_role": _atom_role(parsed["atom_1"]),
                "atom_2_role": _atom_role(parsed["atom_2"]),
            }
            atom_identity = (
                _position_key(parsed["residue_1"]) + (parsed["atom_1"],)
                + _position_key(parsed["residue_2"]) + (parsed["atom_2"], interaction_class)
            )
            contact["contact_id"] = "contact-" + hashlib.sha256(
                json.dumps(atom_identity).encode("utf-8")
            ).hexdigest()[:14]
            if interaction_class == "hydrogen_bond":
                # LigPlot .hhb columns are explicitly Donor then Acceptor.
                contact.update(
                    {
                        "donor": {**parsed["residue_1"], "atom": parsed["atom_1"]},
                        "acceptor": {**parsed["residue_2"], "atom": parsed["atom_2"]},
                        "hydrogen_bond_subtype": _hydrogen_bond_subtype(
                            parsed["atom_1"], parsed["atom_2"]
                        ),
                    }
                )
            contacts.append(contact)
    return contacts, warnings


def _hydrogen_bond_subtype(donor_atom, acceptor_atom):
    donor_main = donor_atom.upper() in BACKBONE_ATOMS
    acceptor_main = acceptor_atom.upper() in BACKBONE_ATOMS
    if donor_main and acceptor_main:
        return "M-M"
    if not donor_main and not acceptor_main:
        return "S-S"
    return "M-S"


def _orient_contact(contact, first_key, second_key):
    oriented = dict(contact)
    key_1 = _position_key(contact["residue_1"])
    key_2 = _position_key(contact["residue_2"])
    if key_1 == first_key and key_2 == second_key:
        return oriented
    if key_1 == second_key and key_2 == first_key:
        oriented["residue_1"], oriented["residue_2"] = contact["residue_2"], contact["residue_1"]
        oriented["atom_1"], oriented["atom_2"] = contact["atom_2"], contact["atom_1"]
        oriented["atom_1_role"], oriented["atom_2_role"] = contact["atom_2_role"], contact["atom_1_role"]
    return oriented


def _contact_sort_key(contact):
    distance = contact.get("distance")
    return (distance is None, distance or 0.0, contact.get("atom_1", ""), contact.get("atom_2", ""))


def _group_contacts(contacts, mutation_key):
    groups = defaultdict(list)
    for contact in contacts:
        key_1 = _position_key(contact["residue_1"])
        key_2 = _position_key(contact["residue_2"])
        if key_1 == key_2:
            # Intra-residue covalent geometry is not a WT-vs-mutant residue
            # interaction and must never be reported as gained/lost.
            continue
        if mutation_key == key_1:
            first_key, second_key = key_1, key_2
        elif mutation_key == key_2:
            first_key, second_key = key_2, key_1
        else:
            first_key, second_key = sorted((key_1, key_2))
        group_key = (first_key, second_key, contact["interaction_class"])
        groups[group_key].append(_orient_contact(contact, first_key, second_key))
    for group in groups.values():
        group.sort(key=_contact_sort_key)
    return groups


def _residue_name(residues, position_key, contacts, side):
    residue = residues.get(position_key)
    if residue and residue.get("residue_name"):
        return residue["residue_name"]
    field = f"residue_{side}"
    for contact in contacts:
        if _position_key(contact[field]) == position_key:
            return contact[field].get("residue_name")
    return None


def _representative_contact(contacts):
    return min(contacts, key=_contact_sort_key) if contacts else None


def _contact_signature(contact):
    return (contact.get("atom_1"), contact.get("atom_2"), round(contact.get("distance") or -1, 3))


def _record_id(group_key):
    raw = json.dumps(group_key, sort_keys=True, default=list).encode("utf-8")
    return "interaction-" + hashlib.sha256(raw).hexdigest()[:14]


def compare_interactions(
    wt_contacts,
    mutant_contacts,
    mutation_position,
    wt_residues=None,
    mutant_residues=None,
):
    """Create the one canonical comparison used by every X-Pro consumer."""
    wt_residues = wt_residues or {}
    mutant_residues = mutant_residues or {}
    mutation_key = _position_key(mutation_position)
    wt_groups = _group_contacts(wt_contacts, mutation_key)
    mutant_groups = _group_contacts(mutant_contacts, mutation_key)
    records = []

    for group_key in sorted(set(wt_groups) | set(mutant_groups)):
        first_key, second_key, interaction_class = group_key
        wt_group = wt_groups.get(group_key, [])
        mutant_group = mutant_groups.get(group_key, [])
        first_position = _position_from_key(first_key)
        second_position = _position_from_key(second_key)
        wt_name_1 = _residue_name(wt_residues, first_key, wt_group, 1)
        wt_name_2 = _residue_name(wt_residues, second_key, wt_group, 2)
        mutant_name_1 = _residue_name(mutant_residues, first_key, mutant_group, 1)
        mutant_name_2 = _residue_name(mutant_residues, second_key, mutant_group, 2)
        wt_representative = _representative_contact(wt_group)
        mutant_representative = _representative_contact(mutant_group)

        if wt_group and mutant_group:
            result = "Retained"
        elif mutant_group:
            result = "Gained"
        else:
            result = "Lost"

        wt_count = len(wt_group)
        mutant_count = len(mutant_group)
        record = {
            "id": _record_id(group_key),
            "position_1": first_position,
            "position_2": second_position,
            "chain_wt": first_position["chain"],
            "residue_number": first_position["residue_number"],
            "insertion_code": first_position["insertion_code"],
            "residue_wt": wt_name_1,
            "residue_mutant": mutant_name_1,
            "partner_chain": second_position["chain"],
            "partner_residue": second_position["residue_number"],
            "partner_insertion_code": second_position["insertion_code"],
            "partner_residue_name": mutant_name_2 or wt_name_2,
            "partner_residue_name_wt": wt_name_2,
            "partner_residue_name_mutant": mutant_name_2,
            "wt_residue_pair": (
                f"{_residue_label(first_position, wt_name_1)}–{_residue_label(second_position, wt_name_2)}"
                if wt_group
                else None
            ),
            "mutant_residue_pair": (
                f"{_residue_label(first_position, mutant_name_1)}–{_residue_label(second_position, mutant_name_2)}"
                if mutant_group
                else None
            ),
            "wt_atom_1": wt_representative.get("atom_1") if wt_representative else None,
            "wt_atom_2": wt_representative.get("atom_2") if wt_representative else None,
            "mutant_atom_1": mutant_representative.get("atom_1") if mutant_representative else None,
            "mutant_atom_2": mutant_representative.get("atom_2") if mutant_representative else None,
            "wt_distance": wt_representative.get("distance") if wt_representative else None,
            "mutant_distance": mutant_representative.get("distance") if mutant_representative else None,
            "interaction_class": interaction_class,
            "interaction_type": INTERACTION_TYPES[interaction_class]["label"],
            "interaction_source": INTERACTION_TYPES[interaction_class]["source"],
            "wt_contact_count": wt_count,
            "mutant_contact_count": mutant_count,
            "delta_contact_count": mutant_count - wt_count,
            "result": result,
            "geometry_or_atoms_changed": bool(
                wt_group
                and mutant_group
                and (
                    wt_count != mutant_count
                    or [_contact_signature(item) for item in wt_group]
                    != [_contact_signature(item) for item in mutant_group]
                )
            ),
            "wt_contacts": wt_group,
            "mutant_contacts": mutant_group,
        }
        records.append(record)
    return records


def _backbone_validation(structure, mutation_key, state):
    residue = structure["residues"].get(mutation_key)
    if not residue:
        return {
            "state": state,
            "residue_found": False,
            "record_types": [],
            "present_backbone_atoms": [],
            "missing_backbone_atoms": ["N", "CA", "C", "O"],
        }
    required = {"N", "CA", "C", "O"}
    present = required & residue["atoms"]
    return {
        "state": state,
        "residue_found": True,
        "record_types": sorted(residue["record_types"]),
        "present_backbone_atoms": sorted(present),
        "missing_backbone_atoms": sorted(required - present),
    }


def build_canonical_comparison(
    wt_pdb,
    mutant_pdb,
    wt_hhb,
    mutant_hhb,
    wt_nnb,
    mutant_nnb,
    mutation_chain,
    mutation_residue,
    expected_mutant_name=None,
):
    """Parse engine outputs and build a fully provenance-labelled dataset."""
    residue_number, insertion_code = split_residue_identifier(mutation_residue)
    mutation_position = _position(mutation_chain, residue_number, insertion_code)
    mutation_key = _position_key(mutation_position)
    wt_structure = parse_pdb_structure(wt_pdb)
    mutant_structure = parse_pdb_structure(mutant_pdb)

    wt_hbonds, wt_hbond_warnings = parse_ligplot_contact_file(wt_hhb, "hydrogen_bond", wt_structure)
    mutant_hbonds, mutant_hbond_warnings = parse_ligplot_contact_file(
        mutant_hhb, "hydrogen_bond", mutant_structure
    )
    wt_nonbonded, wt_nonbonded_warnings = parse_ligplot_contact_file(
        wt_nnb, "non_bonded_contact", wt_structure
    )
    mutant_nonbonded, mutant_nonbonded_warnings = parse_ligplot_contact_file(
        mutant_nnb, "non_bonded_contact", mutant_structure
    )
    records = compare_interactions(
        wt_hbonds + wt_nonbonded,
        mutant_hbonds + mutant_nonbonded,
        mutation_position,
        wt_structure["residues"],
        mutant_structure["residues"],
    )

    wt_residue = wt_structure["residues"].get(mutation_key, {})
    mutant_residue_data = mutant_structure["residues"].get(mutation_key, {})
    warnings = (
        wt_structure["warnings"]
        + mutant_structure["warnings"]
        + wt_hbond_warnings
        + mutant_hbond_warnings
        + wt_nonbonded_warnings
        + mutant_nonbonded_warnings
    )
    wt_validation = _backbone_validation(wt_structure, mutation_key, "WT")
    mutant_validation = _backbone_validation(mutant_structure, mutation_key, "mutant")
    for validation in (wt_validation, mutant_validation):
        if validation["missing_backbone_atoms"]:
            warnings.append(
                f"{validation['state']} mutation site is missing backbone atoms: "
                + ", ".join(validation["missing_backbone_atoms"])
            )
        if validation["residue_found"] and "ATOM" not in validation["record_types"]:
            warnings.append(f"{validation['state']} mutation site is not represented by ATOM records.")
    actual_mutant_name = mutant_residue_data.get("residue_name")
    if expected_mutant_name and actual_mutant_name and actual_mutant_name != expected_mutant_name.upper():
        warnings.append(
            f"Mutant PDB contains {actual_mutant_name} at {_position_text(mutation_position)}, "
            f"but {expected_mutant_name.upper()} was requested."
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "authoritative": True,
        "description": "Canonical WT-vs-mutant local interaction comparison",
        "model_pair": {
            "wt": {
                "filename": Path(wt_pdb).name,
                "sha256": _sha256_file(wt_pdb),
            },
            "mutant": {
                "filename": Path(mutant_pdb).name,
                "sha256": _sha256_file(mutant_pdb),
            },
            "same_pair_used_for_all_records": True,
        },
        "mutation": {
            **mutation_position,
            "residue_wt": wt_residue.get("residue_name"),
            "residue_mutant": actual_mutant_name,
        },
        "scientific_sources": {
            key: value["source"] for key, value in INTERACTION_TYPES.items()
        },
        "distance_method": (
            "Euclidean distance from the corresponding first-model WT or mutant PDB coordinates; "
            "engine-reported distance is retained only when an atom coordinate is unavailable."
        ),
        "summary_contact_rule": (
            "The shortest engine-detected atom contact supplies the summary atom pair and distance; "
            "wt_contacts and mutant_contacts retain every atom-level contact."
        ),
        "comparison_identity": (
            "Structural residue-position pair (chain, residue number, insertion code) plus interaction class; "
            "residue names are labels, while atom identities remain in each contact detail."
        ),
        "atom_contact_identity": (
            "Both structural positions, both actual atom names, and interaction class; each has a stable contact_id."
        ),
        "classification_rule": {
            "Gained": "Residue-position pair and interaction class absent in WT, present in mutant",
            "Lost": "Residue-position pair and interaction class present in WT, absent in mutant",
            "Retained": "Residue-position pair and interaction class present in both structures",
        },
        "backbone_validation": {"wt": wt_validation, "mutant": mutant_validation},
        "records": records,
        "warnings": sorted(set(warnings)),
    }


CSV_FIELDS = [
    "id",
    "chain_wt",
    "residue_number",
    "insertion_code",
    "residue_wt",
    "residue_mutant",
    "partner_chain",
    "partner_residue",
    "partner_insertion_code",
    "partner_residue_name",
    "wt_atom_1",
    "wt_atom_2",
    "mutant_atom_1",
    "mutant_atom_2",
    "wt_distance",
    "mutant_distance",
    "interaction_type",
    "wt_contact_count",
    "mutant_contact_count",
    "delta_contact_count",
    "result",
]


def write_comparison_exports(dataset, json_filename, csv_filename):
    Path(json_filename).write_text(json.dumps(dataset, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    with Path(csv_filename).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(dataset.get("records", []))


def write_text_summary(dataset, filename):
    groups = defaultdict(list)
    for record in dataset.get("records", []):
        groups[record["result"]].append(record)
    lines = [
        "X-PRO CANONICAL WT-VS-MUTANT INTERACTION COMPARISON",
        "The JSON/CSV canonical dataset is authoritative for the interaction table.",
        "",
    ]
    for result in ("Gained", "Lost", "Retained"):
        lines.append(result.upper())
        if not groups[result]:
            lines.append("None")
        for record in groups[result]:
            pair = record.get("mutant_residue_pair") or record.get("wt_residue_pair") or record["id"]
            lines.append(
                f"{pair} | {record['interaction_type']} | contacts "
                f"{record['wt_contact_count']} -> {record['mutant_contact_count']} | "
                f"delta {record['delta_contact_count']:+d}"
            )
        lines.append("")
    if dataset.get("warnings"):
        lines.append("WARNINGS / EXTERNAL-TOOL LIMITATIONS")
        lines.extend(f"- {warning}" for warning in dataset["warnings"])
    Path(filename).write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def create_comparison_outputs(
    output_dir,
    pdb_name,
    wt_pdb,
    mutant_pdb,
    mutation_chain,
    mutation_residue,
    expected_mutant_name=None,
    extra_warnings=None,
    annotate_ligplots=False,
):
    """Build the authoritative canonical dataset and export it once."""
    output_dir = Path(output_dir)
    dataset = build_canonical_comparison(
        wt_pdb=wt_pdb,
        mutant_pdb=mutant_pdb,
        wt_hhb=output_dir / f"altered_{pdb_name}_ligplot.hhb",
        mutant_hhb=output_dir / f"mutated_{pdb_name}_ligplot.hhb",
        wt_nnb=output_dir / f"altered_{pdb_name}_ligplot.nnb",
        mutant_nnb=output_dir / f"mutated_{pdb_name}_ligplot.nnb",
        mutation_chain=mutation_chain,
        mutation_residue=mutation_residue,
        expected_mutant_name=expected_mutant_name,
    )
    if extra_warnings:
        dataset["warnings"] = sorted(set(dataset["warnings"] + list(extra_warnings)))
    if annotate_ligplots:
        dataset["ligplot_annotations"] = annotate_ligplot_pair(
            dataset, output_dir, pdb_name
        )
    write_comparison_exports(
        dataset,
        output_dir / f"interaction_comparison_{pdb_name}.json",
        output_dir / f"interaction_comparison_{pdb_name}.csv",
    )
    write_text_summary(dataset, output_dir / "ligplot_differences.txt")
    return dataset


# Red common-residue markers are annotations, not bounding boxes: these
# PostScript-point limits keep every circle readable and stop a large side
# chain's marker from covering the diagram's chemistry.
LIGPLOT_ANNOTATION_MIN_RADIUS = 24.0
LIGPLOT_ANNOTATION_MAX_RADIUS = 44.0


def _canonical_state_residues(dataset):
    """Return canonical interacting residue identities present in each state."""
    residues = {"wt": set(), "mutant": set()}
    for record in dataset.get("records", []):
        for state in residues:
            if int(record.get(f"{state}_contact_count") or 0) <= 0:
                continue
            names = (
                record.get("residue_wt" if state == "wt" else "residue_mutant"),
                record.get(
                    "partner_residue_name_wt"
                    if state == "wt"
                    else "partner_residue_name_mutant"
                ),
            )
            for position, residue_name in zip(
                (record["position_1"], record["position_2"]), names
            ):
                if residue_name:
                    residues[state].add(
                        _position_key(position) + (str(residue_name).upper(),)
                    )
    return residues


def _parse_ligplot_drawing(filename):
    """Read LigPlot's flattened per-residue atom geometry from its .drw file."""
    residues = {}
    current = None
    expect_residue = False
    reading_atoms = False
    for raw_line in Path(filename).read_text(errors="replace").splitlines():
        line = raw_line.rstrip()
        if line == "#R":
            expect_residue = True
            reading_atoms = False
            continue
        if expect_residue:
            parts = line.split()
            if len(parts) < 4 or not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", parts[0]):
                raise ValueError(f"Invalid LigPlot residue geometry in {filename}: {line!r}")
            mode = int(parts[0][0])
            residue_name = parts[0][1:].upper()
            residue_number, insertion_code = split_residue_identifier(parts[1])
            if len(parts) >= 5 and not re.fullmatch(r"[-+]?\d+(?:\.\d+)?", parts[2]):
                chain, x_index = parts[2].strip(), 3
            else:
                chain, x_index = "", 2
            current = {
                "mode": mode,
                "residue_name": residue_name,
                "position": _position(chain, residue_number, insertion_code),
                "center": (float(parts[x_index]), float(parts[x_index + 1])),
                "atoms": [],
            }
            residues[_position_key(current["position"])] = current
            expect_residue = False
            continue
        if line == "#A":
            reading_atoms = current is not None
            continue
        if line.startswith("#"):
            reading_atoms = False
            continue
        if reading_atoms and line.strip():
            parts = line.split()
            if len(parts) >= 3:
                current["atoms"].append((parts[0], float(parts[1]), float(parts[2])))
    if not residues:
        raise ValueError(f"No LigPlot residue geometry found in {filename}")
    return residues


def _parse_ligplot_postscript(filename):
    """Read exact atom and residue-label coordinates from LigPlot PostScript."""
    text = Path(filename).read_text(errors="replace")
    atom_section_match = re.search(
        r"% Atoms\s*(.*?)(?=\n% Hydrophobic interactions)", text, re.DOTALL
    )
    if not atom_section_match:
        raise ValueError(f"LigPlot atom section missing from {filename}")
    atom_points = [
        (float(match.group(1)), float(match.group(2)))
        for match in re.finditer(
            r"^\s*([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+"
            r"(?:Ligatom_radius|Nligatom_radius)\s+Sphere\s*$",
            atom_section_match.group(1),
            re.MULTILINE,
        )
    ]
    labels = {}
    label_pattern = re.compile(
        r"^\s*([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+moveto\s*$\n"
        r"\(([A-Za-z]{3})(-?\d+)([A-Za-z]?)\(([^()]*)\)\)\s+"
        r"(?:Ligresname_size|Nligresnam_size|Hydrophnam_size)\s+Center",
        re.MULTILINE,
    )
    for match in label_pattern.finditer(text):
        key = (
            match.group(6).strip(),
            int(match.group(4)),
            match.group(5).upper(),
        )
        labels[key] = {
            "point": (float(match.group(1)), float(match.group(2))),
            "text": f"{match.group(3)}{match.group(4)}{match.group(5)}({match.group(6).strip()})",
        }
    return text, atom_points, labels


def _fit_axis_transform(source, target):
    """Fit the scale/offset of one unordered LigPlot coordinate axis."""
    if len(source) != len(target) or len(source) < 2:
        raise ValueError("LigPlot coordinate transform requires at least two atom pairs")
    source_span = max(source) - min(source)
    target_span = max(target) - min(target)
    if source_span <= 1e-12 or target_span <= 1e-12:
        raise ValueError("LigPlot coordinate transform is degenerate")
    scale = target_span / source_span
    offset = min(target) - scale * min(source)
    residual = max(
        min(abs(scale * source_value + offset - target_value) for target_value in target)
        for source_value in source
    )
    return scale, offset, residual


def _ellipse_for_ligplot_residue(residue, transform, label):
    """Bound the exact native residue geometry and return a padded ellipse.

    A residue drawn with its full side chain can span most of the diagram, and
    an ellipse that large would sit on top of the ligand chemistry it is meant
    to annotate. Such an ellipse is therefore clamped and re-centred on the
    residue's own LigPlot label, which is what identifies the residue.
    """
    scale_x, offset_x, scale_y, offset_y = transform

    def map_point(point):
        return scale_x * point[0] + offset_x, scale_y * point[1] + offset_y

    def clamp(radius):
        return min(
            LIGPLOT_ANNOTATION_MAX_RADIUS,
            max(LIGPLOT_ANNOTATION_MIN_RADIUS, radius),
        )

    label_half_width = (
        max(18.0, len(label["text"]) * 4.3 / 2.0) if label else 20.0
    )

    if residue["mode"] == 3:
        center = label["point"] if label else map_point(residue["center"])
        return center[0], center[1], clamp(label_half_width + 7.0), 30.0

    points = [map_point((x, y)) for _, x, y in residue["atoms"]]
    if not points:
        center = label["point"] if label else map_point(residue["center"])
        return center[0], center[1], 30.0, 30.0
    xs, ys = [point[0] for point in points], [point[1] for point in points]
    left, right = min(xs) - 6.0, max(xs) + 6.0
    bottom, top = min(ys) - 6.0, max(ys) + 6.0
    if label:
        left = min(left, label["point"][0] - label_half_width)
        right = max(right, label["point"][0] + label_half_width)
        bottom = min(bottom, label["point"][1] - 7.0)
        top = max(top, label["point"][1] + 7.0)
    padding = 9.0
    natural_x = (right - left) / 2.0 + padding
    natural_y = (top - bottom) / 2.0 + padding
    center_x, center_y = (left + right) / 2.0, (bottom + top) / 2.0
    radius_x, radius_y = clamp(natural_x), clamp(natural_y)
    if label and (
        natural_x > LIGPLOT_ANNOTATION_MAX_RADIUS
        or natural_y > LIGPLOT_ANNOTATION_MAX_RADIUS
    ):
        center_x, center_y = label["point"]
        radius_x = max(radius_x, clamp(label_half_width + 8.0))
    return center_x, center_y, radius_x, radius_y


def _render_postscript_to_png(ps_path, png_path, density=150):
    """Rasterize a PostScript file to PNG.

    Deliberately mirrors X-Pro's own enhanced_ps_to_png() exactly -- same
    tools (ImageMagick ``convert`` first, Ghostscript ``gs`` fallback), same
    settings, and critically the same shell=True invocation style. A
    no-shell subprocess.run([...]) list form relies purely on the calling
    process's inherited PATH to locate the binaries via execvp(); under a
    restricted web-server environment (e.g. PHP exec() spawning www-data)
    that PATH can be minimal enough that the list form raises
    FileNotFoundError while `sh -c "..."` still resolves the same binaries
    via the shell's own built-in fallback search path. enhanced_ps_to_png()
    already works in that environment, so matching it exactly here avoids
    silently failing in a way enhanced_ps_to_png() wouldn't.

    Returns True if a PNG was actually written, False otherwise (never
    raises -- rendering is a display convenience, not part of the canonical
    comparison dataset -- but every failure is printed so it isn't silently
    invisible in job logs).
    """
    ps_path, png_path = str(ps_path), str(png_path)
    try:
        cmd_convert = f"convert -density {density} {ps_path} {png_path}"
        result = subprocess.run(cmd_convert, shell=True, capture_output=True, timeout=30)
        if result.returncode == 0 and Path(png_path).exists():
            return True
        print(f"[LigPlot annotation] ImageMagick convert failed for {ps_path} "
              f"(exit {result.returncode}): {result.stderr.decode(errors='replace').strip()}")
    except Exception as e:
        print(f"[LigPlot annotation] ImageMagick convert not available for {ps_path}: {e}")
    try:
        cmd_gs = (f"gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m "
                  f"-r{density} -sOutputFile={png_path} {ps_path}")
        result = subprocess.run(cmd_gs, shell=True, capture_output=True, timeout=30)
        if result.returncode == 0 and Path(png_path).exists():
            return True
        print(f"[LigPlot annotation] Ghostscript failed for {ps_path} "
              f"(exit {result.returncode}): {result.stderr.decode(errors='replace').strip()}")
    except Exception as e:
        print(f"[LigPlot annotation] Ghostscript not available for {ps_path}: {e}")
    return False


def annotate_ligplot_postscript(dataset, state, postscript, drawing, output_postscript):
    """Add red common-environment ellipses behind native LigPlot chemistry.

    Residue commonality comes only from the canonical comparison. Geometry
    comes only from LigPlot's .drw and PostScript files; no pixels or OCR are
    inspected and no interaction is recalculated.
    """
    if state not in {"wt", "mutant"}:
        raise ValueError(f"Unsupported LigPlot annotation state: {state}")
    residues = _parse_ligplot_drawing(drawing)
    ps_text, ps_atom_points, labels = _parse_ligplot_postscript(postscript)
    drawn_atoms = [
        (x, y)
        for residue in residues.values()
        if residue["mode"] != 3
        for _, x, y in residue["atoms"]
    ]
    if len(drawn_atoms) != len(ps_atom_points):
        raise ValueError(
            f"LigPlot atom geometry mismatch for {postscript}: "
            f"{len(drawn_atoms)} drawing atoms vs {len(ps_atom_points)} PostScript atoms"
        )
    scale_x, offset_x, _ = _fit_axis_transform(
        [point[0] for point in drawn_atoms], [point[0] for point in ps_atom_points]
    )
    scale_y, offset_y, _ = _fit_axis_transform(
        [point[1] for point in drawn_atoms], [point[1] for point in ps_atom_points]
    )
    geometry_error = max(
        min(
            math.hypot(
                scale_x * source_x + offset_x - target_x,
                scale_y * source_y + offset_y - target_y,
            )
            for target_x, target_y in ps_atom_points
        )
        for source_x, source_y in drawn_atoms
    )
    if geometry_error > 0.15:
        raise ValueError(
            f"LigPlot coordinate transform residual is too large ({geometry_error:.3f})"
        )

    state_residues = _canonical_state_residues(dataset)
    mutation = dataset.get("mutation", {})
    mutation_key = (
        (mutation.get("chain") or "").strip(),
        int(mutation.get("residue_number") or 0),
        (mutation.get("insertion_code") or "").strip().upper(),
    )
    common_identities = sorted(
        identity
        for identity in (state_residues["wt"] & state_residues["mutant"])
        if identity[:3] != mutation_key
    )
    common = [identity[:3] for identity in common_identities]
    missing = [
        identity
        for identity in common_identities
        if identity[:3] not in residues
        or residues[identity[:3]]["residue_name"] != identity[3]
    ]
    if missing:
        missing_text = ", ".join(
            f"{identity[0] or '_'}:{identity[3]}{identity[1]}{identity[2]}"
            for identity in missing
        )
        raise ValueError(f"Common canonical residues missing from LigPlot geometry: {missing_text}")

    transform = (scale_x, offset_x, scale_y, offset_y)
    ellipses = [
        (key, _ellipse_for_ligplot_residue(residues[key], transform, labels.get(key)))
        for key in common
    ]
    commands = [
        "% X-Pro canonical common-environment annotations",
        "% Red no-fill ellipses are geometry overlays; native LigPlot chemistry follows unchanged.",
        "gsave",
        "1.000 0.000 0.000 setrgbcolor",
        "1.250 setlinewidth",
        "[] 0 setdash",
    ]
    for key, (center_x, center_y, radius_x, radius_y) in ellipses:
        position_text = f"{key[0] or '_'}:{key[1]}{key[2]}"
        commands.extend([
            f"% X-Pro common residue {position_text}",
            "matrix currentmatrix",
            f"{center_x:.3f} {center_y:.3f} translate",
            f"{radius_x:.3f} {radius_y:.3f} scale",
            "newpath 0 0 1 0 360 arc closepath",
            "setmatrix stroke",
        ])
    commands.extend(["grestore", "% End X-Pro annotations", ""])
    marker = "% Atoms"
    if marker not in ps_text:
        raise ValueError(f"LigPlot interaction marker missing from {postscript}")
    annotated = ps_text.replace(marker, "\n".join(commands) + marker, 1)

    if state == "mutant":
        # This annotated copy is what X-ProResult.php displays as the third
        # ("Differences-Ligplot") panel, while the plain (non-annotated)
        # mutant .ps -- sharing this same native LigPlot run -- is displayed
        # as the second ("Mutant-Ligplot") panel. Both would otherwise show
        # the identical "Mutant-Ligplot" caption LigPlot baked into the
        # diagram (see enhanced_ligplot_interface() in xpro.py, which is what
        # controls that text via the input filename it feeds LigPlot). Swap
        # it here, on this annotated copy only, so the two panels read
        # differently despite coming from one LigPlot run.
        annotated = annotated.replace("Mutant-Ligplot", "Differences-Ligplot")

    Path(output_postscript).write_text(annotated, encoding="utf-8")

    # Render the annotated PostScript to PNG. This is what makes the circles
    # actually visible anywhere -- the PostScript file alone isn't something
    # a browser can display. This writes ONLY the standalone
    # "*_annotated.png"; it must NOT touch the canonical "*_ligplot.png"
    # that enhanced_ps_to_png() already produced from the plain PostScript
    # -- that file is the plain/native diagram shown elsewhere, and
    # overwriting it would put circles on views that are supposed to stay
    # circle-free.
    output_png = Path(output_postscript).with_suffix(".png")
    png_rendered = _render_postscript_to_png(output_postscript, output_png)
    if not png_rendered:
        print(f"[LigPlot annotation] Could not rasterize {output_postscript} to PNG; "
              f"the annotated PostScript was still written, but no image is available.")

    return {
        "source": "canonical interaction records + native LigPlot .drw/PostScript geometry",
        "common_positions": [
            {
                "chain": identity[0],
                "residue_number": identity[1],
                "insertion_code": identity[2],
                "residue_name": identity[3],
            }
            for identity in common_identities
        ],
        "circle_count": len(ellipses),
        "coordinate_transform_max_error": round(geometry_error, 6),
        "native_chemistry_preserved": True,
        "output_postscript": Path(output_postscript).name,
        "output_png": output_png.name if png_rendered else None,
    }


def annotate_ligplot_pair(dataset, output_dir, pdb_name):
    """Annotate both native diagrams from one canonical common-position set.

    Writes output_dir/altered_{pdb_name}_ligplot_annotated.png and
    output_dir/mutated_{pdb_name}_ligplot_annotated.png (plus the matching
    .ps files). The plain output_dir/altered_{pdb_name}_ligplot.png and
    mutated_{pdb_name}_ligplot.png that enhanced_ps_to_png() already wrote
    are left untouched -- callers that want the plain diagram and callers
    that want the circled one read two different files. Call this only
    after enhanced_ps_to_png() has run for this job.
    """
    output_dir = Path(output_dir)
    manifests = {}
    for state, prefix in (("wt", "altered"), ("mutant", "mutated")):
        base = output_dir / f"{prefix}_{pdb_name}_ligplot"
        manifests[state] = annotate_ligplot_postscript(
            dataset,
            state,
            base.with_suffix(".ps"),
            base.with_suffix(".drw"),
            output_dir / f"{prefix}_{pdb_name}_ligplot_annotated.ps",
        )
    if manifests["wt"]["common_positions"] != manifests["mutant"]["common_positions"]:
        raise ValueError("WT and mutant LigPlot common-position annotations disagree")
    return manifests
