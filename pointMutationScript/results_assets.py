"""Reproducible, path-independent result assets for X-Pro downloads."""

import argparse
from pathlib import Path


def _residue_name(pdb_path, chain, residue_number, insertion_code=""):
    """Return the actual residue name at a chain-aware PDB position."""
    path = Path(pdb_path)
    if not path.is_file():
        return None
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[21:22].strip() != str(chain).strip():
            continue
        if line[22:26].strip() != str(residue_number).strip():
            continue
        if line[26:27].strip() != str(insertion_code).strip():
            continue
        return line[17:20].strip().upper() or None
    return None


def generate_pymol_script(
    output_path,
    wt_pdb,
    mutant_pdb,
    chain,
    residue_number,
    insertion_code="",
):
    """Write a PML file that loads the exact analyzed WT and mutant PDBs.

    The PML contains only relative basenames. The web UI downloads the PDBs as
    those same names, so placing the three downloaded files together is enough
    to reproduce the setup without leaking a server filesystem path.
    """
    output_path = Path(output_path)
    wt_pdb = Path(wt_pdb)
    mutant_pdb = Path(mutant_pdb)
    wt_name = _residue_name(wt_pdb, chain, residue_number, insertion_code) or "WT"
    mutant_name = _residue_name(mutant_pdb, chain, residue_number, insertion_code) or "MUT"
    residue_id = f"{residue_number}{insertion_code}".strip()
    chain_value = str(chain).strip()
    chain_clause = f"chain {chain_value} and " if chain_value else "chain '' and "

    output_path.parent.mkdir(parents=True, exist_ok=True)
    script = f"""# X-Pro WT-vs-mutant visualization
# Keep this file beside {wt_pdb.name} and {mutant_pdb.name}.
# These are the downloadable X-Pro structures; see canonical JSON for the exact analysis scope.
reinitialize
load {wt_pdb.name}, wt
load {mutant_pdb.name}, mutant
align mutant, wt
hide everything, all
show cartoon, wt or mutant
color cyan, wt
color lightblue, mutant
select wt_mutation_site, wt and {chain_clause}resi {residue_id}
select mutant_mutation_site, mutant and {chain_clause}resi {residue_id}
show sticks, wt_mutation_site or mutant_mutation_site
color orange, wt_mutation_site
color magenta, mutant_mutation_site
label wt_mutation_site and name CA, "{wt_name}{residue_id}"
label mutant_mutation_site and name CA, "{wt_name}{residue_id} -> {mutant_name}{residue_id}"
set label_color, black
set label_outline_color, white
bg_color white
orient wt_mutation_site or mutant_mutation_site
zoom wt_mutation_site or mutant_mutation_site, 10
"""
    output_path.write_text(script)
    return output_path

def main():
    parser = argparse.ArgumentParser(description="Generate an X-Pro PyMOL script")
    parser.add_argument("output")
    parser.add_argument("wt_pdb")
    parser.add_argument("mutant_pdb")
    parser.add_argument("chain")
    parser.add_argument("residue_number")
    parser.add_argument("--insertion-code", default="")
    args = parser.parse_args()
    generate_pymol_script(
        args.output,
        args.wt_pdb,
        args.mutant_pdb,
        args.chain,
        args.residue_number,
        args.insertion_code,
    )


if __name__ == "__main__":
    main()
