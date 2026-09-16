#!/usr/bin/python3
"""
automate_point_mutation.py

Structure-only stability automation for the X-Pro -> pathogenicity-
prediction bridge (xpro_to_mutxplor.php).

This lives in /var/www/html/pointMutationScript -- deliberately NOT in
/var/www/html/mutXplorScripts, and it deliberately does NOT touch
mutXplorScripts/automate.py. X-Pro only ever hands this script a single
PDB-based mutation with no FASTA, so all of automate.py's FASTA-generation,
sequence-tool, chain-offset and tool-selection machinery simply doesn't
apply here -- this script only ever runs the three structural stability
tools:

    - FoldX             -> outFiles/foldX.out
    - DynaMut2          -> outFiles/dynamut2.out
    - DUET (one submission, three results):
          -> outFiles/duet.out
          -> outFiles/mCSM.out
          -> outFiles/sdm.out
    - HOPE (best-effort reachability check only -- see auto_hope.py's own
      header for why this doesn't produce a real prediction yet):
          -> outFiles/hope_status.json

Usage (same argument shape the PHP bridge invokes it with):
    python3 automate_point_mutation.py <randNum> <chainID> <pdbID>

Expects, inside /var/www/html/fileUpload/<randNum>/:
    - a PDB file already copied there by the PHP bridge (named
      <pdbID>.pdb, or any *.pdb/*.ent file already present)
    - pending_mutation.txt, written by the PHP bridge as
      "<residueNumber>:<targetAA single-letter>"

On completion (success or failure of individual tools) this writes
job_complete.txt into the job directory, which the PHP status page uses to
tell "still running" apart from "finished, but this tool produced no
result".
"""

import sys
import os
import re
import glob
import time
import warnings

# ── Make sure the sibling auto_*.py modules (in this same directory) and
#    the usual site-packages locations are importable regardless of how
#    this script is launched ────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/chandreyee/.local/lib/python3.10/site-packages')
sys.path.insert(0, '/usr/local/lib/python3.10/dist-packages')
sys.path.insert(0, '/usr/lib/python3/dist-packages')

os.environ.setdefault('http_proxy', 'http://sthiyaga%40ibab%2Eac%2Ein:Ibab2025@proxy.ibab.ac.in:3128')
os.environ.setdefault('https_proxy', 'http://sthiyaga%40ibab%2Eac%2Ein:Ibab2025@proxy.ibab.ac.in:3128')
os.environ.setdefault('HTTP_PROXY', 'http://sthiyaga%40ibab%2Eac%2Ein:Ibab2025@proxy.ibab.ac.in:3128')
os.environ.setdefault('HTTPS_PROXY', 'http://sthiyaga%40ibab%2Eac%2Ein:Ibab2025@proxy.ibab.ac.in:3128')

from multiprocessing import Process

warnings.filterwarnings("ignore")

try:
    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1
except ImportError as e:
    print(f"ERROR: Cannot import BioPython: {e}")
    sys.exit(1)

try:
    import mechanicalsoup  # noqa: F401 -- required by auto_duet / auto_dynamut2
except ImportError as e:
    print(f"ERROR: Cannot import mechanicalsoup: {e}")
    sys.exit(1)

# Structure tools this bridge runs -- copies of the MutXplor tool scripts,
# relocated to this same directory (see each file's own header comment).
import auto_duet
import auto_dynamut2
import auto_foldx
import auto_hope

STANDARD_AA_3LETTER = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
}


def normalize_hetatm_standard_residues(pdb_path):
    """
    Some structures flag genuine standard residues as HETATM (e.g. a
    disordered/low-occupancy terminal residue -- observed in 1ATP chain E
    residue 15). Re-flag them as ATOM so both our own resolver below and
    the external prediction servers see a correctly-formed file.
    (Mirrors the equivalent fix in mutXplorScripts/automate.py.)
    """
    try:
        with open(pdb_path, "r") as fh:
            lines = fh.readlines()
        changed = 0
        for i, line in enumerate(lines):
            if line.startswith("HETATM"):
                resname = line[17:20].strip().upper()
                if resname in STANDARD_AA_3LETTER:
                    lines[i] = "ATOM  " + line[6:]
                    changed += 1
        if changed:
            with open(pdb_path, "w") as fh:
                fh.writelines(lines)
            print(f"[NORMALIZE] Converted {changed} HETATM record(s) to ATOM in {pdb_path}")
    except Exception as e:
        print(f"[NORMALIZE] Could not normalize {pdb_path}: {e}")


def resolve_pending_mutation(pending_file, pdb_file, chain_id, out_file="mutationList.txt"):
    """
    Read pending_mutation.txt ("<residueNumber>:<targetAA>"), look up the
    wild-type residue at that PDB author sequence number in chain_id,
    write the standard mutation string (e.g. "V15A") to out_file, and
    return it. Returns None on failure.
    """
    try:
        with open(pending_file) as fh:
            content = fh.read().strip()

        parts = content.split(":")
        if len(parts) != 2:
            print(f"[MUTATION] Invalid pending_mutation.txt format: '{content}'")
            return None

        resi_num = int(parts[0].strip())
        target_aa = parts[1].strip().upper()

        print(f"[MUTATION] Resolving pending mutation: position {resi_num} -> {target_aa} in chain {chain_id}")

        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        structure = parser.get_structure("struct", pdb_file)

        for model in structure:
            if chain_id not in [c.id for c in model]:
                continue
            for residue in model[chain_id]:
                if residue.id[1] != resi_num:
                    continue
                if residue.resname.strip().upper() not in STANDARD_AA_3LETTER:
                    continue
                orig_aa = seq1(residue.resname)
                mutation_str = f"{orig_aa}{resi_num}{target_aa}"
                with open(out_file, "w") as fh:
                    fh.write(mutation_str + "\n")
                print(f"[MUTATION] Resolved and written: {mutation_str}")
                return mutation_str
            break

        print(f"[MUTATION] Residue {resi_num} not found in chain {chain_id} of {pdb_file}")
        return None

    except Exception as e:
        print(f"[MUTATION] Error resolving pending mutation: {e}")
        return None


def main():
    if len(sys.argv) < 4:
        print("Error: Insufficient arguments. Expected: randNum chainID pdbID")
        sys.exit(1)

    rndNum, chainID, pdbID = sys.argv[1], sys.argv[2], sys.argv[3]

    print("Starting X-Pro point-mutation stability automation...")
    print(f"Job ID: {rndNum}, Chain: {chainID}, PDB ID: {pdbID}")

    upload_dir = f"/var/www/html/fileUpload/{rndNum}"
    if not os.path.exists(upload_dir):
        print(f"Error: Directory {upload_dir} does not exist")
        sys.exit(1)

    os.chdir(upload_dir)
    print(f"Working directory: {os.getcwd()}")

    os.makedirs("outFiles", exist_ok=True)
    os.makedirs("logFiles", exist_ok=True)

    pdb_files = glob.glob("*.pdb") + glob.glob("*.ent") + glob.glob("*.PDB") + glob.glob("*.ENT")
    if not pdb_files:
        print(f"Error: No PDB file found in {os.getcwd()}")
        sys.exit(1)
    pdbFile = pdb_files[0]
    print(f"Found PDB file: {pdbFile}")

    normalize_hetatm_standard_residues(pdbFile)

    pending_file = "pending_mutation.txt"
    if not os.path.exists(pending_file):
        print("Error: pending_mutation.txt not found; nothing to run.")
        sys.exit(1)

    mutation_str = resolve_pending_mutation(pending_file, pdbFile, chainID)
    if not mutation_str:
        print("Error: Could not resolve the pending mutation from the PDB.")
        sys.exit(1)
    try:
        os.remove(pending_file)
    except OSError:
        pass

    mutationFile = "mutationList.txt"

    print("** Automation in progress.... **")
    start_time = time.time()

    processes = [
        Process(target=auto_foldx.runFoldx, args=(pdbFile, mutationFile, chainID, rndNum)),
        Process(target=auto_dynamut2.runDynamut, args=(pdbFile, mutationFile, chainID, rndNum)),
        Process(target=auto_duet.runDuet, args=(pdbFile, mutationFile, chainID, rndNum)),
        Process(target=auto_hope.runHope, args=(pdbFile, mutationFile, chainID, rndNum)),
    ]
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    execution_time = time.time() - start_time
    minutes = int(execution_time // 60)
    seconds = int(execution_time % 60)
    print(f"** Automation process completed in {minutes} minutes {seconds} seconds **")

    # Lets the status page distinguish "still running" from "finished, but
    # this particular tool never wrote a result row".
    with open("job_complete.txt", "w") as fh:
        fh.write("done\n")


if __name__ == "__main__":
    main()
