#!/usr/bin/env python3
import sys
import warnings
import requests
import os
import json
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
from Bio.PDB.Polypeptide import is_aa

# Ignore Biopython warnings
warnings.filterwarnings("ignore")

# ----- PROXY SETTINGS WITH AUTHENTICATION -----
proxies = {
    "http": "http://chandreyeenandi2001%40gmail.com:Chandreyee%401234@proxy.ibab.ac.in:3128",
    "https": "http://chandreyeenandi2001%40gmail.com:Chandreyee%401234@proxy.ibab.ac.in:3128"
}

# ----- FUNCTION: Extract Chain and Residue Info -----
def get_chain_ids(struct_file):
    """Parse a PDB or PDBx/mmCIF file to extract chain and residue information."""
    try:
        ext = os.path.splitext(struct_file)[1].lower()

        if ext == ".cif":
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(PERMISSIVE=1, QUIET=True)

        structure = parser.get_structure('PDB_structure', struct_file)
        chains_residue_ids = {}

        for model in structure:
            for chain in model:
                residues = {}
                for residue in chain:
                    if is_aa(residue):
                        # Get residue number as string
                        res_num = str(residue.get_id()[1])
                        # Get residue name (3-letter code)
                        res_name = residue.resname
                        residues[res_num] = res_name
                
                if residues:  # Only add chain if it has residues
                    chains_residue_ids[chain.id] = residues
        
        return chains_residue_ids
    except Exception as e:
        print(f"[ERROR] Failed to parse structure ({ext or 'unknown'} format): {e}", file=sys.stderr)
        return None

# ----- FUNCTION: Convert mmCIF to legacy PDB format -----
def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert a PDBx/mmCIF file to legacy PDB format.

    The downstream mutagenesis pipeline (xpro.py / MODELLER / PyMOL / LigPlot)
    only understands legacy .pdb files, so a .cif upload needs a companion
    .pdb version written alongside it for those tools to consume.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('structure', cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

# ----- FUNCTION: Download and Clean PDB -----
def download_pdb(pdb_id, base_dir, rand_num):
    """Download a PDB file from RCSB, clean it, and save in work directory."""
    output_dir = os.path.join(base_dir, rand_num)
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if this is an uploaded file
    if pdb_id.startswith("UPLOAD:"):
        # Extract the actual filename
        actual_filename = pdb_id.split(":", 1)[1]
        pdb_path = os.path.join(output_dir, actual_filename)
        
        print(f"[INFO] Processing uploaded structure file: {pdb_path}", file=sys.stderr)

        # Only .pdb and .cif (PDBx/mmCIF) uploads are supported
        upload_ext = os.path.splitext(actual_filename)[1].lower()
        if upload_ext not in (".pdb", ".cif"):
            print(f"[ERROR] Unsupported file extension '{upload_ext}'. Only .pdb and .cif are supported.", file=sys.stderr)
            print("404")
            return

        if os.path.exists(pdb_path):
            try:
                # Normalize filename casing while preserving the actual format extension
                # (.pdb stays .pdb, .cif stays .cif) so the correct parser is used later
                base_name = os.path.splitext(actual_filename)[0]
                standard_path = os.path.join(output_dir, f"{base_name}{upload_ext}")

                # If the uploaded filename doesn't match the normalized path, copy it
                if pdb_path != standard_path:
                    import shutil
                    shutil.copy2(pdb_path, standard_path)
                    pdb_path = standard_path
                    print(f"[INFO] Copied to standard path: {pdb_path}", file=sys.stderr)
                
                chains = get_chain_ids(pdb_path)
                if chains:
                    print(f"[INFO] Extracted {len(chains)} chains and residues", file=sys.stderr)
                    for chain_id, residues in chains.items():
                        print(f"  Chain {chain_id}: {len(residues)} residues", file=sys.stderr)

                    # If this was a .cif upload, also write a companion .pdb file
                    # so downstream tools (MODELLER/PyMOL/LigPlot via xpro.py),
                    # which only understand legacy PDB format, can find it.
                    if upload_ext == ".cif":
                        companion_pdb_path = os.path.join(output_dir, f"{base_name}.pdb")
                        try:
                            convert_cif_to_pdb(pdb_path, companion_pdb_path)
                            print(f"[INFO] Converted CIF to companion PDB: {companion_pdb_path}", file=sys.stderr)
                        except Exception as e:
                            print(f"[ERROR] Failed to convert CIF to PDB for downstream analysis: {e}", file=sys.stderr)
                            print("404")
                            return

                    # Print JSON format for PHP to capture
                    json_output = json.dumps(chains)
                    print(json_output)
                else:
                    print("[ERROR] No chains/residues extracted from structure file", file=sys.stderr)
                    print("404")
            except Exception as e:
                print(f"[ERROR] Failed to parse uploaded structure file: {e}", file=sys.stderr)
                print("404")
        else:
            print(f"[ERROR] Uploaded structure file not found at: {pdb_path}", file=sys.stderr)
            print("404")
        
        return
    
    # For PDB ID downloads
    pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    pdb_id = pdb_id.upper().strip()
    
    print(f"[INFO] Downloading PDB: {pdb_id}", file=sys.stderr)
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"[INFO] Downloading from: {pdb_url}", file=sys.stderr)

    try:
        # Try with proxy first
        try:
            response = requests.get(pdb_url, proxies=proxies, timeout=30)
        except:
            # Fallback without proxy
            print(f"[INFO] Proxy failed, trying direct connection...", file=sys.stderr)
            response = requests.get(pdb_url, timeout=30)
        
        if response.status_code == 200:
            # Clean the file while writing
            with open(pdb_path, "w") as file:
                for line in response.text.splitlines():
                    if not (line.startswith("SEQRES") or line.startswith("ANISOU") or line.startswith("HETATM")):
                        if line.startswith("ATOM") or line.startswith("TER") or line.startswith("END"):
                            file.write(line + "\n")

            print(f"[INFO] ✅ Cleaned PDB file saved to: {pdb_path}", file=sys.stderr)

            # Verify that the file contains ATOM lines
            has_atoms = False
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith("ATOM"):
                        has_atoms = True
                        break
            
            if not has_atoms:
                print(f"[WARNING] ⚠️ No ATOM records found in {pdb_id}.pdb", file=sys.stderr)

            # Display chain and residue info
            chains = get_chain_ids(pdb_path)
            if chains:
                print(f"[INFO] Extracted {len(chains)} chains and residues", file=sys.stderr)
                for chain_id, residues in chains.items():
                    residue_list = list(residues.keys())
                    print(f"  Chain {chain_id}: {len(residues)} residues (first few: {residue_list[:5]})", file=sys.stderr)
                
                # Print JSON format for PHP to capture
                json_output = json.dumps(chains)
                print(json_output)
            else:
                print("[ERROR] Failed to extract chain/residue information", file=sys.stderr)
                print("404")
        else:
            print(f"[ERROR] ❌ HTTP {response.status_code}: Could not fetch {pdb_id}.pdb", file=sys.stderr)
            print("404")
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] ❌ Network error: {e}", file=sys.stderr)
        print("404")
    except Exception as e:
        print(f"[ERROR] ❌ Unexpected error: {e}", file=sys.stderr)
        print("404")

# ----- MAIN EXECUTION -----
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 download_pdbFile.py <PDB_ID or UPLOAD:filename> <base_dir> <random_number>")
        print("Example: python3 download_pdbFile.py 1A8O /var/www/html/fileUpload 1234")
        print("Example: python3 download_pdbFile.py UPLOAD:1a8o.pdb /var/www/html/fileUpload 1234")
        print("Example: python3 download_pdbFile.py UPLOAD:1a8o.cif /var/www/html/fileUpload 1234")
        sys.exit(1)

    pdb_id = sys.argv[1]
    base_dir = sys.argv[2]
    rand_num = sys.argv[3]

    download_pdb(pdb_id, base_dir, rand_num)
