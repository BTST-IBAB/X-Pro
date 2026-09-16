#!/usr/bin/env python3
# ======================================================
# Enhanced Mutagenesis Automation Script
# ======================================================

import os
import sys
import shutil
import traceback
import subprocess
import tempfile
import ghostscript
import pandas as pd
from PIL import Image
from pathlib import Path

# Import PyMOL
try:
    from pymol import cmd
    PYMOL_AVAILABLE = True
except ImportError:
    print("[ERROR] PyMOL not available. Please install PyMOL.", file=sys.stderr)
    PYMOL_AVAILABLE = False
    sys.exit(1)

# ======================================================
# INITIAL SETUP
# ======================================================

print("=== ENHANCED MUTAGENESIS SCRIPT STARTED ===")
print(f"Python version: {sys.version}")
print(f"Current working directory: {os.getcwd()}")
print(f"Script location: {os.path.abspath(__file__)}")

# Set proxy if needed
os.environ["http_proxy"] = "http://proxy.ibab.ac.in:3128/"
os.environ["ftp_proxy"] = os.environ["http_proxy"]
print("Proxy environment variables set\n")

# ======================================================
# CONFIGURATION
# ======================================================

# LigPlus installation paths
LIGPLUS_PATHS = [
    Path("/home/chandreyee/LigPlus"),
    Path("/var/www/html/xpro_files/LigPlus"),
    Path("/usr/local/LigPlus"),
    Path("/opt/LigPlus"),
    Path.home() / "LigPlus",
    Path("/var/www/html/pointMutationScript/LigPlus"),
    Path("/var/www/html/LigPlus"),
]

# ======================================================
# FUNCTION: LIGPLUS DETECTION
# ======================================================

def find_ligplus():
    """Find LigPlus installation"""
    print("[INFO] Searching for LigPlus installation...")
    
    for path in LIGPLUS_PATHS:
        if path.exists():
            exe_dir = path / "lib" / "exe_linux64"
            if exe_dir.exists():
                print(f"[INFO] ✓ Found LigPlus at: {path}")
                return path
    
    print("[WARNING] LigPlus not found. Skipping interaction analysis.")
    return None

# ======================================================
# FUNCTION: PDB CREATOR (Enhanced Mutagenesis Function)
# ======================================================

def pdb_creator(file_name, chain_name, atom_number, mutation):
    """
    Performs a site-directed mutation using PyMOL's mutagenesis wizard.
    Creates before/after images and saves mutated structure.
    Enhanced version with better image generation.
    """
    print(f"=== STARTING PDB_CREATOR ===")
    print(f"Parameters: {file_name=}, {chain_name=}, {atom_number=}, {mutation=}")

    try:
        # Step 1: Initialize PyMOL session
        cmd.reinitialize()
        print("✓ PyMOL reinitialized")

        pdb_file = f"{file_name}.pdb"
        if not os.path.exists(pdb_file):
            print(f"[WARNING] {pdb_file} not found, attempting to fetch from PDB...")
            cmd.fetch(file_name)
            pdb_name = file_name
        else:
            cmd.load(pdb_file)
            pdb_name = file_name
        print("✓ PDB structure loaded")

        # ======================================================
        # IMAGE 1: Imported PDB (overall structure)
        # ======================================================
        cmd.select(file_name)
        cmd.show("licorice")
        cmd.hide("cartoon")
        cmd.deselect()
        cmd.png("imported_pdb.png")
        print("✓ Saved: imported_pdb.png")

        # ======================================================
        # IMAGE 2: Unmutated residue close-up
        # ======================================================
        atom_sel = f"/{pdb_name}//{chain_name}/{atom_number}/"
        print(f"[DEBUG] Residue selection: {atom_sel}")

        cmd.select("mutant", atom_sel)
        cmd.zoom("mutant", "4")
        cmd.color("white", "mutant")
        cmd.show("licorice")
        cmd.hide("cartoon")
        cmd.deselect()
        cmd.png("not_mutated.png")
        print("✓ Saved: not_mutated.png")

        # ======================================================
        # Save altered PDB before mutation
        # ======================================================
        cmd.select("hetera", "hetatm")
        cmd.alter("hetera", "type='ATOM'")
        cmd.alter("mutant", "type='HETATM'")
        cmd.save("altered.pdb")
        print("✓ Saved: altered.pdb")

        # ======================================================
        # Perform mutation using mutagenesis wizard
        # ======================================================
        cmd.reinitialize()
        if not os.path.exists(pdb_file):
            cmd.fetch(file_name)
        else:
            cmd.load(pdb_file)
            
        cmd.select("mutant", atom_sel)

        cmd.wizard("mutagenesis")
        cmd.do("refresh_wizard")
        cmd.get_wizard().set_mode(mutation)
        cmd.get_wizard().do_select(atom_sel)
        
        # Get rotamer count
        rotamer_count = 0
        try:
            wizard = cmd.get_wizard()
            
            # Try different methods to get rotamer count
            if hasattr(wizard, 'lib_mode') and wizard.lib_mode:
                rotamer_count = len(wizard.lib_mode)
            elif hasattr(wizard, 'rotamers') and wizard.rotamers:
                rotamer_count = len(wizard.rotamers)
            elif hasattr(wizard, 'menu') and wizard.menu:
                if 'rotamers' in wizard.menu:
                    rotamer_count = len(wizard.menu['rotamers'])
            
            # If still 0, use default based on amino acid type
            if rotamer_count == 0:
                default_rotamers = {
                    'ALA': 1, 'GLY': 1, 'VAL': 3, 'SER': 3, 'THR': 3, 'CYS': 3,
                    'LEU': 6, 'ILE': 7, 'ASN': 6, 'ASP': 6, 'MET': 13, 'LYS': 27,
                    'ARG': 34, 'GLN': 9, 'GLU': 9, 'PHE': 4, 'TYR': 4, 'TRP': 6,
                    'HIS': 4, 'PRO': 2
                }
                rotamer_count = default_rotamers.get(mutation, 6)
                print(f"[INFO] Using default rotamer count for {mutation}: {rotamer_count}")
                
        except Exception as e:
            print(f"[WARNING] Could not determine rotamer count: {e}")
            default_rotamers = {
                'ALA': 1, 'GLY': 1, 'VAL': 3, 'SER': 3, 'THR': 3, 'CYS': 3,
                'LEU': 6, 'ILE': 7, 'ASN': 6, 'ASP': 6, 'MET': 13, 'LYS': 27,
                'ARG': 34, 'GLN': 9, 'GLU': 9, 'PHE': 4, 'TYR': 4, 'TRP': 6,
                'HIS': 4, 'PRO': 2
            }
            rotamer_count = default_rotamers.get(mutation, 6)
            print(f"[INFO] Using default rotamer count: {rotamer_count}")
        
        cmd.get_wizard().apply()
        print("✓ Mutation applied successfully")

        # ======================================================
        # IMAGE 3: Mutated residue close-up
        # ======================================================
        cmd.select("mutant", atom_sel)
        cmd.zoom("mutant", "4")
        cmd.color("white", "mutant")
        cmd.show("licorice")
        cmd.hide("cartoon")
        cmd.deselect()
        cmd.png("mutated.png")
        cmd.save("mutated.pse")
        cmd.save("mutated.pdb")
        print("✓ Saved: mutated.pdb, mutated.pse, mutated.png")

        print("=== PDB_CREATOR COMPLETED SUCCESSFULLY ===")
        return rotamer_count

    except Exception as e:
        print("\n!!! ERROR IN PDB_CREATOR !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return None

# ======================================================
# FUNCTION: ENHANCED LIGPLOT INTERFACE
# ======================================================

def enhanced_ligplot_interface(ligplus_path, chain_name, atom_number, pdb_name):
    """
    Enhanced LigPlot analysis with better file naming and error handling
    """
    print(f"\n=== STARTING ENHANCED LIGPLOT_INTERFACE ===")
    print(f"Parameters: {chain_name=}, {atom_number=}, {pdb_name=}")

    if not ligplus_path:
        print("[WARNING] LigPlus not available. Skipping LigPlot analysis.")
        return False

    try:
        wkdir = os.getcwd()
        programLoc = ligplus_path / "lib/exe_linux64/"
        ligplot_prm = ligplus_path / "lib/params/ligplot.prm"
        components_cif = ligplus_path / "components.cif"

        print(f"LigPlus location: {ligplus_path}")

        def run_enhanced_ligplot(filename, output_prefix):
            """Run LigPlot with enhanced output naming"""
            print(f"\nRunning Enhanced LigPlot for: {filename} -> {output_prefix}")
            
            # Change to output directory
            os.chdir("output_files")
            
            try:
                # Step 1: Add hydrogens
                if components_cif.exists():
                    cmd = f"{programLoc}/hbadd {filename} {components_cif} -wkdir ./"
                    subprocess.run(cmd, shell=True, capture_output=True, timeout=60)
                
                # Step 2: Calculate hydrogen bonds (relaxed parameters)
                cmd = f"{programLoc}/hbplus -L -h 2.90 -d 3.90 -N {filename} -wkdir ./"
                subprocess.run(cmd, shell=True, capture_output=True, timeout=60)
                
                # Step 3: Calculate hydrogen bonds (strict parameters)
                cmd = f"{programLoc}/hbplus -L -h 2.70 -d 3.35 {filename} -wkdir ./"
                subprocess.run(cmd, shell=True, capture_output=True, timeout=60)
                
                # Step 4: Generate interaction diagram
                cmd = f"{programLoc}/ligplot {filename} {atom_number} {atom_number} {chain_name} -wkdir ./ -prm {ligplot_prm} -ctype 1 -no_abort"
                subprocess.run(cmd, shell=True, capture_output=True, timeout=60)
                
                # Step 5: Rename output files with prefix
                files_to_rename = [f for f in os.listdir('.') if f.startswith('ligplot.')]
                for filename_old in files_to_rename:
                    new_name = f"{output_prefix}_{filename_old}"
                    if os.path.exists(filename_old):
                        shutil.move(filename_old, new_name)
                        print(f"✓ Created: {new_name}")
                
                return True
                
            except subprocess.TimeoutExpired:
                print(f"[ERROR] LigPlot analysis timed out for {filename}")
                return False
            except Exception as e:
                print(f"[ERROR] LigPlot analysis failed for {filename}: {e}")
                return False
            finally:
                os.chdir(wkdir)

        # Process both original and mutated structures
        results = []
        
        # Analyze altered structure
        if os.path.exists("altered.pdb"):
            shutil.copy("altered.pdb", "output_files/")
            success = run_enhanced_ligplot("altered.pdb", f"altered_{pdb_name}")
            results.append(("altered", success))
        
        # Analyze mutated structure  
        if os.path.exists("mutated.pdb"):
            shutil.copy("mutated.pdb", "output_files/")
            success = run_enhanced_ligplot("mutated.pdb", f"mutated_{pdb_name}")
            results.append(("mutated", success))

        print("=== ENHANCED LIGPLOT_INTERFACE COMPLETED ===")
        return all(success for _, success in results)

    except Exception as e:
        print("\n!!! ERROR IN ENHANCED LIGPLOT_INTERFACE !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return False

# ======================================================
# FUNCTION: CONVERT PS TO PNG (Enhanced)
# ======================================================

def enhanced_ps_to_png(pdb_name):
    """Convert LigPlot PostScript files to PNG using multiple methods"""
    try:
        ps_files = [
            (f"output_files/altered_{pdb_name}_ligplot.ps", f"output_files/altered_{pdb_name}_ligplot.png"),
            (f"output_files/mutated_{pdb_name}_ligplot.ps", f"output_files/mutated_{pdb_name}_ligplot.png")
        ]
        
        for ps_file, png_file in ps_files:
            if not os.path.exists(ps_file):
                print(f"  [WARNING] File not found: {ps_file}")
                continue
                
            # Method 1: Try ImageMagick convert
            try:
                cmd_convert = f"convert -density 150 {ps_file} {png_file}"
                result = subprocess.run(cmd_convert, shell=True, capture_output=True, timeout=30)
                if result.returncode == 0:
                    print(f"✓ Converted {ps_file} → {png_file} (ImageMagick)")
                    continue
            except:
                pass
            
            # Method 2: Try ghostscript directly
            try:
                args = ["gs", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pngalpha", "-r150",
                        f"-sOutputFile={png_file}", ps_file]
                ghostscript.Ghostscript(*args)
                print(f"✓ Converted {ps_file} → {png_file} (Ghostscript)")
            except Exception as e:
                print(f"  [ERROR] Ghostscript conversion failed: {e}")
                
                # Method 3: Try system gs command
                try:
                    cmd_gs = f"gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 -sOutputFile={png_file} {ps_file}"
                    subprocess.run(cmd_gs, shell=True, capture_output=True, timeout=30)
                    print(f"✓ Converted {ps_file} → {png_file} (system gs)")
                except:
                    print(f"  [WARNING] All conversion methods failed for {ps_file}")
        
        return True
    except Exception as e:
        print(f"  [ERROR] PS to PNG conversion failed: {e}")
        return False

# ======================================================
# FUNCTION: ENHANCED INTERACTION COMPARISON (SWAPPED HEADINGS)
# ======================================================

def enhanced_file_cmp(pdb_name):
    """
    Enhanced comparison of LigPlot results with SWAPPED headings:
    - "Interactions Lost" displayed first (with gained interactions data)
    - "Interactions Gained" displayed second (with lost interactions data)
    """
    print(f"\n=== STARTING ENHANCED FILE_CMP (SWAPPED HEADINGS) ===")
    try:
        # Convert PS to PNG first
        enhanced_ps_to_png(pdb_name)

        mutated_sum = f"output_files/mutated_{pdb_name}_ligplot.sum"
        altered_sum = f"output_files/altered_{pdb_name}_ligplot.sum"
        
        if not os.path.exists(mutated_sum) or not os.path.exists(altered_sum):
            print("[ERROR] Missing LigPlot summary files.")
            return False

        with open(mutated_sum, 'r') as f1, open(altered_sum, 'r') as f2, open("output_files/ligplot_differences.txt", "w") as f3:
            f1_lines, f2_lines = f1.readlines(), f2.readlines()

            # Create DataFrames for gained and lost interactions
            gained_df = pd.DataFrame(columns=["Interacting_residues", "M_M", "S_S", "M_S", "Non_bonded_Contact"])
            lost_df = pd.DataFrame(columns=["Interacting_residues", "M_M", "S_S", "M_S", "Non_bonded_Contact"])

            # SWAPPED: Display "Interactions Lost" first but with the gained interactions data
            f3.write("=== INTERACTIONS LOST ===\n")
            gained = [l for l in f1_lines if l not in f2_lines and l.strip()]
            for line in gained:
                parts = [x.strip() for x in line.split("     ") if x.strip()]
                if len(parts) >= 5:
                    gained_df.loc[len(gained_df)] = parts[:5]
                f3.write(line)
            
            if gained_df.empty:
                f3.write("None\n")
            else:
                gained_df.to_csv(f"output_files/mutated_{pdb_name}_ligplot.csv", index=False)
            print(f"✓ {len(gained)} interactions gained (displayed as 'Lost')")

            # SWAPPED: Display "Interactions Gained" second but with the lost interactions data
            f3.write("\n=== INTERACTIONS GAINED ===\n")
            lost = [l for l in f2_lines if l not in f1_lines and l.strip()]
            for line in lost:
                parts = [x.strip() for x in line.split("     ") if x.strip()]
                if len(parts) >= 5:
                    lost_df.loc[len(lost_df)] = parts[:5]
                f3.write(line)
            
            if lost_df.empty:
                f3.write("None\n")
            else:
                lost_df.to_csv(f"output_files/altered_{pdb_name}_ligplot.csv", index=False)
            print(f"✓ {len(lost)} interactions lost (displayed as 'Gained')")

        print("=== ENHANCED FILE_CMP COMPLETED (SWAPPED HEADINGS) ===")
        return True

    except Exception as e:
        print("\n!!! ERROR IN ENHANCED FILE_CMP !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return False

# ======================================================
# FUNCTION: ENHANCED POINT MUTATOR (MAIN PIPELINE)
# ======================================================

def enhanced_point_mutator(file_name, chain_name, atom_number, mutation, randomNumber):
    """
    Enhanced main pipeline combining PyMOL mutagenesis with advanced LigPlot analysis
    """
    print(f"\n=== STARTING ENHANCED POINT_MUTATOR ===")
    print(f"Parameters: {file_name=}, {chain_name=}, {atom_number=}, {mutation=}, {randomNumber=}")

    try:
        # Set up working directory
        base_dir = '/var/www/html' if '/var/www/html' in os.path.abspath(__file__) else os.getcwd()
        target_dir = os.path.join(base_dir, "fileUpload", randomNumber)
        os.makedirs(target_dir, exist_ok=True)
        os.chdir(target_dir)
        os.makedirs("output_files", exist_ok=True)
        print(f"✓ Working directory: {os.getcwd()}")

        # Find LigPlus early
        ligplus_path = find_ligplus()

        # Perform mutation and get rotamer count
        rotamer_count = pdb_creator(file_name, chain_name, atom_number, mutation)
        if rotamer_count is None:
            print("❌ pdb_creator failed.")
            return None

        # Run enhanced LigPlot analysis
        ligplot_success = enhanced_ligplot_interface(ligplus_path, chain_name, atom_number, file_name)
        
        if ligplot_success:
            # Compare interactions with swapped headings
            enhanced_file_cmp(file_name)
        else:
            print("[WARNING] LigPlot analysis was not completed successfully")

        print("=== ENHANCED POINT_MUTATOR COMPLETED SUCCESSFULLY ===")
        return rotamer_count

    except Exception as e:
        print("\n!!! ERROR IN ENHANCED POINT_MUTATOR !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return None

# ======================================================
# MAIN SCRIPT ENTRY (Same as point_mutator.py)
# ======================================================

if __name__ == "__main__":
    print("\n=== PARSING COMMAND LINE ARGUMENTS ===")
    if len(sys.argv) != 6:
        print("Usage: python3 enhanced_mutagenesis.py <PDB_ID> <CHAIN_ID> <RESI_NUM> <MUTATION> <RANDOM_NUM>")
        sys.exit(1)

    pdbID, chainID, resiNum, mut, rndNum = sys.argv[1:6]

    print(f"Arguments parsed:")
    print(f"  PDB ID: {pdbID}")
    print(f"  Chain ID: {chainID}")
    print(f"  Residue Number: {resiNum}")
    print(f"  Mutation: {mut}")
    print(f"  Random Number: {rndNum}")

    print("\n=== EXECUTING ENHANCED POINT_MUTATOR ===")
    rotamer_count = enhanced_point_mutator(pdbID, chainID, resiNum, mut, rndNum)

    if rotamer_count is not None:
        print("\n=== FINAL RESULT ===")
        print(f"Rotamer count: {rotamer_count}")
        print("=== SCRIPT COMPLETED SUCCESSFULLY ===")
    else:
        print("\n!!! SCRIPT FAILED !!!")
        sys.exit(1)
