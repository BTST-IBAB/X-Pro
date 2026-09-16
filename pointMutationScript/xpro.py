#!/usr/bin/env python3
# ======================================================
# Integrated Mutagenesis Automation Script
# Using direct imports with os.system as requested
# ======================================================

import os
import sys
import shutil
import traceback
import subprocess
import time
import pandas as pd
import json
import urllib.request
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFIO
import interaction_comparison

def write_progress_step(target_dir, step_file):
    """
    Write a progress-step marker file so X-Pro-Progress.php can show the
    person a live, step-by-step status while this script runs in the
    background (it's launched detached, not something PHP waits on
    synchronously). Uses the absolute target_dir rather than relying on the
    current working directory, since this may be called from inside a
    try/except where an earlier os.chdir() could theoretically not have
    taken effect. Failure to write is non-fatal -- the progress page just
    keeps showing the previous step a little longer.
    """
    try:
        with open(os.path.join(target_dir, step_file), "w") as fh:
            fh.write(str(time.time()))
    except Exception as e:
        print(f"[WARNING] Could not write progress marker {step_file}: {e}")


# ======================================================
# DEBUG: Check write permissions
# ======================================================
print("=== CHECKING PERMISSIONS ===")
print(f"Current user: {os.getuid()}")
print(f"Current directory: {os.getcwd()}")
print(f"Directory writable: {os.access(os.getcwd(), os.W_OK)}")
try:
    test_file = "test_write.txt"
    with open(test_file, 'w') as f:
        f.write("test")
    print(f"✓ Can write files in current directory")
    os.remove(test_file)
except Exception as e:
    print(f"✗ Cannot write files: {e}")

# ======================================================
# DIRECT IMPORTS - Using os.system approach as requested
# ======================================================
print("=== IMPORTING MODULES ===")

# Import modeller directly
try:
    import modeller
    from modeller import *
    from modeller.optimizers import MolecularDynamics, ConjugateGradients
    from modeller.automodel import autosched
    MODELLER_AVAILABLE = True
    print("✓ Modeller imported successfully")
except ImportError as e:
    print(f"✗ Modeller import failed: {e}")
    MODELLER_AVAILABLE = False

# Import pymol directly
try:
    # Set headless display
    if "DISPLAY" not in os.environ:
        os.environ["DISPLAY"] = ":99"
    os.environ.setdefault("PYMOL_HIDE_WELCOME", "1")
    
    import pymol
    from pymol import cmd
    # Launch PyMOL in quiet mode
    pymol.finish_launching(['pymol', '-qc'])
    PYMOL_AVAILABLE = True
    print("✓ PyMOL imported successfully")
except ImportError as e:
    print(f"✗ PyMOL import failed: {e}")
    PYMOL_AVAILABLE = False

print(f"PyMOL available: {PYMOL_AVAILABLE}")
print(f"Modeller available: {MODELLER_AVAILABLE}")

# ======================================================
# INITIAL SETUP
# ======================================================

print("=== INTEGRATED MUTAGENESIS SCRIPT STARTED ===")
print(f"Python version: {sys.version}")
print(f"Current working directory: {os.getcwd()}")
print(f"Script location: {os.path.abspath(__file__)}")

# Set proxy
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
# FUNCTION: DOWNLOAD PDB IF NOT EXISTS
# ======================================================

def ensure_pdb_file(file_name):
    """Download PDB file from RCSB if it doesn't exist"""
    pdb_file = f"{file_name}.pdb"
    
    if os.path.exists(pdb_file):
        print(f"✓ PDB file already exists: {pdb_file}")
        return True
    
    print(f"[INFO] PDB file {pdb_file} not found. Attempting to download...")
    
    # Try multiple URLs for downloading
    urls = [
        f"https://files.rcsb.org/download/{file_name.upper()}.pdb",
        f"https://www.rcsb.org/pdb/files/{file_name.upper()}.pdb",
        f"https://files.rcsb.org/view/{file_name.upper()}.pdb"
    ]
    
    for url in urls:
        try:
            print(f"  Trying: {url}")
            urllib.request.urlretrieve(url, pdb_file)
            if os.path.exists(pdb_file) and os.path.getsize(pdb_file) > 100:
                print(f"✓ Successfully downloaded from: {url}")
                return True
        except Exception as e:
            print(f"  Failed: {e}")
            continue
    
    # If all URLs fail, try using PyMOL if available
    if PYMOL_AVAILABLE:
        try:
            print("  Attempting download via PyMOL...")
            cmd.reinitialize()
            cmd.fetch(file_name)
            cmd.save(pdb_file)
            if os.path.exists(pdb_file):
                print(f"✓ Successfully downloaded via PyMOL")
                return True
        except Exception as e:
            print(f"  PyMOL download failed: {e}")
    
    print(f"[ERROR] Could not download PDB file for {file_name}")
    return False

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
# FUNCTION: PYMOL VISUALIZATION
# ======================================================

def pymol_visualization(file_name, chain_name, atom_number, mutation, method_used):
    """
    Generate PyMOL visualization images for both original and mutated structures.
    """
    print(f"=== STARTING PYMOL VISUALIZATION ===")
    print(f"Parameters: {file_name=}, {chain_name=}, {atom_number=}, {mutation=}, {method_used=}")

    if not PYMOL_AVAILABLE:
        print("[ERROR] PyMOL is not available. Cannot generate visualization images.")
        return False

    try:
        # Ensure output_files directory exists
        os.makedirs("output_files", exist_ok=True)

        # Ensure PDB file exists
        if not ensure_pdb_file(file_name):
            print("[ERROR] Cannot proceed with visualization without PDB file.")
            return False

        # ----------------------------------------------------------------
        # IMAGE SET 1: Original / Before-Mutation structure
        # ----------------------------------------------------------------
        cmd.reinitialize()
        print("✓ PyMOL reinitialized")

        pdb_file = f"{file_name}.pdb"
        cmd.load(pdb_file, file_name)
        print("✓ Original PDB structure loaded")

        # White background, cartoon + highlight target residue
        cmd.bg_color("white")
        cmd.hide("everything", "all")
        cmd.show("cartoon", "all")
        cmd.color("lightblue", "all")

        atom_sel = f"/{file_name}//{chain_name}/{atom_number}/"
        cmd.select("target_res", atom_sel)
        cmd.show("sticks", "target_res")
        cmd.color("red", "target_res")
        cmd.deselect()

        # Save full-structure overview
        overview_orig = os.path.abspath(f"output_files/original_{file_name}_overview.png")
        cmd.zoom("all", 5)
        cmd.png(overview_orig, width=800, height=600, ray=0, quiet=1)
        print(f"✓ Saved original overview: {overview_orig}")

        # Save residue close-up
        closeup_orig = os.path.abspath(f"output_files/altered_{file_name}_ligplot.png")
        cmd.zoom("target_res", 6)
        cmd.png(closeup_orig, width=800, height=600, ray=0, quiet=1)
        print(f"✓ Saved before-mutation close-up: {closeup_orig}")

        # ----------------------------------------------------------------
        # IMAGE SET 2: Mutated / After-Mutation structure
        # ----------------------------------------------------------------
        cmd.reinitialize()
        mutated_file = "mutated.pdb"
        if not os.path.exists(mutated_file):
            print(f"[ERROR] Mutated PDB file '{mutated_file}' not found for visualization.")
            return False

        cmd.load(mutated_file, "mutated")
        print("✓ Mutated PDB structure loaded")

        cmd.bg_color("white")
        cmd.hide("everything", "all")
        cmd.show("cartoon", "all")
        cmd.color("palegreen", "all")

        atom_sel_mut = f"/mutated//{chain_name}/{atom_number}/"
        cmd.select("mutant_res", atom_sel_mut)
        cmd.show("sticks", "mutant_res")
        cmd.color("orange", "mutant_res")
        cmd.deselect()

        # Save full-structure overview
        overview_mut = os.path.abspath(f"output_files/mutated_{file_name}_overview.png")
        cmd.zoom("all", 5)
        cmd.png(overview_mut, width=800, height=600, ray=0, quiet=1)
        print(f"✓ Saved mutated overview: {overview_mut}")

        # Save residue close-up
        closeup_mut = os.path.abspath(f"output_files/mutated_{file_name}_ligplot.png")
        cmd.zoom("mutant_res", 6)
        cmd.png(closeup_mut, width=800, height=600, ray=0, quiet=1)
        print(f"✓ Saved after-mutation close-up: {closeup_mut}")

        # Verify images were actually written
        for img_path in [closeup_orig, closeup_mut]:
            if os.path.exists(img_path) and os.path.getsize(img_path) > 0:
                print(f"✓ Image verified: {img_path} ({os.path.getsize(img_path)} bytes)")
            else:
                print(f"[WARNING] Image missing or empty: {img_path}")

        print("=== PYMOL VISUALIZATION COMPLETED SUCCESSFULLY ===")
        return True

    except Exception as e:
        print("\n!!! ERROR IN PYMOL VISUALIZATION !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return False

# ======================================================
# FUNCTION: GENERATE 3D HTML - From xpro_prev(2).py
# ======================================================

def generate_3d_html(file_name, chain_name, atom_number, mutation, randomNumber, method=None):
    """
    Generate interactive 3D HTML files with enhanced roving details
    """
    print(f"=== GENERATING 3D HTML FROM PYMOL ===")
    print(f"Current directory: {os.getcwd()}")
    print(f"Output directory: output_files/")

    # Viewer styling constants (kept local to this function, mirroring the
    # module-level ones used for the same side-by-side viewer elsewhere).
    VIEWER_BACKGROUND_COLOR = "#d9e1e4"
    CARTOON_COLOR = "spectrum"
    FULL_CARTOON_OPACITY = 1.0
    HIGHLIGHT_BACKGROUND_OPACITY = 0.45  # 55% transparent in highlight mode only.
    LOCAL_NEIGHBOR_RADIUS = 6.0

    try:
        output_dir = "output_files"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"✓ Created output directory: {output_dir}")
        
        os.chmod(output_dir, 0o777)
        print(f"✓ Set permissions on output directory")
        
        # Check if PDB files exist
        original_pdb = f"{file_name}.pdb"
        mutated_pdb = "mutated.pdb"
        
        print(f"Looking for original PDB: {original_pdb}")
        if not os.path.exists(original_pdb):
            pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]
            print(f"Available PDB files: {pdb_files}")
            if pdb_files:
                original_pdb = pdb_files[0]
                print(f"Using alternative original PDB: {original_pdb}")
            else:
                print("[ERROR] No original PDB file found")
                return False
        
        print(f"Looking for mutated PDB: {mutated_pdb}")
        if not os.path.exists(mutated_pdb):
            mutated_files = [f for f in os.listdir('.') if 'mutated' in f.lower() and f.endswith('.pdb')]
            print(f"Available mutated files: {mutated_files}")
            if mutated_files:
                mutated_pdb = mutated_files[0]
                print(f"Using alternative mutated PDB: {mutated_pdb}")
            else:
                print("[ERROR] No mutated PDB file found")
                return False
        
        # Read PDB files
        print(f"Reading original PDB: {original_pdb}")
        with open(original_pdb, 'r') as f:
            original_pdb_content = f.read()
        print(f"✓ Read {len(original_pdb_content)} bytes from original PDB")
        
        print(f"Reading mutated PDB: {mutated_pdb}")
        with open(mutated_pdb, 'r') as f:
            mutated_pdb_content = f.read()
        print(f"✓ Read {len(mutated_pdb_content)} bytes from mutated PDB")
        
        # Panel subtitle text: xpro.py doesn't build a model_provenance.json
        # (that belongs to the separate provenance-tracking pipeline), so the
        # WT/mutant descriptor is derived directly from the modelling method
        # already known to the caller instead.
        wt_provenance_label = "Experimental/available PDB structure"
        _method_key = str(method or "").strip().lower()
        if _method_key == "modeller":
            mutant_provenance_label = "Modeller model"
        elif _method_key == "pymol":
            mutant_provenance_label = "PyMOL model"
        else:
            mutant_provenance_label = "Mutant model"
        comparison_files = sorted(
            Path(output_dir).glob("interaction_comparison_*.json")
        )
        canonical_comparison = {}
        if comparison_files:
            try:
                canonical_comparison = json.loads(
                    comparison_files[0].read_text(encoding="utf-8")
                )
            except (OSError, ValueError):
                canonical_comparison = {}

        def canonical_contacts_for(state):
            contacts = []
            contact_key = "mutant_contacts" if state == "mutant" else "wt_contacts"
            for record in canonical_comparison.get("records", []):
                for contact in record.get(contact_key, []):
                    contacts.append({
                        "record_id": record.get("id"),
                        "result": record.get("result"),
                        "interaction_class": record.get("interaction_class"),
                        "distance": contact.get("distance"),
                        "residue_1": contact.get("residue_1"),
                        "atom_1": contact.get("atom_1"),
                        "residue_2": contact.get("residue_2"),
                        "atom_2": contact.get("atom_2"),
                    })
            return contacts

        def residue_name_at_target(pdb_content):
            for line in pdb_content.splitlines():
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                if line[21:22].strip() != str(chain_name).strip():
                    continue
                if line[22:26].strip() != str(atom_number).strip():
                    continue
                return line[17:20].strip().upper()
            return "???"

        # Both labels come from the exact coordinate files embedded in the two
        # viewers; no residue name is assumed from the request.
        wildtype_3let = residue_name_at_target(original_pdb_content)
        mutant_3let = residue_name_at_target(mutated_pdb_content)

        orig_mutation_label = f"{wildtype_3let}{atom_number}"
        # Function to create the synchronized molecular-viewer panels.
        def create_enhanced_html(pdb_content, title, is_mutated=False):
            """Generate a full-structure view plus canonical local highlighting."""
            
            pdb_content_escaped = json.dumps(pdb_content)
            
            res_num = atom_number
            chain = chain_name
            highlight_color = "red" if is_mutated else "orange"
            site_label = f"{mutant_3let}{atom_number}" if is_mutated else orig_mutation_label
            panel_title = "Mutant" if is_mutated else "Original / WT"
            panel_subtitle = (
                f"Chain {chain or '(blank)'} · {mutant_3let}{atom_number}"
                if is_mutated else f"Chain {chain or '(blank)'} · {wildtype_3let}{atom_number}"
            )
            pane_id = "mutant" if is_mutated else "wt"
            provenance_label = mutant_provenance_label if is_mutated else wt_provenance_label
            canonical_contacts = canonical_contacts_for(
                "mutant" if is_mutated else "wt"
            )
            canonical_contacts_escaped = json.dumps(canonical_contacts)
            
            html = f'''<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <meta charset="utf-8">
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        }}
        
        body, html {{
            width: 100%;
            height: 100%;
            overflow: hidden;
            background-color: {VIEWER_BACKGROUND_COLOR};
        }}
        
        /* Dedicated layout rows keep headers, controls and legends off the molecule. */
        .container {{
            width: 100%;
            height: 100%;
            min-height: 0;
            display: grid;
            grid-template-rows: auto minmax(0, 1fr) auto;
            background: {VIEWER_BACKGROUND_COLOR};
        }}

        .model-header {{
            min-height: 58px;
            padding: 9px 14px;
            background: #f7fafb;
            color: #173e4b;
            border-bottom: 1px solid #cad6db;
            display: flex;
            flex-direction: column;
            justify-content: center;
            position: relative;
            z-index: 2;
        }}

        .model-header strong {{
            font-size: 14px;
            line-height: 1.2;
        }}

        .model-header span {{
            margin-top: 3px;
            font-size: 12px;
            color: #52646d;
            font-weight: 600;
        }}

        .model-header small {{
            margin-top: 2px;
            color: #6b7f88;
            font-size: 10px;
        }}

        .viewer-stage {{
            position: relative;
            min-height: 0;
            background: {VIEWER_BACKGROUND_COLOR};
        }}

        #viewer {{
            width: 100%;
            height: 100%;
            position: absolute;
            inset: 0;
        }}

        .viewer-footer {{
            position: relative;
            z-index: 2;
            padding: 8px 10px;
            background: #ffffff;
            color: #314b55;
            border-top: 1px solid #cad6db;
        }}

        .control-bar {{
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            justify-content: center;
            gap: 7px;
        }}

        .btn {{
            background: #ffffff;
            color: #175161;
            border: 1px solid #8da7b0;
            padding: 7px 12px;
            border-radius: 5px;
            font-size: 11px;
            font-weight: 700;
            cursor: pointer;
            transition: background 0.2s ease, border-color 0.2s ease;
            white-space: nowrap;
        }}

        .btn:hover {{
            background: #eef5f7;
            border-color: #4f7d8c;
        }}

        .btn.active {{
            background: #175f78;
            color: white;
            border-color: #4babc9;
        }}

        .btn.highlight.active {{
            background: #a54b18;
            border-color: #e67e22;
        }}

        .viewer-help {{
            display: flex;
            flex-wrap: wrap;
            justify-content: center;
            gap: 6px 14px;
            margin-top: 7px;
            color: #667b83;
            font-size: 10px;
        }}

        .legend {{
            display: none;
            flex-wrap: wrap;
            justify-content: center;
            gap: 5px 12px;
            margin: 7px auto 0;
            padding: 7px 9px;
            max-width: 460px;
            background: #f6f9fa;
            border: 1px solid #ccd8dc;
            border-radius: 5px;
            font-size: 9px;
        }}

        .legend.active {{ display: flex; }}

        .legend-item {{
            display: flex;
            align-items: center;
            min-width: 0;
        }}

        .legend-dash {{
            width: 20px;
            height: 2px;
            background: repeating-linear-gradient(90deg, #52646d, #52646d 5px, transparent 5px, transparent 9px);
            margin-right: 6px;
            flex: 0 0 auto;
        }}

        .legend-line {{
            width: 20px;
            height: 2px;
            background: #27845b;
            margin-right: 6px;
            flex: 0 0 auto;
        }}

        .legend-cross {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
            background: #d63434;
            margin-right: 9px;
            flex: 0 0 auto;
        }}

        .legend-dot {{
            width: 20px;
            height: 4px;
            margin-right: 6px;
            flex: 0 0 auto;
            background: radial-gradient(circle, #52646d 1.5px, transparent 2px) 0 0 / 6px 4px repeat-x;
        }}

        .legend-site {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
            margin-right: 9px;
            flex: 0 0 auto;
            background: {highlight_color};
            box-shadow: 0 0 0 3px rgba(255,255,255,.75);
        }}

        .legend-sphere {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
            margin: 0 10px 0 3px;
            flex: 0 0 auto;
        }}

        .sphere-blue {{
            background: #3366ff;
            box-shadow: 0 0 5px #3366ff;
        }}

        @media (max-width: 520px) {{
            .btn {{ padding: 6px 9px; font-size: 10px; }}
        }}

        /* Loading indicator */
        .loading {{
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            color: #666;
            font-size: 14px;
            z-index: 1;
        }}
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.8.0/3Dmol-min.js"></script>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
</head>
<body>
    <div class="container">
        <header class="model-header">
            <strong>{panel_title}</strong>
            <span>{panel_subtitle}</span>
            <small>{provenance_label}</small>
        </header>

        <div class="viewer-stage">
            <div id="viewer"></div>
        </div>

        <footer class="viewer-footer">
            <div class="control-bar">
                <button class="btn" type="button" onclick="toggleCartoon()" id="btnCartoon">Cartoon</button>
                <button class="btn" type="button" onclick="toggleStick()" id="btnStick">Stick</button>
                <button class="btn highlight" type="button" onclick="toggleMutationHighlight()" id="btnHighlight">Highlight mutation site</button>
            </div>
            <div class="viewer-help">
                <span>Rotate: drag</span><span>Zoom: scroll</span><span>Pan: right drag</span>
            </div>
            <div class="legend" id="legend">
                <div class="legend-item"><div class="legend-line"></div><span>Local atoms/sticks</span></div>
                <div class="legend-item"><div class="legend-cross"></div><span>Oxygen atoms</span></div>
                <div class="legend-item"><div class="legend-sphere sphere-blue"></div><span>Nitrogen atoms</span></div>
                <div class="legend-item"><div class="legend-dash"></div><span>HBPLUS hydrogen bonds</span></div>
                <div class="legend-item"><div class="legend-dot"></div><span>LigPlot non-bonded contacts</span></div>
                <div class="legend-item"><div class="legend-site"></div><span>Mutation site</span></div>
            </div>
        </footer>
    </div>

    <script>
        let viewer;
        let cartoonOn = true;
        let stickOn = false;
        let highlightMode = false;
        let viewMode = "whole";
        
        // Mutation site info
        const mutationInfo = {{
            resi: "{res_num}",
            chain: "{chain}",
            color: "{highlight_color}",
            label: "{site_label}",
            pane: "{pane_id}"
        }};
        const canonicalContacts = {canonical_contacts_escaped};
        let applyingSyncedView = false;
        
        // Button references
        const btnCartoon = document.getElementById('btnCartoon');
        const btnStick = document.getElementById('btnStick');
        const btnHighlight = document.getElementById('btnHighlight');
        const legend = document.getElementById('legend');
        
        $(document).ready(function() {{
            try {{
                viewer = $3Dmol.createViewer($("#viewer"), {{
                    backgroundColor: "{VIEWER_BACKGROUND_COLOR}",
                    antialias: true,
                    fog: false
                }});
                
                // Add PDB data
                let pdbData = {pdb_content_escaped};
                viewer.addModel(pdbData, "pdb");
                
                // Apply initial style
                applyNormalStyle();
                
                // Zoom to structure
                viewer.zoomTo({{}});
                viewer.render();
                
                // Set initial button states
                updateButtonStates();
                
                if (typeof viewer.setViewChangeCallback === "function") {{
                    viewer.setViewChangeCallback(function() {{
                        if (applyingSyncedView || window.parent === window) return;
                        window.parent.postMessage({{
                            type: "xpro-view-change",
                            pane: mutationInfo.pane,
                            view: viewer.getView()
                        }}, window.location.origin);
                    }});
                }}
                if (window.parent !== window) {{
                    window.parent.postMessage({{
                        type: "xpro-viewer-ready",
                        pane: mutationInfo.pane,
                        view: viewer.getView()
                    }}, window.location.origin);
                }}
                console.log("Viewer initialized");
                
            }} catch (e) {{
                console.error("Error:", e);
                $("#viewer").html('<div style="color: red; padding: 20px; text-align: center;">Error loading viewer: ' + e.message + '</div>');
            }}
        }});
        
        window.addEventListener("message", function(event) {{
            if (event.origin !== window.location.origin || !viewer) return;
            const message = event.data || {{}};
            if (message.type === "xpro-apply-view" && Array.isArray(message.view) &&
                    message.pane !== mutationInfo.pane) {{
                applyingSyncedView = true;
                viewer.setView(message.view);
                viewer.render();
                window.setTimeout(function() {{ applyingSyncedView = false; }}, 0);
                return;
            }}
            if (message.type === "xpro-viewer-command") {{
                if (message.command === "whole") {{
                    viewMode = "whole";
                    setHighlightMode(false);
                    viewer.zoomTo({{}});
                    viewer.render();
                }} else if (message.command === "focus") {{
                    viewMode = "focus";
                    setHighlightMode(false);
                    viewer.zoomTo(canonicalEnvironmentSelection());
                    viewer.zoom(0.78);
                    viewer.render();
                }} else if (message.command === "highlight") {{
                    setHighlightMode(true);
                }} else if (message.command === "report-view") {{
                    window.parent.postMessage({{
                        type: "xpro-view-change",
                        pane: mutationInfo.pane,
                        view: viewer.getView()
                    }}, window.location.origin);
                }}
            }}
        }});
        
        function residueSelection(residue) {{
            if (!residue) return null;
            return {{
                chain: residue.chain || "",
                resi: String(residue.residue_number) + String(residue.insertion_code || "")
            }};
        }}

        function residueKey(residue) {{
            const selection = residueSelection(residue);
            return selection ? selection.chain + "|" + selection.resi : "";
        }}

        function isMutationResidue(residue) {{
            const selection = residueSelection(residue);
            return Boolean(selection) &&
                String(selection.chain).trim() === String(mutationInfo.chain).trim() &&
                String(selection.resi) === String(mutationInfo.resi);
        }}

        function collectCanonicalPartners() {{
            const partners = new Map();
            canonicalContacts.forEach(function(contact) {{
                [contact.residue_1, contact.residue_2].forEach(function(residue) {{
                    if (!residue || isMutationResidue(residue)) return;
                    const key = residueKey(residue);
                    if (key && !partners.has(key)) partners.set(key, residue);
                }});
            }});
            return Array.from(partners.values());
        }}

        const canonicalPartners = collectCanonicalPartners();
        const changedPartnerKeys = new Set();
        canonicalContacts.forEach(function(contact) {{
            if (contact.result === "Retained") return;
            [contact.residue_1, contact.residue_2].forEach(function(residue) {{
                if (residue && !isMutationResidue(residue)) {{
                    changedPartnerKeys.add(residueKey(residue));
                }}
            }});
        }});

        function canonicalEnvironmentSelection() {{
            const selections = [{{resi: mutationInfo.resi, chain: mutationInfo.chain}}];
            canonicalPartners.forEach(function(residue) {{
                const selection = residueSelection(residue);
                if (selection) selections.push(selection);
            }});
            if (selections.length > 1) return {{or: selections}};
            return {{byres: true, within: {{
                distance: {LOCAL_NEIGHBOR_RADIUS}, sel: selections[0]
            }}}};
        }}

        function selectedResidueAtoms(residue) {{
            const selection = residueSelection(residue);
            return selection ? viewer.selectedAtoms(selection) : [];
        }}

        function averagePosition(atoms) {{
            if (!atoms.length) return null;
            return {{
                x: atoms.reduce((sum, atom) => sum + atom.x, 0) / atoms.length,
                y: atoms.reduce((sum, atom) => sum + atom.y, 0) / atoms.length,
                z: atoms.reduce((sum, atom) => sum + atom.z, 0) / atoms.length
            }};
        }}

        function mutationAtoms() {{
            return viewer.selectedAtoms({{
                resi: mutationInfo.resi,
                chain: mutationInfo.chain
            }});
        }}

        function mutationCenter() {{
            return averagePosition(mutationAtoms());
        }}

        function selectCanonicalAtom(residue, atomName) {{
            const selection = residueSelection(residue);
            if (!selection || !atomName) return null;
            selection.atom = atomName;
            const atoms = viewer.selectedAtoms(selection);
            return atoms.length ? atoms[0] : null;
        }}

        function pointBetween(start, end, fraction) {{
            return {{
                x: start.x + (end.x - start.x) * fraction,
                y: start.y + (end.y - start.y) * fraction,
                z: start.z + (end.z - start.z) * fraction
            }};
        }}

        function residueDisplayName(residue) {{
            const insertion = String(residue.insertion_code || "");
            const chain = String(residue.chain || "").trim();
            return (chain ? chain + ":" : "") + String(residue.residue_name || "") +
                String(residue.residue_number) + insertion;
        }}

        function addMutationMarker() {{
            const center = mutationCenter();
            if (!center) return;
            viewer.addSphere({{
                center: center,
                radius: highlightMode ? 0.52 : 0.38,
                color: mutationInfo.color,
                opacity: highlightMode ? 0.36 : 0.24
            }});
            viewer.addLabel(mutationInfo.label, {{
                position: center,
                fontColor: "#ffffff",
                backgroundColor: mutationInfo.color,
                fontSize: highlightMode ? 11 : 10,
                backgroundOpacity: 0.90,
                borderThickness: 1,
                borderColor: "#ffffff",
                inFront: true,
                alignment: "bottomLeft",
                screenOffset: {{x: 14, y: -14}}
            }});
        }}

        function addPartnerLabels() {{
            const orderedPartners = canonicalPartners.slice().sort(function(left, right) {{
                return Number(changedPartnerKeys.has(residueKey(right))) -
                    Number(changedPartnerKeys.has(residueKey(left)));
            }});
            let labelledPartners = viewMode === "focus"
                ? orderedPartners.slice(0, 8)
                : orderedPartners.filter(function(residue) {{
                    return changedPartnerKeys.has(residueKey(residue));
                }}).slice(0, 4);
            if (!labelledPartners.length) labelledPartners = orderedPartners.slice(0, 4);
            labelledPartners.forEach(function(residue, index) {{
                const atoms = selectedResidueAtoms(residue);
                if (!atoms.length) return;
                const anchor = atoms.find(function(atom) {{ return atom.atom === "CA"; }}) || atoms[0];
                const angle = index * 2.399963;
                viewer.addLabel(residueDisplayName(residue), {{
                    position: {{x: anchor.x, y: anchor.y, z: anchor.z}},
                    fontColor: "#263840",
                    backgroundColor: "#f7fafb",
                    fontSize: 8,
                    backgroundOpacity: 0.82,
                    borderThickness: 1,
                    borderColor: "#93a5ac",
                    inFront: true,
                    screenOffset: {{
                        x: Math.round(Math.cos(angle) * 16),
                        y: Math.round(Math.sin(angle) * 13)
                    }}
                }});
            }});
        }}

        function addCanonicalAtomLabels() {{
            const seen = new Set();
            let labelCount = 0;
            const ordered = canonicalContacts.slice().sort(function(left, right) {{
                const leftChanged = left.result === "Retained" ? 1 : 0;
                const rightChanged = right.result === "Retained" ? 1 : 0;
                return leftChanged - rightChanged;
            }});
            ordered.forEach(function(contact) {{
                [[contact.residue_1, contact.atom_1], [contact.residue_2, contact.atom_2]]
                    .forEach(function(pair) {{
                        if (labelCount >= 8) return;
                        const key = residueKey(pair[0]) + "|" + String(pair[1] || "");
                        if (!pair[1] || seen.has(key)) return;
                        const atom = selectCanonicalAtom(pair[0], pair[1]);
                        if (!atom) return;
                        seen.add(key);
                        labelCount += 1;
                        viewer.addLabel(String(pair[1]), {{
                            position: {{x: atom.x, y: atom.y, z: atom.z}},
                            fontColor: "#ffffff",
                            backgroundColor: "#2d3b42",
                            fontSize: 6,
                            backgroundOpacity: 0.58,
                            inFront: true,
                            screenOffset: {{
                                x: labelCount % 2 ? 5 : -5,
                                y: labelCount % 3 ? 5 : -5
                            }}
                        }});
                    }});
            }});
        }}

        function drawCanonicalContact(contact) {{
            const start = selectCanonicalAtom(contact.residue_1, contact.atom_1);
            const end = selectCanonicalAtom(contact.residue_2, contact.atom_2);
            if (!start || !end) return false;
            const colour = contact.result === "Gained" ? "#176f91" :
                (contact.result === "Lost" ? "#bd3f32" : "#7d898f");
            if (contact.interaction_class === "hydrogen_bond") {{
                const segments = 10;
                for (let index = 0; index < segments; index += 2) {{
                    viewer.addCylinder({{
                        start: pointBetween(start, end, index / segments),
                        end: pointBetween(start, end, (index + 1) / segments),
                        radius: 0.055,
                        color: colour,
                        opacity: 1.0,
                        fromCap: 1,
                        toCap: 1
                    }});
                }}
            }} else {{
                for (let index = 1; index < 10; index += 1) {{
                    viewer.addSphere({{
                        center: pointBetween(start, end, index / 10),
                        radius: 0.065,
                        color: colour,
                        opacity: 1.0
                    }});
                }}
            }}
            return true;
        }}

        function styleCanonicalPartners() {{
            canonicalPartners.forEach(function(residue) {{
                const selection = residueSelection(residue);
                if (!selection) return;
                viewer.addStyle(selection, {{stick: {{
                    colorscheme: "Jmol", radius: 0.14, opacity: 1.0
                }}}});
            }});
        }}

        function styleMutationResidue() {{
            viewer.addStyle({{resi: mutationInfo.resi, chain: mutationInfo.chain}}, {{
                stick: {{colorscheme: "Jmol", radius: 0.30, opacity: 1.0}}
            }});
        }}

        function styleCanonicalEnvironment() {{
            styleCanonicalPartners();
            styleMutationResidue();
        }}

        function drawCanonicalDetails() {{
            let rendered = 0;
            canonicalContacts.forEach(function(contact) {{
                if (drawCanonicalContact(contact)) rendered += 1;
            }});
            addPartnerLabels();
            if (viewMode === "focus") addCanonicalAtomLabels();
            addMutationMarker();
            console.log("Rendered " + rendered + " canonical interaction contacts");
        }}

        function applyBaseRepresentation(cartoonOpacity) {{
            viewer.setStyle({{}}, {{}});
            if (cartoonOn) {{
                viewer.setStyle({{}}, {{cartoon: {{
                    color: "{CARTOON_COLOR}", opacity: cartoonOpacity
                }}}});
            }}
            if (stickOn) {{
                viewer.addStyle({{}}, {{stick: {{
                    colorscheme: "greenCarbon",
                    radius: highlightMode ? 0.10 : 0.14,
                    opacity: highlightMode ? 0.30 : 1.0
                }}}});
            }}
        }}

        function applyNormalStyle() {{
            viewer.removeAllShapes();
            viewer.removeAllLabels();
            applyBaseRepresentation({FULL_CARTOON_OPACITY});
            styleMutationResidue();
            if (viewMode === "focus") {{
                styleCanonicalPartners();
                drawCanonicalDetails();
            }} else addMutationMarker();
        }}

        function applyHighlightStyle() {{
            viewer.removeAllShapes();
            viewer.removeAllLabels();
            applyBaseRepresentation({HIGHLIGHT_BACKGROUND_OPACITY});
            styleCanonicalEnvironment();
            drawCanonicalDetails();
        }}

        function applyCurrentStyle() {{
            if (highlightMode) applyHighlightStyle();
            else applyNormalStyle();
        }}

        function setHighlightMode(active) {{
            highlightMode = Boolean(active);
            applyCurrentStyle();
            updateButtonStates();
            viewer.render();
        }}

        function toggleMutationHighlight() {{
            setHighlightMode(!highlightMode);
        }}

        function toggleCartoon() {{
            cartoonOn = !cartoonOn;
            applyCurrentStyle();
            updateButtonStates();
            viewer.render();
        }}

        function toggleStick() {{
            stickOn = !stickOn;
            applyCurrentStyle();
            updateButtonStates();
            viewer.render();
        }}

        function updateButtonStates() {{
            btnCartoon.classList.toggle('active', cartoonOn);
            btnStick.classList.toggle('active', stickOn);
            btnHighlight.classList.toggle('active', highlightMode);
            legend.classList.toggle('active', highlightMode || viewMode === "focus");
        }}

        // Keyboard shortcuts
        document.addEventListener('keydown', function(e) {{
            if (e.key === 'h' || e.key === 'H') {{
                toggleMutationHighlight();
            }} else if (e.key === 'c' || e.key === 'C') {{
                toggleCartoon();
            }} else if (e.key === 's' || e.key === 'S') {{
                toggleStick();
            }}
        }});
    </script>
</body>
</html>'''
            return html
        
        # Generate HTML files
        print("Generating enhanced HTML for original structure...")
        original_html = create_enhanced_html(original_pdb_content, f"Original Structure - {file_name}", False)
        
        # Save original viewer
        original_html_path = os.path.join(output_dir, "original_3d.html")
        with open(original_html_path, 'w') as f:
            f.write(original_html)
        print(f"  ✓ Created: {original_html_path}")
        
        print("Generating enhanced HTML for mutated structure...")
        mutated_html = create_enhanced_html(mutated_pdb_content, f"Mutated Structure - {file_name}{mutation}", True)
        
        mutated_html_path = os.path.join(output_dir, "mutated_3d.html")
        with open(mutated_html_path, 'w') as f:
            f.write(mutated_html)
        print(f"  ✓ Created: {mutated_html_path}")
        
        # Create combined viewer (this is what PHP expects)
        combined_html_path = os.path.join(output_dir, "interactive_3d_viewer.html")
        combined_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <title>X-Pro synchronized 3D viewer</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
        * {{ box-sizing: border-box; }}
        html, body {{
            margin: 0; width: 100%; height: 100%; min-height: 0;
            background: {VIEWER_BACKGROUND_COLOR};
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        }}
        .viewer-shell {{
            display: grid; grid-template-rows: auto minmax(0, 1fr);
            width: 100%; height: 100%; min-height: 0;
        }}
        .shared-toolbar {{
            display: flex; flex-wrap: wrap; align-items: center; justify-content: center;
            gap: 8px 14px; padding: 8px 12px; color: #dce6ea;
            background: #111b22; border-bottom: 1px solid #52646d;
        }}
        .shared-toolbar button {{
            color: #eef5f7; background: #175f78; border: 1px solid #4babc9;
            border-radius: 5px; padding: 7px 11px; font-size: 11px;
            font-weight: 700; cursor: pointer;
        }}
        .shared-toolbar label {{ font-size: 11px; font-weight: 700; }}
        .split-view {{
            display: grid; grid-template-columns: repeat(2, minmax(0, 1fr));
            width: 100%; height: 100%; min-height: 0;
        }}
        .split-pane {{
            min-width: 0; min-height: 0; border-right: 1px solid #65747b;
            overflow: hidden;
        }}
        .split-pane:last-child {{ border-right: 0; }}
        iframe {{ display: block; width: 100%; height: 100%; border: 0; background: {VIEWER_BACKGROUND_COLOR}; }}
        @media (max-width: 720px) {{
            html, body {{ height: auto; min-height: 1100px; }}
            .viewer-shell {{ min-height: 1100px; }}
            .split-view {{
                grid-template-columns: 1fr;
                grid-template-rows: repeat(2, minmax(510px, 1fr));
                min-height: 1040px;
            }}
            .split-pane {{ border-right: 0; border-bottom: 1px solid #65747b; }}
            .split-pane:last-child {{ border-bottom: 0; }}
        }}
    </style>
</head>
<body>
    <div class="viewer-shell">
        <nav class="shared-toolbar" aria-label="Shared 3D viewer controls">
            <button type="button" data-command="whole">Whole protein</button>
            <button type="button" data-command="focus">Focus mutation site</button>
            <button type="button" id="resetCommonView">Reset common view</button>
            <label><input type="checkbox" id="syncViews" checked> Synchronize views</label>
        </nav>
        <main class="split-view">
            <section class="split-pane"><iframe id="wtViewer" title="Original wild-type structure" src="original_3d.html"></iframe></section>
            <section class="split-pane"><iframe id="mutantViewer" title="Mutant structure" src="mutated_3d.html"></iframe></section>
        </main>
    </div>
    <script>
        const wtViewer = document.getElementById('wtViewer');
        const mutantViewer = document.getElementById('mutantViewer');
        const syncViews = document.getElementById('syncViews');
        let wtReferenceView = null;

        function sendCommand(command) {{
            [wtViewer, mutantViewer].forEach(function(frame) {{
                if (frame.contentWindow) frame.contentWindow.postMessage({{
                    type: 'xpro-viewer-command', command: command
                }}, window.location.origin);
            }});
        }}
        function requestWtView() {{
            if (wtViewer.contentWindow) wtViewer.contentWindow.postMessage({{
                type: 'xpro-viewer-command', command: 'report-view'
            }}, window.location.origin);
        }}
        document.querySelectorAll('[data-command]').forEach(function(button) {{
            button.addEventListener('click', function() {{
                sendCommand(button.dataset.command);
                window.setTimeout(requestWtView, 120);
            }});
        }});
        document.getElementById('resetCommonView').addEventListener('click', function() {{
            sendCommand('whole');
            window.setTimeout(requestWtView, 120);
        }});

        window.addEventListener('message', function(event) {{
            if (event.origin !== window.location.origin) return;
            const message = event.data || {{}};
            if (message.type === 'xpro-viewer-ready' && Array.isArray(message.view)) {{
                if (message.pane === 'wt') {{
                    wtReferenceView = message.view;
                    if (mutantViewer.contentWindow) mutantViewer.contentWindow.postMessage({{
                        type: 'xpro-apply-view', pane: 'wt', view: wtReferenceView
                    }}, window.location.origin);
                }} else if (wtReferenceView && mutantViewer.contentWindow) {{
                    mutantViewer.contentWindow.postMessage({{
                        type: 'xpro-apply-view', pane: 'wt', view: wtReferenceView
                    }}, window.location.origin);
                }}
                return;
            }}
            if (message.type !== 'xpro-view-change' || !Array.isArray(message.view)) return;
            if (message.pane === 'wt') wtReferenceView = message.view;
            if (!syncViews.checked) return;
            const target = message.pane === 'wt' ? mutantViewer : wtViewer;
            if (target.contentWindow) target.contentWindow.postMessage({{
                type: 'xpro-apply-view', pane: message.pane, view: message.view
            }}, window.location.origin);
        }});
    </script>
</body>
</html>"""
        with open(combined_html_path, 'w') as f:
            f.write(combined_html)
        print(f"  ✓ Created: {combined_html_path}")
        
        # Also save a backup copy
        backup_path = os.path.join(output_dir, "viewer_backup.html")
        shutil.copy(combined_html_path, backup_path)
        print(f"  ✓ Created backup: {backup_path}")
        
        # Set permissions
        for f in os.listdir(output_dir):
            if f.endswith('.html'):
                os.chmod(os.path.join(output_dir, f), 0o644)
        
        # Verify files were created
        print("\n📁 VERIFYING 3D HTML FILES:")
        if os.path.exists(combined_html_path):
            file_size = os.path.getsize(combined_html_path)
            print(f"  ✅ interactive_3d_viewer.html ({file_size} bytes)")
        else:
            print(f"  ❌ interactive_3d_viewer.html not found")
            
        if os.path.exists(original_html_path):
            print(f"  ✅ original_3d.html ({os.path.getsize(original_html_path)} bytes)")
        else:
            print(f"  ❌ original_3d.html not found")
            
        if os.path.exists(mutated_html_path):
            print(f"  ✅ mutated_3d.html ({os.path.getsize(mutated_html_path)} bytes)")
        else:
            print(f"  ❌ mutated_3d.html not found")
        
        return True
        
    except Exception as e:
        print(f"\n!!! ERROR GENERATING 3D HTML: {e}")
        traceback.print_exc()
        return False

# ======================================================
# FUNCTION: CREATE VIEWER DATA
# ======================================================

def create_viewer_data(file_name, chain_name, atom_number, mutation):
    """
    Creates JSON data for 3D visualization
    """
    print("=== CREATING 3D VIEWER DATA ===")
    
    try:
        # Ensure PDB file exists
        if not ensure_pdb_file(file_name):
            print("[ERROR] Cannot create viewer data without PDB file.")
            return False
        
        # Read PDB files to extract atom positions for highlighting
        def extract_pdb_data(pdb_file, is_mutated=False):
            """Extract atom data from PDB file for 3D visualization"""
            atoms = []
            if not os.path.exists(pdb_file):
                print(f"  Warning: {pdb_file} not found")
                return atoms
                
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        try:
                            atom_name = line[12:16].strip()
                            residue_name = line[17:20].strip()
                            chain = line[21:22].strip()
                            res_num = int(line[22:26].strip())
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            element = line[76:78].strip() or atom_name[0]
                            
                            # Determine if this is the mutation site
                            is_target = (res_num == int(atom_number) and chain == chain_name)
                            
                            atoms.append({
                                'atom_name': atom_name,
                                'residue_name': residue_name,
                                'chain': chain,
                                'res_num': res_num,
                                'x': x, 'y': y, 'z': z,
                                'element': element,
                                'is_target': is_target,
                                'is_mutated': is_mutated
                            })
                        except Exception as e:
                            continue
            return atoms
        
        # Extract data from both structures
        original_atoms = extract_pdb_data(f"{file_name}.pdb", False)
        mutated_atoms = extract_pdb_data("mutated.pdb", True)
        
        print(f"  Original structure: {len(original_atoms)} atoms")
        print(f"  Mutated structure: {len(mutated_atoms)} atoms")
        
        # Find target residue position for initial zoom
        target_pos = None
        for atom in original_atoms:
            if atom['is_target']:
                target_pos = [atom['x'], atom['y'], atom['z']]
                break
        
        if not target_pos and mutated_atoms:
            for atom in mutated_atoms:
                if atom['is_target']:
                    target_pos = [atom['x'], atom['y'], atom['z']]
                    break
        
        if not target_pos and original_atoms:
            target_pos = [original_atoms[0]['x'], original_atoms[0]['y'], original_atoms[0]['z']]
        else:
            target_pos = [0, 0, 0]
        
        # Create JSON data for the web viewer
        viewer_data = {
            'pdb_id': file_name,
            'chain': chain_name,
            'residue_number': int(atom_number),
            'wild_type': mutation if mutated_atoms else '',
            'mutation': mutation,
            'target_position': target_pos,
            'structures': {
                'original': {
                    'file': f"{file_name}.pdb",
                    'atoms': original_atoms[:2000]  # Limit for performance
                },
                'mutated': {
                    'file': "mutated.pdb",
                    'atoms': mutated_atoms[:2000]
                }
            }
        }
        
        # Save JSON data for PHP to use
        json_path = os.path.join("output_files", "viewer_data.json")
        with open(json_path, 'w') as f:
            json.dump(viewer_data, f, indent=2)
        
        print("✓ Created viewer data: viewer_data.json")
        return True
        
    except Exception as e:
        print(f"\n!!! ERROR CREATING VIEWER DATA !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return False

# ======================================================
# FUNCTION: PYMOL MUTAGENESIS
# ======================================================

def pymol_mutagenesis(file_name, chain_name, atom_number, mutation):
    """
    Performs a site-directed mutation using PyMOL's mutagenesis wizard.
    """
    print(f"=== STARTING PYMOL MUTAGENESIS ===")
    print(f"Parameters: {file_name=}, {chain_name=}, {atom_number=}, {mutation=}")

    if not PYMOL_AVAILABLE:
        print("[ERROR] PyMOL is not available. Cannot perform PyMOL mutagenesis.")
        return None

    try:
        # Ensure PDB file exists
        if not ensure_pdb_file(file_name):
            print("[ERROR] Cannot proceed with mutagenesis without PDB file.")
            return None

        # Default rotamer counts by amino acid
        default_rotamers = {
            'ALA': 1, 'GLY': 1, 'VAL': 3, 'SER': 3, 'THR': 3, 'CYS': 3,
            'LEU': 6, 'ILE': 7, 'ASN': 6, 'ASP': 6, 'MET': 13, 'LYS': 27,
            'ARG': 34, 'GLN': 9, 'GLU': 9, 'PHE': 4, 'TYR': 4, 'TRP': 6,
            'HIS': 4, 'PRO': 2
        }

        # Step 1: Initialize PyMOL session
        cmd.reinitialize()
        print("✓ PyMOL reinitialized")

        pdb_file = f"{file_name}.pdb"
        cmd.load(pdb_file, file_name)
        print("✓ PDB structure loaded")

        # Save altered PDB before mutation (for LigPlot)
        cmd.select("hetera", "hetatm")
        cmd.alter("hetera", "type='ATOM'")
        # Create selection for the residue to be mutated
        atom_sel = f"/{file_name}//{chain_name}/{atom_number}/"
        cmd.select("mutant", atom_sel)
        cmd.alter("mutant", "type='HETATM'")
        cmd.save("altered.pdb")
        print("✓ Saved: altered.pdb")

        # Perform mutation using mutagenesis wizard
        cmd.reinitialize()
        cmd.load(pdb_file, file_name)
            
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
                rotamer_count = default_rotamers.get(mutation, 6)
                print(f"[INFO] Using default rotamer count for {mutation}: {rotamer_count}")
                
        except Exception as e:
            print(f"[WARNING] Could not determine rotamer count: {e}")
            rotamer_count = default_rotamers.get(mutation, 6)
            print(f"[INFO] Using default rotamer count: {rotamer_count}")
        
        cmd.get_wizard().apply()
        print("✓ Mutation applied successfully")

        # Save mutated structure
        cmd.save("mutated.pdb")
        print("✓ Saved: mutated.pdb")

        # Save rotamer count to file
        os.makedirs("output_files", exist_ok=True)
        with open(f"output_files/rotamer_count_{file_name}.txt", "w") as f:
            f.write(str(rotamer_count))

        print("=== PYMOL MUTAGENESIS COMPLETED SUCCESSFULLY ===")
        return rotamer_count

    except Exception as e:
        print("\n!!! ERROR IN PYMOL MUTAGENESIS !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return None

# ======================================================
# FUNCTION: MODELLER MUTAGENESIS
# ======================================================

def modeller_mutagenesis(file_name, chain_name, atom_number, mutation):
    """
    Performs in-silico point mutation using Modeller
    """
    print(f"=== STARTING MODELLER MUTAGENESIS ===")
    print(f"Parameters: {file_name=}, {chain_name=}, {atom_number=}, {mutation=}")

    if not MODELLER_AVAILABLE:
        print("[ERROR] Modeller is not available. Cannot perform Modeller mutagenesis.")
        return None

    try:
        pdb_file = f"{file_name}.pdb"

        # First, create the altered.pdb (original structure) for LigPlot
        if os.path.exists(pdb_file):
            shutil.copy(pdb_file, "altered.pdb")
            print("✓ Created: altered.pdb (original structure for LigPlot)")

        # Modeller optimization functions
        def optimize(atmsel, sched):
            for step in sched:
                step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
            cg = ConjugateGradients()
            cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)

        def refine(atmsel):
            md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                                   md_return='FINAL')
            init_vel = True
            for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                        (200, 600,
                                         (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
                for temp in temps:
                    md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                                 max_iterations=its, equilibrate=equil)
                    init_vel = False

        def make_restraints(mdl1, aln):
            rsr = mdl1.restraints
            rsr.clear()
            s = Selection(mdl1)
            for typ in ('stereo', 'phi-psi_binormal'):
                rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
            for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
                rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                         spline_dx=0.3, spline_min_points=5, aln=aln,
                         spline_on_site=True)

        log.verbose()

        env = environ(rand_seed=-49837)
        env.io.hetatm = True
        env.edat.dynamic_sphere = False
        env.edat.dynamic_lennard = True
        env.edat.contact_shell = 4.0
        env.edat.update_dynamic = 0.39

        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')

        # Read the original PDB file
        mdl1 = Model(env, file=pdb_file)
        ali = Alignment(env)
        ali.append_model(mdl1, atom_files=pdb_file, align_codes=file_name)

        # Select residue using PDB residue number and chain (e.g. '15:E')
        # This is the correct Modeller approach - do NOT use 0-based integer indexing,
        # which maps to internal Modeller indices and not PDB author residue numbers.
        residue_id = f"{atom_number}:{chain_name}"
        try:
            target_residue = mdl1.residues[residue_id]
        except KeyError:
            # Fallback: search all residues for matching PDB resnum and chain
            target_residue = None
            for res in mdl1.residues:
                if str(res.num) == str(atom_number) and res.chain.name == chain_name:
                    target_residue = res
                    break
            if target_residue is None:
                print(f"[ERROR] Residue {atom_number} in chain {chain_name} not found in PDB.")
                print(f"  Available chains: {[ch.name for ch in mdl1.chains]}")
                print(f"  First 10 residues: {[(r.num, r.chain.name, r.pdb_name) for r in list(mdl1.residues)[:10]]}")
                return None

        print(f"✓ Target residue identified: {target_residue.pdb_name} {target_residue.num}:{chain_name}")

        # Perform mutation
        s = Selection(target_residue)
        s.mutate(residue_type=mutation)
        ali.append_model(mdl1, align_codes=file_name)

        # Generate molecular topology for mutant
        mdl1.clear_topology()
        mdl1.generate_topology(ali[-1])
        mdl1.transfer_xyz(ali)
        mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

        mdl2 = Model(env, file=pdb_file)
        mdl1.res_num_from(mdl2, ali)

        # Write temporary file
        tmp_file = f"{file_name}{mutation}{atom_number}.tmp"
        mdl1.write(file=tmp_file)
        mdl1.read(file=tmp_file)

        # Set up restraints
        make_restraints(mdl1, ali)
        mdl1.env.edat.nonbonded_sel_atoms = 1

        sched = autosched.loop.make_for_model(mdl1)
        # Re-select the target residue after re-reading the tmp file
        try:
            target_residue_new = mdl1.residues[residue_id]
        except KeyError:
            target_residue_new = None
            for res in mdl1.residues:
                if str(res.num) == str(atom_number) and res.chain.name == chain_name:
                    target_residue_new = res
                    break
        if target_residue_new is None:
            print(f"[ERROR] Could not re-select residue {atom_number}:{chain_name} after topology rebuild.")
            return None
        s = Selection(target_residue_new)

        mdl1.restraints.unpick_all()
        mdl1.restraints.pick(s)
        s.energy()
        s.randomize_xyz(deviation=4.0)

        mdl1.env.edat.nonbonded_sel_atoms = 2
        optimize(s, sched)
        mdl1.env.edat.nonbonded_sel_atoms = 1
        optimize(s, sched)
        s.energy()

        # Save mutated structure with consistent naming
        output_name = f"{file_name}{mutation}{atom_number}.pdb"
        mdl1.write(file=output_name)
        
        # Copy to mutated.pdb for consistent workflow
        shutil.copy(output_name, "mutated.pdb")
        print(f"✓ Saved: {output_name} and mutated.pdb")

        # Clean up temporary file
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

        print("=== MODELLER MUTAGENESIS COMPLETED SUCCESSFULLY ===")
        return 1

    except Exception as e:
        print("\n!!! ERROR IN MODELLER MUTAGENESIS !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return None

# ======================================================
# FUNCTION: CONVERT MUTATED PDB TO PDBx/mmCIF
# ======================================================

def convert_mutated_pdb_to_cif(pdb_path="mutated.pdb", cif_path="mutated.cif"):
    """
    Convert the mutated structure (legacy PDB format) into PDBx/mmCIF format,
    so the user can download the mutated structure as .cif as well as .pdb.
    """
    print(f"\n=== CONVERTING {pdb_path} TO PDBx/mmCIF ===")
    try:
        if not os.path.exists(pdb_path):
            print(f"[WARNING] {pdb_path} not found. Skipping CIF conversion.")
            return False

        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        structure = parser.get_structure("mutated_structure", pdb_path)

        io = MMCIFIO()
        io.set_structure(structure)
        io.save(cif_path)

        print(f"✓ Saved: {cif_path}")
        return True
    except Exception as e:
        print(f"[ERROR] Failed to convert {pdb_path} to mmCIF: {e}")
        traceback.print_exc()
        return False

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
        
        # Analyze altered structure (original / wild-type)
        if os.path.exists("altered.pdb"):
            shutil.copy("altered.pdb", "output_files/")
            # LigPlot embeds the exact input filename as a label directly into
            # the rendered diagram, so feed it a same-content copy named
            # "WT-Ligplot.pdb" purely so the image displays "WT-Ligplot"
            # instead of "altered". All other internal naming (altered_*,
            # altered.pdb, CSV/summary files, download buttons) is untouched
            # and still keyed off the "altered_{pdb_name}" output prefix below.
            wildtype_input = "WT-Ligplot.pdb"
            shutil.copy("altered.pdb", os.path.join("output_files", wildtype_input))
            success = run_enhanced_ligplot(wildtype_input, f"altered_{pdb_name}")
            results.append(("altered", success))
        else:
            print("[WARNING] altered.pdb not found for LigPlot analysis")
        
        # Analyze mutated structure  
        if os.path.exists("mutated.pdb"):
            shutil.copy("mutated.pdb", "output_files/")
            # Same trick as above: feed LigPlot a same-content copy named
            # "Mutant-Ligplot.pdb" so the rendered diagram's baked-in caption
            # reads "Mutant-Ligplot" instead of "mutated". The annotated
            # ("Differences-Ligplot") panel reuses this same native LigPlot
            # output -- see annotate_ligplot_postscript() in
            # interaction_comparison.py, which relabels that caption again
            # specifically on the annotated copy, so the two panels display
            # different captions despite sharing one LigPlot run. All other
            # internal naming is still keyed off "mutated_{pdb_name}" below.
            mutant_input = "Mutant-Ligplot.pdb"
            shutil.copy("mutated.pdb", os.path.join("output_files", mutant_input))
            success = run_enhanced_ligplot(mutant_input, f"mutated_{pdb_name}")
            results.append(("mutated", success))
        else:
            print("[WARNING] mutated.pdb not found for LigPlot analysis")

        print("=== ENHANCED LIGPLOT_INTERFACE COMPLETED ===")
        return all(success for _, success in results)

    except Exception as e:
        print("\n!!! ERROR IN ENHANCED LIGPLOT_INTERFACE !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return False

# ======================================================
# FUNCTION: CONVERT PS TO PNG
# ======================================================

def enhanced_ps_to_png(pdb_name):
    """Convert LigPlot PostScript files to PNG using system commands"""
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
                else:
                    print(f"  [WARNING] ImageMagick failed with return code {result.returncode}")
            except Exception as e:
                print(f"  [WARNING] ImageMagick not available: {e}")
            
            # Method 2: Try system gs command (Ghostscript)
            try:
                cmd_gs = f"gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 -sOutputFile={png_file} {ps_file}"
                result = subprocess.run(cmd_gs, shell=True, capture_output=True, timeout=30)
                if result.returncode == 0:
                    print(f"✓ Converted {ps_file} → {png_file} (Ghostscript)")
                else:
                    print(f"  [WARNING] Ghostscript failed with return code {result.returncode}")
            except Exception as e:
                print(f"  [WARNING] Ghostscript not available: {e}")
        
        return True
    except Exception as e:
        print(f"  [ERROR] PS to PNG conversion failed: {e}")
        return False

# ======================================================
# FUNCTION: CORRECTED INTERACTION COMPARISON
# ======================================================

def enhanced_file_cmp(pdb_name):
    """
    CORRECTED comparison of LigPlot results with proper logic and CSV formatting
    """
    print(f"\n=== STARTING CORRECTED FILE_CMP ===")
    try:
        # Convert PS to PNG first
        enhanced_ps_to_png(pdb_name)

        mutated_sum = f"output_files/mutated_{pdb_name}_ligplot.sum"
        altered_sum = f"output_files/altered_{pdb_name}_ligplot.sum"
        
        if not os.path.exists(mutated_sum) or not os.path.exists(altered_sum):
            print("[ERROR] Missing LigPlot summary files.")
            return False

        # Read and parse the summary files
        def parse_sum_file(filename):
            """Parse LigPlot summary file and return DataFrame with proper columns"""
            interactions_data = []
            with open(filename, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        if line.startswith('Residues:'):
                            parts = line.split('] -> [')
                            if len(parts) == 2:
                                residue1_part = parts[0].replace('Residues: [', '').strip()
                                residue2_and_numbers = parts[1].replace(']', '').strip()
                                
                                residue2_parts = residue2_and_numbers.split()
                                if len(residue2_parts) >= 3:
                                    residue2 = ' '.join(residue2_parts[:3])
                                    numbers = residue2_parts[3:]
                                    
                                    m_m = numbers[0] if len(numbers) > 0 else '0'
                                    s_s = numbers[1] if len(numbers) > 1 else '0'
                                    m_s = numbers[2] if len(numbers) > 2 else '0'
                                    non_bonded = numbers[3] if len(numbers) > 3 else '0'
                                    
                                    interaction_data = {
                                        'Interacting_residues': f"[{residue1_part}] -> [{residue2}]",
                                        'Residue1': residue1_part,
                                        'Residue2': residue2,
                                        'M_M': m_m,
                                        'S_S': s_s,
                                        'M_S': m_s,
                                        'Non_bonded_Contact': non_bonded,
                                        'Full_Interaction': line
                                    }
                                    interactions_data.append(interaction_data)
            return pd.DataFrame(interactions_data)

        # Parse both files into DataFrames
        print("Parsing altered summary file...")
        original_df = parse_sum_file(altered_sum)
        print("Parsing mutated summary file...")
        mutated_df = parse_sum_file(mutated_sum)
        
        print(f"Original structure interactions: {len(original_df)}")
        print(f"Mutated structure interactions: {len(mutated_df)}")

        # Find differences
        original_interactions_set = set(original_df['Full_Interaction']) if not original_df.empty else set()
        mutated_interactions_set = set(mutated_df['Full_Interaction']) if not mutated_df.empty else set()
        
        lost_interactions_set = original_interactions_set - mutated_interactions_set
        gained_interactions_set = mutated_interactions_set - original_interactions_set

        lost_df = original_df[original_df['Full_Interaction'].isin(lost_interactions_set)]
        gained_df = mutated_df[mutated_df['Full_Interaction'].isin(gained_interactions_set)]

        # Save CSV files
        if not lost_df.empty:
            lost_df[['Interacting_residues', 'M_M', 'S_S', 'M_S', 'Non_bonded_Contact']].to_csv(
                f"output_files/altered_{pdb_name}_ligplot.csv", index=False)
            print(f"✓ Saved CSV: output_files/altered_{pdb_name}_ligplot.csv")
        
        if not gained_df.empty:
            gained_df[['Interacting_residues', 'M_M', 'S_S', 'M_S', 'Non_bonded_Contact']].to_csv(
                f"output_files/mutated_{pdb_name}_ligplot.csv", index=False)
            print(f"✓ Saved CSV: output_files/mutated_{pdb_name}_ligplot.csv")

        # Save detailed versions
        if not original_df.empty:
            original_df[['Interacting_residues', 'M_M', 'S_S', 'M_S', 'Non_bonded_Contact']].to_csv(
                f"output_files/original_interactions_{pdb_name}.csv", index=False)

        if not mutated_df.empty:
            mutated_df[['Interacting_residues', 'M_M', 'S_S', 'M_S', 'Non_bonded_Contact']].to_csv(
                f"output_files/mutated_interactions_{pdb_name}.csv", index=False)

        if not lost_df.empty:
            lost_df[['Interacting_residues', 'M_M', 'S_S', 'M_S', 'Non_bonded_Contact']].to_csv(
                f"output_files/lost_interactions_{pdb_name}.csv", index=False)

        if not gained_df.empty:
            gained_df[['Interacting_residues', 'M_M', 'S_S', 'M_S', 'Non_bonded_Contact']].to_csv(
                f"output_files/gained_interactions_{pdb_name}.csv", index=False)

        # Write comparison results
        with open("output_files/ligplot_differences.txt", "w") as f3:
            f3.write("=== INTERACTIONS LOST (Present in original, lost after mutation) ===\n")
            if not lost_df.empty:
                for _, row in lost_df.iterrows():
                    f3.write(row['Full_Interaction'] + "\n")
            else:
                f3.write("None\n")

            f3.write("\n=== INTERACTIONS GAINED (New interactions after mutation) ===\n")
            if not gained_df.empty:
                for _, row in gained_df.iterrows():
                    f3.write(row['Full_Interaction'] + "\n")
            else:
                f3.write("None\n")

            f3.write(f"\n=== SUMMARY ===\n")
            f3.write(f"Total interactions in original structure: {len(original_df)}\n")
            f3.write(f"Total interactions in mutated structure: {len(mutated_df)}\n")
            f3.write(f"Interactions lost: {len(lost_df)}\n")
            f3.write(f"Interactions gained: {len(gained_df)}\n")
            f3.write(f"Net change: {len(gained_df) - len(lost_df)}\n")

        print("=== CORRECTED FILE_CMP COMPLETED SUCCESSFULLY ===")
        return True

    except Exception as e:
        print("\n!!! ERROR IN CORRECTED FILE_CMP !!!")
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()
        return False

# ======================================================
# FUNCTION: ENHANCED POINT MUTATOR (MAIN PIPELINE)
# ======================================================

def enhanced_point_mutator(file_name, chain_name, atom_number, mutation, randomNumber, method=None):
    """
    Enhanced main pipeline combining mutagenesis with advanced analysis
    """
    print(f"\n=== STARTING ENHANCED POINT_MUTATOR ===")
    print(f"Parameters: {file_name=}, {chain_name=}, {atom_number=}, {mutation=}, {randomNumber=}, {method=}")

    target_dir = None
    try:
        # Set up working directory
        base_dir = '/var/www/html' if '/var/www/html' in os.path.abspath(__file__) else os.getcwd()
        target_dir = os.path.join(base_dir, "fileUpload", randomNumber)
        os.makedirs(target_dir, exist_ok=True)
        os.chdir(target_dir)
        os.makedirs("output_files", exist_ok=True)
        print(f"✓ Working directory: {os.getcwd()}")

        # Test write permissions
        test_file = os.path.join("output_files", "test_write.txt")
        try:
            with open(test_file, 'w') as f:
                f.write("test")
            print(f"✓ Can write to output_files/")
            os.remove(test_file)
        except Exception as e:
            print(f"✗ Cannot write to output_files/: {e}")

        # Ensure PDB file exists
        if not ensure_pdb_file(file_name):
            print("[ERROR] Cannot proceed without PDB file.")
            return None

        # Find LigPlus early
        ligplus_path = find_ligplus()

        # Choose mutagenesis method
        if method and method.lower() == 'modeller' and MODELLER_AVAILABLE:
            method = "Modeller"
            print("✓ Using Modeller for mutagenesis (user specified)")
        elif method and method.lower() == 'pymol' and PYMOL_AVAILABLE:
            method = "PyMOL"
            print("✓ Using PyMOL for mutagenesis (user specified)")
        elif PYMOL_AVAILABLE:
            method = "PyMOL"
            print("✓ Using PyMOL for mutagenesis (default)")
        elif MODELLER_AVAILABLE:
            method = "Modeller"
            print("✓ Using Modeller for mutagenesis (PyMOL not available)")
        else:
            print("❌ No mutagenesis methods available. Please install PyMOL or Modeller.")
            return None

        # Perform mutation based on chosen method
        if method == "PyMOL":
            result = pymol_mutagenesis(file_name, chain_name, atom_number, mutation)
        else:  # Modeller
            result = modeller_mutagenesis(file_name, chain_name, atom_number, mutation)

        # Generate PyMOL visualization images for BOTH methods
        if result is not None and PYMOL_AVAILABLE:
            visualization_success = pymol_visualization(file_name, chain_name, atom_number, mutation, method)
            if not visualization_success:
                print("[WARNING] PyMOL visualization failed, but mutagenesis was successful")

        if result is None:
            print("❌ Mutagenesis failed.")
            return None

        write_progress_step(target_dir, "progress_1_model.txt")

        # Convert the mutated structure to PDBx/mmCIF format as well,
        # so the user can download it in either format
        convert_mutated_pdb_to_cif("mutated.pdb", "mutated.cif")

        # Create 3D viewer data for web integration
        if result is not None:
            viewer_data_success = create_viewer_data(file_name, chain_name, atom_number, mutation)
            if viewer_data_success:
                print("✓ 3D viewer data created successfully")
            else:
                print("[WARNING] Failed to create 3D viewer data")

        # Run enhanced LigPlot analysis
        ligplot_success = enhanced_ligplot_interface(ligplus_path, chain_name, atom_number, file_name)
        if not ligplot_success:
            print("[WARNING] LigPlot analysis was not completed successfully")

        # Convert LigPlot's .ps diagrams to the actual displayed .png files.
        # This must run regardless of the CSV/comparison logic below: it's
        # what overwrites the earlier PyMOL close-up screenshot (written to
        # the same output_files/altered_{pdb}_ligplot.png /
        # mutated_{pdb}_ligplot.png filenames during 3D visualization) with
        # the real 2D LigPlot interaction diagram. enhanced_file_cmp() used
        # to trigger this as a side effect; calling it directly here keeps
        # that working now that enhanced_file_cmp() itself is no longer
        # called on the normal success path.
        enhanced_ps_to_png(file_name)

        write_progress_step(target_dir, "progress_2_ligplot.txt")

        # Build the canonical WT-vs-mutant interaction comparison (replaces
        # the old enhanced_file_cmp CSV-based Interactions Lost/Gained
        # tables). This reads the same .hhb/.nnb LigPlot output files
        # enhanced_ligplot_interface() just produced, so no changes were
        # needed to the LigPlot invocation itself.
        try:
            print("\n=== BUILDING CANONICAL INTERACTION COMPARISON ===")
            interaction_comparison.create_comparison_outputs(
                output_dir="output_files",
                pdb_name=file_name,
                wt_pdb="altered.pdb",
                mutant_pdb="mutated.pdb",
                mutation_chain=chain_name,
                mutation_residue=atom_number,
                expected_mutant_name=mutation,
                annotate_ligplots=True,
            )
            print("✓ Canonical interaction comparison written to "
                  f"output_files/interaction_comparison_{file_name}.json")
        except Exception as e:
            print(f"[WARNING] Could not build canonical interaction comparison: {e}")
            traceback.print_exc()
            # Fall back to the legacy CSV-based comparison so the page still
            # has something to show rather than nothing at all.
            try:
                enhanced_file_cmp(file_name)
            except Exception as fallback_error:
                print(f"[WARNING] Legacy comparison fallback also failed: {fallback_error}")

        write_progress_step(target_dir, "progress_3_comparison.txt")

        # Build the viewer only after the canonical comparison exists (or has
        # failed over to the legacy fallback above) so its "Focus mutation
        # site" overlay can draw the same HBPLUS/LigPlot contact records
        # shown in output_files/interaction_comparison_{file_name}.json.
        print("\n" + "="*50)
        print("🚀 ATTEMPTING TO GENERATE 3D HTML FILES")
        print("="*50)
        
        sys.stdout.flush()
        
        html_success = generate_3d_html(file_name, chain_name, atom_number, mutation, randomNumber, method)
        
        if html_success:
            print("✅ 3D HTML files generated successfully")
        else:
            print("⚠️ Failed to generate 3D HTML files, but continuing...")

        print("\n=== ENHANCED POINT_MUTATOR COMPLETED SUCCESSFULLY ===")
        write_progress_step(target_dir, "progress_4_complete.txt")
        try:
            marker = os.path.join(target_dir, "processing_started.txt")
            if os.path.exists(marker):
                os.remove(marker)
        except Exception:
            pass
        return result

    except Exception as e:
        print("\n!!! ERROR IN ENHANCED POINT_MUTATOR !!!")
        traceback.print_exc()
        
        # Try to create empty CSV files even on error
        try:
            os.makedirs("output_files", exist_ok=True)
            with open(f"output_files/altered_{file_name}_ligplot.csv", 'w') as f:
                f.write("Interacting_residues,M_M,S_S,M_S,Non_bonded_Contact\n")
            with open(f"output_files/mutated_{file_name}_ligplot.csv", 'w') as f:
                f.write("Interacting_residues,M_M,S_S,M_S,Non_bonded_Contact\n")
            print("✓ Created empty CSV files as final fallback")
        except:
            pass

        # Clean up the "analysis in progress" lock file even on failure, and
        # write a distinct marker the progress page can use to stop polling
        # and show an error instead of waiting forever.
        if target_dir:
            try:
                marker = os.path.join(target_dir, "processing_started.txt")
                if os.path.exists(marker):
                    os.remove(marker)
            except Exception:
                pass
            write_progress_step(target_dir, "progress_failed.txt")

        return None

# ======================================================
# MAIN SCRIPT ENTRY
# ======================================================

if __name__ == "__main__":
    print("\n=== PARSING COMMAND LINE ARGUMENTS ===")
    
    if len(sys.argv) == 6:
        pdbID, chainID, resiNum, mut, rndNum = sys.argv[1:6]
        method = None
    elif len(sys.argv) == 7:
        pdbID, chainID, resiNum, mut, method, rndNum = sys.argv[1:7]
    else:
        print("Usage: python3 xpro.py <PDB_ID> <CHAIN_ID> <RESI_NUM> <MUTATION> <RANDOM_NUM>")
        print("   or: python3 xpro.py <PDB_ID> <CHAIN_ID> <RESI_NUM> <MUTATION> <METHOD> <RANDOM_NUM>")
        print("Example: python3 xpro.py 1ATP E 50 ALA 12345")
        print("Example: python3 xpro.py 1ATP E 50 ALA modeller 12345")
        sys.exit(1)

    print(f"Arguments parsed:")
    print(f"  PDB ID: {pdbID}")
    print(f"  Chain ID: {chainID}")
    print(f"  Residue Number: {resiNum}")
    print(f"  Mutation: {mut}")
    if method:
        print(f"  Method: {method}")
    print(f"  Random Number: {rndNum}")

    print("\n=== EXECUTING ENHANCED POINT_MUTATOR ===")
    rotamer_count = enhanced_point_mutator(pdbID, chainID, resiNum, mut, rndNum, method)

    if rotamer_count is not None:
        print("\n=== FINAL RESULT ===")
        print(f"Mutagenesis completed successfully")
        if isinstance(rotamer_count, int) and rotamer_count > 1:
            print(f"Rotamer count: {rotamer_count}")
        print("\n=== OUTPUT FILES GENERATED ===")
        print("  • mutated.pdb - Mutated structure")
        print("  • altered.pdb - Original structure")
        print("  • output_files/interactive_3d_viewer.html - Combined 3D viewer")
        print("  • output_files/original_3d.html - Original structure viewer")
        print("  • output_files/mutated_3d.html - Mutated structure viewer")
        print("  • output_files/altered_*.csv - Lost interactions")
        print("  • output_files/mutated_*.csv - Gained interactions")
        print("  • output_files/ - LigPlot analysis results")
        print("\n=== SCRIPT COMPLETED SUCCESSFULLY ===")
    else:
        print("\n!!! SCRIPT FAILED !!!")
        sys.exit(1)
