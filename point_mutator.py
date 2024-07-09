# Modules Imported
from pymol import cmd
import subprocess
import os
import argparse
import shutil
import sys
import pandas as pd

os.environ["http_proxy"] = "http://proxy.ibab.ac.in:3128/"
os.environ["ftp_proxy"] = os.environ["http_proxy"]

# Function to do the point mutations in a pdb file using PyMol
def pdb_creator(file_name, chain_name, atom_number, mutation):
    cmd.reinitialize()
    # pymol.finish_launching()
    atom_selected = chain_name + '/' + atom_number + '/'

    # Load molecule
    cmd.load(f"{file_name}.pdb")
    cmd.select(file_name)
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.deselect()
    cmd.save("imported_pdb.png")

    # Image of un-mutated residue
    cmd.select(file_name)
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.select("mutant", atom_selected)
    cmd.zoom("mutant", "4")
    cmd.color("white", "mutant")
    cmd.show("licorice")
    cmd.deselect()
    cmd.save("not_mutated.png")

    # Saving pdb file
    cmd.select("hetera", "hetatm")
    cmd.alter("hetera", "type='ATOM'")
    cmd.alter("mutant", "type='HETATM'")
    cmd.save("altered.pdb")

    # Clearing Screen
    cmd.reinitialize()
    cmd.fetch(file_name)
    cmd.select("mutant", atom_selected)

    # Turn on the wizard
    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")

    # Apply mutation
    cmd.get_wizard().set_mode(mutation)
    cmd.get_wizard().do_select(atom_selected)
    nStates = cmd.count_states("mutation")
    cmd.get_wizard().apply()

    # Saving Mutated pdb
    cmd.select("hatm", "hetatm")
    cmd.alter("hatm", "type='ATOM'")
    cmd.select("new", atom_selected)
    cmd.alter("new", "type='HETATM'")
    cmd.save("mutated.pdb")

    # Image of Mutated Residue
    cmd.select("mutant", atom_selected)
    cmd.zoom("mutant", "4")
    cmd.color("white", "mutant")
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.deselect()
    cmd.save("mutated.png")
    cmd.save("mutated.pse")
    return nStates

    # Quitting PyMol
    # cmd.quit()

# Running the LigPlus Software to find the interactions
def ligplot_interface(chain_name, atom_number):
    wkdir = os.getcwd()
    ligPlotLoc = "/var/www/html/pointMutationScript/LigPlus"
    programLoc = f"{ligPlotLoc}/lib/exe_linux64/"
    ligplot_prm = f"{ligPlotLoc}/lib/params/ligplot.prm"
    
    def ligplot(filename):
        os.system(f"{programLoc}hbadd {filename} {ligPlotLoc}/components.cif -wkdir ./")
        os.system(f"{programLoc}hbplus -L -f hbplus.rc -h 2.90 -d 3.90 -N {filename} -wkdir ./")
        os.system(f"{programLoc}hbplus -L -f hbplus.rc -h 2.70 -d 3.35 {filename} -wkdir ./")
        os.system(f"{programLoc}ligplot {filename} {atom_number} {atom_number} {chain_name} -wkdir ./ -prm {ligplot_prm} -ctype 1 -no_abort")

        # Rename trash files
        files_to_rename = [_ for _ in os.listdir('.') if 'ligplot.' in _[0:8]]

        for i in files_to_rename:
            # print(i, filename)
            subprocess.check_call(['mv', i, filename.split(".")[0] + '_' + i])
            src_path = "./" + filename.split(".")[0] + '_' + i
            dst_path = "./output_files/" + filename.split(".")[0] + '_' + i
            shutil.move(src_path, dst_path)

    # Main function
    def main_run():
        # Get list of pdb files in the directory
        pdb_files = [_ for _ in os.listdir(wkdir) if
                     (_[-4:] == '.pdb') and ('ligplot' not in _)]
        print("The following .pdb files were found: ")
        print(pdb_files)

        # Running the LigPlus Software
        for pdb_file in pdb_files:
            os.chdir(wkdir)
            print("File being opened : ", pdb_file)
            ligplot(filename = pdb_file)

    main_run()

# Comparing the differences in the summary of the interactions of the two PDBs
def file_cmp():
    # reading files
    f1 = open("output_files/mutated_ligplot.sum", "r")
    f2 = open("output_files/altered_ligplot.sum", "r")
    f3 = open("output_files/ligplot_differences.txt", "w")

    f1_data = f1.readlines()
    f2_data = f2.readlines()
    syntax = f1_data[1:5]

    f3.write("Interactions Gained: \n")
    f3.writelines(syntax)

    key = 0

    df = pd.DataFrame(columns=["Interacting_residues", "M_M", "S_S", "M_S", "Non_bonded_Contact"])
    for line1 in f1_data[1:]:
        if line1 not in f2_data:
            key = 1
            f3.write(line1)
            df.loc[len(df.index)] = [i.strip() for i in line1.split("     ")]
    df.to_csv("output_files/mutated_ligplot.csv", index=False)
    df.drop(df.index , inplace=True)
    
    if key == 0:
        f3.write("None \n")

    f3.write("\n")
    f3.write("\n")

    key = 0
    f3.write("Interactions Lost: \n ")
    f3.writelines(syntax)

    for line2 in f2_data[1:]:
        if line2 not in f1_data:
            key = 1
            f3.write(line2)
            df.loc[len(df.index)] = [i.strip() for i in line2.split("     ")]
    df.to_csv("output_files/altered_ligplot.csv", index=False)
    if key == 0:
        f3.write("None \n")

    # closing files
    f1.close()
    f2.close()
    f3.close()
        
def point_mutator(file_name, chain_name, atom_number, mutation, randomNumber):
    os.chdir(os.getcwd()+f"/fileUpload/{randomNumber}/")
    os.makedirs("output_files", 0o777, exist_ok=True)
    rotamerC = pdb_creator(file_name, chain_name, atom_number, mutation)
    ligplot_interface(chain_name, atom_number)
    file_cmp()
    return rotamerC

pdbID = sys.argv[1]
chainID = sys.argv[2]
resiNum = sys.argv[3]
mut = sys.argv[4]
rndNum = sys.argv[5]
rotamerNumber = point_mutator(pdbID, chainID, resiNum, mut, rndNum)
print(rotamerNumber)
# point_mutator("1GAG", "A", "1200", "Cys")