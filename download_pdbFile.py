import sys
import warnings
import requests
from Bio import PDB
from Bio.PDB import PDBParser

warnings.filterwarnings("ignore")
proxies = {'https':'http://proxy.ibab.ac.in:3128', 'http':'http://proxy.ibab.ac.in:3128'}
        
def get_chain_ids(pdb_file):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('PDB_structure', pdb_file)
    chains_residue_ids = {}
    for model in structure:
        for chain in model:
            if chain.id == chain.id:
                residues = [residue.get_id()[1] for residue in chain]
                chains_residue_ids[chain.id] = residues
    return chains_residue_ids

def download_pdb(pdbFileName, randNum):
    pdb_url = f'https://files.rcsb.org/download/{pdbFileName}.pdb'
    response = requests.get(pdb_url, proxies=proxies)
    if response.status_code == 200:
        with open(f'fileUpload/{randNum}/{pdbFileName}.pdb', 'wb') as file:
            file.write(response.content)
        print(get_chain_ids(f'fileUpload/{randNum}/{pdbFileName}.pdb'))
    else:
        print(response.status_code)

pdbID = sys.argv[1]
randNum = sys.argv[2]
download_pdb(pdbID, randNum)