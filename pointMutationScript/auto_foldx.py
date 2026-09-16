#!/usr/bin/python3
"""
auto_foldx.py — FoldX runner for the X-Pro point-mutation stability bridge.

This is a copy of the original MutXplor auto_foldx.py, relocated to
/var/www/html/pointMutationScript (alongside automate_point_mutation.py)
instead of /var/www/html/mutXplorScripts. Only the two paths that pointed
at mutXplorScripts have been updated to point at pointMutationScript:

    - the rotabase.txt this needs is copied from pointMutationScript
    - the foldx binary itself is invoked from pointMutationScript

Everything else (the BuildModel workflow, the output format written to
outFiles/foldX.out) is unchanged, so it stays compatible with the existing
job-status polling logic.
"""

import sys
import os
import re
import shutil
import pandas as pd

def runFoldx(pdbFilename,mutationFilename,fchain,rndNum):
    os.chdir("/var/www/html/fileUpload/" + rndNum)
    shutil.copy2("../../pointMutationScript/rotabase.txt", "rotabase.txt")
    os.chmod("rotabase.txt", 0o777)
    exFolderFile = ['Foldxout', 'individual_list.txt', 'Unrecognized_molecules.txt']
    b = [os.system(f"rm -rf {i}") for i in exFolderFile if os.path.exists(i) == True]
    pdb=open(pdbFilename,"r")
    pdbFile=pdb.readlines()
    for lines in pdbFile:
        if re.match("^ATOM",lines):
            firstRes=int(re.split('\s+',lines)[5])
            break
    with open(mutationFilename,"r") as mutFile:
        muts = mutFile.readlines()
    muts = [x.strip() for x in muts]
    newRes = [x[-1] for x in muts]
    oldRes = [x[0] for x in muts]
    pos = []
    for x in muts:
        position=int(''.join(re.findall('\d',x)))
        if position < firstRes:
            pos.append(0)
        else:
            pos.append(int(''.join(re.findall('\d',x))))
    os.makedirs("logFiles", exist_ok=True)
    os.system("mkdir Foldxout")
    fwrite = open('outFiles/foldX.out', "w")
    fwrite.write('Mutation,Stability,Score\n')
    logFile = open('logFiles/foldX.log', "w")
    for i in range(len(pos)):
        if pos[i] == 0:
            pass
        else:
            #Run FoldX in command-line
            mutFile = open("individual_list.txt","w")
            mutFile.write(oldRes[i]+fchain+str(pos[i])+newRes[i]+";\n")
            mutFile.close()
            res = os.system("/var/www/html/pointMutationScript/foldx --command=BuildModel --pdb="+pdbFilename+" --mutant-file=individual_list.txt --out-pdb=false --screen=false --output-dir=./Foldxout")
            if res == 0:
                file1 = open("Foldxout/Dif_"+pdbFilename.split(".")[0]+".fxout","r")
                lines = file1.readlines()[8:]
                file1.close()
                out = [line.split("\t") for line in lines]
                df = pd.DataFrame(out)
                df.columns = df.iloc[0]
                df = df.drop([0]).reset_index(drop=True)
                mutation_str = oldRes[i]+str(pos[i])+newRes[i]
                score = df['total energy'][0]
                if float(score) < 0:
                    stability = "Stabilizing"
                else:
                    stability = "Destabilizing"
                fwrite.write(f"{mutation_str},{stability},{score}\n")
                logFile.write(f"Mutation: {mutation_str}\nStability: {stability}\nScore: {score}\n\n")
                rmF = [os.remove(f"Foldxout/{i}") for i in os.listdir('Foldxout')]
            else:
                logFile.write(f"FoldX BuildModel failed (exit code {res}) for mutation "
                               f"{oldRes[i]+fchain+str(pos[i])+newRes[i]}.\n\n")
    fwrite.close()
    logFile.close()
    os.system("rm -rf Foldxout molecules outFile individual_list.txt Unrecognized_molecules.txt rotabase.txt")
# runFoldx('pdbfile.pdb','mutationList.txt','A','965')
