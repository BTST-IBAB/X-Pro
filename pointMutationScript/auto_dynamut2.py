#!usr/bin/python3
"""
auto_dynamut2.py — DynaMut2 runner for the X-Pro point-mutation stability
bridge.

Copy of the original MutXplor auto_dynamut2.py, relocated to
/var/www/html/pointMutationScript. It doesn't hard-code anything under
mutXplorScripts, so no internal paths needed to change -- only its location
on disk did.
"""

import re
import mechanicalsoup as msp
import sys
import pandas as pd
import os

def runDynamut(pdbFile,mutationFile,fChain,rndNum):
    os.chdir("/var/www/html/fileUpload/" + rndNum)
    try:
        chain = fChain
        file1 = open(mutationFile,"r")
        muts = file1.readlines()
        file1.close()
        muts = [x.strip() for x in muts]
        newRes = [x[-1] for x in muts]
        oldRes = [x[0] for x in muts]
        pos = []
        for x in muts:
            position = ''.join(re.findall('\d',x))
            pos.append(int(position))
        pdbfile=open(pdbFile,"r")
        pdblines=pdbfile.readlines()
        for lines in pdblines:
            if re.match("^ATOM",lines):
                firstRes=re.split('\s+',lines)[5]
                break
        missing = []
        for lines in pdblines:
            if re.match("REMARK 465",lines):
                missing.append(lines)
        chainA=[]
        for x in missing[7:]:
            if re.split('\s+',x)[3]==chain:
                chainA.append(int(re.split('\s+',x)[4]))
        common=[x for x in pos if x in chainA]
        #Removing missing residues from mutation file
        notCommon=[]
        for x in pos:
            if x in common or x < int(firstRes):
                notCommon.append(0)
            else:
                notCommon.append(x)
        muts=[]
        for i in range(len(pos)):
            if notCommon[i]==0:
                pass
            else:
                muts.append(oldRes[i]+str(notCommon[i])+newRes[i])
        outFile=open("outFiles/dynamut2.out","w")
        outFile.write("Mutation,Stability,ddG\n")
        outFile.close()
        os.makedirs("logFiles", exist_ok=True)
        logFile=open("logFiles/dynamut2.log","w")
        outFile=open("outFiles/dynamut2.out","a")
        proxies = {'https':'http://sthiyaga%40ibab%2Eac%2Ein:Ibab2025@proxy.ibab.ac.in:3128', 'http':'http://sthiyaga%40ibab%2Eac%2Ein:Ibab2025@proxy.ibab.ac.in:3128'}
        url1 = "https://biosig.lab.uq.edu.au/dynamut2/submit_prediction"
        br = msp.StatefulBrowser()
        br.session.proxies = proxies   # THIS IS THE SOLUTION!
        for i in range(len(muts)):
            #Submission of mutations to the web server
            br.open(url1)
            form = br.select_form(nr=0)
            # form.set("pdb_file_single",pdbFile)
            form["pdb_file_single"] = open(pdbFile, "rb")
            br["mutation_single"]=muts[i]
            br["chain_single"] = chain
            br.submit_selected()
            if "Error" in br.page.text:
                print("DynaMut2: Error in submission of a mutation. Check log file for details.")
                logFile.write(br.page.text.split("Submission Error")[2].split("arrow")[0])
            else:
                #Fetching results from the webpage
                res = br.get_url()
                while("Predicted Stability Change" not in br.page.text):
                    br.open(res)
                result=re.split("Predicted Stability Change",br.page.text)[1].split('\n\n')[0]
                score=result.split('\n')[1].split()[0]
                stability=result.split('\n')[2][1:-1]
                outFile.write(muts[i]+","+stability+","+score+"\n")
                logFile.write(f"Mutation: {muts[i]}\nStability: {stability}\nScore: {score}\n\n")
        logFile.close()
        outFile.close()
    except:
        print(str(sys.exc_info()[0])+" occured and DynaMut2 could not run")

# runDynamut("pdbfile.pdb", "mutationList.txt", "A", "775")
