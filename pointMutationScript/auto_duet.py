#!/usr/bin/python3
"""
auto_duet.py — DUET/mCSM/SDM runner for the X-Pro point-mutation stability
bridge.

Copy of the original MutXplor auto_duet.py, relocated to
/var/www/html/pointMutationScript. It doesn't hard-code anything under
mutXplorScripts, so no internal paths needed to change -- only its location
on disk did.
"""

import mechanicalsoup as msp
import sys
import os

def runDuet(pdbFile, mutation, fChain, rndNum):
    os.chdir("/var/www/html/fileUpload/" + rndNum)
    muts = open(mutation, "r").readlines()
    os.makedirs("logFiles", exist_ok=True)

    tools = ['mCSM', 'sdm', 'duet']
    for tool in tools:
        with open(f'outFiles/{tool}.out', "w") as fwrite:
            fwrite.write('Mutation,Stability,Score\n')

    # One log file per tool, since a single DUET submission yields three
    # separate results (DUET, mCSM, SDM) that each get their own row on the
    # status page -- keep their logs equally separate rather than one shared
    # "duet.log" for all three.
    log_handles = {tool: open(f'logFiles/{tool}.log', "w") for tool in tools}

    url1 = "http://biosig.unimelb.edu.au/duet/stability"
    proxies = {'https': 'http://chandreyeenandi2001%40gmail%2Ecom:Chandreyee%401234@proxy.ibab.ac.in:3128', 'http': 'http://chandreyeenandi2001%40gmail%2Ecom:Chandreyee%401234@proxy.ibab.ac.in:3128'}

    try:
        br = msp.StatefulBrowser()
        br.session.proxies = proxies

        for mut in muts:
            mut = mut.strip()
            br.open(url1)
            form = br.select_form(nr=1)
            # form.set("wild", pdbFile)
            form["wild"] = open(pdbFile, "rb")
            br["mutation"] = mut
            br["chain"] = fChain
            res = br.submit_selected()

            content = br.page.find('div', attrs={"class": "well"})
            if content:
                stability = content.findAll("font", attrs={"size": "4"})
                for c, s in enumerate(stability):
                    tool = tools[c]
                    stability_label = s.contents[1].contents[0]
                    score = s.contents[0].lstrip().rstrip().replace(' kcal/mol (', '')
                    with open(f'outFiles/{tool}.out', "a") as fapnd:
                        fapnd.write(mut + ',' + stability_label + ',' + score + '\n')
                    log_handles[tool].write(
                        f"Mutation: {mut}\nStability: {stability_label}\nScore: {score}\n\n"
                    )
            else:
                error = br.page.find("div", attrs={"class": "alert alert-error"})
                if error and "Error" in error.text:
                    # This single submission covers all three tools at once,
                    # so a submission-level error applies to all three logs.
                    for tool in tools:
                        log_handles[tool].write(f"Mutation: {mut}\n{error.text}\n\n")

        for handle in log_handles.values():
            handle.close()

    except Exception as e:
        print(f"{type(e).__name__} occurred and DUET could not run: {e}")
        for handle in log_handles.values():
            if not handle.closed:
                handle.close()

# runDuet("pdbfile.pdb", "mutationList.txt", "A", "775")
