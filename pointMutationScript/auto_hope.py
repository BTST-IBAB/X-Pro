#!/usr/bin/python3
"""
auto_hope.py — HOPE ("Have (y)our Protein Explained") status check for the
X-Pro point-mutation bridge.

IMPORTANT / HONEST LIMITATION:
HOPE's real submission flow is a multi-step web form followed by an
asynchronous, email-delivered report -- unlike DUET/DynaMut2/FoldX, there is
no synchronous "submit and scrape the result page" path to reuse, and no
existing HOPE automation script was available to base this on. So this
module does NOT parse or return a real pathogenicity prediction. It writes
a "NO RESULT" status the PHP bridge displays, always with the message:

    "HOPE server unavailable: the service could not be contacted."

outFiles/hope_status.json is a small JSON object (not a CSV, unlike the
other tools here) since there is no numeric stability score to report:
    {"state": "no_result", "message": "HOPE server unavailable: the service could not be contacted."}

If HOPE submission/result-parsing is implemented later, this should instead
write {"state": "done", "stability": "...", "score": "..."} once a genuine
result is available, and the PHP side already knows how to render that.
"""

import json
import os

HOPE_URL = "https://www3.cmbi.umcn.nl/hope/"


def _write_status(state, message, stability=None, score=None):
    payload = {"state": state, "message": message}
    if stability is not None:
        payload["stability"] = stability
    if score is not None:
        payload["score"] = score
    with open("outFiles/hope_status.json", "w") as fh:
        json.dump(payload, fh)


def runHope(pdbFile, mutationFile, fChain, rndNum):
    os.chdir("/var/www/html/fileUpload/" + rndNum)
    os.makedirs("outFiles", exist_ok=True)
    os.makedirs("logFiles", exist_ok=True)

    try:
        import requests
        resp = requests.get(HOPE_URL, timeout=15)
        resp.raise_for_status()
    except Exception as e:
        with open("logFiles/hope.log", "w") as fh:
            fh.write(f"{type(e).__name__}: {e}\n")

    # Real HOPE submission + result-retrieval isn't implemented yet (see the
    # module docstring above), so this always reports the same "no result"
    # status regardless of whether the reachability check above succeeded.
    _write_status(
        "no_result",
        "HOPE server unavailable: the service could not be contacted.",
    )


# runHope('pdbfile.pdb', 'mutationList.txt', 'A', '965')
