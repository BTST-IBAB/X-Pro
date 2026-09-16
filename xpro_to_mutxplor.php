 <?php
/**
 * xpro_to_mutxplor.php
 *
 * Bridge between X-Pro and the stability predictors (FoldX, DynaMut2, DUET
 * -- the last of which also yields mCSM and SDM from the same submission).
 *
 * X-Pro already has everything these tools need (PDB/PDB ID, chain, residue
 * number, mutated residue), so this script re-uses that data directly. There
 * is no "which tools should run" choice any more -- this bridge always runs
 * exactly these five results (FoldX, DynaMut2, DUET, mCSM, SDM), so the
 * intermediate tool-selection page is gone.
 *
 *   GET (first click from X-ProResult.php, no "view" param):
 *     1. Locates the wild-type structure X-Pro already produced (altered.pdb).
 *     2. Sets up a fresh job directory under /var/www/html/fileUpload/.
 *     3. Copies the wild-type structure into it.
 *     4. Writes pending_mutation.txt ("<residueNumber>:<mutation>") -- this is
 *        the exact mechanism the automation script resolves itself: it looks
 *        up the wild-type residue at that position from the PDB itself, so
 *        we don't need to pass it again.
 *     5. Launches automate_point_mutation.py in the background. Unlike the
 *        full 9-tool MutXplor pipeline (mutXplorScripts/automate.py, driven
 *        by mutXplor.php -> mutXplorResult.php, both untouched by this file),
 *        this is a separate, much smaller script that only ever runs FoldX,
 *        DynaMut2 and DUET. It lives in /var/www/html/pointMutationScript,
 *        not /var/www/html/mutXplorScripts.
 *     6. Redirects to this same script with "view=1" and the new job's
 *        randNum, so a bookmark/back-button/page-reload never re-launches
 *        the job -- it only ever re-renders the status of the existing one.
 *
 *   GET (view=1&randNum=...):
 *     Renders the live status page: one row per result (FoldX, DynaMut2,
 *     DUET, mCSM, SDM), each filling in as that tool's output file appears
 *     under fileUpload/<randNum>/outFiles/. The page polls this same script
 *     (ajax_status=1) every couple of seconds to update itself.
 *
 *   GET (ajax_status=1&randNum=...):
 *     Returns the current per-tool status as JSON for the page above to
 *     poll. No other parameters are required for this branch.
 */

error_reporting(E_ALL);
ini_set('display_errors', 0);
ini_set('log_errors', 1);

function bridge_fail($msg) {
    error_log("xpro_to_mutxplor: " . $msg);
    echo "<!DOCTYPE html><html><head><title>MutXplor</title>";
    echo "<style>";
    echo "body { background-image: url('background.jpg'); background-repeat: no-repeat; background-attachment: fixed; background-size: cover; background-position: center; font-family: 'Verdana', 'sans-serif'; font-size: 12px; color: #222; margin: 0; padding: 0; }";
    echo ".bridge-fail-card { max-width: 480px; margin: 80px auto; background: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); padding: 30px 35px; text-align: center; }";
    echo ".bridge-fail-card h3 { color: #175161; margin-top: 0; }";
    echo ".bridge-fail-card .errorMessage { color: #721c24; background-color: #f8d7da; border: 1px solid #f5c6cb; border-radius: 0.25rem; padding: 0.75rem 1.25rem; text-align: left; }";
    echo ".bridge-fail-card a { color: #175161; }";
    echo "</style>";
    echo "</head>";
    echo "<body>";
    echo "<div class='bridge-fail-card'>";
    echo "<h3>Unable to start pathogenicity prediction</h3>";
    echo "<p class='errorMessage'>" . htmlspecialchars($msg) . "</p>";
    echo "<p><a href='X-Pro.php'>Return to X-Pro</a></p>";
    echo "</div>";
    echo "</body></html>";
    exit;
}

// The result rows this bridge always runs. "file" is the outFiles/*.out
// basename each underlying tool script writes; a single DUET submission
// produces three of these files (duet.out, mCSM.out, sdm.out) at once.
// HOPE is different: it reports through a small JSON status file rather
// than a CSV result row, since automated submission/parsing for HOPE isn't
// implemented -- it only ever reports whether the service could be reached.
$RESULT_TOOLS = [
    'foldx'    => ['label' => 'FoldX',     'file' => 'foldX.out',         'type' => 'csv'],
    'dynamut2' => ['label' => 'DynaMut2',  'file' => 'dynamut2.out',      'type' => 'csv'],
    'duet'     => ['label' => 'DUET',      'file' => 'duet.out',          'type' => 'csv'],
    'mcsm'     => ['label' => 'mCSM',      'file' => 'mCSM.out',          'type' => 'csv'],
    'sdm'      => ['label' => 'SDM',       'file' => 'sdm.out',           'type' => 'csv'],
    'hope'     => ['label' => 'HOPE',      'file' => 'hope_status.json',  'type' => 'hope'],
];

$THREE_TO_ONE = [
    'ALA' => 'A', 'ARG' => 'R', 'ASN' => 'N', 'ASP' => 'D', 'CYS' => 'C',
    'GLN' => 'Q', 'GLU' => 'E', 'GLY' => 'G', 'HIS' => 'H', 'ILE' => 'I',
    'LEU' => 'L', 'LYS' => 'K', 'MET' => 'M', 'PHE' => 'F', 'PRO' => 'P',
    'SER' => 'S', 'THR' => 'T', 'TRP' => 'W', 'TYR' => 'Y', 'VAL' => 'V',
];

/** Read the last data row (if any) of an outFiles/<name>.out CSV. */
function xpro_read_tool_result($upload_dir, $filename) {
    $path = $upload_dir . 'outFiles/' . $filename;
    if (!is_file($path)) {
        return null; // file not written yet
    }
    $lines = @file($path, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
    if ($lines === false || count($lines) < 2) {
        return []; // header only -- no result row yet
    }
    $last = str_getcsv($lines[count($lines) - 1]);
    return [
        'mutation'  => $last[0] ?? '',
        'stability' => $last[1] ?? '',
        'score'     => $last[2] ?? '',
    ];
}

/** Read HOPE's outFiles/hope_status.json, e.g. {"state":"no_result","message":"..."}. */
function xpro_read_hope_status($upload_dir, $filename) {
    $path = $upload_dir . 'outFiles/' . $filename;
    if (!is_file($path)) {
        return null; // not written yet
    }
    $decoded = json_decode((string) file_get_contents($path), true);
    if (!is_array($decoded) || !isset($decoded['state'])) {
        return null;
    }
    return $decoded;
}

/** Compute {job_done, elapsed, tools:{key: {label,state,stability,score,message}}} for a job. */
function xpro_job_status($upload_dir, $result_tools) {
    $job_done = is_file($upload_dir . 'job_complete.txt');
    $started_file = $upload_dir . 'job_started.txt';
    $started = is_file($started_file) ? (int) trim((string) file_get_contents($started_file)) : time();
    $elapsed = max(0, time() - $started);

    $tools = [];
    foreach ($result_tools as $key => $meta) {
        if (($meta['type'] ?? 'csv') === 'hope') {
            $hope = xpro_read_hope_status($upload_dir, $meta['file']);
            if ($hope === null) {
                $tools[$key] = [
                    'label' => $meta['label'],
                    'state' => $job_done ? 'no_result' : 'running',
                    'stability' => null,
                    'score' => null,
                    'message' => $job_done ? 'HOPE server unavailable: the service could not be contacted.' : null,
                ];
            } else {
                $tools[$key] = [
                    'label' => $meta['label'],
                    'state' => $hope['state'] === 'done' ? 'done' : 'no_result',
                    'stability' => $hope['stability'] ?? null,
                    'score' => $hope['score'] ?? null,
                    'message' => $hope['message'] ?? null,
                ];
            }
            continue;
        }

        $data = xpro_read_tool_result($upload_dir, $meta['file']);
        if ($data === null || $data === []) {
            $tools[$key] = [
                'label' => $meta['label'],
                'state' => $job_done ? 'error' : 'running',
                'stability' => null,
                'score' => null,
                'message' => null,
            ];
        } else {
            $tools[$key] = [
                'label' => $meta['label'],
                'state' => 'done',
                'stability' => $data['stability'],
                'score' => $data['score'],
                'message' => null,
            ];
        }
    }
    return ['job_done' => $job_done, 'elapsed' => $elapsed, 'tools' => $tools];
}

/** Look up the residue name (3-letter) at chain:residueNumber in a PDB file. */
function xpro_residue_name($pdb_path, $chain, $residue_number) {
    if (!is_file($pdb_path)) {
        return null;
    }
    $handle = @fopen($pdb_path, 'r');
    if (!$handle) {
        return null;
    }
    $found = null;
    while (($line = fgets($handle)) !== false) {
        $record = trim(substr($line, 0, 6));
        if ($record !== 'ATOM' && $record !== 'HETATM') {
            continue;
        }
        if (trim(substr($line, 21, 1)) !== (string) $chain) {
            continue;
        }
        if (trim(substr($line, 22, 4)) !== (string) $residue_number) {
            continue;
        }
        $found = strtoupper(trim(substr($line, 17, 3)));
        break;
    }
    fclose($handle);
    return $found;
}

// ================================================================================
// BRANCH 1: ajax_status -- JSON polling endpoint, no other params required
// ================================================================================
if (isset($_GET['ajax_status']) && $_GET['ajax_status'] === '1') {
    header('Content-Type: application/json');
    $randNum = isset($_GET['randNum']) ? trim($_GET['randNum']) : '';
    if (!ctype_digit($randNum)) {
        http_response_code(400);
        echo json_encode(['error' => 'Invalid job reference.']);
        exit;
    }
    $upload_dir = "/var/www/html/fileUpload/" . $randNum . "/";
    if (!is_dir($upload_dir)) {
        http_response_code(404);
        echo json_encode(['error' => 'Job not found.']);
        exit;
    }
    echo json_encode(xpro_job_status($upload_dir, $RESULT_TOOLS));
    exit;
}

// ================================================================================
// BRANCH 1b: download_logs -- zips up this job's tool log files and streams them
// ================================================================================
if (isset($_GET['download_logs']) && $_GET['download_logs'] === '1') {
    $randNum = isset($_GET['randNum']) ? trim($_GET['randNum']) : '';
    if (!ctype_digit($randNum)) {
        bridge_fail("Invalid job reference.");
    }
    $upload_dir = "/var/www/html/fileUpload/" . $randNum . "/";
    if (!is_dir($upload_dir)) {
        bridge_fail("This pathogenicity prediction job could not be found. Please re-run it from X-Pro.");
    }
    if (!class_exists('ZipArchive')) {
        bridge_fail("Log download is unavailable on this server (the PHP zip extension is not installed).");
    }

    // Every per-tool log under logFiles/ (duet.log, dynamut2.log, hope.log
    // when HOPE failed, etc.), plus the orchestrator's own combined
    // stdout+stderr capture at the job root, if present.
    $log_files = [];
    $log_dir = $upload_dir . 'logFiles/';
    if (is_dir($log_dir)) {
        foreach (scandir($log_dir) as $entry) {
            $full = $log_dir . $entry;
            if (is_file($full)) {
                $log_files[$entry] = $full;
            }
        }
    }
    if (is_file($upload_dir . 'automate_debug.log')) {
        $log_files['automate_debug.log'] = $upload_dir . 'automate_debug.log';
    }

    if (empty($log_files)) {
        bridge_fail("No log files are available for this job yet.");
    }

    $zip_path = tempnam(sys_get_temp_dir(), 'xpro_logs_') . '.zip';
    $zip = new ZipArchive();
    if ($zip->open($zip_path, ZipArchive::CREATE | ZipArchive::OVERWRITE) !== true) {
        bridge_fail("Could not prepare the log archive. Please try again.");
    }
    foreach ($log_files as $name => $full_path) {
        $zip->addFile($full_path, $name);
    }
    $zip->close();

    $download_name = "xpro_tool_logs_" . $randNum . ".zip";
    header('Content-Type: application/zip');
    header('Content-Disposition: attachment; filename="' . $download_name . '"');
    header('Content-Length: ' . filesize($zip_path));
    readfile($zip_path);
    unlink($zip_path);
    exit;
}

// ================================================================================
// BRANCH 2: view=1 -- render the live status page for an already-launched job
// ================================================================================
if (isset($_GET['view']) && $_GET['view'] === '1') {
    $randNum       = isset($_GET['randNum']) ? trim($_GET['randNum']) : '';
    $pdbID         = isset($_GET['pdbID']) ? trim($_GET['pdbID']) : '';
    $chainID       = isset($_GET['chainID']) ? trim($_GET['chainID']) : '';
    $residueNumber = isset($_GET['residueNumber']) ? trim($_GET['residueNumber']) : '';
    $mutation      = isset($_GET['mutation']) ? strtoupper(trim($_GET['mutation'])) : '';
    $wtOneLetter   = isset($_GET['wt']) ? strtoupper(trim($_GET['wt'])) : '';

    if (!ctype_digit($randNum)) {
        bridge_fail("Invalid job reference.");
    }
    $upload_dir = "/var/www/html/fileUpload/" . $randNum . "/";
    if (!is_dir($upload_dir)) {
        bridge_fail("This pathogenicity prediction job could not be found. Please re-run it from X-Pro.");
    }

    $mutationOneLetter = $THREE_TO_ONE[$mutation] ?? $mutation;
    $mutationTitle = ($wtOneLetter !== '' ? $wtOneLetter : '?') . htmlspecialchars($residueNumber) . $mutationOneLetter;

    $initial = xpro_job_status($upload_dir, $RESULT_TOOLS);
    $totalTools = count($RESULT_TOOLS);
    $initialDoneCount = 0;
    foreach ($initial['tools'] as $t) {
        if ($t['state'] === 'done' || $t['state'] === 'error') {
            $initialDoneCount++;
        }
    }
    ?>
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>Stability Predictions — <?= htmlspecialchars($mutationTitle) ?></title>
        <style>
            body {
                background-image: url('background.jpg');
                background-repeat: no-repeat;
                background-attachment: fixed;
                background-size: cover;
                background-position: center;
                font-family: "Verdana", "sans-serif";
                margin: 0;
                padding: 0;
                color: #222;
            }
            .container {
                max-width: 640px;
                margin: 40px auto 60px;
                padding: 0;
                background: white;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                overflow: hidden;
            }
            .status-head {
                background: linear-gradient(135deg, #175161, #1f728a, #2789a5, #17809E);
                padding: 20px 30px;
            }
            .status-body { padding: 24px 30px 30px; }

            .mut-title { font-size: 26px; font-weight: bold; color: white; margin: 0 0 4px; }
            .mut-subtitle { font-size: 13px; color: #d9edf6; margin: 0; }
            .mut-subtitle a { color: white; text-decoration: underline; cursor: pointer; }

            .mut-desc { font-size: 14px; line-height: 1.6; color: #444; margin: 0 0 22px; }

            #summary-text { font-size: 14px; color: #333; margin: 0 0 10px; }
            #summary-text b { color: #111; }

            .overall-bar-track { background: #e9ecef; border-radius: 4px; height: 6px; overflow: hidden; margin-bottom: 28px; }
            .overall-bar-fill { background: #175161; height: 100%; width: 0%; transition: width .4s ease; }

            .tool-card { border: 1px solid #e3e6e8; border-radius: 8px; overflow: hidden; }
            .tool-row { padding: 16px 18px; border-bottom: 1px solid #eee; }
            .tool-row:last-child { border-bottom: none; }
            .tool-row-head { display: flex; justify-content: space-between; align-items: center; }
            .tool-name { font-weight: bold; font-size: 15px; color: #111; }
            .tool-status { font-size: 11px; font-weight: bold; letter-spacing: .04em; }
            .tool-status.done { color: #175161; }
            .tool-status.done.stabilizing { color: #1a7f37; }
            .tool-status.done.destabilizing { color: #b91c1c; }
            .tool-status.running { color: #8592a0; }
            .tool-status.error { color: #b91c1c; }
            .tool-status.no_result { color: #b45309; }

            .tool-bar-track { background: #e9ecef; border-radius: 3px; height: 5px; margin-top: 10px; overflow: hidden; position: relative; }
            .tool-bar-fill { position: absolute; top: 0; bottom: 0; background: #175161; }
            .tool-bar-fill.indeterminate { width: 40%; animation: xpro-indeterminate 1.4s infinite ease-in-out; }
            @keyframes xpro-indeterminate {
                0%   { left: -40%; }
                100% { left: 100%; }
            }
            .tool-bar-fill.complete { left: 0; width: 100%; }

            .tool-result-row { display: flex; justify-content: space-between; align-items: baseline; margin-top: 10px; font-size: 13px; }
            .tool-result-label { font-weight: bold; }
            .tool-result-label.destabilizing { color: #b91c1c; }
            .tool-result-label.stabilizing { color: #1a7f37; }
            .tool-result-label.unknown { color: #555; }
            .tool-result-value { font-family: "Courier New", monospace; color: #333; }

            .tool-substatus { color: #8a97a0; font-size: 13px; margin-top: 8px; }
            .tool-substatus.no-result-note { color: #175161; }

            .log-download-row { margin-top: 20px; text-align: center; }
            .log-download-link {
                display: inline-block;
                font-size: 13px;
                font-weight: bold;
                color: #175161;
                text-decoration: none;
                border: 1px solid #cbd8dc;
                border-radius: 6px;
                padding: 8px 16px;
            }
            .log-download-link:hover { background: #f0f7f9; }
        </style>
    </head>
    <body>
        <div class="container">
            <div class="status-head">
                <div class="mut-title"><?= htmlspecialchars($mutationTitle) ?></div>
                <div class="mut-subtitle">
                    PDB <?= htmlspecialchars($pdbID) ?> · chain <?= htmlspecialchars($chainID) ?> ·
                    residue <?= htmlspecialchars($residueNumber) ?> ·
                    from <a onclick="history.back(); return false;" href="#">the structural analysis</a>
                </div>
            </div>
            <div class="status-body">

            <div id="summary-text"></div>
            <div class="overall-bar-track"><div id="overall-bar-fill" class="overall-bar-fill"></div></div>

            <div class="tool-card">
                <?php foreach ($RESULT_TOOLS as $key => $meta): $t = $initial['tools'][$key]; ?>
                <div class="tool-row" id="tool-<?= htmlspecialchars($key) ?>">
                    <div class="tool-row-head">
                        <span class="tool-name"><?= htmlspecialchars($meta['label']) ?></span>
                        <span class="tool-status <?= htmlspecialchars($t['state']) ?>"></span>
                    </div>
                    <div class="tool-body"></div>
                </div>
                <?php endforeach; ?>
            </div>

            <div class="log-download-row">
                <a class="log-download-link"
                   href="xpro_to_mutxplor.php?download_logs=1&randNum=<?= urlencode($randNum) ?>">
                    ⬇ Download result files (.zip)
                </a>
            </div>
            </div>
        </div>

        <script>
            var RAND_NUM = <?= json_encode($randNum) ?>;
            var TOTAL_TOOLS = <?= json_encode($totalTools) ?>;
            var STATUS_URL = 'xpro_to_mutxplor.php?ajax_status=1&randNum=' + encodeURIComponent(RAND_NUM);

            function stateLabel(state) {
                if (state === 'done') return 'DONE';
                if (state === 'error') return 'ERROR';
                if (state === 'no_result') return 'NO RESULT';
                return 'RUNNING';
            }

            function stabilityClass(stability) {
                var s = (stability || '').toLowerCase();
                if (s.indexOf('destab') !== -1) return 'destabilizing';
                if (s.indexOf('stab') !== -1) return 'stabilizing';
                return 'unknown';
            }

            function renderTool(key, info) {
                var row = document.getElementById('tool-' + key);
                if (!row) return;
                var statusEl = row.querySelector('.tool-status');
                var body = row.querySelector('.tool-body');

                if (info.state === 'done') {
                    // The score is intentionally left off this page -- it's
                    // only included in the downloadable log .zip. Here we
                    // just show the tool name (left) and its
                    // Stabilizing/Destabilizing verdict (right).
                    var cls = stabilityClass(info.stability);
                    statusEl.textContent = info.stability || 'DONE';
                    statusEl.className = 'tool-status done ' + cls;
                    body.innerHTML = '';
                    return;
                }

                statusEl.textContent = stateLabel(info.state);
                statusEl.className = 'tool-status ' + info.state;

                if (info.state === 'error') {
                    body.innerHTML = '<div class="tool-substatus">No result returned.</div>';
                } else if (info.state === 'no_result') {
                    body.innerHTML = '<div class="tool-substatus no-result-note">' +
                        (info.message || 'No result returned.') + '</div>';
                } else {
                    body.innerHTML =
                        '<div class="tool-bar-track"><div class="tool-bar-fill indeterminate"></div></div>' +
                        '<div class="tool-substatus">Running</div>';
                }
            }

            function renderSummary(doneCount, elapsed) {
                document.getElementById('summary-text').innerHTML =
                    '<b>' + doneCount + '</b> of ' + TOTAL_TOOLS + ' tools have returned results · ' +
                    (TOTAL_TOOLS - doneCount) + ' still running · ' + elapsed + 's elapsed';
                document.getElementById('overall-bar-fill').style.width = (doneCount / TOTAL_TOOLS * 100) + '%';
            }

            function applyStatus(data) {
                var doneCount = 0;
                Object.keys(data.tools).forEach(function (key) {
                    renderTool(key, data.tools[key]);
                    if (data.tools[key].state === 'done' || data.tools[key].state === 'error') {
                        doneCount++;
                    }
                });
                renderSummary(doneCount, data.elapsed);
                return doneCount;
            }

            // Initial render from what the page already knows, before the first poll.
            applyStatus(<?= json_encode($initial) ?>);

            function poll() {
                fetch(STATUS_URL)
                    .then(function (r) { return r.json(); })
                    .then(function (data) {
                        if (data.error) { return; }
                        applyStatus(data);
                        if (!data.job_done) {
                            setTimeout(poll, 2500);
                        }
                    })
                    .catch(function () { setTimeout(poll, 4000); });
            }

            if (!<?= json_encode($initial['job_done']) ?>) {
                setTimeout(poll, 2500);
            }
        </script>
    </body>
    </html>
    <?php
    exit;
}

// ================================================================================
// BRANCH 3: initial click from X-ProResult.php -- validate, launch, redirect
// ================================================================================

// ── Read and validate the parameters passed from X-ProResult.php ──────────────
$pdbID         = isset($_REQUEST['pdbID']) ? trim($_REQUEST['pdbID']) : "";
$chainID       = isset($_REQUEST['chainID']) ? trim($_REQUEST['chainID']) : "";
$residueNumber = isset($_REQUEST['residueNumber']) ? trim($_REQUEST['residueNumber']) : "";
$mutation      = isset($_REQUEST['mutation']) ? strtoupper(trim($_REQUEST['mutation'])) : "";
$xproRandNum   = isset($_REQUEST['xproRandNum']) ? trim($_REQUEST['xproRandNum']) : "";

if ($pdbID === "" || $chainID === "" || $residueNumber === "" || $mutation === "" || $xproRandNum === "") {
    bridge_fail("Missing required parameters. Please re-run the mutation from X-Pro.");
}
if (!ctype_digit($residueNumber)) {
    bridge_fail("Invalid residue number.");
}

// X-Pro sends 3-letter amino acid codes (e.g. "ALA"), not single-letter --
// the pending-mutation resolver expects a single-letter target AA, so
// convert here.
if (!isset($THREE_TO_ONE[$mutation])) {
    bridge_fail("Invalid mutation residue code.");
}
$mutationOneLetter = $THREE_TO_ONE[$mutation];

if (!preg_match('/^[A-Za-z0-9_\-]{1,50}$/', $chainID)) {
    bridge_fail("Invalid chain ID.");
}
if (!preg_match('/^[A-Za-z0-9_\-]{1,50}$/', $pdbID)) {
    bridge_fail("Invalid PDB ID.");
}
if (!ctype_digit($xproRandNum)) {
    bridge_fail("Invalid session reference.");
}

// Confirm the source structure exists before allocating anything.
$xpro_upload_dir = "/var/www/html/fileUpload/" . $xproRandNum . "/";
$source_pdb = $xpro_upload_dir . "altered.pdb";
if (!file_exists($source_pdb)) {
    bridge_fail("The original structure from your X-Pro session could not be found. Please re-run the mutation in X-Pro first.");
}

// Resolve the wild-type residue letter at the mutation site, purely for
// display on the status page's title (e.g. "V15A").
$wtResidueName = xpro_residue_name($source_pdb, $chainID, $residueNumber);
$wtOneLetter = ($wtResidueName !== null && isset($THREE_TO_ONE[$wtResidueName])) ? $THREE_TO_ONE[$wtResidueName] : '';

// ── Allocate a fresh job directory ─────────────────────────────────────────────
// Same 1000-9999 range used elsewhere, checked for collisions since X-Pro's
// own numbering (1001-2000) overlaps this range.
$base_fileUpload = "/var/www/html/fileUpload/";
$randNum = 0;
$upload_dir = "";
for ($attempt = 0; $attempt < 20; $attempt++) {
    $candidate = rand(1000, 9999);
    $candidate_dir = $base_fileUpload . $candidate . "/";
    if (!file_exists($candidate_dir)) {
        $randNum = $candidate;
        $upload_dir = $candidate_dir;
        break;
    }
}
if ($randNum === 0) {
    bridge_fail("Could not allocate a job directory. Please try again.");
}

if (!mkdir($upload_dir, 0777, true)) {
    bridge_fail("Could not create working directory for the pathogenicity prediction job.");
}

// ── Copy the wild-type structure into the new job directory ───────────────────
$dest_pdb = $upload_dir . $pdbID . ".pdb";
if (!copy($source_pdb, $dest_pdb)) {
    bridge_fail("Could not copy the structure into the pathogenicity prediction job.");
}

// ── Write the pending-mutation marker automate_point_mutation.py resolves ─────
file_put_contents($upload_dir . "pending_mutation.txt", $residueNumber . ":" . $mutationOneLetter);

// ── Record the start time for the status page's elapsed-time display ──────────
file_put_contents($upload_dir . "job_started.txt", (string) time());

// ── Launch automate_point_mutation.py in the background ───────────────────────
// This is a separate script from mutXplorScripts/automate.py, kept in its own
// directory (/var/www/html/pointMutationScript) so it and the tool scripts it
// imports (auto_foldx.py, auto_dynamut2.py, auto_duet.py, also relocated
// there) never affect the full 9-tool MutXplor pipeline.
$automate_path = "/var/www/html/pointMutationScript/automate_point_mutation.py";
if (!file_exists($automate_path)) {
    bridge_fail("The stability-prediction automation script was not found on the server.");
}

$debug_log = $upload_dir . "automate_debug.log";
$command = "python3 " . escapeshellarg($automate_path) . " " .
           escapeshellarg((string)$randNum) . " " .
           escapeshellarg($chainID) . " " .
           escapeshellarg($pdbID) . " > " . escapeshellarg($debug_log) .
           " 2>&1 & echo $!;";

$output = [];
exec($command, $output);
$pid = isset($output[0]) ? (int) trim($output[0]) : 0;

if ($pid <= 0) {
    bridge_fail("Could not start the pathogenicity prediction job.");
}

error_log("xpro_to_mutxplor: Started stability-prediction job PID=$pid randNum=$randNum for " .
    "$pdbID $chainID $residueNumber$mutation ($mutationOneLetter) " .
    "(source X-Pro randNum=$xproRandNum)");

// ── Redirect to the live status page for this job ──────────────────────────────
$redirect_params = http_build_query([
    'view'          => '1',
    'randNum'       => $randNum,
    'pdbID'         => $pdbID,
    'chainID'       => $chainID,
    'residueNumber' => $residueNumber,
    'mutation'      => $mutation,
    'wt'            => $wtOneLetter,
]);
header("Location: xpro_to_mutxplor.php?" . $redirect_params);
exit;
