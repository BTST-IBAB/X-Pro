<?php
/**
 * X-Pro-Progress.php
 *
 * Sits between X-Pro.php (the mutation-selection form) and X-ProResult.php
 * (the results page). X-Pro.php's "Analyze" button now redirects here
 * first, with the exact same query parameters it used to send straight to
 * X-ProResult.php.
 *
 * What this page does:
 *   1. If results already exist for this job (same file checks
 *      X-ProResult.php itself uses), skip straight to X-ProResult.php --
 *      there's nothing to wait for.
 *   2. If a job is already running (processing_started.txt present),
 *      just render the progress view again without relaunching anything.
 *   3. Otherwise, launch xpro.py in the background (unlike
 *      X-ProResult.php's own synchronous exec(), which blocks the browser
 *      with no feedback until the whole multi-minute pipeline finishes),
 *      and render a live 4-step progress page:
 *          1. Generating PyMOL/Modeller model
 *          2. Running HBPLUS and LigPlot
 *          3. Building canonical interaction comparison
 *          4. Analysis complete
 *      Each step lights up as xpro.py writes its corresponding
 *      progress_*.txt marker file (see write_progress_step() in xpro.py).
 *      Once step 4 appears, this page redirects the browser on to
 *      X-ProResult.php, which picks up the now-finished results through
 *      its existing "results already exist" branch -- no changes needed
 *      there at all.
 *
 * processing_started.txt is the same lock file X-ProResult.php's own
 * synchronous path already uses, so the two can never launch xpro.py twice
 * for the same job even if someone bypasses this page directly.
 */

// The four pipeline stages, in order, and the marker file xpro.py writes
// for each (see write_progress_step() calls in enhanced_point_mutator()).
$STEPS = [
    1 => ['label' => 'Generating PyMOL/Modeller model',          'file' => 'progress_1_model.txt'],
    2 => ['label' => 'Running HBPLUS and LigPlot',                'file' => 'progress_2_ligplot.txt'],
    3 => ['label' => 'Building canonical interaction comparison', 'file' => 'progress_3_comparison.txt'],
    4 => ['label' => 'Analysis complete',                         'file' => 'progress_4_complete.txt'],
];

/** Highest step number reached, plus whether the job failed, plus elapsed seconds. */
function xpro_progress_status($upload_dir, $steps) {
    $failed = is_file($upload_dir . 'progress_failed.txt');
    $current = 0;
    foreach ($steps as $num => $meta) {
        if (is_file($upload_dir . $meta['file'])) {
            $current = $num;
        }
    }
    $started_file = $upload_dir . 'processing_started.txt';
    if (is_file($started_file)) {
        $started = strtotime((string) file_get_contents($started_file));
        $started = $started !== false ? $started : filemtime($started_file);
    } else {
        $started = filemtime($upload_dir) ?: time();
    }
    $elapsed = max(0, time() - $started);
    return ['step' => $current, 'failed' => $failed, 'elapsed' => $elapsed];
}

// ================================================================================
// BRANCH: ajax_status -- JSON polling endpoint. This MUST be checked and
// exit before session_start()/header.php below ever run: header.php prints
// the full page chrome (nav bar etc.), and if that happens before this
// branch's json_encode() output, the polling fetch() on the client gets a
// response body of [header HTML][JSON] instead of pure JSON. r.json() then
// throws on every single poll, which silently retries in the page's
// catch() handler forever without ever updating the UI -- exactly the
// "stuck at step 1 forever, even though the job actually finished" symptom.
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
    echo json_encode(xpro_progress_status($upload_dir, $STEPS));
    exit;
}

// ================================================================================
// Everything below is the normal full-page render only.
// ================================================================================
session_start();
include "header.php";

function isValidAminoAcid($code) {
    $validAminoAcids = [
        'Ala', 'Arg', 'Asn', 'Asp', 'Cys',
        'Gln', 'Glu', 'Gly', 'His', 'Ile',
        'Leu', 'Lys', 'Met', 'Phe', 'Pro',
        'Ser', 'Thr', 'Trp', 'Tyr', 'Val'
    ];
    return in_array(ucfirst(strtolower($code)), $validAminoAcids);
}

function progress_fail($msg) {
    echo "<!DOCTYPE html><html><head><title>X-Pro</title></head>";
    echo "<body style='font-family:Verdana,sans-serif;text-align:center;margin-top:60px;'>";
    echo "<h3>Unable to start the analysis</h3>";
    echo "<p>" . htmlspecialchars($msg) . "</p>";
    echo "<p><a href='X-Pro.php'>Return to X-Pro</a></p>";
    echo "</body></html>";
    exit;
}

// ================================================================================
// Normal page load: validate params (same as X-ProResult.php expects)
// ================================================================================
if (empty($_GET['store_PDB_ID']) || empty($_GET['ChainID']) ||
    empty($_GET['ResidueNum']) || empty($_GET['Mutation']) ||
    empty($_GET['MutagenesisMethod']) || empty($_GET['storeID'])) {
    progress_fail("PDB ID, Chain ID, Residue Number, Mutation, Method, or Folder ID is missing. Please start again from X-Pro.");
}

$randNum           = $_GET['storeID'];
$pdbID             = $_GET['store_PDB_ID'];
$chainID           = $_GET['ChainID'];
$residueNumber     = $_GET['ResidueNum'];
$wildRes           = $_GET['wildTypeResidue'] ?? '';
$mutagenesisMethod = strtolower($_GET['MutagenesisMethod']);
$mutation          = strtoupper($_GET['Mutation']);

if (!ctype_digit((string) $randNum)) {
    progress_fail("Invalid session reference.");
}

$upload_dir = "/var/www/html/fileUpload/" . $randNum . "/";
if (!is_dir($upload_dir)) {
    mkdir($upload_dir, 0777, true);
}

// Same query string X-ProResult.php itself expects, reused for both the
// "already done, skip straight there" redirect and the final
// "step 4 reached" redirect at the end of the progress page's JS.
$result_query = http_build_query([
    'store_PDB_ID'      => $pdbID,
    'ChainID'           => $chainID,
    'ResidueNum'        => $residueNumber,
    'wildTypeResidue'   => $wildRes,
    'Mutation'          => $mutation,
    'MutagenesisMethod' => $mutagenesisMethod,
    'storeID'           => $randNum,
]);

// ── Already have results? Skip the progress page entirely. ────────────────────
$altered_csv = $upload_dir . "output_files/altered_" . $pdbID . "_ligplot.csv";
$mutated_csv = $upload_dir . "output_files/mutated_" . $pdbID . "_ligplot.csv";
$mutated_pdb = $upload_dir . "mutated.pdb";
$altered_pdb = $upload_dir . "altered.pdb";

$altered_csv_exists = file_exists($altered_csv) || file_exists($upload_dir . "output_files/original_interactions_" . $pdbID . ".csv");
$mutated_csv_exists = file_exists($mutated_csv) || file_exists($upload_dir . "output_files/mutated_interactions_" . $pdbID . ".csv");
$results_exist = $altered_csv_exists && $mutated_csv_exists && file_exists($mutated_pdb) && file_exists($altered_pdb);

if ($results_exist) {
    header("Location: X-ProResult.php?" . $result_query);
    exit;
}

// ── Already failed? Show that instead of relaunching or waiting forever. ──────
$already_failed = is_file($upload_dir . 'progress_failed.txt');

// ── Already running? Don't relaunch -- just render the progress view. ─────────
$already_running = is_file($upload_dir . 'processing_started.txt');

if (!$already_running && !$already_failed) {
    if (!isValidAminoAcid($mutation)) {
        progress_fail(htmlspecialchars($mutation) . " is not a valid amino acid code.");
    }

    file_put_contents($upload_dir . "processing_started.txt", date('Y-m-d H:i:s'));

    // Clear out any stale markers from a previous failed attempt in this
    // same folder before starting a fresh run.
    foreach ($STEPS as $meta) {
        @unlink($upload_dir . $meta['file']);
    }
    @unlink($upload_dir . 'progress_failed.txt');

    $xpro_script_path = "/var/www/html/pointMutationScript/xpro.py";
    if (!file_exists($xpro_script_path)) {
        progress_fail("The mutagenesis automation script was not found on the server.");
    }

    $python_path = "/usr/bin/python3";
    $wrapper_script = $upload_dir . "run_mutagenesis.sh";
    $wrapper_content = "#!/bin/bash\n";
    $wrapper_content .= "cd " . escapeshellarg($upload_dir) . "\n\n";
    $wrapper_content .= "echo \"=== RUNNING XPRO.PY ===\"\n";
    $wrapper_content .= $python_path . " " . escapeshellarg($xpro_script_path) . " " .
                       escapeshellarg($pdbID) . " " .
                       escapeshellarg($chainID) . " " .
                       escapeshellarg($residueNumber) . " " .
                       escapeshellarg($mutation) . " " .
                       escapeshellarg($mutagenesisMethod) . " " .
                       escapeshellarg($randNum) . " 2>&1\n";

    file_put_contents($wrapper_script, $wrapper_content);
    chmod($wrapper_script, 0755);

    // Launched detached, unlike X-ProResult.php's own blocking exec() --
    // this is what lets the browser get this progress page back
    // immediately instead of hanging for the whole pipeline's duration.
    $debug_log = $upload_dir . "xpro_progress_debug.log";
    $command = escapeshellarg($wrapper_script) . " > " . escapeshellarg($debug_log) . " 2>&1 & echo $!;";
    $output = [];
    exec($command, $output);
    $pid = isset($output[0]) ? (int) trim($output[0]) : 0;

    if ($pid <= 0) {
        @unlink($upload_dir . "processing_started.txt");
        progress_fail("Could not start the mutation analysis. Please try again.");
    }
}

$initial = xpro_progress_status($upload_dir, $STEPS);
$mutationTitle = htmlspecialchars(strtoupper($wildRes) . $residueNumber . $mutation);
?>
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1"/>
    <title>Running Analysis | X-Pro</title>
    <style>
        body {
            background-image: url('background.jpg');
            background-repeat: no-repeat;
            background-attachment: fixed;
            background-size: cover;
            background-position: center;
            font-family: "Verdana", "sans-serif";
            font-size: 12px;
            color: #222;
            margin: 0;
            padding: 0;
        }
        .xpp-container { max-width: 700px; margin: 40px auto; padding: 0 20px 60px; }
        .xpp-header {
            background: linear-gradient(135deg, #175161, #1f728a, #2789a5, #17809E);
            color: #fff;
            padding: 20px 24px;
            border-radius: 8px 8px 0 0;
        }
        .xpp-header .title { font-size: 26px; font-weight: bold; margin: 0 0 4px; }
        .xpp-header .subtitle { font-size: 13px; opacity: .9; }
        .xpp-body {
            background: #fff;
            border-radius: 0 0 8px 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.12);
            padding: 24px;
        }
        .xpp-lede { font-size: 13px; color: #555; margin: 0 0 20px; }

        .xpp-summary { font-size: 13px; color: #333; margin: 0 0 8px; }
        .xpp-summary b { color: #111; }
        .xpp-overall-bar-track { background: #e9ecef; border-radius: 4px; height: 6px; overflow: hidden; margin-bottom: 22px; }
        .xpp-overall-bar-fill { background: #175161; height: 100%; width: 0%; transition: width .4s ease; }

        .xpp-step-list { border: 1px solid #e3e6e8; border-radius: 8px; overflow: hidden; }
        .xpp-step { padding: 16px 18px; border-bottom: 1px solid #eee; }
        .xpp-step:last-child { border-bottom: none; }
        .xpp-step-head { display: flex; justify-content: space-between; align-items: center; }
        .xpp-step-label { font-size: 14px; font-weight: bold; color: #444; }
        .xpp-step-status { font-size: 11px; font-weight: bold; letter-spacing: .04em; color: #9aa8ad; }
        .xpp-step.done .xpp-step-label { color: #175161; }
        .xpp-step.done .xpp-step-status { color: #175161; }
        .xpp-step.active .xpp-step-label { color: #111; }
        .xpp-step.active .xpp-step-status { color: #8592a0; }
        .xpp-step.failed .xpp-step-label { color: #b91c1c; }
        .xpp-step.failed .xpp-step-status { color: #b91c1c; }

        .xpp-step-bar-track { background: #e9ecef; border-radius: 3px; height: 5px; margin-top: 10px; overflow: hidden; position: relative; }
        .xpp-step-bar-fill { position: absolute; top: 0; bottom: 0; left: 0; width: 0%; background: #175161; }
        .xpp-step.done .xpp-step-bar-fill { width: 100%; }
        .xpp-step.failed .xpp-step-bar-fill { width: 100%; background: #b91c1c; }
        .xpp-step.active .xpp-step-bar-fill {
            width: 40%;
            animation: xpp-indeterminate 1.4s infinite ease-in-out;
        }
        @keyframes xpp-indeterminate {
            0%   { left: -40%; }
            100% { left: 100%; }
        }
        .xpp-elapsed { text-align: center; font-size: 12px; color: #8a97a0; margin-top: 16px; }
        .xpp-error-box {
            margin-top: 18px;
            background: #fdecea;
            border: 1px solid #f5c6cb;
            color: #721c24;
            border-radius: 6px;
            padding: 12px 14px;
            font-size: 13px;
        }
    </style>
</head>
<body>
<div class="xpp-container">
    <div class="xpp-header">
        <div class="title"><?= $mutationTitle ?></div>
        <div class="subtitle">PDB <?= htmlspecialchars($pdbID) ?> · chain <?= htmlspecialchars($chainID) ?> · residue <?= htmlspecialchars($residueNumber) ?></div>
    </div>
    <div class="xpp-body">
        <p class="xpp-lede">Running the mutagenesis and interaction analysis pipeline. This page will move on to the results automatically once it's finished -- feel free to leave it open.</p>

        <div class="xpp-summary" id="xppSummary"></div>
        <div class="xpp-overall-bar-track"><div class="xpp-overall-bar-fill" id="xppOverallBarFill"></div></div>

        <div class="xpp-step-list" id="xppStepList">
            <?php foreach ($STEPS as $num => $meta): ?>
            <div class="xpp-step" id="xpp-step-<?= $num ?>" data-step="<?= $num ?>">
                <div class="xpp-step-head">
                    <span class="xpp-step-label"><?= htmlspecialchars($meta['label']) ?></span>
                    <span class="xpp-step-status"></span>
                </div>
                <div class="xpp-step-bar-track"><div class="xpp-step-bar-fill"></div></div>
            </div>
            <?php endforeach; ?>
        </div>

        <div class="xpp-elapsed" id="xppElapsed"></div>
        <div class="xpp-error-box" id="xppErrorBox" style="display:none;">
            The analysis could not be completed. Please
            <a href="X-Pro.php">return to X-Pro</a> and try again.
        </div>
    </div>
</div>

<script>
    var RAND_NUM = <?= json_encode($randNum) ?>;
    var TOTAL_STEPS = <?= count($STEPS) ?>;
    var RESULT_URL = 'X-ProResult.php?' + <?= json_encode($result_query) ?>;
    var STATUS_URL = 'X-Pro-Progress.php?ajax_status=1&randNum=' + encodeURIComponent(RAND_NUM);

    function applyStatus(data) {
        var doneCount = 0;
        for (var i = 1; i <= TOTAL_STEPS; i++) {
            var row = document.getElementById('xpp-step-' + i);
            if (!row) continue;
            row.classList.remove('done', 'active', 'failed');
            var statusEl = row.querySelector('.xpp-step-status');

            if (data.failed && i === data.step + 1) {
                row.classList.add('failed');
                statusEl.textContent = 'FAILED';
            } else if (i <= data.step) {
                row.classList.add('done');
                statusEl.textContent = 'DONE';
                doneCount++;
            } else if (i === data.step + 1 && !data.failed) {
                row.classList.add('active');
                statusEl.textContent = 'RUNNING';
            } else {
                statusEl.textContent = '';
            }
        }
        document.getElementById('xppSummary').innerHTML =
            '<b>' + doneCount + '</b> of ' + TOTAL_STEPS + ' steps complete';
        document.getElementById('xppOverallBarFill').style.width = (doneCount / TOTAL_STEPS * 100) + '%';
        document.getElementById('xppElapsed').textContent = data.elapsed + 's elapsed';
        if (data.failed) {
            document.getElementById('xppErrorBox').style.display = 'block';
        }
    }

    function poll() {
        fetch(STATUS_URL)
            .then(function (r) { return r.json(); })
            .then(function (data) {
                if (data.error) { return; }
                applyStatus(data);
                if (data.step >= TOTAL_STEPS) {
                    window.location.href = RESULT_URL;
                    return;
                }
                if (!data.failed) {
                    setTimeout(poll, 2000);
                }
            })
            .catch(function () { setTimeout(poll, 4000); });
    }

    applyStatus(<?= json_encode($initial) ?>);
    if (<?= json_encode($initial['step']) ?> >= TOTAL_STEPS) {
        window.location.href = RESULT_URL;
    } else if (!<?= json_encode($initial['failed']) ?>) {
        setTimeout(poll, 2000);
    }
</script>
</body>
</html>
