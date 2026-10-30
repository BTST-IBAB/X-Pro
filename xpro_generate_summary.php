<?php
/**
 * xpro_generate_summary.php
 * Called via AJAX from X-ProResult.php.
 * Runs xpro_summarise.py with the user's hint and disease context,
 * saves the output to output_files/summary.txt, and returns JSON.
 *
 * POST params:
 *   job_id   — the X-Pro job folder name under /var/www/html/fileUpload/
 *   pdb_id   — PDB ID (e.g. 4FDI)
 *   hint     — "auto" | "benign" | "deleterious"
 *   disease  — free-text disease context (optional)
 */

header('Content-Type: application/json');

// ── Validate inputs ────────────────────────────────────────────────────────
$job_id  = preg_replace('/[^a-zA-Z0-9_\-.]/', '', $_POST['job_id']  ?? '');
$pdb_id  = preg_replace('/[^a-zA-Z0-9]/',     '', $_POST['pdb_id']  ?? '');
$hint    = $_POST['hint']    ?? 'auto';
$disease = trim($_POST['disease'] ?? '');

if (!$job_id || !$pdb_id) {
    echo json_encode(['error' => 'Missing job_id or pdb_id.']);
    exit;
}

$upload_dir    = "/var/www/html/fileUpload/{$job_id}/";
$canonical_json = $upload_dir . "output_files/interaction_comparison_{$pdb_id}.json";
$summary_file  = $upload_dir . "output_files/summary.txt";
$script        = "/var/www/html/pointMutationScript/xpro_summarise.py";

if (!file_exists($canonical_json)) {
    echo json_encode(['error' => "Canonical interaction JSON not found for {$pdb_id}. Run X-Pro first."]);
    exit;
}

if (!file_exists($script)) {
    echo json_encode(['error' => "xpro_summarise.py not found at {$script}."]);
    exit;
}

// ── Build override hint message from user selection ────────────────────────
$hint_map = [
    'benign'      => 'The user believes this mutation is BENIGN. Use this as a strong prior when making your verdict, but only if the structural data supports it.',
    'deleterious' => 'The user believes this mutation is DELETERIOUS. Use this as a strong prior when making your verdict, but only if the structural data supports it.',
    'auto'        => '',  // auto-detect from data — no override
];
$hint_text = $hint_map[$hint] ?? '';

// ── Set up proxy environment ───────────────────────────────────────────────
$proxy_conf = '/etc/prokgs/proxy.conf';
$proxy_url  = '';
if (file_exists($proxy_conf)) {
    foreach (file($proxy_conf) as $line) {
        $line = trim($line);
        if (strpos($line, 'PROXY_URL=') === 0) {
            $proxy_url = substr($line, strlen('PROXY_URL='));
            break;
        }
    }
}
if (!$proxy_url) {
    $proxy_url = getenv('HTTPS_PROXY') ?: getenv('HTTP_PROXY') ?: '';
}

// ── Build command ──────────────────────────────────────────────────────────
$cmd_parts = [
    'env',
    'NO_PROXY=' . escapeshellarg('localhost,127.0.0.1'),
];
if ($proxy_url) {
    $cmd_parts[] = 'HTTPS_PROXY=' . escapeshellarg($proxy_url);
    $cmd_parts[] = 'HTTP_PROXY='  . escapeshellarg($proxy_url);
}
$cmd_parts[] = '/usr/bin/python3';
$cmd_parts[] = escapeshellarg($script);
$cmd_parts[] = escapeshellarg($canonical_json);
$cmd_parts[] = '--auto';
$cmd_parts[] = '--out ' . escapeshellarg($summary_file);
$cmd_parts[] = '--model phi4-mini-local';

if ($disease) {
    $cmd_parts[] = '--disease ' . escapeshellarg($disease);
}

// Pass user hint via environment variable read by xpro_summarise.py
if ($hint_text) {
    $cmd_parts[] = '--user-hint ' . escapeshellarg($hint_text);
}

$cmd = implode(' ', $cmd_parts) . ' 2>&1';

// ── Run ────────────────────────────────────────────────────────────────────
$output      = [];
$return_code = 0;
exec($cmd, $output, $return_code);

// ── Read result ────────────────────────────────────────────────────────────
if ($return_code !== 0 || !file_exists($summary_file)) {
    $log = implode("\n", array_slice($output, -10));
    echo json_encode([
        'error' => "Ollama generation failed (exit {$return_code}). Check that Ollama is running and phi4-mini-local is installed.\n\nLast log:\n{$log}"
    ]);
    exit;
}

$summary_text = trim(file_get_contents($summary_file));
if (!$summary_text) {
    echo json_encode(['error' => 'Ollama returned an empty response. Try again.']);
    exit;
}

$summary_url = str_replace('/var/www/html', '', $summary_file);

echo json_encode([
    'summary'     => $summary_text,
    'summary_url' => $summary_url,
    'hint_used'   => $hint,
]);
