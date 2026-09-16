<?php
/** Presentation helpers for the authoritative Python comparison export. */

function loadCanonicalComparison($uploadDir, $pdbId) {
    $path = $uploadDir . 'output_files/interaction_comparison_' . $pdbId . '.json';
    if (!is_file($path)) {
        error_log('X-Pro: canonical comparison missing: ' . $path);
        return null;
    }
    $decoded = json_decode(file_get_contents($path), true);
    if (!is_array($decoded) || !isset($decoded['records']) || !is_array($decoded['records'])) {
        error_log('X-Pro: invalid canonical comparison JSON: ' . $path);
        return null;
    }
    return $decoded;
}

function xproEsc($value) {
    return htmlspecialchars((string)$value, ENT_QUOTES, 'UTF-8');
}

function xproDisplayValue($value) {
    return ($value === null || $value === '') ? '—' : xproEsc($value);
}

function xproFormatDistance($value) {
    if ($value === null || $value === '') {
        return '—';
    }
    return number_format((float)$value, 2) . ' Å';
}

function xproFormatAtomPair($record, $prefix) {
    $atom1 = $record[$prefix . '_atom_1'] ?? null;
    $atom2 = $record[$prefix . '_atom_2'] ?? null;
    if (!$atom1 || !$atom2) {
        return '—';
    }
    return xproEsc($atom1) . '–' . xproEsc($atom2);
}

function xproContactDetails($contacts) {
    if (!is_array($contacts) || count($contacts) === 0) {
        return '';
    }
    $html = '<details class="contact-details"><summary>' . (count($contacts) === 1 ? 'Contact details' : 'All ' . count($contacts) . ' contacts') . '</summary>';
    $html .= '<div class="contact-detail-scroll"><table><thead><tr>';
    $html .= '<th>Atom pair</th><th>Distance</th><th>Roles/subtype</th>';
    $html .= '</tr></thead><tbody>';
    foreach ($contacts as $contact) {
        $roles = ($contact['atom_1_role'] ?? '') . '–' . ($contact['atom_2_role'] ?? '');
        if (!empty($contact['hydrogen_bond_subtype'])) {
            $roles .= ' (' . $contact['hydrogen_bond_subtype'] . ')';
        }
        if (!empty($contact['donor']) && !empty($contact['acceptor'])) {
            $roles .= '; donor ' . ($contact['donor']['atom'] ?? '?') . ' → acceptor ' . ($contact['acceptor']['atom'] ?? '?');
        }
        $html .= '<tr>';
        $html .= '<td>' . xproEsc($contact['atom_1'] ?? '') . '–' . xproEsc($contact['atom_2'] ?? '') . '</td>';
        $html .= '<td>' . xproFormatDistance($contact['distance'] ?? null) . '</td>';
        $html .= '<td>' . xproEsc($roles) . '</td>';
        $html .= '</tr>';
    }
    $html .= '</tbody></table></div></details>';
    return $html;
}

function xproResiduePairCell($record) {
    $wt = $record['wt_residue_pair'] ?? null;
    $mutant = $record['mutant_residue_pair'] ?? null;
    return '<div><strong>WT:</strong> ' . xproDisplayValue($wt) . '</div>' .
           '<div><strong>Mutant:</strong> ' . xproDisplayValue($mutant) . '</div>';
}

function xproDelta($value) {
    $number = (int)$value;
    return ($number > 0 ? '+' : '') . $number;
}

/**
 * Under-diagram summary for one LigPlot panel (state = 'wt' or 'mutant').
 *
 * The WT panel highlights interactions that disappear in the mutant (Lost);
 * the mutant panel highlights interactions that are new (Gained). Retained
 * interactions are relevant to both panels but are the bulk of the list, so
 * they're tucked into a collapsed <details> beneath the primary list rather
 * than repeating everything the main comparison table below already shows
 * in full.
 */
function xproLigPlotHighlights($comparison, $state) {
    $records = $comparison['records'] ?? [];
    if (!is_array($records) || count($records) === 0) {
        return '<div class="ligplot-highlight-layer"><p class="ligplot-highlight-empty">No canonical interaction records available.</p></div>';
    }

    $primaryResult = ($state === 'wt') ? 'Lost' : 'Gained';
    $swatchClass = ($primaryResult === 'Gained') ? 'highlight-gained' : 'highlight-lost';
    $distanceKey = ($state === 'wt') ? 'wt_distance' : 'mutant_distance';
    $pairKey = ($state === 'wt') ? 'wt_residue_pair' : 'mutant_residue_pair';

    $primary = [];
    $retained = [];
    foreach ($records as $record) {
        $result = (string)($record['result'] ?? '');
        if ($result === $primaryResult) {
            $primary[] = $record;
        } elseif ($result === 'Retained') {
            $retained[] = $record;
        }
    }

    $renderPrimaryItem = function ($record) use ($pairKey, $distanceKey, $swatchClass) {
        $pair = $record[$pairKey] ?? ($record['mutant_residue_pair'] ?? $record['wt_residue_pair'] ?? null);
        return '<li class="ligplot-highlight-item ' . $swatchClass . '">'
            . '<span class="highlight-swatch"></span>'
            . '<span class="highlight-meta"><strong>' . xproDisplayValue($pair) . '</strong> · '
            . xproEsc($record['interaction_type'] ?? '') . ' · '
            . xproFormatDistance($record[$distanceKey] ?? null) . '</span>'
            . '</li>';
    };
    $renderRetainedItem = function ($record) use ($pairKey, $distanceKey) {
        $pair = $record[$pairKey] ?? ($record['mutant_residue_pair'] ?? $record['wt_residue_pair'] ?? null);
        return '<li class="ligplot-highlight-item">'
            . '<span class="highlight-swatch"></span>'
            . '<span class="highlight-meta"><strong>' . xproDisplayValue($pair) . '</strong> · '
            . xproEsc($record['interaction_type'] ?? '') . ' · '
            . xproFormatDistance($record[$distanceKey] ?? null) . '</span>'
            . '</li>';
    };

    $html = '<div class="ligplot-highlight-layer">';
    $html .= '<div class="ligplot-highlight-title">' . xproEsc($primaryResult) . ' interactions</div>';
    if (count($primary) === 0) {
        $html .= '<p class="ligplot-highlight-empty">No ' . strtolower($primaryResult) . ' interactions.</p>';
    } else {
        $html .= '<ul class="ligplot-highlight-list">';
        foreach ($primary as $record) {
            $html .= $renderPrimaryItem($record);
        }
        $html .= '</ul>';
    }

    $html .= '<details class="retained-highlight-details"><summary>'
        . count($retained) . ' retained interaction' . (count($retained) === 1 ? '' : 's') . '</summary>';
    if (count($retained) === 0) {
        $html .= '<p class="ligplot-highlight-empty">No retained interactions.</p>';
    } else {
        $html .= '<ul class="ligplot-highlight-list">';
        foreach ($retained as $record) {
            $html .= $renderRetainedItem($record);
        }
        $html .= '</ul>';
    }
    $html .= '</details>';
    $html .= '</div>';

    return $html;
}

?>
