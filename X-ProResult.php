<?php
    ini_set('display_errors', 1);
    ini_set('display_startup_errors', 1);
    error_reporting(E_ALL);
    
    session_start();
    include "header.php";
    require_once __DIR__ . "/comparison_helpers.php";

    function isValidAminoAcid($code) {
        $validAminoAcids = [
            'Ala', 'Arg', 'Asn', 'Asp', 'Cys',
            'Gln', 'Glu', 'Gly', 'His', 'Ile',
            'Leu', 'Lys', 'Met', 'Phe', 'Pro',
            'Ser', 'Thr', 'Trp', 'Tyr', 'Val'
        ];
        return in_array(ucfirst(strtolower($code)), $validAminoAcids);
    }

    function createTable($upload_dir, $pdbID) {
        // Try the expected names first (from LigPlot)
        $alteredFile = $upload_dir . "output_files/altered_" . $pdbID . "_ligplot.csv";
        $mutatedFile = $upload_dir . "output_files/mutated_" . $pdbID . "_ligplot.csv";
        
        error_log("X-ProResult: Looking for altered file: " . $alteredFile);
        error_log("X-ProResult: Looking for mutated file: " . $mutatedFile);
        error_log("X-ProResult: Altered exists: " . (file_exists($alteredFile) ? 'Yes' : 'No'));
        error_log("X-ProResult: Mutated exists: " . (file_exists($mutatedFile) ? 'Yes' : 'No'));
        
        // If not found, try alternative names
        if (!file_exists($alteredFile)) {
            $alteredFile = $upload_dir . "output_files/original_interactions_" . $pdbID . ".csv";
            error_log("X-ProResult: Trying alternative altered file: " . $alteredFile);
            error_log("X-ProResult: Alternative altered exists: " . (file_exists($alteredFile) ? 'Yes' : 'No'));
        }
        
        if (!file_exists($mutatedFile)) {
            $mutatedFile = $upload_dir . "output_files/mutated_interactions_" . $pdbID . ".csv";
            error_log("X-ProResult: Trying alternative mutated file: " . $mutatedFile);
            error_log("X-ProResult: Alternative mutated exists: " . (file_exists($mutatedFile) ? 'Yes' : 'No'));
        }
        
        // If still not found, try to find any CSV files in output directory
        if (!file_exists($alteredFile) || !file_exists($mutatedFile)) {
            $output_dir = $upload_dir . "output_files/";
            if (is_dir($output_dir)) {
                $csv_files = glob($output_dir . "*.csv");
                error_log("X-ProResult: All CSV files in output directory: " . implode(", ", array_map('basename', $csv_files)));
                
                // Try to find any file with 'altered' or 'original' in name
                foreach ($csv_files as $csv_file) {
                    $basename = basename($csv_file);
                    if ((strpos($basename, 'altered') !== false || strpos($basename, 'original') !== false) && !file_exists($alteredFile)) {
                        $alteredFile = $csv_file;
                        error_log("X-ProResult: Using found altered file: " . $alteredFile);
                    } elseif (strpos($basename, 'mutated') !== false && !file_exists($mutatedFile)) {
                        $mutatedFile = $csv_file;
                        error_log("X-ProResult: Using found mutated file: " . $mutatedFile);
                    }
                }
            }
        }
        
        // If still not found, try looking for any CSV files that might work
        if (!file_exists($alteredFile) || !file_exists($mutatedFile)) {
            $output_dir = $upload_dir . "output_files/";
            if (is_dir($output_dir)) {
                $csv_files = glob($output_dir . "*.csv");
                if (count($csv_files) >= 2) {
                    // Use the first two CSV files as altered and mutated
                    $alteredFile = $csv_files[0];
                    $mutatedFile = $csv_files[1];
                    error_log("X-ProResult: Using first two CSV files: " . basename($alteredFile) . " and " . basename($mutatedFile));
                }
            }
        }
        
        // Final check
        if (!file_exists($alteredFile) || !file_exists($mutatedFile)) {
            error_log("X-ProResult: ERROR - CSV files still not found after all attempts");
            error_log("X-ProResult: Altered file exists: " . (file_exists($alteredFile) ? 'Yes' : 'No'));
            error_log("X-ProResult: Mutated file exists: " . (file_exists($mutatedFile) ? 'Yes' : 'No'));
            return null;
        }
        
        error_log("X-ProResult: Successfully found both CSV files");
        error_log("X-ProResult: Using altered: " . basename($alteredFile));
        error_log("X-ProResult: Using mutated: " . basename($mutatedFile));
        
        $alteredRst = [];
        if (file_exists($alteredFile)) {
            foreach (file($alteredFile) as $line) {
                $alteredRst[] = str_getcsv($line, ",", "\"", "\\");
            }
            if (count($alteredRst) > 0) {
                array_shift($alteredRst); // Remove header row
            }
        }

        $mutatedRst = [];
        if (file_exists($mutatedFile)) {
            foreach (file($mutatedFile) as $line) {
                $mutatedRst[] = str_getcsv($line, ",", "\"", "\\");
            }
            if (count($mutatedRst) > 0) {
                array_shift($mutatedRst); // Remove header row
            }
        }

        return array($alteredRst, $mutatedRst);
    }

    function isScriptSuccessful($output) {
        foreach ($output as $line) {
            if (strpos($line, 'SCRIPT COMPLETED SUCCESSFULLY') !== false ||
                strpos($line, 'INTEGRATED MUTAGENESIS COMPLETED SUCCESSFULLY') !== false ||
                strpos($line, '=== ENHANCED POINT_MUTATOR COMPLETED SUCCESSFULLY ===') !== false ||
                strpos($line, '=== PYMOL MUTAGENESIS COMPLETED SUCCESSFULLY ===') !== false ||
                strpos($line, '=== MODELLER MUTAGENESIS COMPLETED SUCCESSFULLY ===') !== false ||
                strpos($line, 'Mutagenesis completed successfully') !== false ||
                strpos($line, '=== SCRIPT COMPLETED SUCCESSFULLY ===') !== false) {
                return true;
            }
        }
        return false;
    }

    function getRotamerCount($output) {
        foreach ($output as $line) {
            if (preg_match('/Rotamer count: (\d+)/', $line, $matches)) {
                return (int)$matches[1];
            }
            if (preg_match('/Possible rotamers found after mutation: (\d+)/', $line, $matches)) {
                return (int)$matches[1];
            }
            if (preg_match('/rotamer count: (\d+)/i', $line, $matches)) {
                return (int)$matches[1];
            }
        }
        return 12;
    }

    function getMethodUsed($output) {
        foreach ($output as $line) {
            if (strpos($line, 'Using Modeller for mutagenesis') !== false) return 'Modeller';
            if (strpos($line, 'Using PyMOL for mutagenesis') !== false) return 'PyMOL';
            if (strpos($line, 'Using Modeller for mutagenesis (user specified)') !== false) return 'Modeller';
            if (strpos($line, 'Using PyMOL for mutagenesis (user specified)') !== false) return 'PyMOL';
            if (strpos($line, 'Using PyMOL for mutagenesis (default)') !== false) return 'PyMOL';
            if (strpos($line, 'Using Modeller for mutagenesis (PyMOL not available)') !== false) return 'Modeller';
        }
        return 'PyMOL';
    }

    // Ensure original PDB file exists
    function ensureOriginalPDB($upload_dir, $pdbID) {
        $original_pdb = $upload_dir . strtolower($pdbID) . '.pdb';
        
        if (file_exists($original_pdb)) {
            return true;
        }
        
        $altered_pdb = $upload_dir . 'altered.pdb';
        if (file_exists($altered_pdb)) {
            if (copy($altered_pdb, $original_pdb)) {
                error_log("X-ProResult: Copied altered.pdb to " . $original_pdb);
                return true;
            }
        }
        
        return false;
    }

    // Initialize variables
    $pdbID = "";
    $chainID = "";
    $residueNumber = "";
    $wildRes = "";
    $mutation = "";
    $mutagenesisMethod = "";
    $randNum = "";
    $upload_dir = "";
    $rst = null;
    $comparison = null;   // canonical WT-vs-mutant interaction comparison (interaction_comparison.py output)
    $warningVariable = "";
    $abc = 0;
    $analysisComplete = false;
    $methodUsed = "PyMOL";
    $showProgress = false;
    $processingStarted = false;
    
    // Variables for images
    $before_image = null;
    $after_image = null;
    $annotated_mutant_image = null;
    $viewer_file_url = "";
    $viewer_file_exists = false;

    // Define path to xpro.py
    $xpro_script_path = "/var/www/html/pointMutationScript/xpro.py";
    
    if (!file_exists($xpro_script_path)) {
        error_log("X-ProResult: ERROR - xpro.py not found at: " . $xpro_script_path);
    }

    if (!empty($_GET['store_PDB_ID']) && !empty($_GET['ChainID']) && 
        !empty($_GET['ResidueNum']) && !empty($_GET['Mutation']) && 
        !empty($_GET['MutagenesisMethod']) && !empty($_GET['storeID'])){
        
        $randNum = $_GET['storeID'];
        $upload_dir = "/var/www/html/fileUpload/" . $randNum . "/";
        $pdbID = $_GET['store_PDB_ID'];
        $chainID = $_GET['ChainID'];
        $residueNumber = $_GET['ResidueNum'];
        $wildRes = $_GET['wildTypeResidue'];
        $mutagenesisMethod = strtolower($_GET['MutagenesisMethod']);
        $mutation = strtoupper($_GET['Mutation']);

        error_log("X-ProResult: Processing parameters - PDB: $pdbID, Chain: $chainID, Residue: $residueNumber, Mutation: $mutation, Method: $mutagenesisMethod");

        if (!file_exists($upload_dir)) {
            mkdir($upload_dir, 0777, true);
            error_log("X-ProResult: Created upload directory: $upload_dir");
        }

        $altered_csv = $upload_dir . "output_files/altered_" . $pdbID . "_ligplot.csv";
        $mutated_csv = $upload_dir . "output_files/mutated_" . $pdbID . "_ligplot.csv";
        $mutated_pdb = $upload_dir . "mutated.pdb";
        $altered_pdb = $upload_dir . "altered.pdb";
        $interactive_viewer = $upload_dir . "output_files/interactive_3d_viewer.html";
        
        // Check if results already exist (using flexible file checking)
        $altered_csv_exists = file_exists($altered_csv);
        $mutated_csv_exists = file_exists($mutated_csv);
        
        // Also check for alternative names
        if (!$altered_csv_exists) {
            $altered_csv_exists = file_exists($upload_dir . "output_files/original_interactions_" . $pdbID . ".csv");
        }
        if (!$mutated_csv_exists) {
            $mutated_csv_exists = file_exists($upload_dir . "output_files/mutated_interactions_" . $pdbID . ".csv");
        }
        
        $results_exist = $altered_csv_exists && $mutated_csv_exists &&
                        file_exists($mutated_pdb) && file_exists($altered_pdb);
        
        error_log("X-ProResult: Results exist: " . ($results_exist ? 'Yes' : 'No'));
        
        if ($results_exist) {
            // Results already exist
            ensureOriginalPDB($upload_dir, $pdbID);
            
            $methodUsed = $_SESSION["methodUsed"] ?? ucfirst($mutagenesisMethod);
            
            $rotamer_file = $upload_dir . "output_files/rotamer_count_" . $pdbID . ".txt";
            if (file_exists($rotamer_file)) {
                $abc = (int)file_get_contents($rotamer_file);
                error_log("X-ProResult: Rotamer count from file: $abc");
            } else {
                $abc = $_SESSION["possibleRotamer"] ?? 12;
                error_log("X-ProResult: Rotamer count from session: $abc");
            }
            
            $rst = createTable($upload_dir, $pdbID);
            $comparison = loadCanonicalComparison($upload_dir, $pdbID);
            if ($rst !== null || $comparison !== null) {
                $analysisComplete = true;
                error_log("X-ProResult: Table created successfully");
            } else {
                error_log("X-ProResult: Failed to create table");
                $warningVariable = "Warning: CSV files exist but could not be parsed properly.";
            }
            
            // Check for 3D viewer
            if (file_exists($interactive_viewer)) {
                $viewer_file_exists = true;
                $viewer_file_url = str_replace('/var/www/html', '', $interactive_viewer);
                error_log("X-ProResult: Found 3D viewer: $viewer_file_url");
            }
            
            $output_dir = $upload_dir . "output_files/";
            if (is_dir($output_dir)) {
                $before_images = glob($output_dir . "*altered*.png");
                $before_image = !empty($before_images) ? $before_images[0] : null;
                
                $after_images = glob($output_dir . "*mutated*.png");
                $after_images = array_filter($after_images, function($file) {
                    return strpos($file, 'altered') === false;
                });
                $after_image = !empty($after_images) ? reset($after_images) : null;
                
                $annotated_images = glob($output_dir . "*mutated*_ligplot_annotated.png");
                $annotated_mutant_image = !empty($annotated_images) ? $annotated_images[0] : null;
                
                error_log("X-ProResult: Before image: " . ($before_image ? basename($before_image) : 'None'));
                error_log("X-ProResult: After image: " . ($after_image ? basename($after_image) : 'None'));
            }

        } else {
            $processing_file = $upload_dir . "processing_started.txt";
            if (file_exists($processing_file)) {
                $processingStarted = true;
                $showProgress = true;
                $warningVariable = "Analysis is already in progress. Please wait or check back later.";
                error_log("X-ProResult: Analysis already in progress");
            } else {
                if (isValidAminoAcid($_GET['Mutation'])) {
                    error_log("X-ProResult: Starting new analysis");
                    file_put_contents($processing_file, date('Y-m-d H:i:s'));
                    
                    // Use system python directly
                    $python_path = "/usr/bin/python3";
                    
                    // Create wrapper script
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
                    error_log("X-ProResult: Created wrapper script: $wrapper_script");
                    
                    $output = array();
                    $return_var = 0;
                    
                    exec($wrapper_script, $output, $return_var);
                    
                    file_put_contents($upload_dir . "xpro_output.log", implode("\n", $output));
                    error_log("X-ProResult: Script completed with return code: $return_var");
                    
                    if (file_exists($processing_file)) {
                        unlink($processing_file);
                    }
                    
                    $mutated_pdb_exists = file_exists($upload_dir . "mutated.pdb");
                    $altered_pdb_exists = file_exists($upload_dir . "altered.pdb");
                    
                    $script_successful = isScriptSuccessful($output);
                    error_log("X-ProResult: Script successful: " . ($script_successful ? 'Yes' : 'No'));
                    
                    if ($script_successful || ($mutated_pdb_exists && $altered_pdb_exists)) {
                        
                        $methodUsed = getMethodUsed($output);
                        if (empty($methodUsed)) {
                            $methodUsed = ucfirst($mutagenesisMethod);
                        }
                        
                        $abc = getRotamerCount($output);
                        error_log("X-ProResult: Rotamer count: $abc");
                        
                        $rotamer_file = $upload_dir . "output_files/rotamer_count_" . $pdbID . ".txt";
                        if (!file_exists($upload_dir . "output_files")) {
                            mkdir($upload_dir . "output_files", 0777, true);
                        }
                        file_put_contents($rotamer_file, $abc);
                        
                        // Wait a moment for files to be written
                        sleep(2);
                        
                        $rst = createTable($upload_dir, $pdbID);
                        $comparison = loadCanonicalComparison($upload_dir, $pdbID);
                        
                        // Check for 3D viewer
                        if (file_exists($interactive_viewer)) {
                            $viewer_file_exists = true;
                            $viewer_file_url = str_replace('/var/www/html', '', $interactive_viewer);
                            error_log("X-ProResult: Found 3D viewer after analysis");
                        }
                        
                        if ($rst === null && $comparison === null) {
                            $created_files = glob($upload_dir . "*");
                            $output_files = glob($upload_dir . "output_files/*");
                            $all_files = array_merge($created_files, $output_files);
                            $warningVariable = "Analysis completed but CSV files not found. ";
                            $warningVariable .= "Created files: " . implode(", ", array_map('basename', $all_files));
                            error_log("X-ProResult: " . $warningVariable);
                        } else {
                            $analysisComplete = true;
                            error_log("X-ProResult: Analysis complete and CSV files found");
                            
                            ensureOriginalPDB($upload_dir, $pdbID);
                            
                            $output_dir = $upload_dir . "output_files/";
                            if (is_dir($output_dir)) {
                                $before_images = glob($output_dir . "*altered*.png");
                                $before_image = !empty($before_images) ? $before_images[0] : null;
                                
                                $after_images = glob($output_dir . "*mutated*.png");
                                $after_images = array_filter($after_images, function($file) {
                                    return strpos($file, 'altered') === false;
                                });
                                $after_image = !empty($after_images) ? reset($after_images) : null;
                                
                                $annotated_images = glob($output_dir . "*mutated*_ligplot_annotated.png");
                                $annotated_mutant_image = !empty($annotated_images) ? $annotated_images[0] : null;
                                
                                error_log("X-ProResult: Images found - Before: " . ($before_image ? 'Yes' : 'No') . ", After: " . ($after_image ? 'Yes' : 'No'));
                            }
                        }
                    } else {
                        $warningVariable = "Mutation analysis failed. Please check the parameters and try again. ";
                        $warningVariable .= "Return code: " . $return_var . "<br>";
                        
                        $error_lines = [];
                        foreach ($output as $line) {
                            if (strpos($line, 'ERROR') !== false || strpos($line, 'Error') !== false) {
                                $error_lines[] = $line;
                            }
                        }
                        
                        if (!empty($error_lines)) {
                            $warningVariable .= "Error details: " . implode("<br>", array_slice($error_lines, -5));
                        } else {
                            $warningVariable .= "Last output: " . implode("<br>", array_slice($output, -10));
                        }
                        error_log("X-ProResult: " . $warningVariable);
                    }

                    $_SESSION["Mutation"] = $mutation;
                    $_SESSION["possibleRotamer"] = $abc;
                    $_SESSION["methodUsed"] = $methodUsed;
                } else {
                    $warningVariable = $_GET['Mutation'] . " is not a valid amino acid code.";
                    error_log("X-ProResult: " . $warningVariable);
                }
            }
        }
    } else {
        $warningVariable = "PDB ID, Chain ID, Residue Number, Mutation, Method, or Folder ID is missing.";
        error_log("X-ProResult: " . $warningVariable);
    }
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8" />
    <meta http-equiv="X-UA-Compatible" content="IE=edge"/>
    <meta name="viewport" content="width=device-width, initial-scale=1"/>
    <title>X-Pro Results | BioToolSuite</title>

    <!-- Local assets -->
    <link rel="stylesheet" href="bootstrap_4.0.css">
    <link rel="stylesheet" href="headerStyle.css">
    <link rel="stylesheet" href="fontawesome_5.8.1.css">
    <link type="text/css" rel="stylesheet" href="css/dataTables.min.css">
    <link type="text/css" rel="stylesheet" href="css/jquery.dataTables.min.css">
    <link type="text/css" rel="stylesheet" href="css/buttons.dataTables.min.css">

    <script src="jquery_3.2.1.js"></script>
    <script src="bootstrap_4.0.0.js"></script>
    <script src="js/dataTables.min.js"></script>
    <script src="js/dataTables.buttons.min.js"></script>
    <script src="js/buttons.html5.min.js"></script>
    <script src="js/highcharts.js"></script>
    <script src="js/exporting.js"></script>
    <script src="js/export-data.js"></script>

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

        .page-title {
            text-align: center;
            font-size: 24px;
            font-weight: bold;
            margin: 20px 0;
            color: #0b3a4d;
        }

        .errorMessage {
            color: #721c24;
            background-color: #f8d7da;
            border: 1px solid #f5c6cb;
            border-radius: 0.25rem;
            padding: 0.75rem 1.25rem;
            margin: 1rem auto;
            max-width: 1150px;
        }

        .warningMessage {
            color: #856404;
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 0.25rem;
            padding: 0.75rem 1.25rem;
            margin: 1rem auto;
            max-width: 1150px;
        }

        .successMessage {
            color: #155724;
            background-color: #d4edda;
            border: 1px solid #c3e6cb;
            border-radius: 0.25rem;
            padding: 0.75rem 1.25rem;
            margin: 1rem auto;
            max-width: 1150px;
        }

        .headerTable {
            max-width: 1150px;
            margin: 20px auto;
            padding: 0 20px;
        }

        .headSummary {
            background: linear-gradient(135deg, #175161, #1f728a, #2789a5, #17809E);
            padding: 12px 15px;
            border-radius: 8px;
            margin-bottom: 1rem;
            border: none;
        }

        .headSummary strong {
            color: white;
        }

        .headerTable table {
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 6px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .headerTable table tr td {
            padding: 12px 15px;
            border-bottom: 1px solid #e0e0e0;
        }

        .headerTable table tr:last-child td {
            border-bottom: none;
        }

        .headerTable table tr td:first-child {
            font-weight: bold;
            width: 40%;
            background-color: #f8f9fa;
            color: #175161;
        }

        .method-badge {
            display: inline-block;
            background: #175161;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 11px;
            margin-left: 10px;
            font-weight: bold;
        }

        .content {
            max-width: 1150px;
            margin: 20px auto;
            padding: 0 20px;
            overflow: hidden;
        }

        .mutated_PDB, .not_mutated_PDB {
            width: 48%;
            float: left;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }

        .mutated_PDB {
            margin-right: 4%;
        }

        .LigPlotHead {
            position: relative;
            padding: 12px 15px;
            border-radius: 8px 8px 0 0;
        }

        .ligplotImageDownload {
            font-size: 16px;
            float: right;
            margin-top: 1px;
            color: white;
            text-decoration: none;
        }

        .ligplotImageDownload:hover {
            color: #f0f0f0;
        }

        .LigPlot_Image {
            width: 100%;
            height: auto;
            display: block;
            border-radius: 0 0 8px 8px;
        }

        .resTable {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 1rem;
            background: white;
        }

        .resTable th, .resTable td {
            border: 1px solid #dee2e6;
            padding: 8px 4px;
            text-align: center;
            vertical-align: middle;
        }

        .resTable thead th {
            background-color: #f8f9fa;
            font-weight: bold;
        }

        .button-container {
            text-align: center;
            margin: 30px auto;
            clear: both;
            max-width: 1150px;
            padding: 0 20px;
        }

        .downloadButton {
            color: black;
            background: #efefef;
            border-radius: 5px;
            padding: 8px 16px;
            text-decoration: none;
            display: inline-block;
            margin: 0 10px;
            border: 1px solid #ccc;
            cursor: pointer;
            transition: all 0.3s ease;
        }

        .downloadButton:hover {
            background: #e0e0e0;
            transform: translateY(-2px);
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }

        .viewer-main-container {
            max-width: 1150px;
            margin: 20px auto;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
        }
        
        .viewer-header {
            background: linear-gradient(135deg, #175161, #1f728a);
            color: white;
            padding: 15px 20px;
            font-size: 16px;
            font-weight: bold;
        }
        
        .viewer-frame {
            width: 100%;
            height: 600px;
            border: none;
            display: block;
        }

        .progress-container {
            max-width: 1150px;
            margin: 20px auto;
            padding: 15px;
            background: rgba(255, 255, 255, 0.95);
            border-radius: 8px;
            border: 1px solid rgba(0,0,0,0.12);
        }
        
        .progress-bar {
            width: 100%;
            background: #e0e0e0;
            border-radius: 4px;
            height: 20px;
            margin: 10px 0;
        }
        
        .progress-fill {
            height: 100%;
            background: linear-gradient(135deg,#175161,#1f728a);
            border-radius: 4px;
            transition: width 0.5s;
            width: 50%;
        }
        
        .auto-refresh-indicator {
            display: inline-block;
            padding: 2px 8px;
            background: #175161;
            color: white;
            border-radius: 12px;
            font-size: 10px;
            margin-left: 10px;
        }

        @media (max-width: 768px) {
            .mutated_PDB, .not_mutated_PDB {
                width: 100%;
                float: none;
                margin-right: 0;
            }
            
            .button-container .downloadButton {
                display: block;
                margin: 10px auto;
                width: 200px;
            }
        }
        
        .footer {
            background-color: #ffffff;
            margin-top: 30px;
            text-align: center;
            padding: 20px 0;
            width: 100%;
        }
        
        .footerContent {
            color: #ffffff;
            font-size: 11px;
            margin-top: 28px;
            padding: 10px 0;
            max-width: 1150px;
            margin: 0 auto;
        }
        
        .clearfix::after {
            content: "";
            clear: both;
            display: table;
        }
    </style>
</head>
<body>

<div class="page-title">X - Pro Results</div>

<?php if ($showProgress && !$analysisComplete): ?>
    <!-- Progress display for ongoing analysis -->
    <div class="progress-container" id="progressContainer">
        <h4>Analysis in Progress... <span class="auto-refresh-indicator">Auto-refreshing</span></h4>
        <div class="progress-bar">
            <div class="progress-fill" id="progressFill" style="width: 50%"></div>
        </div>
        <div id="progressText" style="margin-top: 10px; font-size: 12px;">
            Processing mutation analysis... This may take a few minutes.
        </div>
        <?php if (!empty($warningVariable)): ?>
        <div class="warningMessage" style="margin-top: 15px;">
            <?= htmlspecialchars($warningVariable) ?>
        </div>
        <?php endif; ?>
    </div>
    <script>
        setTimeout(function() {
            location.reload();
        }, 5000);
    </script>
    
<?php elseif ($analysisComplete && ($rst !== null || $comparison !== null)): ?>
    <!-- Display results when analysis is complete -->
    <?php if (!empty($warningVariable)): ?>
    <div class="warningMessage">
        <strong>Note:</strong> <?= htmlspecialchars($warningVariable); ?>
    </div>
    <?php endif; ?>

    <div class="headerTable">
        <div class="headSummary">
            <strong>Summary</strong>
        </div>
        <table>
            <tr>
                <td><strong>PDB ID:</strong></td>
                <td><?= htmlspecialchars($pdbID); ?></td>
            </tr>
            <tr>
                <td><strong>Chain ID:</strong></td>
                <td><?= htmlspecialchars($chainID); ?></td>
            </tr>
            <tr>
                <td><strong>Residue Number:</strong></td>
                <td><?= htmlspecialchars($residueNumber); ?></td>
            </tr>
            <tr>
                <td><strong>Wild-type Residue:</strong></td>
                <td><?= htmlspecialchars($wildRes); ?></td>
            </tr>
            <tr>
                <td><strong>Mutation:</strong></td>
                <td><?= htmlspecialchars($wildRes . $residueNumber . $mutation); ?></td>
            </tr>
            <tr>
                <td><strong>Mutagenesis Method:</strong></td>
                <td>
                    <span class="method-badge"><?= htmlspecialchars(ucfirst($methodUsed)); ?></span>
                </td>
            </tr>
            <tr>
                <td><strong>Possible rotamers found after mutation:</strong></td>
                <td><?= $abc; ?></td>
            </tr>
        </table>
    </div>
    
    <!-- 3D Viewer Section -->
    <?php if ($viewer_file_exists): ?>
    <div class="viewer-main-container">
        <div class="viewer-header">
            <span>🔬 3D Structure Viewer - Side by Side Comparison</span>
        </div>
        <iframe src="<?= htmlspecialchars($viewer_file_url) ?>" class="viewer-frame" allow="fullscreen"></iframe>
    </div>
    <?php endif; ?>

    <!-- Quick downloads, placed right after the 3D viewer and before the
         LigPlot interaction diagrams. These duplicate three of the buttons
         in the "Centered Buttons" row further down the page (which also
         keeps Download Interaction Summary and Reset), so they use their
         own element ids -- wired up in the script block at the bottom of
         the page -- rather than reusing the originals' ids. -->
    <div class="button-container" style="margin-top: 15px;">
        <?php if (file_exists($upload_dir . "mutated.pdb")): ?>
        <a class="downloadButton" id="mutatedPDBTop">Download Mutated PDB</a>
        <?php endif; ?>
        <?php if (file_exists($upload_dir . "mutated.cif")): ?>
        <a class="downloadButton" id="mutatedCIFTop">Download Mutated CIF</a>
        <?php endif; ?>
        <?php if (file_exists($upload_dir . "altered.pdb")): ?>
        <a class="downloadButton" id="originalPDBTop">Download Original PDB</a>
        <?php endif; ?>
    </div>

    <!-- LigPlot Interaction Diagrams -->
    <style>
        .xpro-ligplot-section { max-width: 1150px; margin: 20px auto; padding: 0 20px; }
        .xpro-ligplot-legend { display: flex; flex-wrap: wrap; gap: 9px; margin-bottom: 14px; }
        .xpro-ligplot-legend span { padding: 4px 8px; border-radius: 4px; font-weight: 800; }
        .xpro-ligplot-legend .gained { background: #d9edf6; color: #145c78; }
        .xpro-ligplot-legend .lost { background: #f8dfdc; color: #8e2d23; }
        .xpro-ligplot-legend .retained { background: #eeeccf; color: #60571e; }
        .xpro-ligplot-circle-note { display: flex; align-items: center; gap: 7px; margin-bottom: 14px; color: #175161; font-size: 13px; }
        .xpro-ligplot-circle-note .circle-swatch { width: 12px; height: 12px; border: 2px solid #d32f2f; border-radius: 50%; flex: none; }
        .xpro-ligplot-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 18px; }
        .xpro-ligplot-card { border: 1px solid #d7e1e5; border-radius: 7px; overflow: hidden; background: #fff; }
        .xpro-ligplot-head { padding: 11px 13px; background: #f5f8f9; font-weight: 800; display: flex; justify-content: space-between; align-items: flex-start; gap: 10px; color: #175161; }
        .xpro-ligplot-head a { color: #175161; }
        .xpro-ligplot-title { font-weight: 800; color: #175161; }
        .xpro-ligplot-subtitle { font-weight: 400; color: #6b7c83; font-size: 12px; margin-top: 2px; }
        .xpro-ligplot-image { display: block; width: 100%; height: auto; background: #fff; }
        .ligplot-highlight-layer { padding: 12px; }
        .ligplot-highlight-title { font-weight: 800; color: #175161; }
        .ligplot-highlight-list { list-style: none; padding: 0; margin: 7px 0; }
        .ligplot-highlight-item { display: flex; gap: 8px; border-top: 1px solid #edf1f2; padding: 7px 0; }
        .highlight-swatch { width: 12px; height: 12px; border-radius: 2px; flex: none; margin-top: 4px; background: #716626; }
        .highlight-gained .highlight-swatch { background: #176f91; }
        .highlight-lost .highlight-swatch { background: #bd3f32; }
        .highlight-meta { color: #667b83; }
        .ligplot-highlight-empty { color: #657980; }
        .retained-highlight-details summary { cursor: pointer; font-weight: 700; color: #175161; margin-top: 6px; }
        @media (max-width: 900px) { .xpro-ligplot-grid { grid-template-columns: 1fr; } }
    </style>
    <div class="xpro-ligplot-section">
        <div class="headSummary">
            <strong>Ligplot diagrams</strong>
        </div>
        <div class="xpro-ligplot-circle-note">
            Residues common to the Wild type and mutant are represented in red circles
        </div>
        <div class="xpro-ligplot-grid">
            <div class="xpro-ligplot-card">
                <div class="xpro-ligplot-head">
                    <div>
                        <div class="xpro-ligplot-title">Wild type LigPlot</div>
                        <div class="xpro-ligplot-subtitle"><?= xproEsc($chainID . ':' . $residueNumber . ' ' . strtoupper($wildRes)) ?></div>
                    </div>
                    <?php if ($before_image): ?>
                    <a href="<?= str_replace('/var/www/html', '', $before_image); ?>" download="<?= $pdbID ?>_Normal_Interaction.png">Download WT LigPlot</a>
                    <?php endif; ?>
                </div>
                <?php if ($before_image && file_exists($before_image)): ?>
                <img class="xpro-ligplot-image" src="<?= str_replace('/var/www/html', '', $before_image); ?>" alt="Wild-type native LigPlot diagram">
                <?php else: ?>
                <div style="padding: 40px; text-align: center; background: #f8f9fa; min-height: 300px; display: flex; align-items: center; justify-content: center;">
                    <em>Interaction diagram not available</em>
                </div>
                <?php endif; ?>
            </div>
            <div class="xpro-ligplot-card">
                <div class="xpro-ligplot-head">
                    <div>
                        <div class="xpro-ligplot-title">Mutant LigPlot</div>
                        <div class="xpro-ligplot-subtitle"><?= xproEsc($chainID . ':' . $residueNumber . ' ' . strtoupper($wildRes) . ' → ' . strtoupper($mutation)) ?> · native diagram</div>
                    </div>
                    <?php if ($after_image): ?>
                    <a href="<?= str_replace('/var/www/html', '', $after_image); ?>" download="<?= $pdbID ?>_Mutated_Interaction.png">Download Mutant LigPlot</a>
                    <?php endif; ?>
                </div>
                <?php if ($after_image && file_exists($after_image)): ?>
                <img class="xpro-ligplot-image" src="<?= str_replace('/var/www/html', '', $after_image); ?>" alt="Mutant native LigPlot diagram">
                <?php else: ?>
                <div style="padding: 40px; text-align: center; background: #f8f9fa; min-height: 300px; display: flex; align-items: center; justify-content: center;">
                    <em>Interaction diagram not available</em>
                </div>
                <?php endif; ?>
            </div>
            <div class="xpro-ligplot-card">
                <div class="xpro-ligplot-head">
                    <div>
                        <div class="xpro-ligplot-title">Mutant LigPlot</div>
                        <div class="xpro-ligplot-subtitle"><?= xproEsc($chainID . ':' . $residueNumber . ' ' . strtoupper($wildRes) . ' → ' . strtoupper($mutation)) ?> · common residues circled</div>
                    </div>
                    <?php if ($annotated_mutant_image): ?>
                    <a href="<?= str_replace('/var/www/html', '', $annotated_mutant_image); ?>" download="<?= $pdbID ?>_Mutated_Interaction_Annotated.png">Download annotated</a>
                    <?php endif; ?>
                </div>
                <?php if ($annotated_mutant_image && file_exists($annotated_mutant_image)): ?>
                <img class="xpro-ligplot-image" src="<?= str_replace('/var/www/html', '', $annotated_mutant_image); ?>" alt="Mutant LigPlot diagram with residues common to the WT and mutant interaction environments circled in red">
                <?php else: ?>
                <div style="padding: 40px; text-align: center; background: #f8f9fa; min-height: 300px; display: flex; align-items: center; justify-content: center;">
                    <em>Interaction diagram not available</em>
                </div>
                <?php endif; ?>
            </div>
        </div>
    </div>

    <!-- Interaction changes between WT and mutant (canonical comparison) -->
    <style>
        .xpro-cmp-toolbar { display: flex; flex-wrap: wrap; gap: 12px; margin: 12px 0; }
        .xpro-cmp-toolbar label { font-size: 12px; font-weight: 700; display: block; color: #175161; }
        .xpro-cmp-toolbar select, .xpro-cmp-toolbar input { padding: 7px; border: 1px solid #afc0c6; border-radius: 4px; }
        .xpro-cmp-table-wrap { overflow: auto; }
        .xpro-cmp-table { border-collapse: collapse; width: 100%; font-size: 12px; }
        .xpro-cmp-table th, .xpro-cmp-table td { border: 1px solid #d8e1e4; padding: 8px; vertical-align: top; }
        .xpro-cmp-table th { background: #175161; color: #fff; cursor: pointer; white-space: nowrap; }
        .xpro-cmp-table .number { text-align: right; white-space: nowrap; }
        .xpro-cmp-legend { display: flex; flex-wrap: wrap; gap: 9px; margin-bottom: 14px; }
        .xpro-cmp-legend span, .xpro-result-badge { padding: 4px 8px; border-radius: 4px; font-weight: 800; }
        .xpro-cmp-legend .gained, .xpro-result-gained { background: #d9edf6; color: #145c78; }
        .xpro-cmp-legend .lost, .xpro-result-lost { background: #f8dfdc; color: #8e2d23; }
        .xpro-cmp-legend .retained, .xpro-result-retained { background: #eeeccf; color: #60571e; }
        .xpro-result-cell { display: flex; align-items: center; gap: 6px; flex-wrap: wrap; }
        .xpro-result-gain-lost { background: #e5edf1; color: #244b5a; }
        .xpro-result-direction { font-size: 10px; font-weight: 800; }
        .contact-details summary { cursor: pointer; margin-top: 5px; }
        .contact-detail-scroll { overflow: auto; }
        .contact-detail-scroll table { font-size: 10px; width: 100%; border-collapse: collapse; }
        .contact-detail-scroll th, .contact-detail-scroll td { border: 1px solid #d8e1e4; padding: 4px; }
        .xpro-cmp-warning { background: #fff3e0; border-left: 4px solid #e67e22; padding: 12px 14px; border-radius: 6px; margin: 13px 0; }
    </style>
    <div style="max-width: 1150px; margin: 20px auto; padding: 0 20px;">
        <div class="headSummary">
            <strong>Interaction changes between WT and mutant</strong>
        </div>
        <?php if ($comparison === null): ?>
            <p style="text-align:center; padding: 20px;">Interaction comparison data is not available for this run.</p>
        <?php else: ?>
            <?php if (!empty($comparison['warnings'])): ?>
            <div class="xpro-cmp-warning">
                <strong>Scientific or external-tool limitations</strong>
                <ul>
                    <?php foreach ($comparison['warnings'] as $warning): ?>
                        <li><?= xproEsc($warning) ?></li>
                    <?php endforeach; ?>
                </ul>
            </div>
            <?php endif; ?>

            <div class="xpro-cmp-legend">
                <span class="gained">Gained</span>
                <span class="lost">Lost</span>
                <span class="retained">Retained</span>
            </div>

            <div class="xpro-cmp-toolbar">
                <div>
                    <label for="xproResultFilter">Result</label>
                    <select id="xproResultFilter">
                        <option value="">All</option>
                        <option value="Gain/Lost">Gain/Lost</option>
                        <option value="Retained">Retained</option>
                    </select>
                </div>
                <div>
                    <label for="xproTypeFilter">Interaction type</label>
                    <select id="xproTypeFilter">
                        <option value="">All</option>
                        <?php foreach (array_values(array_unique(array_column($comparison['records'], 'interaction_type'))) as $type): ?>
                            <option><?= xproEsc($type) ?></option>
                        <?php endforeach; ?>
                    </select>
                </div>
                <div>
                    <label for="xproResidueFilter">Residue</label>
                    <input id="xproResidueFilter" type="search" placeholder="Residue number" inputmode="numeric" pattern="[0-9-]*">
                </div>
            </div>

            <div class="xpro-cmp-table-wrap">
                <table class="xpro-cmp-table" id="xproComparisonTable">
                    <thead>
                        <tr>
                            <th>Residue pair</th>
                            <th>WT atom pair</th>
                            <th>Mutant atom pair</th>
                            <th>WT distance</th>
                            <th>Mutant distance</th>
                            <th>Interaction type</th>
                            <th>WT contacts</th>
                            <th>Mutant contacts</th>
                            <th>Δ contacts</th>
                            <th>Result</th>
                        </tr>
                    </thead>
                    <tbody>
                    <?php if (empty($comparison['records'])): ?>
                        <tr><td colspan="10" style="text-align:center;">No interactions were reported by the configured analysis engine.</td></tr>
                    <?php else: foreach ($comparison['records'] as $record):
                        $resultDirection = (string)($record['result'] ?? '');
                        $resultCategory = $resultDirection === 'Retained' ? 'Retained' : 'Gain/Lost';
                        $directionClass = strtolower($resultDirection);
                        $residueSearch = (string)($record['residue_number'] ?? '') . ' ' . (string)($record['partner_residue'] ?? '');
                    ?>
                        <tr data-result="<?= xproEsc($resultCategory) ?>" data-direction="<?= xproEsc($resultDirection) ?>"
                            data-type="<?= xproEsc($record['interaction_type']) ?>" data-residue="<?= xproEsc(strtolower($residueSearch)) ?>"
                            data-record-id="<?= xproEsc($record['id'] ?? '') ?>">
                            <td><?= xproResiduePairCell($record) ?></td>
                            <td><?= xproFormatAtomPair($record, 'wt') ?><?= xproContactDetails($record['wt_contacts'] ?? []) ?></td>
                            <td><?= xproFormatAtomPair($record, 'mutant') ?><?= xproContactDetails($record['mutant_contacts'] ?? []) ?></td>
                            <td class="number"><?= xproFormatDistance($record['wt_distance'] ?? null) ?></td>
                            <td class="number"><?= xproFormatDistance($record['mutant_distance'] ?? null) ?></td>
                            <td><?= xproEsc($record['interaction_type']) ?></td>
                            <td class="number"><?= (int)$record['wt_contact_count'] ?></td>
                            <td class="number"><?= (int)$record['mutant_contact_count'] ?></td>
                            <td class="number"><?= xproEsc(xproDelta($record['delta_contact_count'])) ?></td>
                            <td>
                                <div class="xpro-result-cell">
                                    <span class="xpro-result-badge xpro-result-<?= $resultCategory === 'Retained' ? 'retained' : 'gain-lost' ?>"><?= xproEsc($resultCategory) ?></span>
                                    <?php if ($resultCategory === 'Gain/Lost'): ?>
                                        <span class="xpro-result-direction xpro-result-<?= xproEsc($directionClass) ?>"><?= xproEsc($resultDirection) ?></span>
                                    <?php endif; ?>
                                </div>
                            </td>
                        </tr>
                    <?php endforeach; endif; ?>
                    </tbody>
                </table>
            </div>

            <div class="actions" style="margin-top: 14px;">
                <?php if (file_exists($upload_dir . "output_files/interaction_comparison_" . $pdbID . ".json")): ?>
                <a class="downloadButton" href="<?= str_replace('/var/www/html', '', $upload_dir . "output_files/interaction_comparison_" . $pdbID . ".json") ?>">Download interaction JSON</a>
                <?php endif; ?>
                <?php if (file_exists($upload_dir . "output_files/interaction_comparison_" . $pdbID . ".csv")): ?>
                <a class="downloadButton" href="<?= str_replace('/var/www/html', '', $upload_dir . "output_files/interaction_comparison_" . $pdbID . ".csv") ?>">Download interaction CSV</a>
                <?php endif; ?>
            </div>

            <script>
            (function(){
                const table = document.getElementById('xproComparisonTable');
                if (!table) return;
                const body = table.tBodies[0];
                const rf = document.getElementById('xproResultFilter');
                const tf = document.getElementById('xproTypeFilter');
                const q = document.getElementById('xproResidueFilter');
                function filterRows() {
                    Array.from(body.rows).forEach(row => {
                        if (!row.dataset.result) return;
                        row.style.display = (!rf.value || row.dataset.result === rf.value) &&
                                             (!tf.value || row.dataset.type === tf.value) &&
                                             (!q.value || row.dataset.residue.split(' ').includes(q.value.trim()))
                                             ? '' : 'none';
                    });
                }
                [rf, tf, q].forEach(control => control.addEventListener(control === q ? 'input' : 'change', filterRows));
                Array.from(table.tHead.rows[0].cells).forEach((head, column) => {
                    let asc = true;
                    head.addEventListener('click', () => {
                        const rows = Array.from(body.rows).filter(row => row.dataset.result);
                        rows.sort((a, b) => {
                            const av = a.cells[column].innerText.trim(), bv = b.cells[column].innerText.trim();
                            const an = parseFloat(av.replace(/[^0-9.+-]/g, '')), bn = parseFloat(bv.replace(/[^0-9.+-]/g, ''));
                            const d = !Number.isNaN(an) && !Number.isNaN(bn) ? an - bn : av.localeCompare(bv);
                            return asc ? d : -d;
                        });
                        rows.forEach(row => body.appendChild(row));
                        asc = !asc;
                    });
                });
            })();
            </script>
        <?php endif; ?>
    </div>

    <!-- Primary actions: start over, or move on to pathogenicity predictions -->
    <style>
        .primary-actions { display: flex; justify-content: center; align-items: center; gap: 18px; flex-wrap: wrap; margin: 30px auto; }
        .primary-actions a {
            display: inline-block;
            background: linear-gradient(135deg, #175161, #1f728a);
            color: #ffffff !important;
            font-weight: bold;
            font-size: 15px;
            padding: 14px 30px;
            border-radius: 8px;
            text-decoration: none;
            box-shadow: 0 3px 10px rgba(23, 81, 97, 0.35);
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }
        .primary-actions a:hover,
        .primary-actions a:visited,
        .primary-actions a:active {
            color: #ffffff !important;
        }
        .primary-actions a:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 18px rgba(23, 81, 97, 0.45);
        }
        @media (max-width: 600px) {
            .primary-actions { flex-direction: column; }
            .primary-actions a { width: 80%; text-align: center; }
        }
    </style>
    <div class="primary-actions">
        <a id="resetXpro" onclick="window.location='X-Pro.php';">Start new analysis</a>
        <a id="pathogenicPrediction"
           href="xpro_to_mutxplor.php?pdbID=<?= urlencode($pdbID) ?>&chainID=<?= urlencode($chainID) ?>&residueNumber=<?= urlencode($residueNumber) ?>&mutation=<?= urlencode($mutation) ?>&xproRandNum=<?= urlencode($randNum) ?>"
           target="_blank">
            For pathogenic predictions, click Here
        </a>
    </div>
    
<?php elseif (!$analysisComplete && !empty($warningVariable)): ?>
    <div class="errorMessage">
        <strong>Error:</strong> <?= htmlspecialchars($warningVariable); ?>
    </div>

    <div class="button-container">
        <a class="downloadButton" onclick="window.location='X-Pro.php';">Return to X-Pro</a>
    </div>
<?php endif; ?>

<script>
// Initialize download links
var upload_dir = <?= json_encode(str_replace('/var/www/html', '', $upload_dir)); ?>;
var pdbID = <?= json_encode($pdbID); ?>;
var wildRes = <?= json_encode($wildRes); ?>;
var mutRes = <?= json_encode($mutation); ?>;
var mutResN = <?= json_encode($residueNumber); ?>;

if (document.getElementById("mutatedPDB")) {
    document.getElementById("mutatedPDB").href = upload_dir + 'mutated.pdb';
    document.getElementById("mutatedPDB").download = pdbID + '_' + wildRes + mutResN + mutRes + '_mutated.pdb';
}

if (document.getElementById("mutatedCIF")) {
    document.getElementById("mutatedCIF").href = upload_dir + 'mutated.cif';
    document.getElementById("mutatedCIF").download = pdbID + '_' + wildRes + mutResN + mutRes + '_mutated.cif';
}

if (document.getElementById("originalPDB")) {
    document.getElementById("originalPDB").href = upload_dir + 'altered.pdb';
    document.getElementById("originalPDB").download = pdbID + '_original.pdb';
}

if (document.getElementById("mutatedPDBTop")) {
    document.getElementById("mutatedPDBTop").href = upload_dir + 'mutated.pdb';
    document.getElementById("mutatedPDBTop").download = pdbID + '_' + wildRes + mutResN + mutRes + '_mutated.pdb';
}

if (document.getElementById("mutatedCIFTop")) {
    document.getElementById("mutatedCIFTop").href = upload_dir + 'mutated.cif';
    document.getElementById("mutatedCIFTop").download = pdbID + '_' + wildRes + mutResN + mutRes + '_mutated.cif';
}

if (document.getElementById("originalPDBTop")) {
    document.getElementById("originalPDBTop").href = upload_dir + 'altered.pdb';
    document.getElementById("originalPDBTop").download = pdbID + '_original.pdb';
}

if (document.getElementById("interactionSummary")) {
    document.getElementById("interactionSummary").href = upload_dir + 'output_files/ligplot_differences.txt';
    document.getElementById("interactionSummary").download = pdbID + '_' + wildRes + mutResN + mutRes + '_interaction_differences.txt';
}

console.log("X-ProResult: Page loaded with upload_dir =", upload_dir);
console.log("X-ProResult: PDB ID =", pdbID);
</script>

</body>
</html>

<?php
    include "footer.php";
?>
