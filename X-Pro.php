<?php
    session_start();
    include_once "header.php";

    ini_set('display_errors', 1);
    error_reporting(E_ALL);

    $chains_residue_IDs = "";
    $pdbID = "";
    $randNum = "";
    $showMutationForm = false;
    $errorMessage = "";

    // Create Random Folder 
    function makeDirectory($upload_dir){
        if (!file_exists($upload_dir)) { 
            return mkdir($upload_dir, 0777, true);
        }else{
            $randNum = rand(1001, 2000);
            $upload_dir = "./fileUpload/".$randNum."/";
            mkdir($upload_dir, 0777, true);
            return array($randNum, $upload_dir);
        }
    }

    // Handle form submission
    if ($_SERVER['REQUEST_METHOD'] === 'POST' && !empty($_POST)) {
        $randNum = rand(1001, 2000);
        $upload_dir = "/var/www/html/fileUpload/" . $randNum . "/";
        
        // Create Random Folder 
        if (!file_exists($upload_dir)) {
            mkdir($upload_dir, 0777, true);
        }

        if (!empty($_POST["PDBID"])) {
            $pdbID = strtoupper(trim($_POST["PDBID"]));
            
            $command = "python3 /var/www/html/pointMutationScript/download_pdbFile.py $pdbID /var/www/html/fileUpload $randNum 2>&1";
            $output = array();
            exec($command, $output, $return_var);
            
            error_log("PDB ID command output: " . implode("\n", $output));
            
            $chains_residue_IDs = end($output);
            
            if ($chains_residue_IDs != "404" && !empty($chains_residue_IDs)) {
                $showMutationForm = true;
                $_SESSION["pdbID"] = $pdbID;
                $_SESSION["randNum"] = $randNum;
                $_SESSION["upload_dir"] = $upload_dir;
                $_SESSION["chains_residue_IDs"] = $chains_residue_IDs;
            } else {
                $errorMessage = "Failed to retrieve PDB file. Please check the PDB ID and try again.";
            }
            
        } else if (!empty($_FILES["pdbFile"]['name'])) {
            $pdbFile = $_FILES["pdbFile"]['name'];
            $pdbID = pathinfo($pdbFile, PATHINFO_FILENAME);
            
            $filePath = $upload_dir . $pdbFile;
            
            if (move_uploaded_file($_FILES["pdbFile"]['tmp_name'], $filePath)) {
                $command = "python3 /var/www/html/pointMutationScript/download_pdbFile.py UPLOAD:$pdbFile /var/www/html/fileUpload $randNum 2>&1";
                $output = array();
                exec($command, $output, $return_var);
                
                error_log("Upload command output: " . implode("\n", $output));
                error_log("Return var: " . $return_var);
                
                $chains_residue_IDs = end($output);
                
                if ($chains_residue_IDs != "404" && !empty($chains_residue_IDs)) {
                    $showMutationForm = true;
                    $_SESSION["pdbID"] = $pdbID;
                    $_SESSION["randNum"] = $randNum;
                    $_SESSION["upload_dir"] = $upload_dir;
                    $_SESSION["chains_residue_IDs"] = $chains_residue_IDs;
                } else {
                    $errorMessage = "Failed to process PDB file. Please ensure it's a valid PDB file.";
                    $errorMessage .= "<br>Debug info: " . htmlspecialchars(implode("<br>", $output));
                }
            } else {
                $errorMessage = "Failed to upload file. Please try again.";
            }
        }
    }
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8" />
    <meta http-equiv="X-UA-Compatible" content="IE=edge"/>
    <meta name="viewport" content="width=device-width, initial-scale=1"/>
    <title>X-Pro | BioToolSuite</title>

    <!-- Local assets -->
    <link rel="stylesheet" href="bootstrap_4.0.css">
    <link rel="stylesheet" href="headerStyle.css">
    <link rel="stylesheet" href="select2_min.css">
    <link rel="stylesheet" href="fontawesome_5.8.1.css">

    <script src="jquery_3.2.1.js"></script>
    <script src="bootstrap_4.0.0.js"></script>
    <script src="select2_min.js"></script>

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
            padding: 20px 0;
        }

        .my-3 { margin-top: 1rem; margin-bottom: 1rem; }
        #mdTable { margin-top: 12px; }

        .container-two { 
            display:flex; 
            gap:28px; 
            align-items:flex-start; 
            margin: 20px auto;
            max-width: 1150px;
            padding: 0 20px;
        }
        .desc { 
            width:56%; 
            text-align:justify; 
            padding: 18px 22px; 
            font-size:14px; 
            line-height:1.8;
            background: rgba(255, 255, 255, 0.95);
            border-radius: 6px;
            border: 1px solid rgba(0,0,0,0.12);
        }
        .main-form {
            width:36%;
            padding: 18px;
            border-radius: 6px;
            background: rgba(255, 255, 255, 0.98);
            border: 1px solid rgba(0,0,0,0.12);
            box-shadow: 1px 2px 3px rgba(12,12,12,0.14);
        }

        .strucID, .aaResi, .mutant, select, input[type="text"] {
            padding: 8px; 
            border-radius: 5px; 
            border: 2px solid #ddd; 
            width: 100%; 
            font-size:13px;
            font-family: "Verdana", "sans-serif";
        }

        .file-upload-box {
            border: 2px dashed #175161;
            border-radius: 8px;
            padding: 20px;
            text-align: center;
            background: #f8fafb;
            cursor: pointer;
            transition: all 0.3s;
            margin-bottom: 15px;
        }
        .file-upload-box:hover {
            background: #e3f2fd;
            border-color: #1f728a;
        }
        .file-upload-box.has-file {
            border-color: #4caf50;
            background: #e8f5e9;
        }
        .upload-icon {
            font-size: 36px;
            color: #175161;
            margin-bottom: 10px;
        }
        .file-name-display {
            font-size: 14px;
            color: #333;
            font-weight: 600;
            margin-top: 10px;
        }

        .btn_smt { 
            font-size:12px; 
            padding:8px 12px; 
            border:1px solid #999; 
            background:#efefef; 
            border-radius:3px; 
            cursor:pointer;
            font-family: "Verdana", "sans-serif";
        }
        .analyze-btn { 
            background: linear-gradient(135deg,#175161,#1f728a); 
            color:#fff; 
            border:none; 
            padding:10px 14px; 
            border-radius:6px; 
            cursor:pointer;
            font-family: "Verdana", "sans-serif";
            font-size: 12px;
        }
        .analyze-btn:hover { 
            background: linear-gradient(135deg,#1f728a,#175161); 
        }

        .form-section {
            margin-top: 20px;
            padding-top: 15px;
            border-top: 2px solid #ddd;
            filter: blur(4px);
            opacity: 0.7;
            cursor: not-allowed;
            pointer-events: none;
            transition: filter 0.5s ease, opacity 0.45s ease;
        }
        .form-section.enabled {
            filter: blur(0px);
            opacity: 1;
            pointer-events: auto;
            cursor: default;
        }

        .form-group {
            margin-bottom: 15px;
        }

        .form-group label {
            display: block;
            margin-bottom: 5px;
            font-weight: 600;
            color: #0b3a4d;
            font-size: 12px;
        }

        .help-icon {
            display:inline-block;
            background:#feb22a;
            color:#fff;
            border-radius:50%;
            width:14px;
            height:14px;
            text-align:center;
            line-height:14px;
            font-size:10px;
            margin-left:6px;
            cursor: help;
        }

        .button-group {
            display: flex;
            gap: 8px;
            margin-top: 15px;
        }

        .button-group button {
            flex: 1;
        }

        .errorMessage {
            color: #721c24;
            background-color: #f8d7da;
            border: 1px solid #f5c6cb;
            border-radius: 0.25rem;
            padding: 0.75rem 1.25rem;
            margin-bottom: 1rem;
            font-size: 12px;
            max-width: 1150px;
            margin: 0 auto 20px auto;
        }

        .successBorder {
            border: 2px solid green !important;
            color: green !important;
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

        .page-title {
            text-align: center;
            font-size: 24px;
            font-weight: bold;
            margin: 20px 0;
            color: #0b3a4d;
        }

        @media (max-width:1000px) {
            .container-two { 
                flex-direction:column; 
                padding: 0 20px;
            }
            .desc, .main-form { width:100%; }
        }

        .divider {
            text-align: center;
            margin: 15px 0;
            position: relative;
            font-size: 12px;
            color: #666;
        }
        .divider::before {
            content: "";
            position: absolute;
            top: 50%;
            left: 0;
            right: 0;
            height: 1px;
            background-color: #ddd;
        }
        .divider span {
            background-color: rgba(255, 255, 255, 0.98);
            padding: 0 10px;
            position: relative;
        }

        .method-info {
            font-size: 10px;
            color: #666;
            margin-top: 3px;
            font-style: italic;
        }
        
        /* Progress bar styles */
        .progress-container {
            margin: 20px 0;
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
        }
        
        .progress-steps {
            list-style: none;
            padding: 0;
            margin: 10px 0 0 0;
            font-size: 11px;
        }
        
        .progress-steps li {
            padding: 3px 0;
        }
        
        .progress-steps li.completed {
            color: #4caf50;
        }
    </style>
</head>
<body>

    <?php if (!empty($errorMessage)): ?>
        <div class="errorMessage">
            <strong>Error:</strong> <?= htmlspecialchars($errorMessage) ?>
        </div>
    <?php endif; ?>

    <div class="page-title">X - Pro</div>

    <div class="container-two">
        <!-- Description -->
        <div class="desc">
            <p><strong>X-Pro</strong> is a computational tool designed to unravel the potential rotameric configurations that result from a point mutation within a protein structure. This tool harnesses the powerful functionalities of two readily accessible software programs: <strong>PyMOL and Modeller</strong>. It further employs <strong>LigPlot</strong> to explore changes in interactions between the mutated amino acid and its neighbouring residues. To understand the potential pathogenicity of the mutation, <strong>FoldX, DynaMut2, DUET, mCSM, SDM, and HOPE</strong> are used to assess mutation-induced changes in protein stability and structural properties, providing insights into the possible pathogenic effects of the mutation.</p>
            <p>If you find this tool to be useful and if the result lead to publication, please cite the following: <strong>Aman Vishwakarma, Anurag N, Nadimpally Sai Tharun Goud, Chandreyee Nandi, Sudha Srinivasan, S. Thiyagarajan.</strong></p>
            <p>This work is funded by <strong>Indian Council of Medical Research (ICMR)</strong>.</p>
            <div style="margin-top: 20px; padding-top: 15px; border-top: 1px solid #ecf0f1;">
            </div>
        </div>

        <!-- Form -->
        <div class="main-form">
            <form action="" enctype="multipart/form-data" method="post">
                <div style="font-weight:700; margin-bottom:10px;">
                    Step 1: Provide Protein Structure
                    <span class="help-icon" title="Provide PDB ID or upload PDB file">?</span>
                </div>
                
                <div class="form-group">
                    <label for="PDBID">PDB ID:</label>
                    <input type="text" name="PDBID" id="PDBID" class="strucID" placeholder="e.g., 1A8O" 
                           value="<?= htmlspecialchars($pdbID) ?>" <?= $showMutationForm ? 'disabled' : '' ?>>
                </div>
                
                <div class="divider">
                    <span>OR</span>
                </div>
                
                <div class="form-group">
                    <label>Upload PDB/PDBx File:</label>
                    <div class="file-upload-box" id="uploadBox" onclick="document.getElementById('file-input-pdb').click()">
                        <div class="upload-icon">
                            <i class="fa fa-cloud-upload"></i>
                        </div>
                        <div style="font-size:14px; font-weight:600; color:#175161;">Click to Upload .pdb/.cif File</div>
                        <div style="font-size:11px; color:#666; margin-top:5px;">or drag and drop</div>
                        <div class="file-name-display" id="file-name-pdb">No file chosen</div>
                        <input id="file-input-pdb" name="pdbFile" type="file" accept=".pdb,.cif" style="display:none;" 
                               onchange="showFileName(this)" <?= $showMutationForm ? 'disabled' : '' ?>>
                    </div>
                </div>
                
                <div class="button-group">
                    <button type="submit" class="analyze-btn">
                        <i class="fa fa-search"></i> Validate and prepare input
                    </button>
                    <button type="button" class="btn_smt" onclick="runExample()">
                        <i class="fa fa-play"></i> Run Example
                    </button>
                    <button type="button" class="btn_smt" onclick="window.location='X-Pro.php';">
                        <i class="fa fa-refresh"></i> Reset
                    </button>
                </div>

                <!-- Step 2: Always visible but blurred by default -->
                <div class="form-section <?= $showMutationForm ? 'enabled' : '' ?>" id="mutationSection">
                    <div style="font-weight:700; margin-bottom:10px;">
                        Step 2: Select Mutation Parameters
                        <span class="help-icon" title="Select chain, residue and mutation type">?</span>
                    </div>
                    
                    <div class="form-group">
                        <label for="chain-ID">Chain ID:</label>
                        <select class="form-control" name="ChainID" id="chain-ID" onchange="updateResidues()" <?= !$showMutationForm ? 'disabled' : '' ?>>
                            <option value="">Select Chain</option>
                        </select>
                    </div>

                    <div class="form-group">
                        <label for="residue">Residue Number:</label>
                        <select class="form-control" name="ResidueNum" id="residue" onchange="updateWildType()" <?= !$showMutationForm ? 'disabled' : '' ?>>
                            <option value="">Select Residue</option>
                        </select>
                    </div>

                    <div class="form-group">
                        <label>Wild-type Residue:</label>
                        <div style="padding: 8px; background: #f8f9fa; border-radius: 4px;">
                            <strong id="wildType_Label" style="font-size: 14px; color: #0b3a4d;"></strong>
                            <input type="hidden" id="wildType">
                        </div>
                    </div>

                    <div class="form-group">
                        <label for="mutation">Mutation:</label>
                        <select class="form-control" name="Mutation" id="mutation" <?= !$showMutationForm ? 'disabled' : '' ?>>
                            <option value="">Select Mutation</option>
                            <option value="Ala">Alanine (Ala)</option>
                            <option value="Arg">Arginine (Arg)</option>
                            <option value="Asn">Asparagine (Asn)</option>
                            <option value="Asp">Aspartic acid (Asp)</option>
                            <option value="Cys">Cysteine (Cys)</option>
                            <option value="Gln">Glutamine (Gln)</option>
                            <option value="Glu">Glutamic acid (Glu)</option>
                            <option value="Gly">Glycine (Gly)</option>
                            <option value="His">Histidine (His)</option>
                            <option value="Ile">Isoleucine (Ile)</option>
                            <option value="Leu">Leucine (Leu)</option>
                            <option value="Lys">Lysine (Lys)</option>
                            <option value="Met">Methionine (Met)</option>
                            <option value="Phe">Phenylalanine (Phe)</option>
                            <option value="Pro">Proline (Pro)</option>
                            <option value="Ser">Serine (Ser)</option>
                            <option value="Thr">Threonine (Thr)</option>
                            <option value="Trp">Tryptophan (Trp)</option>
                            <option value="Tyr">Tyrosine (Tyr)</option>
                            <option value="Val">Valine (Val)</option>
                        </select>
                    </div>

                    <!-- Mutagenesis Method Selection -->
                    <div class="form-group">
                        <label for="mutagenesis-method">Mutagenesis Method:</label>
                        <select class="form-control" name="MutagenesisMethod" id="mutagenesis-method" <?= !$showMutationForm ? 'disabled' : '' ?>>
                            <option value="">Select Method</option>
                            <option value="pymol">1. PyMOL</option>
                            <option value="modeller">2. Modeller (recommended)</option>
                        </select>
                    </div>
                    
                    <div class="button-group">
                        <button type="button" class="analyze-btn" onclick="submitToResults()" <?= !$showMutationForm ? 'disabled' : '' ?>>
                            <i class="fa fa-flask"></i> Analyze
                        </button>
                        <button type="button" class="btn_smt" onclick="window.location='X-Pro.php';">
                            <i class="fa fa-refresh"></i> Reset
                        </button>
                    </div>

                    <input type="hidden" id="storePDB_ID" value="<?= htmlspecialchars($pdbID) ?>">
                    <input type="hidden" id="folderID" value="<?= htmlspecialchars($randNum) ?>">
                </div>
            </form>
        </div>
    </div>

    <!-- footer -->
    <div class="footer">
        <div class="footerContent">
            <small>| © 2025, Institute of Bioinformatics and Applied Biotechnology (IBAB) |<br>
            Biotech park, GN Ramachandran Rd, Electronic city phase 1, Electronic city, Bengaluru - 560100<br>
            Phone: +91 8028528900</small>
        </div>
    </div>

<script>
    // Store chain data globally
    var chainsData = null;

    const aminoAcidsMap = {
        'ALA': 'Alanine (Ala)', 'CYS': 'Cysteine (Cys)', 'ASP': 'Aspartic acid (Asp)',
        'GLU': 'Glutamic acid (Glu)', 'PHE': 'Phenylalanine (Phe)', 'GLY': 'Glycine (Gly)',
        'HIS': 'Histidine (His)', 'ILE': 'Isoleucine (Ile)', 'LYS': 'Lysine (Lys)',
        'LEU': 'Leucine (Leu)', 'MET': 'Methionine (Met)', 'ASN': 'Asparagine (Asn)',
        'PRO': 'Proline (Pro)', 'GLN': 'Glutamine (Gln)', 'ARG': 'Arginine (Arg)',
        'SER': 'Serine (Ser)', 'THR': 'Threonine (Thr)', 'VAL': 'Valine (Val)',
        'TRP': 'Tryptophan (Trp)', 'TYR': 'Tyrosine (Tyr)'
    };

    function showFileName(input) {
        const fileNameBox = document.getElementById('file-name-pdb');
        const uploadBox = document.getElementById('uploadBox');
        if (input.files.length > 0) {
            fileNameBox.innerText = input.files[0].name;
            fileNameBox.style.color = "#333";
            uploadBox.classList.add('has-file');
            uploadBox.querySelector('.upload-icon').innerHTML = '<i class="fa fa-check-circle" style="color:#4caf50;"></i>';
        } else {
            fileNameBox.innerText = 'No file chosen';
            fileNameBox.style.color = "#666";
            uploadBox.classList.remove('has-file');
            uploadBox.querySelector('.upload-icon').innerHTML = '<i class="fa fa-cloud-upload"></i>';
        }
    }

    function runExample() {
        window.location.href = 'X-Pro_Example.php';
    }

    // Drag and drop functionality
    document.getElementById('uploadBox').addEventListener('dragover', function(e) {
        e.preventDefault();
        this.style.borderColor = '#4caf50';
        this.style.background = '#e8f5e9';
    });

    document.getElementById('uploadBox').addEventListener('dragleave', function(e) {
        this.style.borderColor = '#175161';
        this.style.background = '#f8fafb';
    });

    document.getElementById('uploadBox').addEventListener('drop', function(e) {
        e.preventDefault();
        var file = e.dataTransfer.files[0];
        if (file && (file.name.toLowerCase().endsWith('.pdb') || file.name.toLowerCase().endsWith('.cif'))) {
            document.getElementById('file-input-pdb').files = e.dataTransfer.files;
            showFileName(document.getElementById('file-input-pdb'));
        } else {
            alert('Please upload a .pdb or .cif file');
        }
        this.style.borderColor = '#175161';
        this.style.background = '#f8fafb';
    });

    function updateResidues() {
        const chainSelect = document.getElementById('chain-ID');
        const residueSelect = document.getElementById('residue');
        const selectedChain = chainSelect.value;
        
        residueSelect.innerHTML = '<option value="">Select Residue</option>';
        
        if (selectedChain && chainsData && chainsData[selectedChain]) {
            const residues = chainsData[selectedChain];
            Object.keys(residues).forEach(function(resNum) {
                const option = document.createElement('option');
                option.value = resNum;
                option.text = resNum + ' (' + residues[resNum] + ')';
                residueSelect.appendChild(option);
            });
            
            if (residueSelect.options.length > 1) {
                residueSelect.selectedIndex = 1;
                updateWildType();
            }
        }
    }

    function updateWildType() {
        const chainSelect = document.getElementById('chain-ID');
        const residueSelect = document.getElementById('residue');
        const wildTypeLabel = document.getElementById('wildType_Label');
        const wildTypeInput = document.getElementById('wildType');
        
        const selectedChain = chainSelect.value;
        const selectedResidue = residueSelect.value;
        
        if (selectedChain && selectedResidue && chainsData) {
            const wildTypeAA = chainsData[selectedChain][selectedResidue];
            const fullName = aminoAcidsMap[wildTypeAA] || wildTypeAA;
            wildTypeLabel.textContent = fullName;
            wildTypeInput.value = wildTypeAA;
        }
    }

    function submitToResults() {
        // Validate all fields are selected
        const chainID = document.getElementById('chain-ID').value;
        const residueNum = document.getElementById('residue').value;
        const mutation = document.getElementById('mutation').value;
        const mutagenesisMethod = document.getElementById('mutagenesis-method').value;
        const pdbID = document.getElementById('storePDB_ID').value;
        const folderID = document.getElementById('folderID').value;
        const wildType = document.getElementById('wildType').value;

        if (!chainID || !residueNum || !mutation || !mutagenesisMethod) {
            alert('Please select all mutation parameters before analyzing.');
            return;
        }

        if (!pdbID || !folderID || !wildType) {
            alert('Missing required parameters. Please submit the PDB again.');
            return;
        }

        const params = new URLSearchParams();
        params.append('store_PDB_ID', pdbID);
        params.append('ChainID', chainID);
        params.append('ResidueNum', residueNum);
        params.append('wildTypeResidue', wildType);
        params.append('Mutation', mutation);
        params.append('MutagenesisMethod', mutagenesisMethod);
        params.append('storeID', folderID);
        
        window.location.href = 'X-Pro-Progress.php?' + params.toString();
    }

    <?php if ($showMutationForm && !empty($chains_residue_IDs)): ?>
    // Initialize dropdowns
    (function() {
        console.log("Initializing X-Pro...");
        
        try {
            chainsData = <?= $chains_residue_IDs ?>;
            console.log("Chains data loaded:", chainsData);
            
            const chainSelect = document.getElementById('chain-ID');
            const chains = Object.keys(chainsData);
            
            console.log("Found chains:", chains);
            
            chains.forEach(function(chainId) {
                const option = document.createElement('option');
                option.value = chainId;
                option.text = 'Chain ' + chainId;
                chainSelect.appendChild(option);
            });
            
            console.log("Chain dropdown populated with", chainSelect.options.length - 1, "chains");
            
            if (chains.length > 0) {
                chainSelect.value = chains[0];
                updateResidues();
                console.log("Auto-selected first chain and populated residues");
            }
            
            console.log("Initialization complete!");
            
        } catch(e) {
            console.error("Error initializing:", e);
            alert("Error loading chain/residue data: " + e.message);
        }
    })();
    <?php endif; ?>
</script>
</body>
</html>
