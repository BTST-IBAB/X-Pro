<?php
    // ini_set('display_errors', 1);
    // ini_set('display_startup_errors', 1);
    // error_reporting(E_ALL);

    error_reporting(0);

    $randNum = "X-Pro";
    $upload_dir = "ExampleFolder/$randNum/";
    $pdbID = "1BBT";
    $chainID = 1;
    $residueNumber = 60;
    $mutation = "Met";

    foreach (file($upload_dir."output_files/altered_ligplot.csv") as $line) {
        $alteredRst[] = str_getcsv($line);
    }
    array_shift($alteredRst);

    foreach (file($upload_dir."output_files/mutated_ligplot.csv") as $line) {
        $mutatedRst[] = str_getcsv($line);
    }
    array_shift($mutatedRst);

    $rst = array($alteredRst, $mutatedRst);

    include "header.php";
?>


<script src="js/dataTables.min.js"></script>
<script src="js/dataTables.buttons.min.js"></script>
<script src="js/buttons.html5.min.js"></script>
<script src="js/highcharts.js"></script>
<script src="js/exporting.js"></script>
<script src="js/export-data.js"></script>

<link type="text/css" rel="stylesheet" href="css/dataTables.min.css">
<link type="text/css" rel="stylesheet" href="css/jquery.dataTables.min.css">
<link type="text/css" rel="stylesheet" href="css/buttons.dataTables.min.css">
<link rel="stylesheet" href="css/fontawesome_5.8.1.css">

<style>
    .headerTable{
        padding: 0px 15px;
    }

    .headerTable table tr td:first-child strong{
        margin-right: 1rem;
    }

    .headerTable table tr td:first-child{
        padding-bottom: 0.4rem;
    }

    .content {
        overflow: auto;
    }

    .left, .right {
        padding: 1em;
    }

    .left, #alteredTable {
        float: left;
    }

    .right, #mutatedTable{
        float: right;
    }

    .headSummary {
        border-width: 2px;
        border-style: solid;
        background-image: linear-gradient(to right, #175161, #1f728a, #2789a5, #17809E);
        padding: 2px 6px;
        border-radius: 5px;
        margin-bottom: .5rem;
    }

    .imagedesign{
        border: 2px solid #000;
        border-radius: 5px;
        background: black;
    }

    .resTable thead th{
        text-align: center;
        vertical-align: middle;
    }

    .resTable tbody td{
        text-align-last: center;
    }

    .column{
        display: inline-block;
        /* float: left; */
        padding: 12px;
        width: 49.8%;
        background: white;
    }

    .downloadButton{
        color: black;
        background: #efefef;
        border-radius: 5px;
    }

    .LigPlotHead{
        position: relative;
    }

    .LigPlot_Image{
        margin: -6rem 0rem 0rem 0rem
    }

    .xProResultTable{
        margin-top: -5.5rem;
    }

    .ligplotImageDownload{
        font-size: 16px;
        float: right;
        margin-top: 1px;
        color: black;
    }
</style>

<!-- Body -->
<p></p>
<div class="mb-3">
    <h2 align="center" style="font-size: 18px;"><strong>X - Pro</strong></h2>
</div>

<?php
    if (!empty($rst)){
        ?>
            <div class="headerTable">
                <div class="headSummary">
                    <strong style="color: white;">Summary</strong>
                </div>
                <table>
                    <tr>
                        <td><strong> PDB ID: </strong></td>
                        <td><?= $pdbID; ?></td>
                    </tr>
                    <tr>
                        <td><strong> Chain ID: </strong></td>
                        <td><?= $chainID; ?></td>
                    </tr>
                    <tr>
                        <td><strong> Residue Number: </strong></td>
                        <td><?= $residueNumber; ?></td>
                    </tr>
                    <tr>
                        <td><strong> Mutation: </strong></td>
                        <td><?= ucfirst(strtolower($mutation)); ?></td>
                    </tr>
                    <tr>
                        <td><strong> Possible rotamers found after mutation: </strong></td>
                        <td>12</td>
                    </tr>
                </table>
            </div>
            <div class="content">
                <div class="mutated_PDB left">
                    <div class="headSummary">
                        <strong style="color: white;">Before Mutation</strong>
                    </div>
                    <div>
                        <img width="535" height="430" class="imagedesign" id="notMutatedImage">
                    </div>
                </div>

                <div class="not_mutated_PDB right">
                    <div class="headSummary">
                        <strong style="color: white;">After Mutation</strong>
                    </div>
                    <div>
                        <img width="535" height="430" class="imagedesign" id="mutatedImage">
                    </div>
                </div>
            </div>

            <div class="content">
                <div class="mutated_PDB left">
                    <div class="headSummary LigPlotHead">
                        <strong style="color: white;">Before Mutation</strong>
                        <a id="LigPlot_Image_UnMutated"><i class="fa ligplotImageDownload">&#xf019;</i></a>
                    </div>
                    <div>
                        <img width="535" height="730" class="LigPlot_Image" id="LigPlot_Unmutated">
                    </div>
                </div>

                <div class="not_mutated_PDB right">
                    <div class="headSummary LigPlotHead">
                        <strong style="color: white;">After Mutation</strong>
                        <a id="LigPlot_Image_Mutated"><i class="fa ligplotImageDownload">&#xf019;</i></a>
                    </div>
                    <div>
                        <img width="535" height="730" class="LigPlot_Image" id="LigPlot_mutated">
                    </div>
                </div>
            </div>

            <div class="xProResultTable">
                <?php
                    $tblID = array("alteredTable", "mutatedTable");
                    $table_main_header = array("Interactions Gained", "Interactions Lost");
                    $c = 0;
                    foreach($rst as $key){
                        ?>
                            <div class="column">
                                <div class="headSummary">
                                    <strong style="color: white;"><?= $table_main_header[$c]; ?></strong>
                                </div>
                                <table border="1" id="<?= $tblID[$c]; ?>" class="resTable table table-striped">
                                    <thead>
                                        <tr>
                                            <th rowspan="2" style="width: 12rem;">Interacting residue</th>
                                            <th colspan="3">Hydrogen bonds</th>
                                            <th rowspan="2" style="width: 7rem;">Non-bonded contact</th>
                                        </tr>
                                        <tr>
                                            <th>M-M</th>
                                            <th>S-S</th>
                                            <th>M-S</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <?php
                                            foreach($key as $val){
                                                echo "<tr>";
                                                foreach ($val as $indiVal){
                                                    echo "<td>$indiVal</td>";
                                                }
                                                echo "<tr>";
                                            }
                                        ?>
                                    </tbody>
                                </table>
                            </div>
                        <?php
                        $c++;
                    }
                ?>
            </div>
            
            <div class="donwRstBtn" align="center">
                <a class="btn btn_smt resetVal downloadButton mr-5" id="mutatedPDB">Mutated PDB</a>
                <a class="btn btn_smt resetVal downloadButton mr-5" id="interactionSummary">Interaction Summary</a>
                <a class="btn btn_smt resetVal downloadButton" id="resetXpro" onclick="window.location='X-Pro.php';" >Reset</a>
            </div>
        <?php
    } else {
        ?>
            <div class="alert alert-danger" role="alert">
                <label class="m-0"><?= $warningVariable; ?></label>
            </div>

            <div class="donwRstBtn" align="center">
                <a class="btn btn_smt resetVal downloadButton" id="resetXpro" onclick="window.location='X-Pro.php';" >Reset</a>
            </div>
        <?php
    }
?>
<p></p>     <!-- Create some bottom margin in the end -->

<script>
    upload_dir = <?= json_encode($upload_dir); ?>;
    pdbID = <?= json_encode($_GET['store_PDB_ID']); ?>;
    document.getElementById("mutatedImage").src = `${upload_dir}/mutated.png`;
    document.getElementById("notMutatedImage").src = `${upload_dir}/not_mutated.png`;

    document.getElementById("mutatedPDB").href = `${upload_dir}/mutated.pdb`;
    document.getElementById("mutatedPDB").download = `${pdbID}_Mutated_PDB.pdb`;

    document.getElementById("LigPlot_mutated").src = `${upload_dir}/output_files/mutated_ligplot.png`;
    document.getElementById("LigPlot_Image_Mutated").href = `${upload_dir}/output_files/mutated_ligplot.png`;
    document.getElementById("LigPlot_Image_Mutated").download = `${pdbID}_Mutated_Interaction.png`;

    document.getElementById("LigPlot_Unmutated").src = `${upload_dir}/output_files/altered_ligplot.png`;
    document.getElementById("LigPlot_Image_UnMutated").href = `${upload_dir}/output_files/altered_ligplot.png`;
    document.getElementById("LigPlot_Image_UnMutated").download = `${pdbID}_Normal_Interaction.png`;

    document.getElementById("interactionSummary").href = `${upload_dir}/output_files/ligplot_differences.txt`;
    document.getElementById("interactionSummary").download = `${pdbID}_differences.txt`;
</script>


<?php
    include "footer.php";

    // Remove directory which contains all files created by ligPlot and PyMol
    // rrmdir($upload_dir);
?>