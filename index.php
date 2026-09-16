
<?php
    include 'connect.php';
    include_once "header.php";

    function getUserIpAddr(){
        if(!empty($_SERVER['HTTP_CLIENT_IP'])){
            $ip = $_SERVER['HTTP_CLIENT_IP'];
        }elseif(!empty($_SERVER['HTTP_X_FORWARDED_FOR'])){
            $ip = $_SERVER['HTTP_X_FORWARDED_FOR'];
        }else{
            $ip = $_SERVER['REMOTE_ADDR'];
        }
        return $ip;
    }
    $ip = getUserIpAddr();
    
    $qry = "SELECT * FROM `ipdb` WHERE `ip` = '$ip'";
    $result = mysqli_query($conn,$qry);
    $num = mysqli_num_rows($result);
    if ($num == 0){
        $qry3 = "INSERT INTO `ipdb`(`ip`) VALUES ('$ip')";
        mysqli_query($conn,$qry3);
        $qry1 = "SELECT * FROM `counter` WHERE `id` = 0";
        $result1 = mysqli_query($conn,$qry1);
        $row1 = mysqli_fetch_array($result1, MYSQLI_ASSOC);
        $count = $row1['visitors'];
        $count = $count + 1;
        $qry2 = "UPDATE `counter` SET `visitors`='$count' WHERE `id`=0";
        $result2=mysqli_query($conn,$qry2);
    }else{
        $qry1 = "SELECT * FROM `counter` WHERE `id` = 0";
        $result1 = mysqli_query($conn,$qry1);
        $row1 = mysqli_fetch_array($result1, MYSQLI_ASSOC);
        $count = $row1['visitors'];
    }
    $numlength = strlen((string)$count);
?>


<main class="main">
    <article>
        <h4><strong> About BioToolSuite </strong></h4>
        <p>BioToolSuite is an open source web-based application, developed by Prof. S. Thiyagarajan's group, that
            integrates comprehensive tools, MetaPredictor and AnnoDUF. This platform provide valuable insights into
            protein function and phenotypic changes.
        </p>
        <!-- <p>
            <a href="predictor.php">Meta-Predictor</a> utilizes supervised machine
            learning, specifically logistic regression with 10 feature tools, to predict the phenotypic effects of
            single amino acid mutations based on genotype. The tool also provides the probability of accuracy for
            each prediction, enhancing its reliability. 
        </p> -->
        <p>
            <a href="annoduf.php">AnnoDUF</a> leverages psi-blast and data mining techniques to annotate proteins with
            domains of unknown function (DUFs). AnnoDUF identifies putative hits and provides valuable information about
            these elusive proteins domain.
        </p>
        <!-- <p>
            <a href="X-Pro.php">X-Pro</a> is a computational tool that uses PyMOL and LigPlot to explore potential
            rotameric configuration resulting from a point mutation within a protein structure. It starts by using
            PyMOL's "mutagenesis" wizard to substitute amino acids and then employs LigPlot to study interaction with
            neighbouring amino acids.
        </p> -->
        <p>
            Together, these tools offer a comprehensive solution for protein analysis and functional annotation.
        </p>
    </article>
    <br>
    <aside class="mx-4 p-1">
        <div class="border border-dark" id="heading">
            <div class="border border-dark p-2" align="center"><strong>Visitor Count</strong></div>
        </div>
        <div class="border border-dark bg-info">
            <div class="border border-dark p-1" align="center"><?= $count; ?></div>
        </div>
    </aside>
</main>

<?php
    include_once "footer.php";
?>