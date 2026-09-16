<?php
    $servername = "localhost";
    $username   = "YOUR_DB_USERNAME";
    $password   = "YOUR_DB_PASSWORD";
    $dbname     = "YOUR_DB_NAME";

    $conn = mysqli_connect($servername, $username, $password, $dbname);
    if (mysqli_connect_error()) {
        echo "Failed to connect to MySQL: " . mysqli_connect_error();
    }
?>
