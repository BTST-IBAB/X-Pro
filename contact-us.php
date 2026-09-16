<?php
    include_once "header.php";
    // header("Set-Cookie: cross-site-cookie=whatever; SameSite=None; Secure");
?>
<link type="text/css" rel="stylesheet" href="css/fontawesome_5.8.1.css">

<style>
    .contactLogo{
        font-size: 1.3rem;
        margin-right: 0.5rem;
    }

    .contactText{
        font-size: 14px;
    }
</style>

<h6 align="center" class="my-3" style="font-size: 18px;"><strong>Contact US</strong></h6>
<div class="mb-3 d-flex">
    <div class="mr-auto ml-5">
        <i class="fa fa-university contactLogo" aria-hidden="true"></i>
        <label class="contactText"><strong>Institute of Bioinformatics and Applied Biotechnology</strong></label><br>

        <i class="fas fa-map-marker-alt contactLogo" style="color: red;"></i>
        <label style="display: inline-grid; margin-left: 5px;" class="contactText">Biotech Park, GN Ramachandran
            Rd, Electronics City Phase 1
            <br>Electronic City, Bengaluru, Karnataka , India
        </label><br>

        <i class="fa fa-phone fa-rotate-90 contactLogo" aria-hidden="true" style="color: green;"></i>
        <label class="contactText"> 080 2852 8124 / 72088 90633</label><br>

        <i class="fa fa-envelope contactLogo" aria-hidden="true" style="color: #0066ff;"></i>
        <label class="contactText">st.bts@ibab.ac.in</label>
    </div>
    <div>
        <iframe style="border: 2px solid black;"
            src="https://www.google.com/maps/embed?pb=!1m18!1m12!1m3!1d3889.9871427503513!2d77.6564913141962!3d12.844106790939566!2m3!1f0!2f0!3f0!3m2!1i1024!2i768!4f13.1!3m3!1m2!1s0x3bae6c8368d47aa3%3A0x8b7ab32388f6f0ad!2sIBAB!5e0!3m2!1sen!2sin!4v1672641160059!5m2!1sen!2sin"
            width="500" height="400" style="border:0;" loading="lazy"
            referrerpolicy="no-referrer-when-downgrade">
        </iframe>
    </div>
</div>

<?php
    include_once "footer.php";
?>