please note: Windows users will not currently be able to install the package due to some compilation issues with the Indigo dependency. We are currently working with Indigo to right this and should hopefully be able to support the Windows OS soon.

# camb - chemistry aware model builder
======

camb - chemistry aware model builder is an R package that can be used for the rapid generation of quantitative predictive models in the area of medicinal chemistry (QSAR, QSPR, QSAM, PCM). It is aimed at both advanced and beginner R users.
Its capabilities include the standardisation of representation of chemical structures, computation of 905 two-dimensional and 14 fingerprint type descriptors for small molecules, 8 types of amino acid descriptors, 13 whole protein sequence descriptors, filter methods for feature selection, generation of predictive models (R package caret), as well as techniques to ensemble these models (R package caretEnsemble).
Results can be visualised through high-quality, customisable plots (R package ggplot2).

This is the root folder which holds the package folder as well as other folders which contain examples of package use.

Two tutorials concerning the application of camb in the context of QSPR and Proteochemometrics are available in the examples folder.

Coding is done with the Google's R style guide: http://google-styleguide.googlecode.com/svn/trunk/Rguide.xml#functiondefinition

# INSTALLATION:

NOTE: If the installation instructions here no longer work for you (because of updates to OSes and R), you can use the below method (tested 11/7/2017 with Ubuntu 16.04.2 x64) to setup a droplet on Digital Ocean which you can access through R Studio Server through the browser.

--- Digital Ocean R Studio Server Method ---

Step 1 (allocating swap space) is only required if you choose the smallest droplet ($5/m). If you chose
the $20/m droplet then you can skip this step and the installation goes much faster.

# step 0 - create a Digital Ocean droplet
- Register with Digital Ocean https://m.do.co/c/ae523dc7d5e4
- Create a droplet with the Ubuntu 16.04.2 x64 operating system (add an ssh key for easy access)
- ssh into your droplet and continue to the next step

# step 1 - allocate swap space
- sudo fallocate -l 2G /swapfile
- sudo chmod 600 /swapfile
- sudo mkswap /swapfile
- sudo swapon /swapfile
- echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab

# step 2 - install requirements for RStudio Server and camb
- sudo apt-get update
- sudo apt-get install -y r-base gdebi-core cmake libcurl4-gnutls-dev
- libssl-dev libfreetype6-dev libfontconfig1-dev
- wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
- sudo gdebi rstudio-server-1.0.143-amd64.deb
- sudo adduser david

# step 3 - install camb
- connect through the browser to http://your_droplet_IP:8787 (the R Studio server)
- > install.packages('devtools')
- > library(devtools)
- > install_github("cambDI/camb/camb")

--- End Digital Ocean R Studio Server Method ---

--- Original installation instructions (may be prone to errors as dependencies have changed) ---

camb can be installed the by typing: library(devtools); install_github("cambDI/camb/camb")

Additionally, one can download the zip file and type in the command line (after unzipping the file): R CMD install camb

General requirements:

1. Make sure that cmake (version >= 2.8) is installed.
2. The following python packages are required to calculate Morgan fingerprints: argparse, numpy, rdkit, operator, sys and os.
3. The rJava package is required to calculate PaDEL descriptors.

Please read the following operating system specific instructions if the above method fails for you at some point.

1. OSX: you'll need to have the Xcode installed. If you don't already, this install should be triggered automatically during the installation although you may need to restart the installation: install_github("cambDI/camb/camb"). If it does not install automatically, you can install it from here: https://developer.apple.com/xcode/downloads/ relevant error: 'xcode-select: note: no developer tools were found at '/Applications/Xcode.app', requesting install. Choose an option in the dialog to download the command line developer tools.'

2. UNIX: In principle, you should not encounter any problem. Please send us a mail if you face any problem when installing camb.

3. CentOs (UNIX): make sure the package rJava can be installed (in case you want to use the function GeneratePadelDescriptors). CentOs users have reported issues in this regard. 

Please email us with any other problems you might encounter! Thanks!

# Tips for SDF files

Line breaks should not appear in the fields of any SDF file processed. 
For instance (line 57 refers to the field -property- name, i.e. Description, whereas lines 58 and 59 refer to the field value, i.e. the description for a given molecule):

- WRONG:

57 Description

58 Weak adenosine receptor antagonist; weak phosphodiesterase 

59 inhibitor; diuretic; smooth muscle relaxant

- Correct:

57 Description

58 Weak adenosine receptor antagonist; weak phosphodiesterase inhibitor; diuretic; smooth muscle relaxant


