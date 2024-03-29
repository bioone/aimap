Installation Instructions
*************************

Copyright (C) 2018-2023 Agriculture and Biology

Quickly install in ubuntu 22.04
==================

Step 1: Python Dependencies
---------------------------

The aimap requires Python3 as well as the following Python3 modules.

biopython
pandas
gffutils

These can be installed with pip3 using:

`pip3 install -r requirements.txt`

#Note that Python package dependencies should automatically be installed.

Step 2: Trim_galore Installation
--------------------------------

### Install cutadapt

`pip3 install --user --upgrade cutadapt`

### Install fastqc

`sudo apt-get install fastqc`

### Check that cutadapt is installed
`cutadapt --version`

### Check that FastQC is installed
`fastqc -v`

### Install Trim Galore
```
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
```
### Add Trim Galore path to $PATH environment variable
```
vi /etc/environment
add "export PATH=$PATH:/path/to/Trim Galore:$PATH" 
source /etc/environment
```

Step 3: Other Dependencies
---------------------------

### Install BWA

`sudo apt-get install bwa`

### Install samtools
`sudo apt-get install samtools`

### Install bcftools
`sudo apt-get install bcftools`

Step 4: Install aimap
---------------------

### Install aimap

`pip3 install aimap`

` chmod a+x /usr/local/bin/Adenosine_to_inosine.py `

### test

`Adenosine_to_inosine.py -h`

Finish!
