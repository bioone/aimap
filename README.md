# README.md (aimap)
## Overview
`aimap` is a Python3 module that provides support for finding the modification from Adenosine to inosine (A to I). This is a one-step automated A to I modification analysis pipeline (AIMAP) for bacteria that works on illumine RNA-seq raw reads. 
The pipeline integrates reads quality control, adaptor removal, reference genome mapping, coverage pileup and synonymous-nonsynonymous annotation.

`aimap` installs one scripts into the `$PATH`:

* `aimap.py` that enables command-line A to I analysis.

## Installation

The easiest way to install `aimap` is to use `pip3`:

```
pip3 install aimap
```

## DEPENDENCIES

Note that Python package dependencies should automatically be installed

### Python packages dependencies

* **Biopython** <https://biopython.org/>
* **pandas** <https://pandas.pydata.org/>
* **gffutils** <https://github.com/daler/gffutils>

#### Other tools

* **trim_galore** <http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>
* **BWA** <https://github.com/lh3/bwa>
* **Samtools** <http://www.htslib.org/>

**Note that you can find the full install instruction in install file.**

## Running `aimap`

The `aimap.py` script - installed as part of this package - enables straightforward ANI analysis at the command-line, and uses the `aimap` module behind the scenes.

You can get a summary of available command-line options with `aimap.py -h`

```
$ aimap.py -h
usage: aimap.py [-h] [--version] -o OUTDIRNAME -g GENOMENAME -a ANNOTATION
                --outfile_name OUTFILE_NAME [-l LENGTH] [-f FILENAME]
                [-f1 FILENAME1] [-f2 FILENAME2] [-m {single,paired}]
                [-t THREADS] [--logfile LOGFILE]
                [--editing_level EDITING_LEVEL] [--coverage COVERAGE]

[…]
```

## Licensing

Unless otherwise indicated, all code is subject to the following agreement:

    (c) Agriculture and Biology 2018, 2019
    Author: Sai Wang

    Contact: skyws@outlook.com

    Address: 
    Sai Wang,
    Agriculture and Biology,
    Errol Road,
    Invergowrie,
    Dundee,
    DD6 9LH,
    ShangHai,
    CN

The MIT License

