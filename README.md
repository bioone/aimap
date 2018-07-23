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

#### other tools

