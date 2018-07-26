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
For example you may type:
```
aimap.py -o /path/to/result -g /path/to/genome/file -a /path/to/annotation/file --outfile_name sample_name -f1 /path/to/r1.fastq -f2 /path/to/r2.fastq -m paired -t 4
```

## Result

At last, you will get a result file, it will be separated by “\t”. You can see the test_result for more information.

## Licensing

Unless otherwise indicated, all code is subject to the following agreement:

    (c) Agriculture and Biology 2018, 2019
    Author: Sai Wang

    Contact: skyws@outlook.com

    Address: 
    Sai Wang,
    Agriculture and Biology,
    Shanghai Jiao Tong University,
    ShangHai,
    CN

The MIT License

Copyright (c) 2018-2019 Agriculture and Biology

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
