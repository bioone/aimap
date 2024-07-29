# README.md (aimap)
## Overview
`aimap` is a Python3 module that provides support for finding the modification from Adenosine to inosine (A to I). This is a one-step automated A to I modification analysis pipeline (AIMAP) for bacteria that works on illumine RNA-seq raw reads. 
The pipeline integrates reads quality control, adaptor removal, reference genome mapping, coverage pileup and synonymous-nonsynonymous annotation.

`aimap` installs one scripts into the `$PATH`:

* `Adenosine_to_inosine.py` that enables command-line A to I analysis.

## Installation

### Recommended installation method

The easiest way to install `aimap` is to use ` conda`:

```
conda create -n aimap_env python=3.10
conda activate aimap_env
conda install bioconda::fastp
conda install -c bioconda bcftools=1.16
conda install -c bioconda samtools=1.16
pip install aimap
```

Need to give execute file permission

```
chmod +x /path/to/aimap_env/bin/Adenosine_to_inosine.py
```

### Manual install 

The easiest way to install `aimap` is to use `pip3`:

```
pip3 install aimap
```

## DEPENDENCIES

Note that Python package dependencies should automatically be installed

### Python packages dependencies

* **Biopython** <https://biopython.org/>
* **pandas** <https://pandas.pydata.org/>

#### Other tools

* **fastp** <[OpenGene/fastp: An ultra-fast all-in-one FASTQ preprocessor (QC/adapters/trimming/filtering/splitting/merging...) (github.com)](https://github.com/OpenGene/fastp)>
* **BWA** <https://github.com/lh3/bwa>
* **Samtools** <http://www.htslib.org/>

**Note that you can find the full install instruction in install file.**

## Running `aimap`

The `Adenosine_to_inosine.py` script - installed as part of this package - enables straightforward ANI analysis at the command-line, and uses the `aimap` module behind the scenes.

You can get a summary of available command-line options with `Adenosine_to_inosine.py -h`

```
$ Adenosine_to_inosine.py -h
usage: Adenosine_to_inosine.py [-h] [--version] -o OUTDIRNAME -g GENOMENAME -a ANNOTATION
                               --outfile_name OUTFILE_NAME [-l LENGTH] [-f FILENAME]
                               [-f1 FILENAME1] [-f2 FILENAME2] [-m {single,paired}] [-t THREADS]
                               [--logfile LOGFILE]

[…]
```
For example you may type:
```
Adenosine_to_inosine.py -o /path/to/result -g /path/to/genome/file -a /path/to/annotation/file --outfile_name sample_name -f1 /path/to/r1.fastq -f2 /path/to/r2.fastq -m paired -t 4
```

## Result

At last, you will get a result file, it will be separated by “\t”. You can see the `*_A_I_result.tsv` for more information.

| Accession | Position  | Old_base  | New_base  | Raw read depth | Coverage  | Edit_level | snp_coverage | snp_f_coverage | snp_r_coverage | Gene_biotype | locus_tag | Gene_name | Gene_strand | gene_start | gene_end  | Amino acid_change | codon_num |
| --------- | --------- | --------- | --------- | -------------- | --------- | ---------- | ------------ | -------------- | -------------- | ------------ | --------- | --------- | ----------- | ---------- | --------- | ----------------- | --------- |
| (Example) | (Example) | (Example) | (Example) | (Example)      | (Example) | (Example)  | (Example)    | (Example)      | (Example)      | (Example)    | (Example) | (Example) | (Example)   | (Example)  | (Example) | (Example)         | (Example) |
| (Example) | (Example) | (Example) | (Example) | (Example)      | (Example) | (Example)  | (Example)    | (Example)      | (Example)      | (Example)    | (Example) | (Example) | (Example)   | (Example)  | (Example) | (Example)         | (Example) |
| (Example) | (Example) | (Example) | (Example) | (Example)      | (Example) | (Example)  | (Example)    | (Example)      | (Example)      | (Example)    | (Example) | (Example) | (Example)   | (Example)  | (Example) | (Example)         | (Example) |

Below are specific descriptions for each column：

| Column Name       | Description                                                  |
| ----------------- | ------------------------------------------------------------ |
| Accession         | The unique identifier assigned to a genetic sequence in a public database, such as NCBI. |
| Position          | The specific location of a nucleotide within the genetic sequence, often given as a numerical coordinate. |
| Old_base          | The original nucleotide base at a specific position before any changes or mutations. |
| New_base          | The mutated or altered nucleotide base at a specific position after a change has occurred. |
| Raw read depth    | The total number of sequencing reads that cover a particular genomic position. |
| Coverage          | The average number of sequencing reads that align to a specific region of the genome, indicating how well that region is covered by the sequencing data. |
| Edit_level        | The percentage of reads supporting a specific edit (mutation) at a given position relative to the total read depth. |
| snp_coverage      | The depth of coverage specifically at the site of a single nucleotide polymorphism (SNP). |
| snp_f_coverage    | The depth of forward reads covering the SNP position.        |
| snp_r_coverage    | The depth of reverse reads covering the SNP position.        |
| Gene_biotype      | The classification of a gene based on its predicted biological function or properties. |
| locus_tag         | A unique identifier assigned to a gene within a genome, often used when the gene has not been given a standard name. |
| Gene_name         | The common or standard name of a gene, which may be based on its function, phenotype, or historical nomenclature. |
| Gene_strand       | The DNA strand (plus or minus) on which the gene is located. |
| gene_start        | The starting position of a gene on the DNA sequence.         |
| gene_end          | The ending position of a gene on the DNA sequence.           |
| Amino acid_change | The specific change in the amino acid sequence resulting from a nucleotide mutation. |
| codon_num         | The position of the codon within the gene sequence that is affected by a mutation, often used to indicate where in the protein sequence the amino acid change occurs. |



## Licensing

Unless otherwise indicated, all code is subject to the following agreement:

    (c) Agriculture and Biology 2018, 2024
    Author: Sai Wang
    
    Contact: skyws@outlook.com
    
    Address: 
    Sai Wang,
    Agriculture and Biology,
    Shanghai Jiao Tong University,
    ShangHai,
    CN

The MIT License

Copyright (c) 2018-2024 Agriculture and Biology

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
