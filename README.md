vcf2align.py
====

A simple Python script to generate a multiple sequence alignment (FASTA, NEXUS, PHYLIP, etc.) from a VCF/BCF file.

For those occasions when a program requires an alignment in plain text, or you simply are old-school and want to use BioEdit.

If a reference genome is indicated, the alignment will include the (invariant) positions in the region that are not in the VCF for all samples, as well as an additional sequence representing the reference sample.

The script does not attempt to reinvent the wheel. It depends on

- pysam (https://pysam.readthedocs.io) to parse and extract records from the VCF file
- biopython (https://biopython.org/) to manipulate and export the sequence alignments


## Usage

```bash
python vcf2align.py [-f format (default: fasta)] [-ref reference] [-o output] VCF region
```


* The region coordinates must be indicated like in samtools/bcftools (1-based, inclusive), i.e. `contig:start-end` (for example, `chr22:4561540-4561550`)

* To allow direct access to the required genomic region, both the reference and the VCF should be indexed (with `samtools faidx` and `tabix`, respectively).

* The script is only intended for diploids.

* Heterozygous positions or missing genotypes are denoted with the standard IUPAC ambiguity codes.

* In the case of positions with alleles that span several bases (indels), the script will insert dashes (`-`) to those positions to keep all sequences in the alignment at the same length. Therefore, sequence lengths may be shorter or longer than expected based on the reference coordinates. In this case, the script will throw a warning to alert the user, indicating the potentially problematic positions.


## Installation

1. Clone this repository and change into the folder


2. Create a clean environment and install the requirements

```bash
# using python venv
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# OR using conda
conda create -n env python=3
conda activate env
conda install --file requirements.txt
```

3. See the command line options by running

```bash
python vcf2align.py -h
```
