# Assignment 4. Small pipeline


Two scripts (BASH/Python) that, given an input FASTA/FASTQ file, implement the following pipeline:

1. Counts the abundance of each possible 3-mers (histogram)

2. Splits the input in FASTA/FASTQ into files of only 1 sequence each

3. Hard-trims (from the right) all sequences from all files 20nt

4. Using the genome mapping tool BWA and the reference genome of the Saccharomyces cerevisiae (any strain will do), aligns each of the files producing its corresponding SAM file.

5. Merges all SAM files ignoring headers

6. Sorts the SAM file by chromosome and position

7. Computes how many reads have been aligned

The output folder will contain the histogram and the alignments sorted.

## Executing BASH script

`bash assigment4.sh <input name-file.fasta/q>`

## Python

`python3 assigment.py <input name-file.fasta/q>`

