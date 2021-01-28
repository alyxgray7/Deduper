# Deduper

## The Problem:

Before sequencing, DNA fragments are amplified and adapter-ligated using PCR. PCR sometimes creates duplicate reads that are present after sequencing. PCR duplicates in RNA-seq analysis will distort the results because duplicates are not actual biological reads, but artifacts of the sequence process. It's important to control for the PCR/sequencing reactions during analysis so variation in expressions are only biologically relevant.

After sequencing, QC/trimming, and demultiplexing, reads are aligned to a reference and outputted into a SAM file. These duplicates can be identified in a SAM file using these pieces of information:

        1. UMI sequence. Duplicates stemming from the same PCR product will share the same UMI.
                a. 96 unique UMIs (8 nucleotides, column 1:end)
        2. Location of the aligned read. Duplicates will have the same chromosome and starting coordinate position.
                a. Chromosome.
                b. Position start.
                        i. Must account for soft clipping (found in CIGAR string).
                        ii. If "S" is set, subtract from the left-most consumed base position.
        3. Strand specificity. Duplicates will arise from the same strand reading direction (forward or reverse).

PCR duplicates identified in the reads by these characteristics need to be removed into a separate file. Test data is from single-end sequencing and UMIs have been placed in QNAME of the header.

**This algorithm is designed for single-end data,** with 96 UMIs. UMI information will be in the QNAME, like so: ```NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT```. Any reads with UMI errors are removed into a separate file.

This Python code:
- Is Python 3 compatible code
- Includes the following argparse options
    - ```-f```, ```--file```: required arg, absolute file path
    - ```-p```, ```--paired```: optional arg, designates file is paired end (however, this is ONLY single-end compatable)
    - ```-u```, ```--umi```: required arg, designates file containing the list of UMIs
    - ```-h```, ```--help```: optional arg, prints a USEFUL help message 
        - If the script is not capable of dealing with a particular option (ex: no paired-end functionality), the script wil print an error message and quit
- Will output the first read encountered if duplicates are found
- Will output a properly formatted SAM file with “_deduped” appended to the filename


