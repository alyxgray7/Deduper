#!/usr/bin/env python

### Import methods ### 
import argparse
import gzip
import re

### ARGPARSE ###
def get_args():
    """
    Will set a command line option to run arguments.
    """
    # Parser main
    parser = argparse.ArgumentParser(description = "Remove PCR duplicates from a SAM file. \
        SAM file must be previously sorted by ascending chromosome using tools such as \
        SAMtools. Sorted SAM files should be unzipped and contain reads from a single-end \
        sequencing experiment.")
    
    # Parser arguments
    parser.add_argument('-d', '--directory', help = "Desired directory for output files to \
        write in. ", required = True)
    parser.add_argument('-f','--SAMfile', help = "Absolute/path/to/<SAMfile.sam>. The file \
        should be presorted by chromosome location using a program such as SAMtools.", \
        required = True)
    parser.add_argument('-e','--which_end', help = "Argument to specify if <SAMfile.sam> \
        contains single or paired-end reads. Use 1 for single and 2 for paired. \
        WARNING: this script does no account for paired-end reads.", required = True, \
        type = int, nargs = 1)
    parser.add_argument('-u', '--UMIfile', help = "Absolute/path/to/<UMIfile.txt> \
        containing a list of UMIs used during the sequencing run.", required = True)
    #parser.add_argument('-h','--help', help = "Print help statement here.")
    
    return parser.parse_args()
args = get_args()

### MAKE UMI DICTIONARY ###
def make_UMI_dict(file):
    """
    This function will read in a file containing the UMIs used for sequencing. It will 
    output a dictionary containing as many elements as UMIs used.

    Example input:      Example output:
    CTGTTCAC            {CTGTTCAC:CTGTTCAC, ACGACTTG:ACGACTTG, ATCCATGG:ATCCATGG}
    ACGACTTG
    ATCCATGG
    """

    dictionary = {}
    with open(args.UMIfile, 'r') as umi_fh:
        LN = 0
        for umi in umi_fh:
            LN += 1
            umi = umi.strip()
            dictionary.setdefault(umi,0)
    return dictionary

### GET FILE NAME FOR OUTFILE NAMES ###
def get_filename(string):
    """
    This function will get the file basename from the argparse file argument and return
    the name as a string.
    """

    for i in args.SAMfile:
        file_name = args.SAMfile.split("/")[-1]
    return file_name
file_name = get_filename(args.SAMfile)

### FILES ###
sam_file = args.SAMfile
error_file = args.directory + file_name + "_err"
duplicates_file = args.directory + file_name + "_duplicates" 
original_file = args.directory + file_name + "_deduped"
output_txt = args.directory + "outputs.txt"

### MAIN ###
def dedupe(file):

    ### COUNTERS ###
    LN = 0          # Line Number
    RN = 0          # Read Number
    CN = 0          # Chromosome Number
    mapped = 0      # Mapped reads
    soft_clip = 0   # Soft-clipped reads
    duplicate = 0   # PCR duplicates
    original = 0    # Original reads
    error_UMI = 0   # UMIs containing errors

    ### OTHERS ###
    UMI_dict = make_UMI_dict(args.UMIfile)  # Dictionary; keys = known UMIs (created from UMI_file)
                                            #             values = none 
    QNAME_record_dict = {}                  # Dictionary; keys = QNAME
                                            #             values = Num of reads for this QNAME
    read_info = {}                          # Dictionary; keys = read characteristics (UMI, POS, strand_direction)
                                            #             values = none
    
    # Read the lines of the SAM_file
    for line in sam_fh:
        LN += 1

        # Write off header lines to all out-files
        if line.startswith("@") == True:
            err_fh.write(line) # UMI error file
            dup_fh.write(line) # Duplicates_file
            og_fh.write(line) # Originals_file

        # Loop through lines that are actually reads
        else:
            RN += 1

            # Split the line by tabs and assign variable name to bitwise flag
            read = line.split("\t")
            FLAG = int(read[1])

            # Check FLAG to see if the read is mapped
            if ((FLAG & 4) != 4):
                mapped += 1
                QNAME = read[2]

                # Check if read sequence is forward = "+"; reverse = "-")
                if ((FLAG & 16) == 16):
                    strand_DIR = "+"
                else:
                    strand_DIR = "-"

                # Check if current chromosome matches the chromosome counter
                if QNAME == CN:
                    # Increment count of reads for that QNAME
                    QNAME_record_dict[QNAME] += 1
                    UMI = read[0].split(":")[7]

                    # Check if UMI contains errors
                    if UMI in UMI_dict.keys():
                        POS = int(read[3])
                        CIGAR = read[5]
                        check_CIGAR = bool(re.search("(\dS)", CIGAR))

                        # Soft-clipping is not present in CIGAR
                        if check_CIGAR == False:

                            # Save read characteristics as tuple
                            info1 = (UMI, POS, strand_DIR)

                            # Read info is original
                            if info1 not in read_info:
                                original += 1
                                read_info.setdefault(info1,0)
                                og_fh.write(line) # Write to originals_file

                            # Read info is duplicated
                            else:
                                duplicate += 1
                                dup_fh.write(line) # Write to duplicates_file
                        
                        # Soft-clipping is present in CIGAR
                        else:
                            soft_clip += 1
                            clip = int(re.split("(\S)", CIGAR, 1)[1])
                            new_POS = POS - clip
                            info2 = (UMI, new_POS, strand_DIR)

                            # Read info is original
                            if info2 not in read_info:
                                original += 1
                                read_info.setdefault(info2,0)
                                og_fh.write(line) # Originals file

                            # Read info is duplicated
                            else:
                                duplicate += 1
                                dup_fh.write(line) # Duplicates file
                        
                    # UMI contains errors
                    else:
                        error_UMI += 1
                        err_fh.write(line) # UMI error file
                
                # Current chromosome does not match chromosome counter
                else:
                    QNAME_record_dict.setdefault(QNAME,1)
                    CN = QNAME # Change chromosome number to loop through next chromosome
                    read_info.clear() # Clear dictionary of read information for the next chromosome set
    
    ### Print return statement in output file ###
    with open("outputs.txt", "w") as output_fh:
        output_fh.write("Read type" + "\t" + "Num. reads" + "\t" + "% /total" + "\n")
        output_fh.write("Total reads" + "\t" + str(RN) + "\t" + str(RN/RN * 100) + "\n")
        percent_mapped = round((((mapped-duplicate-error_UMI)/RN) * 100), 2)
        output_fh.write("Deduplicated" + "\t" + str(mapped-duplicate-error_UMI) + "\t" + str(percent_mapped) + "\n")
        percent_duplicate = round(((duplicate/mapped)*100), 2)
        output_fh.write("PCR Duplicates" + "\t" + str(duplicate) + "\t" + str(percent_duplicate) + "\n")
        percent_UMIerror = round(((error_UMI/mapped) * 100), 2)
        output_fh.write("Errors" + "\t" + str(error_UMI) + "\t" + str(percent_UMIerror) + "\n\n")


        # Sorted table of QNAME and number of reads
        output_fh.write("QNAME" + "\t" + "Number of reads" + "\n")
        for item in sorted(QNAME_record_dict):
            output_fh.write(str(item) + "\t" + str(QNAME_record_dict[item]) + "\n")

    return RN, mapped, original, duplicate, QNAME_record_dict

while args.which_end[0] == 1:
    # Check if the file is gzipped
    try:
        with gzip.open(args.SAMfile, "tr") as sam_fh, \
        open(error_file, "w") as err_fh, \
        open(duplicates_file, "w") as dup_fh, \
        open(original_file, "w") as og_fh:
            dedupe(sam_fh)
        print("Deduping with a compressed SAMfile.")
        break

    # File is not gzipped
    except OSError as error:
        with open(args.SAMfile, "r") as sam_fh, \
        open(error_file, "w") as err_fh, \
        open(duplicates_file, "w") as dup_fh, \
        open(original_file, "w") as og_fh:
            dedupe(sam_fh)
        print("Deduping with an uncompressed SAMfile")
        break

while args.which_end[0] == 2:
    #if args.paired_end == True:
    print("This option can't be completed.")
    print("Current state of the script does not account for paired-end reads, only single-end reads.")
    print("Please use [-e 1] argument.")
    break

while args.which_end[0] >= 3:
    print("Please enter [-e 1] for single-end reads or [-e 2] for paired-end reads.")
    break
