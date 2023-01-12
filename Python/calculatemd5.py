#!/usr/bin/env python

"""
This script takes two directory paths as input and prints the full path of the files in both directories
which have the same basenames. It looks for matching basenames in both directories and then calculates the MD5.
If the MD5 checksums match, it prints a message to the screen. If they do not match, it prints a message to the screen.
It also writes the results to a file.

Usage:
    python calculatemd5.py --dir1 /path/to/directory1 --dir2 /path/to/directory2 --output /path/to/outputfile.csv

"""

# std libraries
import os
import sys
import time
import csv
import subprocess
import argparse
import pandas as pd
import numpy as np

# Create the parser
parser = argparse.ArgumentParser()

parser.add_argument('--dir1', type = str, required = True, help = 'Input directory 1')
parser.add_argument('--dir2', type = str, required = True, help = 'Input directory 2')
parser.add_argument('--output', type = str, required = False, help = 'Output file destination. Default is current directory.')

args = parser.parse_args()


# Check both directories for validity and emptiness
def check_dirs(directory1, directory2):
    """
    Check if a directory exists and is not empty

    Parameters
    ----------
    directory : str
        Path to the directory to check

    Returns
    -------
    None
        
    """
    try:
        if not os.path.isdir(directory1) or not os.path.isdir(directory2):
            raise Exception(f"{directory1} or {directory2} is not a valid directory.")
        
        if not os.listdir(directory1) or not os.listdir(directory2):
            raise Exception(f"{directory1} or {directory2} is empty.")

    except Exception as e:
        print(e)
        sys.exit(1) 


# Calculate the MD5 checksums for a file
def calculate_MD5(file):
    """
    Calculate the MD5 checksums of a file

    Parameters
    ----------
    file : str
        Path to the file to calculate the checksums for

    Returns
    -------
    md5 : str
        The MD5 checksum of the file
        
    """
    
    # Creating the shell script command
    md5_cmd = "md5 " # md5 command
    md5_cmd += f" {file} " # input file
    md5_cmd += f"| " # pipe
    md5_cmd += "cut -d ' ' -f 4" # pull out checksums only
    
    # Running the command
    p = subprocess.run(md5_cmd, 
                       shell = True, 
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       universal_newlines=True)

    md5 = p.stdout.rstrip('\n') # remove trailing newline
    return md5



# Calculate the MD5 checksums for the files in the directories
def findmatches(directory1, directory2, output_file = "./checksum_comparison.csv"):
    """
    This function takes two directory paths as input and prints the full path of the files in both directories
    which have the same basenames. It looks for matching basenames in both directories and then calculates the MD5.
    If the MD5 checksums match, it prints a message to the screen. If they do not match, it prints a message to the screen. 
    It also writes the results to a file.   

    Parameters
    ----------
    directory1 : str
        Path to the first directory
    directory2 : str
        Path to the second directory
    output_file : str
        Path to the output file

    Returns
    -------
    csv file containing the results
    
        
    """
    # Start the timer
    starttime = time.time()

    # Create a dictionary to store the results in
    out_dict = {"File":[], "Dir1":[], "Dir2":[], "MD5Chksum_1":[], "MD5Chksum_2":[], "Match":[]}

    # Create a dictionary to store the basename of the file and the full path of the file in directory1
    files1 = {}
    for root, dirs, files in os.walk(directory1):
        for file in files:
            files1[os.path.basename(file)] = os.path.join(root, file)

    # Create a dictionary to store the basename of the file and the full path of the file in directory2        
    files2 = {}
    for root, dirs, files in os.walk(directory2):
        for file in files:
            files2[os.path.basename(file)] = os.path.join(root, file)
    
    # Find the intersection of the two dictionaries containing files from both directories
    matching_files = set(files1.keys()).intersection(set(files2.keys()))
    if not matching_files:
        print("No matching files were found.")
    else:
        print("Matching files were found.\n")
        print("Calculating MD5 checksums...\n")

        # chksum_count = 0
        for file in matching_files:
            # print(f"Evaluating File: {file}\n")
            # print(f"{file} in {directory1} : {files1[file]}")
            # print(f"{file} in {directory2} : {files2[file]}")
            chksum1 = calculate_MD5(files1[file])
            chksum2 = calculate_MD5(files2[file])
            match = chksum1 == chksum2
        
            if match == True:
                print(f'Checksums match: {file}\n')
                out_dict["File"].append(file)
                out_dict["Dir1"].append(files1[file])
                out_dict["Dir2"].append(files2[file])
                out_dict["MD5Chksum_1"].append(chksum1)
                out_dict["MD5Chksum_2"].append(chksum2)
                out_dict["Match"].append(match)
                # chksum_count += 1

            else:
                print(f'Checksums do not match: {file}\n')
                out_dict["File"].append(file)
                out_dict["Dir1"].append(files1[file])
                out_dict["Dir2"].append(files2[file])
                out_dict["MD5Chksum_1"].append(chksum1)
                out_dict["MD5Chksum_2"].append(chksum2)
                out_dict["Match"].append(match)
    
    # Write the results to a file
    df = pd.DataFrame(out_dict)

    match_count = df["Match"].value_counts()[True]
    # Percent of matching checksums
    percent = match_count/len(matching_files)
    print(f"Percent of matching checksums: {percent:.2%}")

    print("Printing results to file...")

    endtime = time.time()

    print(f"Elapsed time: {endtime - starttime:.4f} seconds")

    return df.to_csv(output_file, index=False) 





def main():

    # Check if the directories are valid
    check_dirs(args.dir1,  args.dir2)

    print('Directories are valid and not empty.\n')

    # Print the input directories
    print(f"Directory 1: {args.dir1}")
    print(f"Directory 2: {args.dir2}\n")

    # Find the matching files
    # Check if the output file is specified
    if args.output:
        print(f"Output file: {args.output}\n")
        findmatches(args.dir1, args.dir2, args.output)
    else:
        print("No output path specified defaulting to ./checksum_comparison.csv.\n")
        findmatches(args.dir1, args.dir2)

if __name__ == "__main__":
    main()

