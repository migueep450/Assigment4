#!/usr/bin/env python3
#We create a workflows in which all the sequences has been analyzed from a fasta/q file. 
#Each sequence is going to be aligned with the yeast genome
#IMPUT: a FASTA/Q file
#OUTPUT: a SAM file with the aligned sequences from the fasta/q file
#FUNCTIONS(2): FROM 15 to 68

import matplotlib.pyplot as plt
from collections import Counter
import shutil
import glob
import os
import sys

def create_folder(name_folder):
    '''
    Create a new folder in which is going to contain some files.
    '''
    try:
        # Create a new direcctory
        os.mkdir(name_folder)
    except FileExistsError:
        # If the direcctory already exits print.
        print(f"The directory {name_folder} already exists.")
	
	
def trim_seq_right(seq: str, trim_right: int, qualities: str = '') -> str:
    '''
    Trims the input sequence `seq` (and the qualities if they are introduced) `trim_right` nucleotides at the enf of
    the sequence
    :param seq: string encoding a DNA sequence
    :param trim_right: length of the trim that is going to be removed at the end of the sequence.
    :param qualities: string containing sequence quality (Optional)
    :return:
        processed_seq and processed_qualities: sequence and qualities trimmed
    '''
    # Split the sequence
    processed_seq = seq[:-trim_right]
    # Split the quality string, as long as it is not empty
    if qualities != '':
        assert len(seq) == len(qualities), 'Something went wrong while trimming qualities and the nucleotides'
        qualities = qualities[:-trim_right]
        qualities = f"+\n{qualities}\n"
    return f"{processed_seq}\n{qualities}"

def plot_hist(dict):
    '''
    Represent the frequencies of trimers, which are included in the sequences from the imput file, in a graph. These frequencies and the
    trimers are introduced with a dicctionary. 
    '''
    plt.figure(figsize=(20, 5))
    #Paint frequency bars for each trimer.
    plt.bar(*zip(*dict.items()))
    #Set label.
    plt.xlabel("Trimers")
    plt.ylabel("Frequency")
    #Put the name for each bar rotated 45º.
    plt.xticks(rotation=45)
    #Save the image in png format
    plt.savefig("trimers_hist.png")
    #Show the graphic created
    plt.show()

def search_trimmers(seq):
    '''
    Return a list in which contains all the possible trimmers from the imput sequence.
    '''
    return [seq[i:i+3] for i in range(len(seq)-2)]

input_file = sys.argv[1]
extension = input_file[-5:]
split_name = "split"
split_dir = "split_seq"
align_dir = "alignments"
#Create the path for split files.
save_split = os.path.join(split_dir, split_name)
nbases_trim = 20
ref_genome = "ref/genome_yeast.fasta"
num_align = 0


# Check if the input_file is filled.
assert input_file is not None, 'Please, introduce the name of the input file'
# Check if the name of the input file is a string.
assert not input_file.isdigit(), 'Please, do not introduce numbers as input file'
# Check if the input file is a fasta or fastq file.
assert input_file[-5:].lower() == 'fasta' or input_file[-5:].lower() == 'fastq', 'The input file must be .fasta or .fastq'

create_folder(split_dir)
create_folder(align_dir)

trimers=[]
num_seq = 0

with open(input_file, 'rt') as f:
    if extension.lower()=='fasta':
        while True:
            tag = f.readline().strip()
            if not tag: break
            seq = f.readline().strip().upper()
            # Store all the trimers that exist in the sequence
            seq_trimers = search_trimmers(seq)
            # Add all the element of the seq_trimers list to the trimers list.
            trimers.extend(seq_trimers)
            # Generate a new file with its tag and seq.
            with open(f'{save_split}{num_seq}.{extension}', 'wt') as f_split:
                f_split.write(f'{tag}\n{seq}\n')
            # Sequences'counter to put it in the name of the new file.
            num_seq += 1
            
    else:
        while True:
            tag = f.readline().strip()
            if not tag: break
            seq = f.readline().strip().upper()
            f.readline()  # Ignore quality tag '+'
            qualities = f.readline().strip()
            # Store all the trimers that exist in the sequence
            seq_trimers = search_trimmers(seq)
            # Add all the element of the seq_trimers list to the trimers list.
            trimers.extend(seq_trimers)
            # Generate a new file with its tag, seq and its qualities
            with open(f'{save_split}{num_seq}.{extension}', 'wt') as f_split:
                f_split.write(f'{tag}\n{seq}\n+\n{qualities}\n')
            # Sequences'counter to put it in the name of the new file.
            num_seq += 1

# Create a dictionary with all the trimers that we found in the sequences'file 
# and count how many of each trimer there are
hist = dict(Counter(sorted(trimers)))
plot_hist(hist)

# Extract each file which has only a sequence, from the first directory which was created at the beggining.
# Create new files in which the sequences and the possible quality's lines are trimmed
for split in glob.glob(f'{save_split}*.{extension}'):
    with open(split, 'rt') as f:
        with open(f'{split[:-6]}trimmed.{extension}', 'wt') as t:
            if extension.lower()=='fasta':
                tag = f.readline().strip()
                if not tag: break
                seq = f.readline().strip().upper()
                #Will be need it to run the function 'trim_seq_rigth'
                qualities = ''
            else:
                tag = f.readline().strip()
                if not tag: break
                seq = f.readline().strip().upper()
                f.readline()  # Ignore '+'
                qualities = f.readline().strip()
            # Trimmed the sequence and the quality line. Then, store them with the name of the sequence in a new file
            t.write(f'{tag}\n{trim_seq_right(seq,nbases_trim,qualities=qualities)}')

# Put the reference genome in the bwa application to align it with ours splited sequences.
os.system(f"bwa index {ref_genome}")
# Generate SAM files with the alignment of each sequence(from the splited sequence files) and the reference genome.
for trim_split in glob.glob(f'{save_split}*trimmed.{extension}'):
    # Store all the new files in the directory aling_dir
    os.system(f"bwa mem {ref_genome} {trim_split} > '{align_dir}{trim_split[9:-6]}.sam'")

# Remove the folder and its contents, in which was all the trimmed sequences files
shutil.rmtree(split_dir)

# Open all SAM files in aling_dir
for align in glob.glob(f'{align_dir}/split*trimmed.sam'):
    # Create a new file in which is going to contain all the alignment from the SAM files
    with open(f'{align_dir}/merge.sam', 'at+') as m:
        with open(align) as a:
            while True:
                line = a.readline().strip()
                if not line: break
                # Ignore the headers
                elif line[0] == '@': continue
                else:
                    # Count how many reads has been aligned
                    num_align += 1
                    # Add the alingment line to the new merge.sam file
                    m.write(f'{line}\n')

# Sort all the alignments in the merged sam's file by chromosome and position
os.system(f"sort -k4 -n {align_dir}/merge.sam | sort -k3 > sorted_aligment.sam")

print(f"{num_align} reads have been aligned with the yeast genome")

# Remove the directory with all its contens, in which was all the SAM's files
shutil.rmtree(align_dir)
