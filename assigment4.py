#!/usr/bin/env python3

import matplotlib.pyplot as plt
from collections import Counter
import shutil
import glob
import os


#TODO posar input file com argv[1] un cop haguem acabat el programa
input_file = "small_real_sample.fastq"
extension = input_file[-5:]
split_name = "split"
split_dir = "split_seq"
align_dir = "alignments"
save_split = os.path.join(split_dir, split_name)
nbases_trim = 20
ref_genome = "ref/genome_yeast.fasta"
num_align = 0
# Per si es vol posar input file com un argument
assert input_file is not None, 'Please, introduce the name of the input file'
assert not input_file.isdigit(), 'Please, do not introduce numbers as input file'
assert input_file[-5:].lower() == 'fasta' or input_file[-5:].lower() == 'fastq', 'The input file must be .fasta or .fastq'

try:
    os.mkdir(split_dir)
except FileExistsError:
    print(f"The directory {split_dir} already exists.")
try:
    os.mkdir(align_dir)
except FileExistsError:
    print(f"The directory {align_dir} already exists.")


def trim_seq_right(seq: str, trim_right: int, qualities: str = '') -> str and str and str:
    """
    Trims the input sequence `seq` (and the qualities if they are introduced) `trim_right` nucleotides at the enf of
    the sequence
    :param seq: string encoding a DNA sequence
    :param trim_right: length of the trim that is going to be removed at the end of the sequence.
    :param qualities: string containing sequence quality (Optional)
    :return:
        processed_seq: sequence trimmed
        qualities: qualities trimmed in the same way as `seq`
    """
    assert len(seq) == len(qualities)
    processed_seq = seq[:-trim_right]
    if qualities != '':
        processed_qualities = qualities[:-trim_right]
        assert len(processed_seq) == len(processed_qualities), 'Something went wrong while trimming qualities and the nucleotides'
        processed_qualities = f"+\n{processed_qualities}\n"
    return f"{processed_seq}\n{processed_qualities}"


def plot_hist(dict):
    #TODO posar una opcio per guardar la imatge
    plt.figure(figsize=(20, 5))
    plt.bar(*zip(*dict.items()))
    plt.xlabel("Trimers")
    plt.ylabel("Frequency")
    plt.xticks(rotation=45)
    plt.savefig("trimers_hist.png")
    plt.show()


trimers=[]
num_seq = 0
with open(input_file, 'rt') as f:
    if extension.lower()=='fasta':
        while True:
            tag = f.readline().strip()
            if not tag: break
            seq = f.readline().strip().upper()
            # Count trimers in seq
            seq_trimers = [seq[i:i+3] for i in range(len(seq)-2)] # a una funcio
            trimers.extend(seq_trimers)
            # Generate new file with one seq
            with open(f'{save_split}{num_seq}.{extension}', 'wt') as f_split:
                f_split.write(f'{tag}\n{seq}\n')
            num_seq += 1
    else:
        while True:
            tag = f.readline().strip()
            if not tag: break
            seq = f.readline().strip().upper()
            f.readline()  # Ignore '+'
            qualities = f.readline().strip()
            seq_trimers = [seq[i:i+3] for i in range(len(seq)-2)] # a una funcio
            trimers.extend(seq_trimers)
            with open(f'{save_split}{num_seq}.{extension}', 'wt') as f_split:
                f_split.write(f'{tag}\n{seq}\n+\n{qualities}\n')
            num_seq += 1

hist = dict(Counter(sorted(trimers)))

plot_hist(hist)

for split in glob.glob(f'{save_split}*.{extension}'):
    with open(split, 'rt') as f:
        with open(f'{split[:-6]}trimmed.{extension}', 'wt') as t:
            if extension.lower()=='fasta':
                tag = f.readline().strip()
                if not tag: break
                seq = f.readline().strip().upper()
                qualities = ''
            else:
                tag = f.readline().strip()
                if not tag: break
                seq = f.readline().strip().upper()
                f.readline()  # Ignore '+'
                qualities = f.readline().strip()
            t.write(f'{tag}\n{trim_seq_right(seq,nbases_trim,qualities=qualities)}')

os.system(f"bwa index {ref_genome}")
for trim_split in glob.glob(f'{save_split}*trimmed.{extension}'):
    os.system(f"bwa mem {ref_genome} {trim_split} > '{align_dir}{trim_split[9:-6]}.sam'")

shutil.rmtree(split_dir)


for align in glob.glob(f'{align_dir}/split*trimmed.sam'):
    with open(f'{align_dir}/merge.sam', 'at+') as m:
        with open(align) as a:
            while True:
                line = a.readline().strip()
                if not line:
                    break
                elif line[0] == '@':
                    continue
                else:
                    num_align += 1
                    m.write(f'{line}\n')

os.system(f"sort -k4 -n {align_dir}/merge.sam | sort -k3 > sorted_aligment.sam")

print(f"{num_align} reads have been aligned with the yeast genome")

shutil.rmtree(align_dir)