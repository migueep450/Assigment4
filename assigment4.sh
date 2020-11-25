#!/bin/bash
# We create a workflows in which all the sequences has been analyzed from a fasta/q file. 
# Each sequence is going to be aligned with the yeast genome
# INPUT: a FASTA/Q file
# OUTPUT: a SAM file with the aligned sequences from the fasta/q file
# TO EXECUTE THE SCRIPT: bash assigment4.sh <input name-file.fasta/q>
# FUNCTIONS(4): From line 9 to 54

get_uniq_trimers () {
	#Concatenate the same file three times. In the second time, the first base of every sequence is going to be removed.
	#In the third time, the first and the second base of every sequence is going to be removed. 
	#Then, it extracts all the possible trimers of the file.
	cat $1 <(sed "s/^[ACGTNacgtn]//g" $1) <(sed "s/^[ACGTNacgtn][ACGTNacgtn]//g" $1) \
	| grep -v '^>' | fold -3 | awk '/[AGTCNagtcn]{3}/{print}' | sort | uniq -c
}

plot_hist () {
	#Plot all the trimers which are obteined from the input file.
	awk '{$1=sprintf("%-*s", $1, ""); gsub(" ", "=", $1); printf("%-10s%s\n", $2, $1)}' $1
}

fastq_to_fasta () {
	#The fastq file is converted into a fasta file
	#Expected input:
	#	$1 Filename of the fastq file
	awk 'NR%4 == 1 {print ">" substr($0, 2)} NR%4 == 2 {print}' $1
}

split_fastx_file () {
	#Separete all the records into differente files. One record per new file.
	#Expected input:
	#	$1 is the file that is going to be splitted and 
	#	$2 is EXT of the file
	if [ "$2" == "fasta" -o "$2" == "FASTA" ]; then
		# Each record has two lines
		NLINES=2
	else
		# Each record has four lines
		NLINES=4
	fi
	SPLIT_FILENAME="split_seq/split%d.$2"
	# Create the new files from the input file
	awk -v lines=$NLINES -v fmt=$SPLIT_FILENAME '{print>sprintf(fmt,1+int((NR-1)/lines))}' $1
}


FILENAME=$1
# Get FILENAME extension
EXT="${FILENAME##*.}"
NTRIMBASES=20
GENOMEFILE="ref/genome_yeast.fasta"

# It will store some files to execute the histogram counting, but it'll be remove afterwards
mkdir dummy 
# It will contain the splitted data before and after trimming
mkdir split_seq 
# Create a folder in which are going to be stored the generated SAM files
mkdir aligments
# This folder will contain the histogram and the aligments sorted by chromosome and position
mkdir output


if [ "$EXT" != "$FILENAME" ]; then
	# Check if the imput file is fasta or fastq.
	if [ "$EXT" == "fastq" -o "$EXT" == "FASTQ" ]; then
		echo "You introduced a FASTQ file."
		# Convert the fastq file into fasta file and then the fasta is stored in the dummy directory
		fastq_to_fasta $FILENAME > dummy/fastq_to_fasta.fasta
		DATA="dummy/fastq_to_fasta.fasta"
		DATA_EXT='FASTA'
	elif [ "$EXT" == "fasta" -o "$EXT" == "FASTA" ]; then
		echo "You introduced a FASTA file."
		DATA=$FILENAME
		DATA_EXT='FASTA'
	else
		echo "Invalid format! Ensure that your file has extension FASTA or FASTQ"
		rm -r dummy
		exit 1
	fi
	# Get all the trimers that have the sequences from the input file
	get_uniq_trimers $DATA > dummy/hist_data.txt
	# Plot the abundancy trimers which are obtained from the input file
	plot_hist dummy/hist_data.txt > output/hist_plot.txt
	echo "The trimers histogram have been saved in the folder as trimers_hist.txt"
	else
		echo "The file has no extension or has not been introduced"
	 	rm -r dummy
	 	exit 1
fi

rm -r dummy

# Separete all the records into different files
split_fastx_file $FILENAME $EXT

echo "$NTRIMBASES are going to be trimmed at the end of each splitted sequence."
echo "All the files are going to be saved in split_seq/splitXtrimmed.fastx"
# Iterate over all the splited files which have been created
for split in split_seq/split*.$EXT; do
	awk -v bases=$NTRIMBASES 'NR%2==1 {print $0}; NR%2==0 {print substr($0,1,length-bases)}' $split > "${split%.*}trimmed.$EXT"
done

echo "After splitting, we will align each split with the yeast reference genome"
# Initialize the bwa application with the reference genome
bwa index $GENOMEFILE
# Itereate in all the splitted and trimmed files.
for split in $(ls split_seq | grep -E 'split[0-9]+trimmed.'$EXT); do
	# Align all the sequences from the files and the reference genome, creating a SAM file with each alignment 
	bwa mem $GENOMEFILE split_seq/$split > "aligments/${split%.*}.sam"
done

rm -r split_seq

echo "Merging all the sam files"
samtools merge -f aligments/merge_out.sam aligments/split*.sam

echo "Removing headers"
grep -v '^@' aligments/merge_out.sam > aligments/clean_merge.sam

echo "Sorting the merged sam file by position and chromosome accession number"
sort -k4 -n aligments/clean_merge.sam | sort -k3 > output/sorted_aligment.sam

# How many sequence has been aligned against the reference genome
NALIGN=$(samtools view -c aligments/merge_out.sam) 
echo "$NALIGN reads have been aligned to the yeast genome"
rm -r aligments
