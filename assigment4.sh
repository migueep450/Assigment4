#!/bin/bash

# This script aligns the input file with the yeast genome and returns merged sam files

FILENAME=$1
EXT="${FILENAME##*.}"
NTRIMBASES=20
GENOMEFILE="ref/genome_yeast.fasta"

mkdir dummy # It will store some files to execute the histogram counting, but it'll be remove afterwards
mkdir split_seq # It will contain the splitted data before and after trimming
mkdir aligments

get_uniq_trimers () {
	cat $1 <(sed "s/^[ACGTNacgtn]//g" $1) <(sed "s/^[ACGTNacgtn][ACGTNacgtn]//g" $1) \
	| grep -v '^>' | fold -3 | awk '/[AGTCNagtcn]{3}/{print}' | sort | uniq -c
}

plot_hist () {
	awk '{$1=sprintf("%-*s", $1, ""); gsub(" ", "=", $1); printf("%-10s%s\n", $2, $1)}' $1
}

fastq_to_fasta () {
	echo '$PATH'
	awk 'NR%4 == 1 {print ">" substr($0, 2)} NR%4 == 2 {print}' split_seq/$1
}


split_fastx_file () {
	# $1 is the file that is going to be splitted and $2 is EXT of the file
	if [ "$2" == "fasta" -o "$2" == "FASTA" ]; then
		NLINES=2
	else
		NLINES=4
	fi
	SPLIT_FILENAME="split_seq/split%d.$2"
	awk -v lines=$NLINES -v fmt=$SPLIT_FILENAME '{print>sprintf(fmt,1+int((NR-1)/lines))}' $1
}


if [ "$EXT" != "$FILENAME" ]; then
	if [ "$EXT" == "fastq" -o "$EXT" == "FASTQ" ]; then
		echo "You introduced a FASTQ file."
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
	get_uniq_trimers $DATA > dummy/hist_data.txt
	plot_hist dummy/hist_data.txt > hist_plot.txt
	echo "The trimers histogram have been saved in the folder as trimers_hist.txt"
	else
		echo "The file has no extension or has not been introduced"
	 	rm -r dummy
	 	exit 1
fi

rm -r dummy

split_fastx_file $FILENAME $EXT

echo "$NTRIMBASES are going to be trimmed at the end of each splitted sequence."
echo "All the files are going to be saved in split_seq/splitXtrimmed.fastx"
for split in split_seq/split*.$EXT; do
	awk -v bases=$NTRIMBASES 'NR%2==1 {print $0}; NR%2==0 {print substr($0,1,length-bases)}' $split > "${split%.*}trimmed.$EXT"
done

echo "After splitting, we will align each split with the yeast reference genome"
bwa index $GENOMEFILE
for split in $(ls split_seq | grep -E 'split[0-9]+trimmed.'$EXT); do
	bwa mem $GENOMEFILE split_seq/$split > "aligments/${split%.*}.sam"
done
rm -r split_seq

echo "Merging all the sam files"
samtools merge -f aligments/merge_out.sam aligments/split*.sam
echo "Removing comments"
grep -v '^@' aligments/merge_out.sam > aligments/clean_merge.sam

echo "Sorting the merged sam file by position and chromosome accession number"
sort -k4 -n aligments/clean_merge.sam | sort -k3 > sorted_aligment.sam

NALIGN=$(samtools view -c aligments/merge_out.sam) 
echo "$NALIGN reads have been aligned to the yeast genome"
rm -r aligments