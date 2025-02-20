#!/bin/bash

handle_error() {
    echo "Error: $1"
    exit 1
}

check_sorted_bam() {
    for bam_file in $1/*.bam; do
        if [[ "$bam_file" == *sorted.bam ]]; then
            index_file="${bam_file}.bai"
            if [[ ! -f "$index_file" ]]; then
                handle_error "Missing index file for: $bam_file"
            fi
        else 
            handle_error "$bam_file is not sorted"
        fi
    done
    return
}

#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh

bam_folder="$1"
proc="$2"
output_dir="$3"

check_sorted_bam "$bam_folder"

if [ -e "$output_dir/input_fred/nreads.txt" ]; then
    rm $output_dir/input_fred/nreads.txt
fi
#conda activate samtools
for bam in $bam_folder/*sorted.bam; do
    nreads=$(samtools view -@ proc -c $bam)
    echo "Processing $bam ..."
    echo "$bam $nreads" >> $output_dir/input_fred/nreads.txt
done

#conda deactivate