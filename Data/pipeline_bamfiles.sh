#!/bin/bash

handle_error() {
    echo "Error: $1"
    exit 1
}

is_compressed_by_content() {
    if file "$1" | grep -qE 'gzip compressed|bzip2 compressed|XZ compressed data|Zip archive'; then
        return 0
    else
        return 1
    fi
}

decompress(){ #takes folder as input
    #check if files are compressed
    for file in $(ls "$1"); do
        if is_compressed_by_content "$1/$file" ; then
            echo "Decompressing $1/$file."
            pigz -d -p $proc "$1/$file"
        fi
    done
    return
}

reads_checks(){
    R1_array=($(ls "$1"*_R1* 2> /dev/null | grep -v '_unpaired' | sed 's/_R1*//'))
    R2_array=($(ls "$1"*_R2* 2> /dev/null | grep -v '_unpaired' | sed 's/_R2*//'))
    # Check if both lists are the same
    if [[ "${#R1_array[@]}" -ne "${#R2_array[@]}" ]]; then
        handle_error "The number of R1 and R2 files does not match."
    fi
    # Compare each element of R1_array and R2_array
    for sample in "${R1_array[@]}"; do
        if [[ ! " ${R2_array[@]} " =~ " $sample " ]]; then
            handle_error "Sample '$sample' has an R1 file but no corresponding R2 file."
        fi
    done
    for sample in "${R2_array[@]}"; do
        if [[ ! " ${R1_array[@]} " =~ " $sample " ]]; then
            handle_error "Sample '$sample' has an R2 file but no corresponding R1 file."
        fi
    done
    echo "All files are properly paired."
    return 
}

unpaired_reads_checks (){
    if [ $2 == 'True' ]; then
        R1s_paired=($(ls "$1"*_R1* 2> /dev/null | grep -v '_unpaired'))
        R1s_unpaired=($(ls "$1"*_R1* 2> /dev/null | grep '_unpaired'))
        R2s_unpaired=($(ls "$1"*_R2* 2> /dev/null | grep '_unpaired')) 
        if [[ "${#R1s_unpaired[@]}" == "${#R1s_paired[@]}"  && "${#R2s_unpaired[@]}" == "${#R1s_paired[@]}" ]]; then
            echo "Unpaired reads were correctly provided"
        else
            handle_error "Unpaired reads not provided for all the samples."
        fi
    fi
    return
}


#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh

#input folders
concatenate_genomes="$1"
input_reads="$2"
output_dir="$3"
proc="$4"
unpaired="$5"

dir=$(pwd)

if [ ! -d "$input_reads" ]; then
    handle_error "Directory '$input_reads' does not exist."
fi

#create input_fred directory
if [ ! -d "$output_dir/input_fred" ]; then
    mkdir "$output_dir/input_fred"
fi 

#create bowtie2 directory
if [ ! -d "$output_dir/input_fred/bowtie2" ]; then
    mkdir "$output_dir/input_fred/bowtie2"
fi 

#check if files are compressed
decompress "$input_reads"

#check reads files
reads_checks "$input_reads"

# Check for unpaired reads in the folder
unpaired_reads_checks "$input_reads" "$unpaired"

#conda activate bowtie2
echo Building index from provided genomes...
bowtie2-build $concatenate_genomes $output_dir/input_fred/bowtie2/index -p $proc -q

echo Alignment step...
cd $input_reads
for f in $(ls *"_R1"* 2> /dev/null | grep -v '_unpaired'); do
    sample=${f%_R1*}
    r=${f%_R1*}_R2${f#*_R1}
    if [ $unpaired == 'False' ]; then
        echo "Aligning $sample with paired reads only..." 
        bowtie2 -q -1 $f -2 $r -x $output_dir/input_fred/bowtie2/index -p $proc -S $output_dir/input_fred/bowtie2/$sample.sam 
    else
        R1s_unpaired=($(ls *_R1* 2> /dev/null | grep '_unpaired'))
        #R2s_unpaired=($(ls "$input_reads"*_R2* 2> /dev/null | grep '_unpaired'))
        for r1_un in "${R1s_unpaired[@]}"; do
            if [[ "$r1_un" == $sample* ]]; then
                uf=$r1_un #${f%_R1*}_R2${f#*_R1}
                ur=${r1_un%_R1*}_R2${r1_un#*_R1}
            fi
        done
        echo "Aligning $sample ..." 
        bowtie2 -q -1 $f -2 $r -U $uf,$ur -x $output_dir/input_fred/bowtie2/index -p $proc -S $output_dir/input_fred/bowtie2/$sample.sam
    fi
done 

cd $dir

#conda deactivate 

#conda activate samtools 
echo "Creating sorted.bam files..."
for sam in $(ls $output_dir/input_fred/bowtie2/*sam); do
    sample=${sam%.s*}
    echo "Processing $sample ..."
    samtools view -@ $proc $sam -b > $sample.bam
    samtools sort -@ $proc $sample.bam > $sample.sorted.bam 
    samtools index -@ $proc $sample.sorted.bam > $sample.index.bam 
done 

#conda deactivate

rm $output_dir/input_fred/bowtie2/index*
rm $output_dir/input_fred/bowtie2/*sam
for file in $output_dir/input_fred/bowtie2/*bam; do
    # Exclude files ending with sorted.bam
    if [[ "$file" != *sorted.bam ]]; then
        echo "Deleting: $file"  # For testing, you can remove this echo when you're sure
        rm "$file"
    fi
done 
