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

#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh

#input folders
input_genomes="$1"
output_dir="$2"
extension_genomes="$3"
proc="$4"

# Check if the directory exists
if [ ! -d "$input_genomes" ]; then
    handle_error "Directory '$input_genomes' does not exist."
fi

#create input_fred directory
if [ ! -d "$output_dir/input_fred" ]; then
    echo "Directory $output_dir/input_fred created successfully"
    mkdir "$output_dir/input_fred"
fi 

#check if files are compressed
decompress "$input_genomes"

#generate concatenated fasta file
files=("$input_genomes"/*"$extension_genomes")
if [ ! -e "${files[0]}" ]; then
    handle_error "No files with extension '$extension_genomes' found in '$input_genomes'."
fi

echo "Creating all_genomes.fa..."
cat "$input_genomes"/*"$extension_genomes" > "$output_dir/input_fred/all_genomes.fa"  || handle_error "Failed to concatenate files."

#create info.txt
echo "Creating binning file info.txt..."

for file in "$input_genomes"/*; do
    if [[ "$file" != *"$extension_genomes" ]]; then
        handle_error "Found file without the extension '$extension_genomes': $file"
    fi
done
python ./inputwriter.py "$input_genomes" "$output_dir/input_fred/info.txt"

if [ $(grep '>' "$output_dir/input_fred/all_genomes.fa" | wc -l) -ne $(wc -l < "$output_dir/input_fred/info.txt") ]; then
    handle_error "Number of detected contings does not match contings in '$input_genomes'."
fi