#!/bin/bash

handle_error() {
    echo "Error: $1"
    exit 1
}

#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh

input_genomes="$1"
output_dir="$2"
extension_genomes="$3"
proc="$4"
database="$5"
sensmode="$6"

dir=$(pwd)

# Check if the directory exists
if [ ! -d "$input_genomes" ]; then
    handle_error "Directory '$input_genomes' does not exist."
fi

#create input_fred directory
if [ ! -d "$output_dir/input_fred" ]; then
    mkdir "$output_dir/input_fred"
    echo "Directory $output_dir/input_fred created successfully"
fi 

files=("$input_genomes"/*"$extension_genomes")
if [ ! -e "${files[0]}" ]; then
    handle_error "No files with extension '$extension_genomes' found in '$input_genomes'."
fi

#conda activate eggnog-mapper_2.1.10

if [ ! -d "$output_dir/input_fred/eggnog_annotations" ]; then
    mkdir  "$output_dir/input_fred/eggnog_annotations"
fi 

cd $input_genomes/

for f in $(ls *); do 
    echo "Annotating $f..."
    name=${f%$extension_genomes*}
    emapper.py -m diamond --itype genome --genepred prodigal  --cpu $proc -i $f -o $name --output_dir $output_dir/input_fred/eggnog_annotations --sensmode $sensmode --data_dir $database --override 
done

rm -r $input_genomes/emappertmp_*

cd $dir
#conda deactivate
