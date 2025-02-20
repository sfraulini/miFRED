# miFRED
**Microbial phenotype-based functional redundancy calculator from metagenomics data**: community-level functional redundancy (FREDc) and redundancy of 86 metabolic and ecological phenotypes (FREDs)

## **Overview**
miFRED is a fast, user-friendly pipeline to calculate functional redundancy (FRED) using as functional units 86 phenotypes predicted by MICROPHERRET machine learning models from KEGG orthologs (KO). 

The pipeline processes metagenomic samples by taking as input:
- metagenomic reads from each sample
- individual fasta files representative of each genome in the samples.
Designed for multiple-sample analysis, it enables direct comparisons across microbial communities in a single run. 

miFRED can automatically generate any missing input file using built-in software if needed. Initial input generation steps include:
- genome functional annotation using eggNOG-mapper
- alignment of metagenomic reads to the genomes using Bowtie2
- sorted bam files generation with Samtools.
eggNOG-mapper annotation and alignment files can be provided directly as input. 

The core process uses MICROPHERRET function-specific ML models for phenotype prediction and calculates FRED.  For the acetoclastic methanogenesis phenotype, the refined version of the model previously provided by the authors was used. The obtained functional profiles are scanned to generate a Jaccard distance matrix, representing the functional overlap within genome pairs. By default, the Jaccard distance between genomes without associated functions is set to 1. 
Alignment results are processed to calculate relative abundances, which, along with the Jaccard distance matrix and predicted functional profiles, are used to compute FREDc and FREDs metrics. 

miFRED calculates the following metrics:
- **FREDc**, based on Ricotta *et al.*
- **FREDs "number"**: number of genomes performing each function in each sample
- **FREDs "proportion"**: the fraction of community carrying out each function
- **FREDs "relative abundance sum"**: total relative abundance of these genomes

Outputs include several files storing the calculated FRED metrics, summary statistics and distribution plots. Intermediate results are also provided, including MICROPHERRET functional predictions, Jaccard distances between organisms and relative abundance values.

## **Installation**
1. Download the miFRED folder from this repository:
   
     ```git clone https://github.com/sfraulini/miFRED/```
   
2. Unzip training set files inside the Data folder
   
3. Download MICROPHERRET folder from github and move it inside miFRED/Data/:
   
     ```git clone https://github.com/BizzoTL/MICROPHERRET/```
  
4. miFRED requires an apposite conda environment, which can be generated as follow using the miFRED.yml file:

     ```conda env create -f miFRED_env.yml```

## **Requirements**
miFRED requires a Linux system and at least 6 CPUs due to the computational demands of the MICROPHERRET training process.

## **Usage**
Below are different command examples depending on the inputs provided by the user. The fewer inputs the user provides, the more miFRED will generate automatically.

- For complete input generation, i.e. user provides only the FASTA files folder and metagenomic reads:

   ```python3 miFRED_core.py -g GENOMES_FOLDER -r READS_FOLDER -u {True,False} -A eggnog_annotation -o output_folder```

- If alignment of metagenomics reads against the provided genomes was already performed by the user:

   ```python3 miFRED_core.py -g GENOMES_FOLDER -B BAM_FILES -A eggnog_annotation -o output_folder```

- To compute FRED from KO for comparison analysis::

   ```python3 miFRED_core.py -g GENOMES_FOLDER -B BAM_FILES -A eggnog_annotation -o output_folder --KO```

- For help type the following or look at the Manual:

   ```python3 miFRED_core.py --help```

## **Manual**
### **Inputs**
??????
- Several parameters can be set to control the calculation, including the minimum fraction of genome with coverage higher than 0 (breadth of coverage) and the minimum relative abundance required to consider a genome as present in a sample; both these parameters are aimed at excluding spurious associations.

Users can define specific phenotypes to be included in the calculation, or by default the 86 phenotypes from high-performance models are considered. Annotation files are processed to extract KEGG Orthologs (KO), used by the models for the predictions.
????????

### **Outputs**
The following folders are generated in the output folder specified by the user.

#### input_fred folder
It stores the results of miFRED's input generation step procedure. The list can change depending on which files were already provided by the user.
- all_genomes.fa: fasta file obtained by concatenating the files provided by the user. It is the reference used in the alignment procedure
- info.txt: ?
- nreads.txt: ?
- bowtie2 folder: sorted and indexed bam files obtained by aligning reads against all_genomes.fa with bowtie2 and processing resulting files with Samtools
- eggnog_annotations folder: contains eggNOG-mapper results
- annotation_matrix.csv: ?
- KO_for_calculation.csv: ?
  
#### output_micropherret folder
If MCROPHERRET models are used, it stores the results of the ML predictions:
- predict_functions.csv: matrix with predicted functions per genome
- predict_sum.csv: number of genomes predicted to perform each function
  
#### output_fred folder
It stores results of FRED calculation procedure:
- fredc.csv file: store results of FREDc calculation for each analysed sample. Additional metrics like alpha diversity (Gini-Simpson index, GSI), a non-normalised version of FREDc(FREDc_tian) and Rao’s entropy for functional diversity are included too.
- freds.csv: store results of FREDs calculation for each function for each analysed sample
- FREDs_statistic.csv: FREDs statistics for each analysed function and sample

Intermediate results are also provided:
- jaccard_distances.csv: storing the pairwise functional diversity within each genome pair
- used_relative_abundances.csv: relative abundances computed by miFRED and used for FRED calculation
- normalised_relative_abundances.csv: relative abundances computed by miFRED, normalised by the number of mapped reads

FREDc and FREDs distribution plot are provided too to aid the analysis. 
