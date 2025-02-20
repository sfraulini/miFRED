# miFRED
**Microbial phenotype-based functional redundancy calculator from metagenomics data**: community-level functional redundancy (FREDc) and redundancy of 86 metabolic and ecological phenotypes (FREDs)

miFRED is a fast, user-friendly pipeline to calculate functional redundancy (FRED) using as functional units 86 phenotypes predicted by MICROPHERRET machine learning models from KEGG orthologs (KO). 

**Overview**
The pipeline processes metagenomic samples by taking as input:
- metagenomic reads from each sample
- individual fasta files representative of each genome in the samples.
Designed for multiple-sample analysis, it enables direct comparisons across microbial communities in a single run. 

miFRED can automatically generate any missing input file using built-in software if needed. Initial input generation steps include:
- genome functional annotation using eggNOG-mapper
- alignment of metagenomic reads to the genomes using Bowtie2
- sorted bam files generation with Samtools.
eggNOG-mapper annotation and alignment files can be provided directly as input. 

The core process employs MICROPHERRET function-specific ML models for phenotype prediction and calculates FRED.  For the acetoclastic methanogenesis phenotype, the refined version of the model previously provided by the authors was used. The obtained functional profiles are then scanned to generate a Jaccard distance matrix, representing the functional overlap within genome pairs. By default, the Jaccard distance between genomes without associated functions is set to 1. 
Alignment results are processed to calculate relative abundances, which, along with the Jaccard distance matrix and predicted functional profiles, are used to compute FREDc and FREDs metrics. 

miFRED calculates the following metrics:
- FREDc, based on Ricotta *et al.*
- FREDs "number": number of genomes performing each function in each sample
- FREDs "proportion": the fraction of community carrying out each function
- FREDs "relative abundance sum": total relative abundance of these genomes

Outputs include several files storing the calculated FRED metrics, summary statistics and distribution plots. Intermediate results are also provided, including MICROPHERRET functional predictions, Jaccard distances between organisms and relative abundance values.

## **Installation**
1. Download the miFRED folder from this repository:
   
     ```git clone https://github.com/sfraulini/miFRED/```
   
2. Unzip training set files inside the Data folder
   
3. Download MICROPHERRET folder from github and move it inside miFRED/Data/:
   
     ```git clone https://github.com/BizzoTL/MICROPHERRET/```
  
4. miFRED requires an apposite conda environment, which can be generated as follow using the miFRED.yml file:

     ```conda env create -f miFRED_env.yml```

   

Additional metrics like alpha diversity (Gini-Simpson index, GSI), a non-normalised version of FREDc [19] (FREDc_tian) and Rao’s entropy for functional diversity are also calculated [47]. Several parameters can be set to control the calculation, including the minimum fraction of genome with coverage higher than 0 (breadth of coverage) and the minimum relative abundance required to consider a genome as present in a sample; both these parameters are aimed at excluding spurious associations.

Users can define specific phenotypes to be included in the calculation, or by default the 86 phenotypes from high-performance models are considered. Annotation files are processed to extract KEGG Orthologs (KO), used by the models for the predictions.
