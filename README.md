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

![Alt Text](Data/miFREDpipeline.png)

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

Examples of miFRED input and output files are available here: [Examples](https://mega.nz/folder/2JYhwCYJ#R94ItGZ8_L25P65fMBcalQ)
