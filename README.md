# miFRED
**Microbial phenotype-based functional redundancy calculator from metagenomics data**: community-level functional redundancy (FREDc) and redundancy of 86 metabolic and ecological phenotypes (FREDs)

miFRED is a fast, user-friendly pipeline to calculate functional redundancy (FRED) using as functional units 86 phenotypes predicted by MICROPHERRET machine learning models from KEGG orthologs. miFRED simplifies the analysis by requiring only genomes or metagenomes fasta files, functional annotations as inputs and NGS read alignments as inputs. Designed for multiple-sample analysis, it enables direct comparisons across microbial communities in a single run. 
miFRED calculates both FREDc and various FREDs values. FREDc, based on Ricotta et al.,  assesses redundancy across samples using relative abundances of species and functional overlaps. FREDs metrics include the number of genomes performing each function in each sample (“number”), the fraction of community carrying out each function (“proportion”) and total relative abundance of these genomes (“relative abundance sum”). 
Outputs include several files storing the calculated FRED metrics, summary statistics and distribution plots. Intermediate results are also provided, including MICROPHERRET functional predictions, Jaccard distances between organisms and relative abundance values.


## **Installation**
1. Download the miFRED folder from this repository:
   
     ```git clone https://github.com/sfraulini/miFRED/```
   
2. Unzip training set files inside the Data folder
   
3. Download MICROPHERRET folder from github and move it inside miFRED/Data/:
   
     ```git clone https://github.com/BizzoTL/MICROPHERRET/```
  
4. miFRED requires an apposite conda environment, which can be generated as follow using the miFRED.yml file:

     ```conda env create -f miFRED_env.yml```


