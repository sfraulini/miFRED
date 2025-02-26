#!/usr/bin/env python

import os
import argparse
import subprocess
import math
import pandas as pd
import pysam
from pickle import load
import warnings
import re
import matplotlib.pyplot as plt
import seaborn as sns

from multiprocessing import Pool
from tqdm import tqdm

warnings.filterwarnings("ignore", category=FutureWarning) 

#########################################################################################################################################

def parse_args():
    parser = argparse.ArgumentParser(
                    prog='miFRED',
                    description='Microbial functional redundancy calculator from metagenomics data: community-level functional redundancy (FREDc) and redundancy of 86 metabolic and ecological phenotypes (FREDs)',
                    epilog='Developed by Metabioinfomics Lab, University of Padova')
    group_genomes = parser.add_mutually_exclusive_group(required = True)
    group_genomes.add_argument('-g', '--genomes_folder', help = 'Directory where genomes fastq files are stored') 
    group_genomes.add_argument('--all_genomes', help = 'Multifasta file obtained concatenating all single fastq files')
    parser.add_argument('--binning_file', help = '.txt file with each line listing a scaffold and the corresponding genome/bin name, tab-seperated')
    group = parser.add_mutually_exclusive_group(required = True)
    group.add_argument('-r', '--reads_folder', help = 'Directory where metagenomic reads fastq files are stored')
    group.add_argument('-B', '--bam_files', help = 'Directory where sample-specific sorted.bam files and indexes are stored')
    parser.add_argument('-u', '--unpaired_reads', choices=['True', 'False'], help = 'True if fasta file for unpaired reads are also provided in reads_folder, False otherwise (default : True)')
    parser.add_argument('-o', '--output_folder', required = True, help = 'Output directory')
    parser.add_argument('-x', '--genomes_extension', default = '.fa', help = 'Genome files extension')
    parser.add_argument('-p', '--processors', default = 5, type = int, help = 'threads (default: 5)')
    group_annotations =  parser.add_mutually_exclusive_group(required = True)
    group_annotations.add_argument('-A', '--eggnog_annotation', help = 'Directory where genomes eggNOG .annotations files are stored. .csv file obtained from parsing .annotations files can be directly provided, with genomes as rows, KO as columns and KO counts as values. First column must be named "Genomes')
    group_annotations.add_argument('-db', '--eggnog_database', help = 'Directory where eggNOG-mapper database is stored')
    parser.add_argument('-sm', '--eggnog_sensmode', default = 'sensitive', choices = ['default', 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive'], help = 'eggNOG-mapper Diamond search option: either default, fast, mid-sensitive, sensitive, more-sensitive, very-sensitive or ultra-sensitive. (default: sensitive)')
    parser.add_argument('-f', '--functions_list', help = '.txt file containing MICROPHERRET functions to be considered for the calculation, one per line. (default: 86 functions whose models were accurate on test set, stored in functions.txt)')
    parser.add_argument('-s', '--training_sets', default = 'training_sets', help = 'Folder containing the dataset.csv and dataset_acetoclastic_methanogenesis.csv files to be used in the training. (default: ./training_sets/ )')
    parser.add_argument('-c', '--covered_genome_fraction', default= 0.1, help = 'Genomes with a fraction of covered bases lower than this are reported as having zero coverage. (default: 0.10)')
    parser.add_argument('-t', '--relative_abundance_threshold', default= 0, help = 'Mininum relative abundance threshold to be considered as present in the sample. (default: 0)')
    group_ko_vs_pred = parser.add_mutually_exclusive_group()
    group_ko_vs_pred.add_argument('-m', '--micropherret_predictions', help = 'csv file containing MICROPHERRET predictions for all the genomes, the first column with genomes names must be unnamed.')
    group_ko_vs_pred.add_argument('-k', '--KO', action='store_true', help = 'perform calcalution based on KO')
    args = parser.parse_args()
    
    # Post-parsing validation
    if args.reads_folder and not args.unpaired_reads:
        raise parser.error("Argument -u/--unpaired_reads is required when -r/reads_folder is specified.")
    
    if args.all_genomes and not args.binning_file:
        raise parser.error("Argument --binning_file is required when --all_genomes is specified.")

    return args

########################################################################################################################################
#Generate annotation matrix - Python version of MICROPHERRET Jupyter script get_annotation_matrix.ipynb

def get_file(path):
    data_list = []
    data =  {}
    comments = []
    with open(path,'r') as file:
        for line in file:
            line = line.strip().split('\t')
            if line[0].startswith('#'): comments.append(line[0])
            data_list.append(line)
    if len(comments) != len(set(comments)):
        indexes = []
        for i in range(len(data_list)):
            if comments[1] in data_list[i]:
                indexes.append(i)
        index_to_divide = indexes[-1] -1
    else:
        index_to_divide = 0
    for line in data_list[index_to_divide:]:
        if line[0].startswith('#'):continue
        data[line[0]] = line[1:]
    return data

def get_kos(path):
    ind_b = path.rfind('/')
    ind = path.index('.emapper') 
    genome = path[ind_b+1:ind]
    annotated_file = get_file(path)
    ko_list = []
    for query in annotated_file.keys():
        ko = annotated_file[query][10]
        if ko != '-' and ',' not in ko:
            ko_list.append(ko)
        elif ',' in ko:
            ko_list += ko.split(',')
    ko_list = [k[3:] for k in ko_list]
    genome_data_ko = {}
    for ko in ko_list:
        genome_data_ko[ko] = ko_list.count(ko)
    return genome, genome_data_ko, len(set(ko_list)) 

#########################################################################################################################################
#Launching MICROPHERRET - Portion of predict_functions.py script

def get_validation_set(to_validate, training_set):
    training_kos = training_set.columns
    validation_kos = to_validate.columns
    if len(set(training_kos).intersection(set(validation_kos))) == 0:
        print('No common KOs between provided ones and training')
        return
    else:
        common = list(set(training_kos).intersection(set(validation_kos)))
        print('{} common KOs'.format(len(common)))
        common_table = to_validate[common]
        #remove orthologs in validation not in the training
        to_remove = set(validation_kos) - set(training_kos)
        print('{} KOs present in user set but not in training set will be removed'.format(len(to_remove)))
        missing = list(set(training_kos) - set(validation_kos))
        print('{} KOs missing from the users et will be add to train the classifiers'.format(len(missing)))
        missed = pd.DataFrame(0, index = to_validate.index, columns= missing)
        to_submit = common_table.merge(missed, left_index = True, right_index = True)
        to_submit = to_submit[list(training_kos)] 
        print('Shape of training dataset: {}, Shape of user dataset: {}'.format(training_set.shape, to_submit.shape))
        if list(to_submit.columns) != list(training_set.columns): 
            print('Error columns ML input matrix --- functions error')
            return
    return to_submit

def read_functions(classes_file):
    classes = []
    with open(classes_file, 'r') as file:
        for line in file:
            if line == ' ':
                continue
            classes.append(line[:-1])
    return classes

def validate(classes, ko_validation, /, function_validation = 0):
    results_per_class = {}
    scores = {}
    all_functions_list = []
    with open('./functions.txt', 'r') as file:
        for line in file:
            if line == ' ': continue
            all_functions_list.append(line[:-1])
    for c in tqdm(range(len(classes)), desc="Predicting functions..."):
        if classes[c] not in all_functions_list:
            print('ERROR: '+classes[c]+' not in the available functions!')
            return
        if classes[c] == 'acetoclastic_methanogenesis': 
            model = load(open('./MICROPHERRET/saved_models/new_model_acetoclastic_methanogenesis_2.sav', 'rb'))
            scaler = load(open('./MICROPHERRET/saved_models/new_scaler_acetoclastic_methanogenesis_2.sav', 'rb'))
            to_validate_norm = scaler.transform(ko_aceto)
            pred = model.predict(to_validate_norm)
        else:
            model = load(open('./MICROPHERRET/saved_models/model_'+classes[c]+'.sav', 'rb'))
            scaler = load(open('./MICROPHERRET/saved_models/scaler_'+classes[c]+'.sav', 'rb'))

            to_validate_norm = scaler.transform(ko_validation)
            pred = model.predict(to_validate_norm)

        results_per_class[classes[c]] = pred

        if type(function_validation) != int:
            scores[classes[c]] = [matthews_corrcoef(function_validation[classes[c]], pred), f1_score(function_validation[classes[c]], pred, zero_division=1), 
                        confusion_matrix(function_validation[classes[c]], pred), accuracy_score(function_validation[classes[c]], pred),
                        hamming_loss(function_validation[classes[c]], pred), zero_one_loss(function_validation[classes[c]], pred),  function_validation[classes[c]].sum()]
            results_df = pd.DataFrame(results_per_class, index = validation_set.index)
    results_df = pd.DataFrame(results_per_class, index = validation_set.index)
    results_df.to_csv(os.path.join(output_folder, 'output_micropherret/predict_functions.csv'))
    sums = results_df.sum()
    sums.to_csv(os.path.join(output_folder,'output_micropherret/predict_sum.csv'))
    return results_df, scores

#########################################################################################################################################
#FRED calculation
#########################################################################################################################################
#Process functions
def index_function(functions):
    #take as input txt file with the list of functions that were predicted, one per line
    print('Processing functions...')
    function_index = {}
    with open(functions, 'r') as file:
        i = 0
        for line in file:
            function_index[line[:-1]] = i
            i+=1
    print('Done')
    return function_index
####################################################################
#Process genomic information
def read_binning_file(file):
    # generation of genome_listofcontig and contig_genome dictionary from input.txt
    f= open(file) 
    genome_contigs={}
    contig_genome={}
    for line in f:
        a=line.strip().split('\t')
        contig_genome[a[0]]=a[1]
        try:
            genome_contigs[a[1]].append(a[0])
        except KeyError:
            genome_contigs[a[1]]=[a[0]]
    f.close()
    return genome_contigs

def read_fasta(fasta):
    # calculating contig lengths from fasta file and saving them into a dictionary
    print('Reading fasta file...') 
    f=open(fasta)
    contigs_lengths_d={}
    contig_name=''
    for line in f:
        if line[0]=='>':
            header=line[1:].strip().split()
            contig_name=header[0]
            contigs_lengths_d[contig_name]=0
        else:
            contigs_lengths_d[contig_name]+=len(line.strip())
    f.close()
    print('Done')
    return contigs_lengths_d

def index_genome(genome_contig):
    genome_index={}
    i=0    
    for genome in genome_contig:
        genome_index[genome]=i
        i+=1
    return(genome_index)

def get_reads_number(reads_file):
    per_samples = {}
    with open(reads_file, 'r') as file:
        for line in file:
            line = line.strip().split(' ')
            per_samples[line[0]] = int(line[1])
    return per_samples 
####################################################################
#Calculate Jaccard distance for FRED
def jaccard_distance(gen_function_1or0,genome_index):
    not_functions = set()
    gen_function_1or0 = gen_function_1or0.T[list(genome_index)].to_dict()
    genome_jds={}
    for gen1 in gen_function_1or0:
        genome_jds[gen1]=[0]*len(genome_index)
        function_gen1=gen_function_1or0[gen1]
        for gen2 in gen_function_1or0:
            if gen1==gen2:
                genome_jds[gen1][genome_index[gen2]]=0
                continue
            function_gen2=gen_function_1or0[gen2]
            intersection=0
            union=0   
            for i in function_gen1:
                if function_gen1[i]== 1 and function_gen2[i]==1:
                    intersection+=1
                elif function_gen1[i] == 0 and function_gen2[i] == 0:
                    continue
                union +=1
            if union == 0:
                genome_jds[gen1][genome_index[gen2]] = 1
                not_functions.add(gen1)
                not_functions.add(gen2)
            else:
                genome_jds[gen1][genome_index[gen2]]= 1- (intersection/union) 
    return(genome_jds)
####################################################################
#Coverage calculation and intermediate steps
def get_genome_length(genome_contigs, contigs_lengths_d):
    # calculation of each genome length for normalization purposes
    genomes=[]     
    total_length_allgenomes= 0 
    genome_length={} 
    for genome in genome_contigs:
        genomes.append(genome)
        genome_length[genome]=0
        for contig in genome_contigs[genome]:
            genome_length[genome]+=contigs_lengths_d[contig]
        total_length_allgenomes+=genome_length[genome]
    return genomes, total_length_allgenomes, genome_length

def get_coverage_with_threshold(genome_contigs, contigs_lengths_d, bam, reads_number_sample, frac = 0, thr_min = 0):
    print(f'Calculating coverage for {bam}...')
    genomes, total_length_allgenomes, genome_length = get_genome_length(genome_contigs, contigs_lengths_d) 
    sam=pysam.AlignmentFile(bam,"rb")
    coverages = {}
    mapped = 0
    sum_coverages = 0
    mapped_per_genome = {}
    covered_fractions = {}
    rel_abs = {}
    rel_abs_norm = {}

    for g in genomes:
        g_size = genome_length[g]
        matches = 0
        reads_mapped_to_genome = 0
        covered_genome = 0
        for contig in genome_contigs[g]:
            covered_bases = set()
            try:
                for line in sam.fetch(contig=contig):
                    if line.flag & 0x4: 
                        continue # Check bit 2 in the flag
                    else:
                        reads_mapped_to_genome += 1
                        cig = re.findall('\d+(?=M)', line.cigarstring)       
                        cigars = sum([int(i) for i in cig])
                        matches += cigars
                        for pos in range(line.reference_start, line.reference_end):
                            covered_bases.add(pos)
            except ValueError:
                pass
            covered_genome+= len(covered_bases)
        mapped_per_genome[g] = reads_mapped_to_genome
        covered_fraction = covered_genome/g_size
        covered_fractions[g] = covered_fraction

        if covered_fraction < frac:
            coverage = 0
            sum_coverages += coverage
            coverages[g] = coverage
        else:
            coverage = round(matches/g_size, 8)
            sum_coverages += coverage
            coverages[g] = coverage
            mapped += reads_mapped_to_genome

    for g in coverages:
        if sum_coverages != 0:
            rel_abs_norm[g] = coverages[g]/sum_coverages * (mapped/reads_number_sample) 
            rel_abs[g] = coverages[g]/sum_coverages 
        else:
            rel_abs[g] = 0 
            rel_abs_norm[g] = 0  

    if thr_min != 0:
        for g in rel_abs:
            if rel_abs[g] < thr_min:
                rel_abs[g] = 0
    print('Done')

    return pd.Series(rel_abs_norm), pd.Series(rel_abs)

####################################################################
#FREDs calculation
def single_functions(function_index,genome_index,gen_function_1or0,gen_ra,genome_jd,total_jd):
    gen_function_1or0 = gen_function_1or0.T.to_dict()
    topr={'Function': [], 'Number':[], 'Proportion':[], 'Relative abundance sum':[],'tShannon index':[]} 
    genomes_present = [gen for gen in gen_function_1or0 if gen_ra[gen] > 0]

    for function in function_index:
        topr['Function'].append(function)
        if len(genomes_present) == 0:
            for measure in topr.keys():
                if measure != 'Function':
                    topr[measure].append(0)
        else:
            gen_with_function=[] 
            abundances=[]  
            abundances_all = []
            for gen in genomes_present:
                abundances_all.append(gen_ra[gen])
                gen_funct_profile=gen_function_1or0[gen]
                if gen_funct_profile[function]==1:
                    gen_with_function.append(gen)
                    abundances.append(gen_ra[gen])
            topr['Number'].append(len(gen_with_function))
            topr['Proportion'].append(len(gen_with_function)/len(genomes_present))

            total_abundance=sum(abundances)
            topr['Relative abundance sum'].append(total_abundance)
            num_shan=0
            for ab in abundances: 
                num_shan+= -1*ab/total_abundance*math.log(ab/total_abundance)
            if len(abundances)<=1:
                topr['tShannon index'].append(0)
            else:
                topr['tShannon index'].append(num_shan/math.log(len(genomes_present)))

    return topr 

#FREDc calculation
def all_functions(gen_function_1or0,genome_jds,gen_ra,genome_index,function_index):
    topr={} 
    value_tian=0
    num_ricotta=0
    den_ricotta=0
    rao = 0
    gsi = 0
    sum_gen_ra = pd.Series(gen_ra).sum()

    for gen_i in genome_jds:
        gen_i_jds=genome_jds[gen_i]
        gen_i_ra=gen_ra[gen_i]
        den_ricotta+=gen_i_ra*(sum_gen_ra-gen_i_ra)
        for gen_j in genome_jds:
            if gen_i==gen_j:
                continue
            num_ricotta+=gen_i_ra*gen_ra[gen_j]*gen_i_jds[genome_index[gen_j]]
            value_tian+=(1-gen_i_jds[genome_index[gen_j]])*gen_i_ra*gen_ra[gen_j]
            rao += gen_i_jds[genome_index[gen_j]]*gen_i_ra*gen_ra[gen_j]
            gsi += gen_i_ra*gen_ra[gen_j]
    topr['GSI'] = gsi
    topr['Rao'] = rao
    topr['FREDc_tian'] = value_tian
    if value_tian == 0:
        topr['FREDc'] = 0
    else:
        topr['FREDc'] = 1-(num_ricotta/den_ricotta)
    return topr

####################################################################
#Generate matrix for FREDs plot
def get_matrix(data, metric):
    for_pca_sumra = {}
    for_pca_sumra_red = {}
    for f in set(data.index):
        for_pca_sumra[f] = {}
        for_pca_sumra_red[f] = {}
        for s in set(data['Sample']):
            for_pca_sumra[f][s] = data.loc[f].loc[data.loc[f]['Sample'] == s][metric].item()
            if data.loc[f].loc[data.loc[f]['Sample'] == s]['Number'].item() <= 1: 
                for_pca_sumra_red[f][s] = 0
            else:
                for_pca_sumra_red[f][s] = data.loc[f].loc[data.loc[f]['Sample'] == s][metric].item()
    for_pca_sumra = pd.DataFrame(for_pca_sumra)
    for_pca_sumra_red = pd.DataFrame(for_pca_sumra_red)
    return for_pca_sumra, for_pca_sumra_red

####################################################################
def bams_analysis(bam, contigs_lengths_d, genome_contigs, function_index, gen_ko_1or0, genome_index, genome_jds, total_jd, reads_number_sample, output_dir,thr = 0.1, thr_min = 0):
    ind = bam.index('.sorted')
    bam_name = bam[:ind] #occhio all'output folder
    name_to_save = bam_name[bam_name.rfind('/'):]
    if name_to_save.startswith('/'): 
        name_to_save = name_to_save[1:]
    rel_ab_norm, rel_ab =get_coverage_with_threshold(genome_contigs, contigs_lengths_d, bam, reads_number_sample, thr, thr_min)
    rel_ab_norm = pd.Series(rel_ab_norm, name = name_to_save).sort_values()
    rel_ab_norm.to_csv(os.path.join(output_dir, name_to_save + '_normalised_ra.csv'))
    rel_ab = pd.Series(rel_ab, name = name_to_save).sort_values()
    rel_ab.to_csv(os.path.join(output_dir, name_to_save + '_used_ra.csv'))

    print(f'Calculatiing FREDs for {bam}...')
    freds = single_functions(function_index,genome_index,gen_ko_1or0,rel_ab,genome_jds,total_jd)
    freds = pd.DataFrame(freds).set_index('Function')
    freds['Sample'] = pd.Series(name_to_save, index = freds.index)
    freds.to_csv(os.path.join(output_dir, name_to_save + '_sfred.csv'))
    print('Done')

    print(f'Calculatiing FREDc for {bam}...')
    fredc = all_functions(gen_ko_1or0,genome_jds,rel_ab,genome_index,function_index)
    fredc = pd.Series(fredc, name = name_to_save)
    fredc.to_csv(os.path.join(output_dir, name_to_save + '_cfred.csv'))
    print('Done')
    return

def launch_analysis():
    pool = Pool(int(args.processors))
    pool.starmap(bams_analysis, [(bams_files[i], contigs_lengths_d, genome_contigs, function_index, predicted_functions, genome_index, genome_jds, total_jd, genome_reads[bams_files[i]], output_folder + '/output_fred', float(args.covered_genome_fraction), float(thr_min)) for i in range(len(bams_files))])

####################################################################
#Statistics FREDs
#########################################################################################################################################
if __name__ == '__main__':

    script_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(script_dir)

    args = parse_args()
    
    if os.path.exists(args.output_folder) == False:
        try:
            os.mkdir(args.output_folder)
            print(f"Directory '{args.output_folder}' created successfully.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{args.output_folder}'.")
    else:
        print(f"Directory '{args.output_folder}' already exists.")
    
    if os.path.exists(os.path.join(args.output_folder, 'output_micropherret')) == False:
        try:
            os.mkdir(os.path.join(args.output_folder, 'output_micropherret'))
            print(f"Directory '{os.path.join(args.output_folder, 'output_micropherret')}' created successfully.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{os.path.join(args.output_folder, 'output_micropherret')}'.")
    else:
        print(f"Directory '{os.path.join(args.output_folder, 'output_micropherret')}' already exists.")

    if os.path.exists(os.path.join(args.output_folder, 'output_fred')) == False:
        try:
            os.mkdir(os.path.join(args.output_folder, 'output_fred'))
            print(f"Directory '{os.path.join(args.output_folder, 'output_fred')}' created successfully.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{os.path.join(args.output_folder, 'output_fred')}'.")
    else:
        print(f"Directory '{os.path.join(args.output_folder, 'output_fred')}' already exists.")

    print('\n')
    
    output_folder = os.path.abspath(args.output_folder)

    print('Input generation steps')
    if not args.all_genomes:
        if args.genomes_folder.endswith('/'):
            genomes_folder = args.genomes_folder[:-1]
        else:
            genomes_folder = args.genomes_folder
    
        if not args.genomes_extension.startswith('.'):
            genomes_extension = '.'+args.genomes_extension
        else:
            genomes_extension = args.genomes_extension
        
        genomes_folder = os.path.abspath(genomes_folder)

        subprocess.check_call(['./input_generation.sh', genomes_folder, output_folder, genomes_extension, str(args.processors)])
        all_fasta = os.path.join(output_folder, 'input_fred/all_genomes.fa') 
        binning_info = os.path.join(output_folder, 'input_fred/info.txt') 
    else:
        all_fasta = args.all_genomes
        binning_info = args.binning_file

    if args.eggnog_annotation:
        annotation_folder = args.eggnog_annotation
    else:
        print('Annotation with eggNOG-mapper...')
        subprocess.check_call(['./pipeline_eggnog.sh', genomes_folder, output_folder, genomes_extension, str(args.processors), args.eggnog_database, args.eggnog_sensmode]) 
        annotation_folder = os.path.join(output_folder, 'input_fred/eggnog_annotations')
        print('Done')
    
    if annotation_folder.endswith('.csv'): 
        ko_df = pd.read_csv(annotation_folder).set_index('Genomes')
    else:
        ###get annotation matrix
        print('Getting annotation matrix from eggNOG-mapper annotations...')
        files = [os.path.join(annotation_folder,i) for i in os.listdir(annotation_folder)]
        data_ko = {}
        genomes_list = []
        ko_number = {}
        for f in files:
            if f.endswith('.annotations'):
                print(f'Processing {f} file...')
                genome, data_ko[genome], ko_number[genome] = get_kos(f)
                genomes_list.append(genome)
        df = pd.DataFrame(data_ko).T
        df.fillna(0, inplace = True)
        ko_df = pd.DataFrame(data = data_ko).T
    user_dataset = ko_df.fillna(0)
    user_dataset.to_csv(os.path.join(output_folder, 'input_fred/annotation_matrix.csv'))
    print('Done')
    

    if args.bam_files:
        bam_folder = os.path.abspath(args.bam_files)
    elif args.reads_folder:
        reads_folder = os.path.abspath(args.reads_folder)
        print('Alignment with Bowtie2...')
        subprocess.check_call(['./pipeline_bamfiles.sh', all_fasta, reads_folder, output_folder, str(args.processors), args.unpaired_reads])
        bam_folder = os.path.join(output_folder, 'input_fred/bowtie2')

    bam_folder = os.path.abspath(bam_folder)
    if bam_folder.endswith('/'):
        bam_folder=bam_folder[:-1]
    print('Getting number of aligned reads per sample...')
    subprocess.check_call(['./get_nreads.sh', bam_folder, str(args.processors), output_folder])
    nreads = os.path.join(output_folder, 'input_fred/nreads.txt')
    print('Done')

    print('Processing inputs for FRED calculation...')
    print('Processing genomic information...')
    genome_contigs = read_binning_file(binning_info)
    contigs_lengths_d = read_fasta(all_fasta)
    genome_index=index_genome(genome_contigs)
    genome_reads = get_reads_number(nreads)
    print('Done\n')

    print(f"Input generation is complete: files are stored in{os.path.join(output_folder, 'input_fred')}\n")

    ##### MICROPHERRET part
    if args.functions_list:
        functions = read_functions(args.functions_list)
        print(functions)
        function_file = args.functions_list
    else:
        functions = read_functions('./functions.txt')
        function_file = './functions.txt'
    
    if args.micropherret_predictions:
        print('miFRED will use provided functions for calculation')
        predicted_functions = pd.read_csv(args.micropherret_predictions).set_index('Unnamed: 0')[functions].loc[list(genome_contigs.keys())]
        function_index= index_function(function_file)
    elif args.KO:
        print('miFRED will use KO annotations for calculation')
        predicted_functions = ko_df.loc[list(genome_contigs.keys())]
        predicted_functions[predicted_functions > 1] = 1
        predicted_functions.to_csv(os.path.join(output_folder,'input_fred/KO_for_calculation.csv'))
        print(f"Used KOs for calculation stored in {os.path.join(output_folder, 'input_fred/KO_for_calculation.csv')}")
        function_index = {j:i for i,j in enumerate(list(predicted_functions.columns))}
    else:
        print('ML prediction using FAPROTAX-derived training set')
        print('Loading training set files and generating input for MICROPHERRET...')
        function_index= index_function(function_file)
        training_dataset = pd.read_csv(os.path.join(args.training_sets, 'dataset.csv')).set_index('Genome').drop('Species', axis = 1)
        training_acetoclastic = pd.read_csv(os.path.join(args.training_sets, 'dataset_acetoclastic_methanogenesis.csv')).set_index('Unnamed: 0').drop('acetoclastic_methanogenesis', axis = 1)

        print('Adjusting annotation matrix to MICROPHERRET training set...')
        validation_set= get_validation_set(user_dataset, training_dataset)
        print('Adjusting annotation matrix to MICROPHERRET acetoclastic methanogenesis updated training set...')
        val_aceto = get_validation_set(user_dataset, training_acetoclastic)
        print('Done')
        predicted_functions = validate(functions, validation_set, val_aceto)[0]
        print(f"Predicted functions stored in {os.path.join(output_folder, 'output_micropherret/predict_functions.csv')}")
        print(f"Number of genomes able to perfomed functions stored in {os.path.join(output_folder, 'output_micropherret/predict_sum.csv')}")
        print('Done\n')

    ####
    #pairwise jaccard distance for each pair of genome
    print('Calculating jaccard distances for FRED...')
    genome_jds=jaccard_distance(predicted_functions,genome_index) 
    genome_jds_df = pd.DataFrame(genome_jds)
    genome_jds_df['genome_name'] = list(genome_index.keys())
    genome_jds_df = genome_jds_df.set_index('genome_name')
    genome_jds_df.to_csv(os.path.join(output_folder,'output_fred/jaccard_distances.csv'))
    total_jd=0 
    for gen in genome_jds:
        total_jd+=sum(genome_jds[gen])
    total_jd=total_jd/(len(genome_contigs)*(len(genome_contigs)-1))
    print('Done\n')

    print('Processing alignment files...')
    bams_files = [os.path.join(bam_folder, b) for b in os.listdir(bam_folder) if b.endswith('sorted.bam')]

    if args.relative_abundance_threshold:
        thr_min = args.relative_abundance_threshold
    else:
        thr_min = 0
    
    launch_analysis()

    #################
    used_ras = pd.concat([pd.read_csv(os.path.join(output_folder,'output_fred/'+i)).set_index('Unnamed: 0') for i in os.listdir(os.path.join(output_folder,'output_fred/')) if i.endswith('_used_ra.csv')], axis = 1).rename_axis( 'Genome')
    used_ras.to_csv(os.path.join(output_folder,'output_fred/'+'used_relative_abundance.csv'))
    for f in os.listdir(os.path.join(output_folder,'output_fred/')):
        if f.endswith('_used_ra.csv'): os.remove(os.path.join(output_folder,'output_fred/'+f))

    norm_ras = pd.concat([pd.read_csv(os.path.join(output_folder,'output_fred/'+i)).set_index('Unnamed: 0') for i in os.listdir(os.path.join(output_folder,'output_fred/')) if i.endswith('_normalised_ra.csv')], axis = 1).rename_axis('Genome')
    norm_ras.to_csv(os.path.join(output_folder,'output_fred/'+'normalised_relative_abundance.csv'))
    for f in os.listdir(os.path.join(output_folder,'output_fred/')):
        if f.endswith('_normalised_ra.csv'): os.remove(os.path.join(output_folder,'output_fred/'+f))

    fredc = pd.concat([pd.read_csv(os.path.join(output_folder,'output_fred/'+i)).set_index('Unnamed: 0') for i in os.listdir(os.path.join(output_folder,'output_fred/')) if i.endswith('cfred.csv')], axis = 1).T.rename_axis('Sample')
    fredc.to_csv(os.path.join(output_folder,'output_fred/'+'fredc.csv'))
    for f in os.listdir(os.path.join(output_folder,'output_fred/')):
        if f.endswith('cfred.csv'): os.remove(os.path.join(output_folder,'output_fred/'+f))

    freds = pd.concat([pd.read_csv(os.path.join(output_folder,'output_fred/'+i)).set_index('Function') for i in os.listdir(os.path.join(output_folder,'output_fred/')) if i.endswith('sfred.csv')]).sort_index()
    freds.to_csv(os.path.join(output_folder,'output_fred/'+'freds.csv'))
    for f in os.listdir(os.path.join(output_folder,'output_fred/')):
        if f.endswith('sfred.csv'): os.remove(os.path.join(output_folder,'output_fred/'+f))
    
    #Statistics FREDs

    FREDs_stat = {}
    for f in set(freds.drop('Sample', axis = 1).index):
        FREDs_stat[f] = freds.loc[f].describe().drop(['count', '25%','50%','75%']).rename_axis('Statistic')
    FREDs_stat = pd.concat(FREDs_stat)
    FREDs_stat.to_csv(os.path.join(output_folder,'output_fred/'+'FREDs_statistic.csv'))

    #Plots FREDs

    for_pca_sumra, for_pca_sumra_red = get_matrix(freds, 'Relative abundance sum')
    for_pca_numb, for_pca_numb_red = get_matrix(freds, 'Number')
    for_pca_prop, for_pca_prop_red = get_matrix(freds, 'Proportion')

    f, ax = plt.subplots(figsize = (8,18))
    sns.heatmap(data = for_pca_prop_red.T, cbar_kws={"shrink": 0.3},linewidths=0.8,linecolor="white", cmap = sns.color_palette("flare", as_cmap=True))
    plt.savefig(os.path.join(output_folder,'output_fred/' + 'FREDs_proportion.png'), dpi = 600, bbox_inches = 'tight')

    f, ax = plt.subplots(figsize = (8,18))
    sns.heatmap(data = for_pca_numb_red.T, cbar_kws={"shrink": 0.3},linewidths=0.8,linecolor="white", cmap = sns.color_palette("flare", as_cmap=True))
    plt.savefig(os.path.join(output_folder,'output_fred/' + 'FREDs_number.png'), dpi = 600, bbox_inches = 'tight')

    f, ax = plt.subplots(figsize = (8,18))
    sns.heatmap(data = for_pca_sumra_red.T, cbar_kws={"shrink": 0.3},linewidths=0.8,linecolor="white", cmap = sns.color_palette("flare", as_cmap=True))
    plt.savefig(os.path.join(output_folder,'output_fred/' + 'FREDs_relab.png'), dpi = 600, bbox_inches = 'tight')

    # Plots FREDc

    f, ax = plt.subplots(figsize=(7, 6))
    sns.boxplot(fredc, y="FREDc", whis=[0, 100], width=.6, ax = ax)
    sns.stripplot(fredc, y="FREDc", size=4, color=".3", ax = ax)
    plt.savefig(os.path.join(output_folder,'output_fred/' + 'FREDc_boxplot.png'), dpi = 600)

    f, ax1 = plt.subplots(figsize=(7, 6))
    sns.histplot(fredc, x="FREDc", bins = 20,kde = True, stat = 'count', color = 'purple', ax = ax1)
    plt.savefig(os.path.join(output_folder,'output_fred/' + 'FREDc_histplot.png'), dpi = 600)

    print(f"FRED calculation is complete, output files and plots are stored in {os.path.join(output_folder,'output_fred/')}")
    print('Enjoy your analysis!')
    

    

    



    




