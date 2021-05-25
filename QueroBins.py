#! /usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqUtils import GC
import pandas as pd
import argparse
import subprocess
import re
import glob
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("--coassembly_table", help="Coassembly instructions .tsv table", type=str)
parser.add_argument("--max_memory", help="Maximum memory to use (GB)", default=500, type=int)
parser.add_argument("--abundance_table", help="Flag to run the abundance calculation modules", default=False, type=bool)
parser.add_argument("--abundance_rpkm", help="Flag to calculate abundance as RPKM", default=False, type=bool)
parser.add_argument("--metagenomes_dir", help="Directories containing metagenome fastq files to be used for abundance calculation",  nargs="+", type=str)
parser.add_argument("--metagenomes_extension", help="Extension of the fastq files in metagenomes_dir to be used for abundance calculation", default="fastq", type=str)
parser.add_argument("--assemble", help="Flag to run the assembly module", default=False, type=bool)
parser.add_argument("--min_scaffold_length", help="Minimum scaffold length", default=2000, type=int)
parser.add_argument("--assemblies_files", help="Fasta files of assembled contigs", nargs="+", type=str)
parser.add_argument("--assembly_format", help="Assembly files format", default="fasta", type=str)
parser.add_argument("--scaffold_info_table", help=".tsv file to write scaffold info", default="Scaffolds_Info.tsv", type=str)
parser.add_argument("--bin_info_table", help=".tsv file to write Bin info", default="Bin_Info.tsv", type=str)
parser.add_argument("--tax_aware_binning", help="Flag to perform taxonomy aware binning module", default=False, type=bool)
parser.add_argument("--cat_db", help="Path to cat db", default="/mnt/lustre/repos/bio/databases/public/cat/20210107/2021-01-07_CAT_database", type=str)
parser.add_argument("--cat_taxonomy", help="Path to cat taxonomy", default="/mnt/lustre/repos/bio/databases/public/cat/20210107/2021-01-07_taxonomy", type=str)
parser.add_argument("--make_bins", help="Flag to run the binning module", default=False, type=bool)
parser.add_argument("--binning_method", help="Binning method to use", default='metabat', type=str)
parser.add_argument("--metabat_preset", help="Metabat binning preset to use", default='metabat', type=str)
parser.add_argument("--call_checkm", help="Flag to run CheckM", default=False, type=bool)
parser.add_argument("--call_GTDBtk", help="Flag to run GTDBtk", default=False, type=bool)
parser.add_argument("--threads", help="The number of threads to be used", default=1, type=int)
parser.add_argument("--parse_only", help="Flag to skip running any programs and only parse their output", default=False, type=bool)
args = parser.parse_args()


def central():
    #Perform the assemblies if specified by the user
    coassembly_df = pd.DataFrame({'DUMMY' : []})
    assemblies_list = args.assemblies_files
    if ((args.coassembly_table) and (args.assemble == True)):
        print(f"Reading coassembly info from {args.coassembly_table}")
        coassembly_df = pd.read_csv(args.coassembly_table, sep="\t",index_col=0,header=0)
        assemblies_list = assembly_module(coassembly_df = coassembly_df, threads = args.threads)
    #Index assembled seq info
    (seq_info,assemblies_list) = index_assemblies(assemblies_list = assemblies_list, rename_seqs = True ,  min_length = args.min_scaffold_length)
    if (args.tax_aware_binning):
        (seq_info,assemblies_list) = annotation_module(assemblies_list = assemblies_list, seq_info = seq_info)
    if (args.make_bins):
        (seq_info,bin_info) = binning_module(assemblies_list = assemblies_list, threads = args.threads, method = args.binning_method, min_length = args.min_scaffold_length, seq_info = seq_info)
    if (args.call_checkm):
        (bin_info) = checkm_module(bin_info = bin_info)
    #Convert the 2d dictionary info_dict into a pandas dataframe and print it to output_dataframe_file in .tsv format
    info_dataframe = pd.DataFrame.from_dict(seq_info)
    info_dataframe.index.name = 'Scaffold'
    info_dataframe.to_csv(args.scaffold_info_table,sep="\t",na_rep='NA')
    
    info_dataframe = pd.DataFrame.from_dict(bin_info)
    info_dataframe.index.name = 'Bin'
    info_dataframe.to_csv(args.bin_info_table,sep="\t",na_rep='NA')

def checkm_module(bin_info = defaultdict(dict)):
    command = f"checkm lineage_wf --tab_table --file CheckM_Bin_Info.tsv --threads {args.threads} --pplacer_threads {args.threads} --extension fa . CheckM_Results"
    if (args.parse_only == False):
        print(f"Running: {command}")
        subprocess.call(command, shell=True)
    checkm_df =  pd.read_csv("CheckM_Bin_Info.tsv", sep="\t",index_col=0,header=0)
    for i,row in checkm_df.iterrows():
        bin_info['Marker_Lineage'][i] = row['Marker lineage']
        bin_info['Completeness'][i] = row['Completeness']
        bin_info['Contamination'][i] = row['Contamination']
        bin_info['Strain_Heterogeneity'][i] = row['Strain heterogeneity']
    return(bin_info)

def annotation_module(assemblies_list = [], seq_info = defaultdict(dict)):
    for assembly_file in assemblies_list:
        assembly_file_prefix = get_prefix(assembly_file,"fasta")
        command = f"CAT contigs -c {assembly_file} -d {args.cat_db} -t {args.cat_taxonomy} --nproc {args.threads}"
        if (args.parse_only == False):
            print(f"Running: {command}")
            subprocess.call(command, shell=True)
    
    return(seq_info)

def index_bins(prefix="Bin",extension="fa",format="fasta",seq_info = defaultdict(dict)):
    bin_info = defaultdict(dict)
    bin_files = glob.glob(f"{prefix}*{extension}")
    print ('Indexing bins')
    for bin_file in bin_files:
        print(f'Processing {bin_file}')
        bin_file_prefix = get_prefix(bin_file,"fa")
        bin_info['Sequences'][bin_file_prefix] = 0
        bin_info['Base_Pairs'][bin_file_prefix] = 0
        bin_info['File'][bin_file_prefix] = bin_file
        #Iterate over genomic sequences in the file. Collect basic Info
        for seqobj in SeqIO.parse(bin_file, "fasta"):
            bin_info['Sequences'][bin_file_prefix] += 1
            bin_info['Base_Pairs'][bin_file_prefix] += len(seqobj.seq)
            seq_info['Bin'][seqobj.id] = bin_file_prefix
            seq_info['Bin_File'][seqobj.id] = bin_file
                
    return(seq_info,bin_info)
    
def binning_module(assemblies_list = [], threads = 1, method = 'metabat', seq_info = defaultdict(dict), min_length = 2000):
    if (method == "metabat"):
        for assembly_file in assemblies_list:
            assembly_file_prefix = get_prefix(assembly_file,"fasta")
            command = f"metabat --inFile {assembly_file} --outFile Bin_{assembly_file_prefix} --minContig {min_length} --numThreads {threads}"
            if (args.parse_only == False):
                print(f"Running: {command}")
                subprocess.call(command, shell=True)
    
    (seq_info,bin_info) = index_bins(prefix="Bin",extension="fa",format="fasta",seq_info = defaultdict(dict))
    return(seq_info,bin_info)
            
def assembly_module(coassembly_df = False, threads = 1, max_memory = 500):
    assemblies_list = []
    if (not coassembly_df.empty):
        for group in set(coassembly_df['Group']):
            print(f"Processing samples from group {group}")
            group_df = coassembly_df[coassembly_df['Group'] == group]
            samples_list = group_df.index
            print(f"Samples {samples_list}")
            r1_files = list(group_df['R1'])
            r1_files = ' -1 '.join(r1_files)
            r2_files = list(group_df['R2'])
            r2_files = ' -2 '.join(r2_files)
            #print(r1_files,r2_files)
            command = f"spades.py --meta --memory {max_memory} --threads {threads} -1 {r1_files} -2 {r2_files} -o Assembly_Group_{group}"
            if (args.parse_only == False):
                print(f"Running: {command}")
                subprocess.call(command, shell=True)
    return assemblies_list

def index_assemblies(assemblies_list = [], rename_seqs = True, min_length = 2000):
    seq_info = defaultdict(dict)
    filtered_assemblies_list = []
    print ('Indexing scaffolds')
    for assembly_file in assemblies_list:
        print(f'Processing {assembly_file}')
        assembly_file_prefix = get_prefix(assembly_file,"fasta")
        out_seq_file = "Filtered_Renamed_"+assembly_file_prefix+'.fasta'
        filtered_assemblies_list.append(out_seq_file)
        with open(out_seq_file, 'w', newline='') as OUT:
            seq_counter = 0
            #Iterate over genomic sequences in the file. Collect basic Info
            for seqobj in SeqIO.parse(assembly_file, "fasta"):
                seq_counter += 1
                #Rename sequences if specified in the function call
                if (rename_seqs == True):
                    new_id = assembly_file_prefix+'_Scaffold_'+str(seq_counter)
                    seq_info['Original_ID'][new_id] = seqobj.id
                    seqobj.id = new_id    
                    #Do not allow duplicated IDs
                    if (seqobj.id in seq_info['Description']):
                         raise Exception(f'Duplicated ID: {seqobj.id} in {assembly_file}')
                seq_info['Description'][seqobj.id] = seqobj.description
                seq_info['GC'][seqobj.id] = round(GC(seqobj.seq),2)
                seq_info['Length'][seqobj.id] = len(seqobj.seq)
                seq_info['Original_File'][seqobj.id] = assembly_file
                if (seq_info['Length'][seqobj.id] >= min_length):
                    SeqIO.write(seqobj, OUT, "fasta")
            
    return (seq_info,filtered_assemblies_list)
    
def get_prefix(file,format):
    prefix_file = re.sub(f'.{format}','',file)
    prefix_file = re.sub('(.)+/','',prefix_file)
    return prefix_file
    
central()