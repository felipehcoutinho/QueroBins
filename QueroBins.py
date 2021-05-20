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
parser.add_argument("--metabat2", help="Flag to perform sequence binning through Metabat2", default=False, type=bool)
parser.add_argument("--abundance_table", help="Flag to run the abundance calculation modules", default=False, type=bool)
parser.add_argument("--abundance_rpkm", help="Flag to calculate abundance as RPKM", default=False, type=bool)
parser.add_argument("--metagenomes_dir", help="Directories containing metagenome fastq files to be used for abundance calculation",  nargs="+", type=str)
parser.add_argument("--metagenomes_extension", help="Extension of the fastq files in metagenomes_dir to be used for abundance calculation", default="fastq", type=str)
parser.add_argument("--assemble", help="Flag to run the assembly module", default=False, type=bool)
parser.add_argument("--samples_dir", help="Directory containing fastq files to be used for assembly", type=str)
parser.add_argument("--samples_extension", help="Extension of the fastq files in samples_dir", default='fastq', type=str)
parser.add_argument("--threads", help="The number of threads to be used", default=1, type=int)
parser.add_argument("--info_output", help="The output table file generated by the index module", default="Seq_Info.tsv", type =str)
parser.add_argument("--parse_only", help="Flag to skip running any programs and only parse their output", default=False, type=bool)
args = parser.parse_args()


def central():
    coassembly_df = pd.DataFrame({'DUMMY' : []})
    if (args.coassembly_table):
        print(f"Reading coassembly info from {args.coassembly_table}")
        coassembly_df = pd.read_csv(args.coassembly_table, sep="\t",index_col=0,header=0)
    seq_info_dict = assembly_module(coassembly_df = coassembly_df, threads = args.threads)

def assembly_module(coassembly_df = [], threads = 1, max_memory = 500):
    assembly_list = []
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
            print(f"Running: {command}")
            subprocess.call(command, shell=True)
    #seq_info = index_assemblies(assembly_list = assembly_list)
    seq_info = []
    return seq_info

central()