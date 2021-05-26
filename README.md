# QueroBins
Automated Analysis of Metagenome Assembled Genomes of Bacteria and Archaea

## Dependencies
- [MetaBat2](https://bitbucket.org/berkeleylab/metabat/src/master/)
- [CheckM](https://ecogenomics.github.io/CheckM/)
- [GTDB-tk](https://github.com/Ecogenomics/GTDBTk)
- [CAT/BAT](https://github.com/dutilh/CAT)
- [SPAdes](https://github.com/ablab/spades)

## Usage
### Assembling with 48 threads and using a maximum of 300 GB of memory
`python3 QueroBins.py assemble True --coassembly_table My_Coassembly_Table.tsv --threads 48 --max_memory 300`

My_Coassembly_Table.tsv specifies how samples should be coassembled together. It must be a .tsv format file with 4 columns and aheader in the following header: Sample, Group, R1, R2.
The contents of the table should specify Sample_ID, Group_ID, path to the R1 file, path to the R2 file. To have a sample be assembled individually (i.e. no coassembly) simply assign it to a group with no other samples.

This command will generate, for each group:
- Assembly_Group_[Group_ID]

### Processing pre-assembled scaffold files and filtering scaffolds shorter than 2000 bp

`python3 QueroBins.py --assemblies_files Assembly_1.fasta Assembly_2.fasta Assembly_3.fasta... --min_scaffold_length 2000`

This command will simply filter out scaffolds shorter than the user defined value and also rename the scaffolds from all assemblies to avoid repeated IDs. Generated files:
- Filtered_Renamed_Assembly_1.fasta, Filtered_Renamed_Assembly_2.fasta, Filtered_Renamed_Assembly_3.fasta...
- Scaffolds_Info.tsv A table with basic infomartion about all the processed scaffolds

### Binning scaffolds with MetaBat2

`python3 QueroBins.py --assemblies Filtered_Renamed_Assembly_1.fasta Filtered_Renamed_Assembly_2.fasta Filtered_Renamed_Assembly_3.fasta... --make_bins True --binning_method metabat --threads 24`

### Taxonomy aware binning with CAT and MetaBat2

`python3 QueroBins.py --assemblies Filtered_Renamed_Assembly_1.fasta Filtered_Renamed_Assembly_2.fasta Filtered_Renamed_Assembly_3.fasta... --make_bins True  --tax_aware_binning True --binning_method metabat --threads 24`

The command above will frist classify contigs taxonomically with CAT, split the assemblies according to contig assigned phylum and perform the binning individiaully on the split files
