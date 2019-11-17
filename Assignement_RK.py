from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import GC
import pandas as pd
import matplotlib.pyplot as plt

## 1 ##
def get_sequences_from_file(fasta_fn):
    sequence_data_dict = {}
    for record in SeqIO.parse(fasta_fn, "fasta"):
        description = record.description.split()
        species_name = description[1] + " " + description[2]
        sequence_data_dict[species_name] = record.seq
    return(sequence_data_dict)

#seq = get_sequences_from_file("penguins_cytb.fasta")

## 2 ##
def translate_function(string_nucleotides):
    seek = ""
    aa_seq_string = ""
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
    for x in range(len(string_nucleotides)):
        seek += str(string_nucleotides[x])
        if (x+1) % 3 == 0:
            if seek in mito_table.stop_codons:
                break
            aa_seq_string += mito_table.forward_table[seek]
            seek = ""
    return(aa_seq_string)

## 3 ##
def DNA_translate_biopython(Input_seq):
    output_seq = Input_seq.translate()
    return output_seq

"""aa_seq = []
aa_biopython_seq = []
aa_GC_content = []
for value in seq.values():
    seq = value
    aa_seq.append(translate_function(seq))
    aa_biopython_seq.append(DNA_translate_biopython(seq))
print(aa_seq)
print(aa_biopython_seq)"""



## 4 ##
def Mol_wt_aa(aa_seq):
    analyzed_seq = []
    for list_item in aa_seq:
        analyzed_seq.append(ProteinAnalysis(list_item).molecular_weight())
    return analyzed_seq
#print(Mol_wt_aa(aa_seq))

## 5 ##
def GC_content(Input_seq):
    output_seq = GC(Input_seq)
    return output_seq

"""get_seq_for_GC = get_sequences_from_file("penguins_cytb.fasta")
seq_list_GC = []
for x in get_seq_for_GC.values():
    seq_GC = x
    seq_list_GC.append(GC_content(seq_GC))
print(seq_list_GC)"""

#%%%%%%%%%%%%%%#
###   MAIN   ###
#%%%%%%%%%%%%%%#

cytb_seqs = get_sequences_from_file("penguins_cytb.fasta")

penguins_df = pd.read_csv("penguins_mass.csv")
species_list = list(penguins_df.species)

## 6 ##
penguins_df['molecular_weight'] = 'NaN'
penguins_df['GC_content'] = 'NaN'

## 7 ##
aa_seq = []
GC_list = []
for key, value in cytb_seqs.items():
    aa_seq.append(translate_function(value))
    GC_list.append(GC_content(value))
penguins_df['molecular_weight'] = Mol_wt_aa(aa_seq)
penguins_df['GC_content'] = GC_list

## 8 ##
## Plot a bar-chart of the mass with the x-axes labeled with species names.
## *Q1* What is the smallest penguin species?
#Eudyptula minor
## *Q2* What else is interesting about this species?
#Average lifespan of this species is 6.5 years
penguins_df.plot(kind='bar', x='species', y='mass')
plt.show()

## 9 ##
## Plot a visualization of the molecular weight (y-axis) as a function of GC-content (x-axis).
penguins_df.plot(kind='scatter', x='molecular_weight', y='GC_content', linestyle='-', color='red')
plt.show()

## 10 ##
penguins_df.to_csv('penguins_mass_cytb.csv')




