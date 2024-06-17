from Bio import SeqIO
from Bio.Seq import Seq
import re

infile = open("PF00186_35.fasta", 'r')

sequences = []
glycine_pos_count = 0
mutant_sequence_count_1 = 0
human_sequence_count_2 = 0
gly_no_pro = 0
comp_no_pro = 0
human_trp_pos_count = 0
human_trp_pos_count_yes = 0
taxonomy_list = open("taxonomy_list.fasta", 'w')

for record in SeqIO.parse(infile, "fasta"):
    sequence = str(record.seq)
    sequences.append(sequence)
    glycine_pos = sequence[65].lower()
    proline_pos_2 = sequence[94].lower()
    proline_pos_1 = sequence[88].lower()
    human_insertion = sequence[273:282].lower()
    human_known_comp_1 = sequence[274].lower()
    human_known_comp_2 = sequence[273].lower()
    human_known_comp_3 = sequence[280].lower()
    human_known_comp_4 = sequence[281].lower()
    human_trp_pos = sequence[482].lower()
    human_trp_pos_region = sequence[470:495].lower()
    methionine_pos = sequence[77].lower()
    proline_region = sequence[80:98].lower()
    glycine_region = sequence[63:67].lower()
    if re.search(r'pp|p-*p', proline_region):
        mutant_sequence_count_1 += 1
        if glycine_pos != 'g':
            print(f"{record.description}")
            print(proline_region)
            print(" " + glycine_region)
            print(human_insertion)
            glycine_pos_count += 1
        if human_known_comp_2 != 'p':
            print(f"{record.description}")
            print(proline_region)
            print(" " + glycine_region)
            print(human_insertion)
            human_sequence_count_2 += 1
        if human_trp_pos != 'w':
            human_trp_pos_count += 1
            print(human_trp_pos_region)
    if not re.search(r'pp|p-*p', proline_region):
        if glycine_pos == 'g':
            gly_no_pro += 1
        if human_known_comp_2 == 'p':
            comp_no_pro += 1
        
    
            
infile.close()
print(glycine_pos_count)
print(mutant_sequence_count_1)
print(human_sequence_count_2)
print(gly_no_pro)
print(comp_no_pro)
print(human_trp_pos_count)
print(human_trp_pos_count_yes)