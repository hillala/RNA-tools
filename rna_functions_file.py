#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

@author: hila hillay@gmail.com
"""
import numpy as np
import random

codontable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}


# translate dna sequance to protein sequance in 3 reading frames without errors 
def translate_dna(sequence, codontable):
    Num=len(sequence) %3
    sequence.replace('T', 'U')
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-Num, 3)]
    protein_sequence = ''.join([codontable[codon] for codon in codons])
    return protein_sequence

# chack if codon are readable 
def InList(letter,codon,codontable):        
    if letter in codon:
        return 'X'            
    else:
        return codontable[codon]

# translate dna sequance to protein sequance in 3 reading frames with errors     
def translate_dna_with_err(sequence, codontable):
    #seq=seq_record.seq.tostring()
    Num=len(sequence) %3
    sequence.replace('T', 'U')
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-Num, 3)]
    protein_sequence=[InList('N',codon,codontable) for codon in codons]
    protein_sequence = ''.join(protein_sequence)
    return protein_sequence

# genereate opposite streand from dna sequance
def OppositeStreand(region):
    anti_data=np.array(list(region))
    ind=np.array([i for i,x in enumerate(region) if x == 'T'])
    if len(ind)>0:
        anti_data[ind]='A'
    ind=np.array([i for i,x in enumerate(region) if x == 'A'])
    if len(ind)>0:
        anti_data[ind]='T'
    ind=np.array([i for i,x in enumerate(region) if x == 'G'])
    if len(ind)>0:
        anti_data[ind]='C'
    ind=np.array([i for i,x in enumerate(region) if x == 'C'])
    if len(ind)>0:
        anti_data[ind]='G'
    anti_data=list(anti_data)
    anti_data=''.join(anti_data)
    return anti_data[::-1]

# return random nucli
def RandomNucleotide():
    r=random.randint(1,4)
    if r==1:
        return 'A'
    if r==2:
        return 'T'
    if r==3:
        return 'G'
    if r==4:
        return 'C'

#  3 frame translation on reverse streand
def FindPeptedesFramRev(region,codontable):
    anti_region=OppositeStreand(region)
    protein_sequence1=translate_dna_with_err(anti_region, codontable)
    protein_sequence2=translate_dna_with_err(anti_region[1:], codontable)
    protein_sequence3=translate_dna_with_err(anti_region[2:], codontable)
    return protein_sequence1,protein_sequence2,protein_sequence3

#  3 frame translation on direct streand       
def FindPeptedesFram(region,codontable):
    protein_sequence1=translate_dna_with_err(region, codontable)
    protein_sequence2=translate_dna_with_err(region[1:], codontable)
    protein_sequence3=translate_dna_with_err(region[2:], codontable)
    return protein_sequence1,protein_sequence2,protein_sequence3

#  tranverse a direct opposite dna sequance to reverse opposite streand   
def DNAFramTransRev(region):
    anti_region=OppositeStreand(region)
    anti_region[::-1]
    return anti_region        

#  3 frame translation on an oposite reverse streand    
def FindPeptedesFramTransRev(region,codontable):
    anti_region=OppositeStreand(region)
    anti_region[::-1]
    protein_sequence1=translate_dna_with_err(anti_region, codontable)
    protein_sequence2=translate_dna_with_err(anti_region[1:], codontable)
    protein_sequence3=translate_dna_with_err(anti_region[2:], codontable)
    return protein_sequence1,protein_sequence2,protein_sequence3 


# write a 6 frame translation into a file: frame index 1-3, peptide sequance, chromosome name, start point index
def WritePeptideFile(frame,direction,peptide,chromosome,initial_location,fileID):
    fileID.write('%s_%s_%s_%s\n' %   (chromosome,initial_location,'frame'+frame\
                                      ,direction))
    fileID.write('%s\n' % (peptide))   
