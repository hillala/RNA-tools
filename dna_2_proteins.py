#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

@author: hila hillay@gmail.com 
"""

from __future__ import division, print_function
import numpy as np  
import sys
import pickle
#from rna_2_proteins import translate_dna
sys.path.append("/../rna_function_file")
import rna_functions_file as rna

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

################################################
# This code translate reference dna to a 6 frame translation  
# within a specific nucleutide distance from aa annotated gene
# a referance dna and a text file that contains location of annotated regions is used 
##################################################
 
# a txt file with chromosomes name  
Namefile='/.../ChrName.txt'
chr_name= [name.rstrip('\n') for name in open(Namefile)]
# number of chromosomes
chr_num=sum(chr_name)

# file with transcribed region list + index
fAnn1='/.../Annotation_'
fAnn2='_Sense_Location.txt'

# referance dna
chr_f1='/.../Mus_musculus.GRCm38.dna.chromosome.'
chr_f2='.fa'

# seperate output file for each chromosome 
f_out1='/.../peptides_dir_'
f_out2='.pckl'

# a constant nuclei range to be translated
Vrange=[200,250,300,350,400,450,500]

for range1 in Vrange:
    peptides=[None]*chr_num
    for k in range(chr_num):
        name=chr_name[k]
        chrom=chr_name[k][3:]
        with open(chr_f1+chrom+chr_f2, 'r') as myfile:
                data=myfile.read().replace('\n', '')    
        myfile.close()
        G = [[int(item) for item in line.split()] for line in open(fAnn1+name+fAnn2)]
        p1,p2=np.array(zip(*G),dtype=np.int)
        pep=[rna.FindPeptedes(data[p2[j]:p2[j]+range1],codontable)\
        for j in range(len(p2))]
        peptides[k]  = [val for sublist in pep for val in sublist]
    f = open(f_out1+str(range1)+f_out2, 'wb')
    pickle.dump(peptides, f)
    f.close()
    print(range1)

#########################
## in case of reverse streand

# file with transcribed region list + index
fAnn1='/.../Annotation_'
fAnn2='_AntiSense_Location.txt'

# referance dna
chr_f1='/.../\
Mus_musculus.GRCm38.dna.chromosome.'
chr_f2='.fa'

# seperate output file for each chromosome
f_out1='/.../peptides_rev_'
f_out2='.pckl'

# a constant nuclei range to be translated
Vrange=[200,250,300,350,400,450,500]

for range1 in Vrange:
    peptides=[None]*chr_num
    for k in range(chr_num):
        name=chr_name[k]
        chrom=chr_name[k][3:]
        with open(chr_f1+chrom+chr_f2, 'r') as myfile:
                data=myfile.read().replace('\n', '')    
        myfile.close()
        G = [[int(item) for item in line.split()] for line in open(fAnn1+name+fAnn2)]
        p1,p2=np.array(zip(*G),dtype=np.int)
        pep=[rna.FindPeptedesRev(data[p1[j]-range1:p1[j]],codontable)\
        for j in range(len(p1))]
        peptides[k]  = [val for sublist in pep for val in sublist]
    f = open(f_out1+str(range1)+f_out2, 'wb')
    pickle.dump(peptides, f)
    f.close()




################################################
# save peptides to txt file: chr name, start position,streand
# reading frame, streand, sequance
##################################################

vec=['200','250','300','350','400','450','500']

# previouse saved files
f_name1='/.../peptides_'
f_name2='.pckl'

# annotation file
fAnn1='/.../Annotation_'
fAnn2='Sense_Location.txt'

# output file
out_file1='/.../peptides_'
out_file2='.txt'

direction=['dir','rev']
for d in direction:
    for v in vec:
        f = open(f_name1+d+'_'+v+f_name2, 'rb')
        obj = pickle.load(f)
        f.close()
        fileID=open(out_file1+d+'_'+v+out_file2,"w")
        for k in range(chr_num):
            name=chr_name[k]
            pept=obj[k]
            if d=='dir':
                G = [[int(item) for item in line.split()] for line in open(fAnn1+name+'_'+fAnn2)]
            else:
                G = [[int(item) for item in line.split()] for line in open(fAnn1+name+'_Anti'+fAnn2)]
            p1,p2=np.array(zip(*G),dtype=np.int)
            for i in range(len(p1)):
                for j in range(3):
                    fileID.write('%s\t %s\t %s\t %s\t\n' % \
                         (name,p1[i],'frame '+str(j+1),d))
                    fileID.write('%s\n' % \
                         (pept[i*3+j]))
            
        fileID.close()


