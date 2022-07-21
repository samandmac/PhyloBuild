#!/usr/bin/env python
# coding: utf-8

# In[22]:


#Making sure that the file input's are in the correct format

from Bio import SeqIO
import os
import re
import sys

#This allows user inputted genome and gene directories.
genomePath = sys.argv[1]
genePath = sys.argv[2]

#Function for determining whether a file is FASTA format or not
def is_fasta(filename):
    with open(filename,"r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta) # False when fasta is empty - i.e. wasn't a fasta file

#Checks each directory containing user input - Genomes, genes, template genomes, and template genes.
for file in os.listdir(genomePath):
    if (file != "Template_Genomes"):
        filepath = genomePath + "/" + file 
        if (is_fasta(filepath) == False):
            #If a file is determined to not be in FASTA format, then the script does a sys.exit().
            sys.exit("Genomes of interest not in FASTA format: see " + filepath)
print("Input: Genomes of interest are in FASTA format")

for file in os.listdir(genePath):
    if (file != "Template_Genes"):
        filepath = genePath + "/" + file 
        if (is_fasta(filepath) == False):
            sys.exit("Genes of interest not in FASTA format: see: " + filepath)
print("Input: Genes of interest are in FASTA format")
            
for file in os.listdir(genomePath):
    if (file == "Template_Genomes"):
        filepath = genomePath + "/" + file 
        for file in os.listdir(filepath):
            file = genomePath + "/" + "Template_Genomes/" + file
            if (is_fasta(file) == False):
                sys.exit("Template Genomes are not in FASTA format: see " + file)
print("Input: Template Genomes are in FASTA format")
                
for file in os.listdir(genePath):
    if (file == "Template_Genes"):
        filepath = genePath + "/" + file 
        for file in os.listdir(filepath):
            file = genePath + "/" + "Template_Genes/" + file
            if (is_fasta(file) == False):
                sys.exit("Template Genes are not in FASTA format: see " + file)
print("Input: Template Genes are in FASTA format")


# In[ ]:




