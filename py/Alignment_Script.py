#!/usr/bin/env python
# coding: utf-8

# In[7]:


import os
import re
import sys
import subprocess
from Bio import SeqIO

#We make genes an empty list
genes = []


working_directory = sys.argv[1]
os.mkdir(working_directory + "/Alignment_Files")
output_directory = sys.argv[2]

#I GOT RID OF THE GENE NAMES PARAGRAPH, SEE IPYNB SCRIPT TO SEE WHAT THAT LOOKED LIKE IF EVERYTHING GOES WRONG

#Now we need to extract the sequence matches from the BlastResults file and concatenate them based on gene type

#For each file in the directory BlastResults (which is each strains sequence matches to the tree genes)
strains = []
for file in os.listdir(working_directory + "/BlastResults/"):
    #print(file)
    strain = re.sub(".fasta","",file)
    strains.append(strain)
    #We create a new dict called sequences
    sequences = {}
    #We then open the file to view all the lines
    for lines in open(working_directory + "/BlastResults/" + file):
        #We strip the newline
        line = lines.rstrip("\n")
        #Indicate those gene names beginning with > as the key value
        if re.search(">", line):
            key = line
            #print(key)
        else:
            #And those that do not (i.e., sequence of the gene) as the value in the dict
            value = line
            #print(value)
            #And combine both to a dict i.e. Gene: Sequence
            sequences[key] = value
    #print(sequences)
    
    #For each gene found in the sequence for a specific file
    for genes_present in sequences.keys():
        file_name = re.sub(">","",genes_present)
        #We open a file called the gene name (which is appendable, so we don't overwrite)
        with open(working_directory+ "/Alignment_Files/" + file_name + ".tsv", 'a') as nfile:
            #And save the sequence for each gene, in FASTA format, using the strain name as the header for each sequence.
            nfile.write(">" + file + "\n" + sequences[genes_present] + "\n")

#Come back to below, this is a method to run bash commands through python (but it doesn't work on my windows computer)
bashCommand = "for x in ../Alignment_Files/*.tsv; do mafft $x > ../Alignment_Files/${x##*/}_aligned.csv; done"
process = subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
output, error = process.communicate()

#This is how to generate dashes for each gene type based on aligned sequences
#This is for the strains that had no matching sequence to a particular gene.
#Set up a dict to store the gene name and the corresponding dash sequence
each_file_dashed = {}            
#Open each file, as long as it's a csv file (will be changed for final script)
for file in os.listdir(working_directory + "/Alignment_Files/"):
    if file.endswith(".csv"):
        for lines in open(working_directory + "/Alignment_Files/" + file):
            line = lines.rstrip("\n")
            #If a header, then restart dashes amount
            if line.startswith(">"):
                #print(line)
                dashes = ""
                key = file
            else:
                dashes = dashes + (len(line) * "-")
        #print(dashes)
        file = re.sub(".tsv_aligned.csv","",file)
        each_file_dashed[file] = dashes

#d = SeqIO.to_dict(SeqIO.parse("../py/>chuA.tsv_aligned.csv", 'fasta'))

#print(d)


for strain in strains:
    sequence = ""
    for file in os.listdir(working_directory + "/Alignment_Files/"):
        if file.endswith(".csv"):
            d = SeqIO.to_dict(SeqIO.parse(working_directory + '/Alignment_Files/' + file, 'fasta'))
            #strain = strain + ".fasta"
            try:
                header = re.sub(".tsv_aligned.csv","",file)
                sequence = sequence + str(d[strain+'.fasta'].seq)
                #with open(strain+"_final.tsv", 'a') as nfile:
                    #And save the sequence for each gene, in FASTA format, using the strain name as the header for each sequence.
                    #nfile.write(header + "\n" + str(d[strain+'.fasta'].seq) + "\n")
                    
                    #print(strain)
                    #print("success")
                    #print(header)
                    #print(str(d[strain+'.fasta'].seq))
            except:
                #print(strain)
                header = re.sub(".tsv_aligned.csv","",file)
                sequence = sequence + each_file_dashed[header]
                #print("failure")
                #print(header)
                #with open(strain+"_final.tsv", 'a') as nfile:
                    #And save the sequence for each gene, in FASTA format, using the strain name as the header for each sequence.
                    #nfile.write(header + "\n" + each_file_dashed[header] + "\n")
                #print(strain)
                #print(each_file_dashed[">aceK_arpA"])
    #print(strain + "\n" + sequence)
    with open(working_directory + "/final_alignment.tsv", 'a') as nfile:
        nfile.write(">" + strain + "\n" + sequence + "\n")
        

