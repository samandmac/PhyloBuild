#!/usr/bin/env python
# coding: utf-8

# In[92]:


#import modules needed
import os
import re
import sys

#Establish a list of lists to save values to, descriptors match values to be held
file_names = []

#This allows user inputted genome and gene directories.
blast_directory = sys.argv[1]
#blast_directory = "/home/sam/Documents/Bioinformatics/Glasgow_Andy_Patricia_Tree/Version1/BlastResults"


for file_name in os.listdir(blast_directory):
    #print(file_name)
    file_names.append(file_name)
    
#print(range(len(file_names)))

for new in range(len(file_names)):
    #print("Running through " + file_names[new] + " for reformatting")
    #Establish lists to append to per iteration
    new_info = dict()
    genes = []
    sequences = []
    percentages = []
    
    
    #Open a new file to save output to
    new_file = open(blast_directory + "/new_" + file_names[new], "w")
    
    
    #print(file_names[new])
    
    #For each line in the file
    for lines in open(blast_directory + "/" + file_names[new]):
        #line = lines.rstrip("\n").rsplit("\t")
        #print(lines)
        #If the line has a ">" then it must be the gene name
        if re.search(">", lines):
            #So we append it to the genes list
            genes.append(lines.rstrip("\n"))
        #Else if the length of the line is greater than 10, then it must be the sequence, not the percentage
        elif range(len(lines) > 10):
            #print(lines)
            #So we can append it to the sequences list
            sequences.append(lines.rstrip("\n"))
        #Otherwise it's the percentage
        else:
            #So we can append it to the percentages list
            percentages.append(lines.rstrip("\n"))

    #print(genes)
    #print(sequences)
    #print(percentages)

    #For each index in genes
    for x in range(len(genes)):
        #print(x)
        #For each index between x + 1 and the max
        for y in range(x+1,len(genes)):
            #print(y)
            #If the gene name of x equals the gene name of y
            if (genes[x] == genes[y]):
                if (genes.count(genes[x]) > 2):
                    indecesthree = [i for i, e in enumerate(genes) if e == genes[x]]
                    first = indecesthree[0]
                    second = indecesthree[1]
                    third = indecesthree[2]
                    if percentages[first] > percentages[second]:
                        if percentages[first] > percentages[third]:
                            new_info[genes[x]] = sequences[first]
                        else:
                            new_info[genes[x]] = sequences[third]
                    else:
                        if percentages[second] > percentages[third]:
                            new_info[genes[x]] = sequences[second]
                        else:
                            new_info[genes[x]] = sequences[third]

                #print(genes[x])
                #print(genes[y])
                #print(percentages[x])
                #print(percentages[y])
                #If the percentage of x is less than y
                #print(genes[x])
                elif (percentages[x] > percentages[y]):
                    #Append the lower sequence to the new_sequences
                    new_info[genes[x]] = sequences[x]
                    #print(new_genes[x])
                #Otherwise, y is less than x (or equal, which wouldn't matter which we chose)
                else:
                    #Append the sequence of y to the new_sequences
                    new_info[genes[x]] = sequences[y]

    #print(new_genes)

    #For each index in the genes
    for x in range(len(genes)):
        #If the gene name isn't in the new list, because it didn't have a match...
        if genes[x] not in new_info:
            #print(genes[x])
            #We then add it to the new_genes list
            new_info[genes[x]] = sequences[x]
        #if new_info.count(genes[x]) > 1:
         #   print("error")
        


    #print(new_genes)
    #And then we iteratively add the contents of both lists to the respective genome blast file match.
    #for x in range(len(new_genes)):
    #    new_file.write(new_genes[x] + "\n" + new_sequences[x] + "\n")
    #print(new_info)
    for x in new_info:
        new_file.write(x + "\n" + new_info[x] + "\n")

