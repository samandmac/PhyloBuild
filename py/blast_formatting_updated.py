#!/usr/bin/env python
# coding: utf-8

# In[10]:


#import modules needed
import os
import re
import sys

#So here we would hope to be doing a correction to the results of the blast so that we're just getting one 
#gene match back per queried gene.

#Establish a list to store the file names in
file_names = []

#This allows user inputted genome and gene directories.
blast_directory = sys.argv[1]
#blast_directory = "/media/sam/Expansion/Bioinformatics/Glasgow_Andy_Patricia_Tree/Version1/BlastResults"

#For each file in the blast directory, add to the file name list
for file_name in os.listdir(blast_directory):
    #print(file_name)
    file_names.append(file_name)
    
#print(range(len(file_names)))

#For each file (strain)
for new in range(len(file_names)):
    print("Running through " + file_names[new] + " for reformatting")
    
    #Establish lists to append to per iteration
    new_info = dict()
    genes = []
    sequences = []
    percentages = []
    
    
    #Open a new file to save output to - this is just a new file per strain
    new_file = open(blast_directory + "/new_" + file_names[new], "w")
    
    
    #print(file_names[new])
    
    #Open each file in directory and iterate through the lines
    for lines in open(blast_directory + "/" + file_names[new]):
        
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

    #Now we're checking to see if the next gene in the sequence matches the gene we're currently at
    #For each index in genes
    for x in range(len(genes)):

        #Iterate through x + 1 to the last gene
        for y in range(x+1,len(genes)):
            #print(y)
            #If the gene name of x equals the gene name of y
            if (genes[x] == genes[y]):
                #If the number of genes of the same name is above 2
                if (genes.count(genes[x]) > 2):
                    #Then grab the index number of each (we make the assumption 3 matches will be enough)
                    indecesthree = [i for i, e in enumerate(genes) if e == genes[x]]
                    #Store each gene index
                    first = indecesthree[0]
                    second = indecesthree[1]
                    third = indecesthree[2]
                    
                    #Now compare the percentage coverage of each match and save the one with the highest match as 
                    #the winning sequence match to use
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

               
                #Now if we only have two of the same gene
                elif (percentages[x] > percentages[y]):
                    #Append the higher match sequence to the new_sequences
                    new_info[genes[x]] = sequences[x]
               
                #Otherwise, y is less than x (or equal, which wouldn't matter which we chose)
                else:
                    #Append the sequence of y to the new_sequences
                    new_info[genes[x]] = sequences[y]


    #For each index in the genes
    for x in range(len(genes)):
        #If the gene name isn't in the new list, because it didn't have a match...
        if genes[x] not in new_info:
            #print(genes[x])
            #We then add it to the new_genes list
            new_info[genes[x]] = sequences[x]


    
    #And then we iteratively add the contents of both lists to the respective genome blast file match.
    #for x in range(len(new_genes)):
    #    new_file.write(new_genes[x] + "\n" + new_sequences[x] + "\n")
    #print(new_info)
    for x in new_info:
        new_file.write(x + "\n" + new_info[x] + "\n")

