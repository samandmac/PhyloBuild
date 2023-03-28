#!/usr/bin/env python
# coding: utf-8

# In[51]:


import os
import re
import sys

#This is a new method to generate a sample sheet to separate out the visualisation and making a phylogenetic tree
new_phylogeny_list = open("sample_sheet.txt", "w")

working_directory = sys.argv[1]
group_option = sys.argv[2]
oi_option = sys.argv[3]
name_option = sys.argv[4]
output_dir = sys.argv[5]

new_phylogeny_list.write("strain" + "\t" + "name_for_plot" + "\t" + "group" + "\t" + "of_interest" + "\n")

#Path will need to be changed depending on the variable, this needs to be given user input (check lessons)
for line in open(output_dir +"/listGenomes.txt"):
#for line in open("Documents/Latest/VM_files/3_clades_no_labels/phylogeny_list.txt"): 
    fields = line.rstrip("\n").split("\t") #split into fields
    strain = fields[0]
    name = "NA"
    group = "group"
    oi = "no"
    
    
    if group_option == "yes":
        for group_line in open(working_directory + "/group_list.txt"):
            group_fields = group_line.rstrip("\n").split("\t")
            group_strain = group_fields[0]
            grouping= group_fields[1]
            if group_strain == strain:
                group = grouping
                
    if oi_option == "yes":
        for oi_line in open(output_dir + "/genomesAdded.txt"):
            oi_fields = oi_line.rstrip("\n").split("\t")
            #print(oi_fields)
            oi_strain = oi_fields[0]
            if oi_strain == strain:
                oi = oi_option        
                
    if name_option == "yes":
        for name_line in open(working_directory + "/new_genome_names.txt"):
            name_fields = name_line.rstrip("\n").split("\t")
            #print(oi_fields)
            name_strain = name_fields[0]
            new_name = name_fields[1]
            if name_strain == strain:
                name = new_name        
    
    #print(strain + "\t" + name + "\t" + group + "\t" + oi + "\t" + "\n")
    
    new_phylogeny_list.write(strain + "\t" + name + "\t" + group + "\t" + oi + "\t" + "\n")

