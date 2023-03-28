#!/usr/bin/env python
# coding: utf-8

# In[53]:


import os
import re
import sys

path = sys.argv[1] 

new_phylogeny_list = open(path + "/new_phylogeny_list.txt", "w")
new_genome_list = open(path + "/new_genome_list.txt","w")

#new_phylogeny_list = open("new_phylogeny_list.txt", "w")
#new_genome_list = open("/new_genome_list.txt","w")

#Path will need to be changed depending on the variable, this needs to be given user input (check lessons)
for line in open(path + "/phylogeny_list.txt"):
#for line in open("Documents/Latest/VM_files/3_clades_no_labels/phylogeny_list.txt"): 
    fields = line.rstrip("\n").split("\t") #split into fields
    clade = fields[0]
    name = fields[1]
    if clade not in ('I','II','III','IV','V'):
        new_phylogeny_list.write(name + "\t" + clade + "\n")
        new_genome_list.write(name + "\n")
    else:
        print("Removed " + name + " as it was designated a clade " + clade + " phylogroup.")
