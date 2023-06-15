#!/usr/bin/env python
# coding: utf-8

# In[86]:


import os
import re
import sys

#Below is the method for splitting the ID files into the file_name, and what we want to replace the bit that matches it with
# - the new_name
path = sys.argv[1]
rename_genomes = sys.argv[2]
working_directory = sys.argv[3]
grouping = sys.argv[4]
rename_genes = sys.argv[5]
genes_oi = sys.argv[6]
genomes_oi = sys.argv[7]

new_file = open(path+ "/newRscript.r", "w")

for line in open(working_directory+"/generate_plot.r"): #take file line by line
    
    fields = line.rstrip("\n")
    
    if re.search('args\[1]', fields):
        fields = re.sub("args\[1],\"","\""+ path ,fields)
        
    if re.search("args\[2]", fields):
        fields = re.sub("args\[2]","\""+ rename_genomes +"\"", fields)
        
    if re.search("args\[3]", fields):
        fields = re.sub("args\[3]","\""+ grouping +"\"", fields)
    if re.search("args\[4]", fields):
    	fields = re.sub("args\[4]","\""+ rename_genes + "\"", fields)
    if re.search("args\[5]", fields):
    	fields = re.sub("args\[5]","\""+ genes_oi + "\"", fields)
    if re.search("args\[6]", fields):
    	fields = re.sub("args\[6]","\""+ genomes_oi + "\"", fields)    
    	   
    new_file.write(fields + "\n")
