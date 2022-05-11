#!/usr/bin/env python
# coding: utf-8

# In[86]:


import os
import re
import sys

#Below is the method for splitting the ID files into the file_name, and what we want to replace the bit that matches it with
# - the new_name
path = sys.argv[1]
clades = sys.argv[2]
labels = sys.argv[3]
working_directory = sys.argv[4]

new_file = open(path+ "/newRscript.r", "w")

for line in open(working_directory+"/generatePlot3.r"): #take file line by line
    
    fields = line.rstrip("\n")
    
    if re.search('args\[1]', fields):
        fields = re.sub("args\[1],\"","\""+path,fields)
        
    if re.search("args\[2]", fields):
        fields = re.sub("args\[2]","\""+clades+"\"",fields)
        
    if re.search("args\[3]", fields):
        fields = re.sub("args\[3]","\""+ labels+"\"",fields)
    
    new_file.write(fields + "\n")

